#![allow(non_snake_case)]

use curve25519_dalek::{
    ristretto::{CompressedRistretto, RistrettoPoint},
    scalar::Scalar,
    traits::VartimeMultiscalarMul,
};

use crate::inner_product_proof::gramschmidt::*;
use crate::inner_product_proof::no_zk;
use crate::transcript::TranscriptProtocol;
use merlin::Transcript;
use std::iter;

/*
Notes:
we take a vector from the kernel of witness_prime_prime; r_prime.
Therefore we do not introduce errors as <r_prime, witness_prime_prime> from def of kernel
AND
We then take the next random value r_prime_prime from the kernel of witness_prime intersect r_prime.
Therefore <r_prime_prime, w_prime> = 0 = <r_prime, r_prime_prime>

This allows us to mask the witness without introducing errors into the protcol.
Similar to bulletproof, we will send these masked values into the IPA_NO_ZK which is just a prove of knowledge.
*/

struct AlmZK {
    C_r: CompressedRistretto,
    NoZK: no_zk::NoZK,
}

fn create(
    transcript: &mut Transcript,
    mut G_Vec: Vec<RistrettoPoint>,
    mut H_Vec: Vec<RistrettoPoint>,
    Q: &RistrettoPoint,
    C_w: RistrettoPoint,
    a_Vec: Vec<Scalar>,
    b_Vec: Vec<Scalar>,
    t: Scalar,
) -> AlmZK {
    let r_prime = sample_gram_schmidt(&b_Vec);

    let r_prime_prime = sample_gram_schmidt_twice(&a_Vec, &r_prime);

    let C_r = RistrettoPoint::vartime_multiscalar_mul(
        r_prime.iter().chain(r_prime_prime.iter()),
        G_Vec.iter().chain(H_Vec.iter()),
    );

    transcript.append_message(b"C_r", C_r.compress().as_bytes());

    let beta = transcript.challenge_scalar(b"beta");

    let a_prime_Vec: Vec<Scalar> = a_Vec
        .iter()
        .zip(r_prime.iter())
        .map(|(a, r)| beta * a + r)
        .collect();

    let b_prime_Vec: Vec<Scalar> = b_Vec
        .iter()
        .zip(r_prime_prime.iter())
        .map(|(b, r)| beta * b + r)
        .collect();

    let P = ((beta * C_w) + C_r + (beta * beta * t * Q)).compress();

    transcript.append_message(b"P", P.as_bytes());

    let no_zk_proof = no_zk::create(transcript, G_Vec, H_Vec, Q, a_prime_Vec, b_prime_Vec);

    AlmZK {
        C_r: C_r.compress(),
        NoZK: no_zk_proof,
    }
}

impl AlmZK {
    fn verify(
        &self,
        transcript: &mut Transcript,
        G_Vec: &[RistrettoPoint],
        H_Vec: &[RistrettoPoint],
        Q: &RistrettoPoint,
        n: usize,
        C_w: RistrettoPoint,
        t: Scalar,
    ) -> bool {
        let C_r = self.C_r.decompress().unwrap();

        transcript.append_message(b"C_r", self.C_r.as_bytes());
        let beta = transcript.challenge_scalar(b"beta");

        let beta_sq_t = beta * beta * t;

        let P = (beta * C_w) + C_r + (beta_sq_t * Q);

        transcript.append_message(b"P", P.compress().as_bytes());

        self.NoZK
            .verify(transcript, G_Vec, H_Vec, Q, n, P, beta_sq_t)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math_utils::*;
    use sha3::Sha3_512;
    #[test]
    fn test_create_almzk_proof() {
        let mut rng = rand::thread_rng();

        let n = 4;

        let (a, b, t) = helper_dot_product(n);

        let G: Vec<RistrettoPoint> = (0..n).map(|_| RistrettoPoint::random(&mut rng)).collect();
        let H: Vec<RistrettoPoint> = (0..n).map(|_| RistrettoPoint::random(&mut rng)).collect();

        let mut transcript = Transcript::new(b"ip_alm_zk");

        let Q = RistrettoPoint::hash_from_bytes::<Sha3_512>(b"test point");

        let C_w = RistrettoPoint::vartime_multiscalar_mul(
            a.iter().chain(b.iter()),
            G.iter().chain(H.iter()),
        );

        let proof = super::create(&mut transcript, G.clone(), H.clone(), &Q, C_w, a, b, t);

        transcript = Transcript::new(b"ip_alm_zk");

        assert!(proof.verify(&mut transcript, &G, &H, &Q, n, C_w, t));
    }
}
