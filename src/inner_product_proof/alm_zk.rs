#![allow(non_snake_case)]
use crate::inner_product_proof::{gramschmidt::*, no_zk};
use crate::transcript::TranscriptProtocol;
use curve25519_dalek::{
    ristretto::{CompressedRistretto, RistrettoPoint},
    scalar::Scalar,
    traits::VartimeMultiscalarMul,
};
use merlin::Transcript;

#[derive(Clone)]
pub struct AlmZK {
    C_r: CompressedRistretto,
    NoZK: no_zk::NoZK,
}

pub fn create(
    transcript: &mut Transcript,
    G_Vec: Vec<RistrettoPoint>,
    H_Vec: Vec<RistrettoPoint>,
    Q: &RistrettoPoint,
    C_w: RistrettoPoint,
    a_Vec: Vec<Scalar>,
    b_Vec: Vec<Scalar>,
    t: Scalar,
) -> AlmZK {
    //1. Compute r' and r''
    //
    let r_prime = sample_gram_schmidt(&b_Vec);

    let r_prime_prime = sample_gram_schmidt_twice(&a_Vec, &r_prime);

    //2. Compute Commitment (C_r) to r' and r''
    //
    let C_r = RistrettoPoint::vartime_multiscalar_mul(
        r_prime.iter().chain(r_prime_prime.iter()),
        G_Vec.iter().chain(H_Vec.iter()),
    );
    transcript.append_message(b"C_r", C_r.compress().as_bytes());

    //3. Compute challenge based on C_r
    //
    let beta = transcript.challenge_scalar(b"beta");

    // 4. Compute masked version of `a` and `b` using r' and r''
    //
    let masked_a_Vec: Vec<Scalar> = a_Vec
        .iter()
        .zip(r_prime.iter())
        .map(|(a, r)| beta * a + r)
        .collect();
    let masked_b_Vec: Vec<Scalar> = b_Vec
        .iter()
        .zip(r_prime_prime.iter())
        .map(|(b, r)| beta * b + r)
        .collect();

    // Since we are no longer proving knowledge of `a` and `b`; it is now masked `a` and masked `b`.
    // We need to modify the commitment `P` and the result of the inner product argument `t` to match this
    //
    // 5. Compute modified `t` value
    //
    let beta_sq_t = beta * beta * t;

    //6.  Compute modified commitment
    //
    let P = ((beta * C_w) + C_r + (beta_sq_t * Q)).compress();
    transcript.append_message(b"P", P.as_bytes());

    //6. Call IPA_NoZK as a sub-protocol
    //
    let no_zk_proof = no_zk::create(transcript, G_Vec, H_Vec, Q, masked_a_Vec, masked_b_Vec);
    AlmZK {
        C_r: C_r.compress(),
        NoZK: no_zk_proof,
    }
}

impl AlmZK {
    pub fn verify(
        &self,
        transcript: &mut Transcript,
        G_Vec: &[RistrettoPoint],
        H_Vec: &[RistrettoPoint],
        Q: &RistrettoPoint,
        n: usize,
        C_w: RistrettoPoint,
        t: Scalar,
    ) -> bool {
        //1. Decompress the commitment
        //
        let C_r = self.C_r.decompress().unwrap();
        transcript.append_message(b"C_r", self.C_r.as_bytes());

        //3. Compute challenge based on C_r
        //
        let beta = transcript.challenge_scalar(b"beta");

        // 4. Compute modified `t` value
        //
        let beta_sq_t = beta * beta * t;

        // 5. Compute modified commitment
        //
        let P = (beta * C_w) + C_r + (beta_sq_t * Q);
        transcript.append_message(b"P", P.compress().as_bytes());

        //6. Call IPA_NoZk as a sub-protocol
        //
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
