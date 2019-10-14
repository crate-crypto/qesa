#![allow(non_snake_case)]
use crate::math_utils::inner_product;
use crate::transcript::TranscriptProtocol;
use curve25519_dalek::{
    ristretto::{CompressedRistretto, RistrettoPoint},
    scalar::Scalar,
    traits::{IsIdentity, VartimeMultiscalarMul},
};
use merlin::Transcript;
use std::iter;
/// NoZK is an optimisation over the bulletproofs IPA.
#[derive(Clone)]
pub struct NoZK {
    // From the literature this would be u_{-1}
    pub(crate) L_vec: Vec<CompressedRistretto>,
    // From the literature this would be u_{+1}
    pub(crate) R_vec: Vec<CompressedRistretto>,
    // From the literature, this would be w'
    pub(crate) a: Scalar,
    // From the literature, this would be w''
    pub(crate) b: Scalar,
}

// a and b are the witness prime and witness prime prime respectively
// G_Vec, H_Vec and Q is g prime, g prime prime and Q of the crs respectively
// This implementation will intentionally try to mirror the dalek-cryptography implementation in design choices and variable naming
// making it easier to draw comparisons between the two and provide benchmarks
pub fn create(
    transcript: &mut Transcript,
    mut G_Vec: Vec<RistrettoPoint>,
    mut H_Vec: Vec<RistrettoPoint>,
    Q: &RistrettoPoint,
    mut a_vec: Vec<Scalar>,
    mut b_vec: Vec<Scalar>,
) -> NoZK {
    let mut a = &mut a_vec[..];
    let mut b = &mut b_vec[..];
    let mut G = &mut G_Vec[..];
    let mut H = &mut H_Vec[..];

    let mut n = G.len();

    // All of the input vectors must have the same length.
    assert_eq!(G.len(), n);
    assert_eq!(H.len(), n);
    assert_eq!(a.len(), n);
    assert_eq!(b.len(), n);

    // All of the input vectors must have a length that is a power of two.
    assert!(n.is_power_of_two());

    transcript.append_u64(b"n", n as u64);

    let lg_n = n.next_power_of_two().trailing_zeros() as usize;
    let mut L_vec: Vec<CompressedRistretto> = Vec::with_capacity(lg_n);
    let mut R_vec: Vec<CompressedRistretto> = Vec::with_capacity(lg_n);

    let alpha = transcript.challenge_scalar(b"alpha");
    let Q = alpha.invert() * Q;

    while n != 1 {
        n = n / 2;

        let (a_L, a_R) = a.split_at_mut(n);
        let (b_L, b_R) = b.split_at_mut(n);
        let (G_L, G_R) = G.split_at_mut(n);
        let (H_L, H_R) = H.split_at_mut(n);

        let c_R = inner_product(a_R, b_L);
        let c_L = inner_product(a_L, b_R);

        let L = RistrettoPoint::vartime_multiscalar_mul(
            a_L.iter().chain(b_R.iter()).chain(iter::once(&c_L)),
            G_R.iter().chain(H_L.iter()).chain(iter::once(&Q)),
        )
        .compress();

        let R = RistrettoPoint::vartime_multiscalar_mul(
            a_R.iter().chain(b_L.iter()).chain(iter::once(&c_R)),
            G_L.iter().chain(H_R.iter().chain(iter::once(&Q))),
        )
        .compress();

        L_vec.push(L);
        R_vec.push(R);

        transcript.append_message(b"L", L.as_bytes());
        transcript.append_message(b"R", R.as_bytes());

        let x_i: Scalar = transcript.challenge_scalar(b"x_i");
        for i in 0..n {
            a_L[i] = a_L[i] * x_i + a_R[i];
            b_L[i] = b_L[i] + x_i * b_R[i];
            G_L[i] = &G_L[i] + x_i * &G_R[i];
            H_L[i] = x_i * &H_L[i] + &H_R[i];
        }

        a = a_L;
        b = b_L;
        G = G_L;
        H = H_L;
    }

    NoZK {
        L_vec: L_vec,
        R_vec: R_vec,
        a: a[0],
        b: b[0],
    }
}
impl NoZK {
    pub fn verify(
        &self,
        transcript: &mut Transcript,
        G_Vec: &[RistrettoPoint],
        H_Vec: &[RistrettoPoint],
        Q: &RistrettoPoint,
        mut n: usize,
        mut P: RistrettoPoint,
        t: Scalar,
    ) -> bool {
        let mut G = G_Vec.to_owned();
        let mut H = H_Vec.to_owned();

        let Ls: Vec<RistrettoPoint> = self.L_vec.iter().map(|L| L.decompress().unwrap()).collect();
        let Rs: Vec<RistrettoPoint> = self.R_vec.iter().map(|R| R.decompress().unwrap()).collect();
        assert_eq!(n, 1 << Ls.len());

        transcript.append_u64(b"n", n as u64);

        let alpha = transcript.challenge_scalar(b"alpha");
        let Q = alpha.invert() * Q;

        P = P - (alpha - Scalar::one()) * t * Q;

        let challenges = generate_challenges(self, transcript);

        for ((L, R), challenge) in Ls.iter().zip(Rs.iter()).zip(challenges.iter()) {
            let challenge_sq = challenge * challenge;
            P = challenge_sq * L + challenge * P + R;
        }

        for challenge in challenges.iter() {
            n = n / 2;

            let (G_L, G_R) = G.split_at_mut(n);
            let (H_L, H_R) = H.split_at_mut(n);

            G = vector_vector_add(G_L, &mut vector_scalar_mul(G_R, challenge));
            H = vector_vector_add(&mut vector_scalar_mul(H_L, challenge), H_R);
        }

        let exp_P = G[0] * self.a + H[0] * self.b + (self.a * self.b) * Q;

        exp_P == P
    }
    pub fn verify_multiexp(
        &self,
        transcript: &mut Transcript,
        G_Vec: &[RistrettoPoint],
        H_Vec: &[RistrettoPoint],
        Q: &RistrettoPoint,
        n: usize,
        P: RistrettoPoint,
        t: Scalar,
    ) -> bool {

        // decode L,R
        let Ls: Vec<RistrettoPoint> = self.L_vec.iter().map(|L| L.decompress().unwrap()).collect();
        let Rs: Vec<RistrettoPoint> = self.R_vec.iter().map(|R| R.decompress().unwrap()).collect();
        assert_eq!(n, 1 << Ls.len());
        let mut n = 1 << Ls.len();

        // challenge data
        transcript.append_u64(b"n", n as u64);
        let alpha = transcript.challenge_scalar(b"alpha");
        let challenges = generate_challenges(self, transcript);

        let logn = challenges.len();

        // {g_i},{h_i}
        let mut g_i: Vec<Scalar> = Vec::with_capacity(n);
        let mut h_i: Vec<Scalar> = Vec::with_capacity(n);
        for x in 0..n {
            let mut i: usize = 1;
            let mut j: usize = 0;
            let mut g = self.a;
            let mut h = self.b;
            while i < n {
                if i & x != 0 {
                    g *= challenges[logn-j-1];
                }
                else {
                    h *= challenges[logn-j-1];
                }
                i <<= 1;
                j += 1;
            }
            g_i.push(g);
            h_i.push(h);
        }

        // {l_j},{r_j}
        let mut l_j: Vec<Scalar> = Vec::with_capacity(n);
        let mut r_j: Vec<Scalar> = Vec::with_capacity(n);
        let mut p = Scalar::one();
        for i in 0..logn {
            let mut l = -challenges[i]*challenges[i];
            let mut r = -Scalar::one();
            for j in (i+1)..logn {
                l *= challenges[j];
                r *= challenges[j];
            }
            l_j.push(l);
            r_j.push(r);
            p *= challenges[i];
        }

        // return value goes here
        let q = alpha.invert()*((alpha - Scalar::one())*t*p + self.a*self.b);


        let R = RistrettoPoint::vartime_multiscalar_mul(
            g_i.iter().chain(h_i.iter()).chain(l_j.iter()).chain(r_j.iter()).chain(iter::once(&q)).chain(iter::once(&-p)),
              G_Vec.iter().chain( H_Vec.iter()).chain( Ls.iter()).chain( Rs.iter()).chain(iter::once( Q)).chain(iter::once(&P))
        );

        R.is_identity()
    }
}

fn generate_challenges(proof: &NoZK, transcript: &mut Transcript) -> Vec<Scalar> {
    let mut challenges: Vec<Scalar> = Vec::new();

    for (L, R) in proof.L_vec.iter().zip(proof.R_vec.iter()) {
        transcript.append_message(b"L", L.as_bytes());
        transcript.append_message(b"R", R.as_bytes());

        let x_i = transcript.challenge_scalar(b"x_i");
        challenges.push(x_i);
    }

    challenges
}

fn vector_scalar_mul(vec_p: &mut [RistrettoPoint], challenge: &Scalar) -> Vec<RistrettoPoint> {
    vec_p.iter().map(|p| p * challenge).collect()
}
fn vector_vector_add(a: &mut [RistrettoPoint], b: &mut [RistrettoPoint]) -> Vec<RistrettoPoint> {
    a.iter().zip(b.iter()).map(|(a, b)| a + b).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math_utils::inner_product;
    use sha3::Sha3_512;
    use std::iter;
    #[test]
    fn test_create_nozk_proof() {
        let n = 4;

        let mut rng = rand::thread_rng();

        let a: Vec<Scalar> = (0..n).map(|_| Scalar::random(&mut rng)).collect();
        let b: Vec<Scalar> = (0..n).map(|_| Scalar::random(&mut rng)).collect();

        let t = inner_product(&a, &b);

        let G: Vec<RistrettoPoint> = (0..n).map(|_| RistrettoPoint::random(&mut rng)).collect();
        let H: Vec<RistrettoPoint> = (0..n).map(|_| RistrettoPoint::random(&mut rng)).collect();
        let Q = RistrettoPoint::hash_from_bytes::<Sha3_512>(b"test point");

        let mut prover_transcript = Transcript::new(b"ip_no_zk");

        let P = RistrettoPoint::vartime_multiscalar_mul(
            a.iter().chain(b.iter()).chain(iter::once(&t)),
            G.iter().chain(H.iter()).chain(iter::once(&Q)),
        );

        // We add the compressed point to the transcript, because we need some non-trivial input to generate alpha
        // If this is not done, then the prover always will be able to predict what the first challenge will be
        prover_transcript.append_message(b"P", P.compress().as_bytes());

        let proof = create(&mut prover_transcript, G.clone(), H.clone(), &Q, a, b);

        let mut verifier_transcript = Transcript::new(b"ip_no_zk");
        verifier_transcript.append_message(b"P", P.compress().as_bytes());

        assert!(proof.verify(&mut verifier_transcript, &G, &H, &Q, n, P, t));
    }

    extern crate test;
    use super::*;
    use test::Bencher;

    #[bench]
    fn bench_verify_multiexp(bnch: &mut Bencher) {
        let n = 256;

        let mut rng = rand::thread_rng();

        let a: Vec<Scalar> = (0..n).map(|_| Scalar::random(&mut rng)).collect();
        let b: Vec<Scalar> = (0..n).map(|_| Scalar::random(&mut rng)).collect();

        let t = inner_product(&a, &b);

        let G: Vec<RistrettoPoint> = (0..n).map(|_| RistrettoPoint::random(&mut rng)).collect();
        let H: Vec<RistrettoPoint> = (0..n).map(|_| RistrettoPoint::random(&mut rng)).collect();
        let Q = RistrettoPoint::hash_from_bytes::<Sha3_512>(b"test point");

        let mut prover_transcript = Transcript::new(b"ip_no_zk");

        let P = RistrettoPoint::vartime_multiscalar_mul(
            a.iter().chain(b.iter()).chain(iter::once(&t)),
            G.iter().chain(H.iter()).chain(iter::once(&Q)),
        );

        // We add the compressed point to the transcript, because we need some non-trivial input to generate alpha
        // If this is not done, then the prover always will be able to predict what the first challenge will be
        prover_transcript.append_message(b"P", P.compress().as_bytes());

        let proof = create(&mut prover_transcript, G.clone(), H.clone(), &Q, a, b);

        let mut verifier_transcript = Transcript::new(b"ip_no_zk");
        verifier_transcript.append_message(b"P", P.compress().as_bytes());

        bnch.iter(|| proof.verify_multiexp(&mut verifier_transcript, &G, &H, &Q, n, P, t));
    }
    #[bench]
    fn bench_verify(bnch: &mut Bencher) {
        let n = 256;

        let mut rng = rand::thread_rng();

        let a: Vec<Scalar> = (0..n).map(|_| Scalar::random(&mut rng)).collect();
        let b: Vec<Scalar> = (0..n).map(|_| Scalar::random(&mut rng)).collect();

        let t = inner_product(&a, &b);

        let G: Vec<RistrettoPoint> = (0..n).map(|_| RistrettoPoint::random(&mut rng)).collect();
        let H: Vec<RistrettoPoint> = (0..n).map(|_| RistrettoPoint::random(&mut rng)).collect();
        let Q = RistrettoPoint::hash_from_bytes::<Sha3_512>(b"test point");

        let mut prover_transcript = Transcript::new(b"ip_no_zk");

        let P = RistrettoPoint::vartime_multiscalar_mul(
            a.iter().chain(b.iter()).chain(iter::once(&t)),
            G.iter().chain(H.iter()).chain(iter::once(&Q)),
        );

        // We add the compressed point to the transcript, because we need some non-trivial input to generate alpha
        // If this is not done, then the prover always will be able to predict what the first challenge will be
        prover_transcript.append_message(b"P", P.compress().as_bytes());

        let proof = create(&mut prover_transcript, G.clone(), H.clone(), &Q, a, b);

        let mut verifier_transcript = Transcript::new(b"ip_no_zk");
        verifier_transcript.append_message(b"P", P.compress().as_bytes());

        bnch.iter(|| proof.verify(&mut verifier_transcript, &G, &H, &Q, n, P, t));
    }
}
