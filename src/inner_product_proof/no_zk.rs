#![allow(non_snake_case)]
use crate::transcript::TranscriptProtocol;
use curve25519_dalek::{
    ristretto::{CompressedRistretto, RistrettoPoint},
    scalar::Scalar,
    traits::VartimeMultiscalarMul,
};
use merlin::Transcript;
use std::iter;

/// NoZK is an optimisation over the bulletproofs IPA.
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
pub fn create(
    transcript: &mut Transcript,
    Q: &RistrettoPoint,
    C: RistrettoPoint,
    mut G_Vec: Vec<RistrettoPoint>,
    mut H_Vec: Vec<RistrettoPoint>,
    mut a_vec: Vec<Scalar>,
    mut b_vec: Vec<Scalar>,
) -> NoZK {
    let mut a = &mut a_vec[..];
    let mut b = &mut b_vec[..];
    let mut G = &mut G_Vec[..];
    let mut H = &mut H_Vec[..];

    let mut L_vec: Vec<CompressedRistretto> = Vec::new();
    let mut R_vec: Vec<CompressedRistretto> = Vec::new();

    let mut n = G.len();

    // transcript.append_u64(b"should be a power of two", n as u64);
    // We add the compressed point to the transcript, because we need to generate alpha
    // If this is not done, then the prover always will be able to predict what the first challenge will be 
    transcript.append_message(b"c", C.compress().as_bytes());

    let alpha = transcript.challenge_scalar(b"alpha");
    let Q = alpha.invert() * Q;

    while n != 1 {
        n = n / 2;

        let (a_L, a_R) = a.split_at_mut(n);
        let (b_L, b_R) = b.split_at_mut(n);
        let (G_L, G_R) = G.split_at_mut(n);
        let (H_L, H_R) = H.split_at_mut(n);

        let v_minus_one = inner_product(a_R, b_L);
        let v_plus_one = inner_product(a_L, b_R);

        let L = RistrettoPoint::vartime_multiscalar_mul(
            a_L.iter().chain(b_R.iter()).chain(iter::once(&v_plus_one)),
            G_R.iter().chain(H_L.iter()).chain(iter::once(&Q)),
        )
        .compress();

        let R = RistrettoPoint::vartime_multiscalar_mul(
            a_R.iter().chain(b_L.iter()).chain(iter::once(&v_minus_one)),
            G_L.iter().chain(H_R.iter().chain(iter::once(&Q))),
        )
        .compress();

        L_vec.push(L);
        R_vec.push(R);

        transcript.append_message(b"u_minus_one", L.as_bytes());
        transcript.append_message(b"u_plus_one", R.as_bytes());

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

pub fn verify(
    transcript: &mut Transcript,
    proof: NoZK,
    C: RistrettoPoint,
    Q: &RistrettoPoint,
    t: Scalar,
    G_Vec: &[RistrettoPoint],
    H_Vec: &[RistrettoPoint],
) {

    let mut G = G_Vec.to_owned();
    let mut H = H_Vec.to_owned();

    let Ls: Vec<RistrettoPoint> = proof
        .L_vec
        .iter()
        .map(|L| L.decompress().unwrap())
        .collect();
    let Rs: Vec<RistrettoPoint> = proof
        .R_vec
        .iter()
        .map(|R| R.decompress().unwrap())
        .collect();

    let mut n = 1 << Ls.len();

    transcript.append_message(b"c", C.compress().as_bytes());

    let alpha = transcript.challenge_scalar(b"alpha");

    let Q = alpha.invert() * Q;

    let mut C = C - (alpha - Scalar::one()) * t * Q;

    let challenges = generate_challenges(&proof, transcript);

    for ((L, R), challenge) in Ls.iter().zip(Rs.iter()).zip(challenges.iter()) {
        let challenge_sq = challenge * challenge;
        C = challenge_sq * L + challenge * C + R;
    }

    for challenge in challenges.iter() {
        n = n / 2;

        let (G_L, G_R) = G.split_at_mut(n);
        let (H_L, H_R) = H.split_at_mut(n);

        G = vector_vector_add(G_L, &mut vector_scalar_mul(G_R, challenge));
        H = vector_vector_add(&mut vector_scalar_mul(H_L, challenge), H_R);
    }

    let exp_C = G[0] * proof.a + H[0] * proof.b + (proof.a * proof.b) * Q;

    if exp_C != C {
        panic!("proof failed. expected Commitment {:?} got {:?}", exp_C.compress().as_bytes(),C.compress().as_bytes());
    }
}

fn generate_challenges(proof: &NoZK, transcript: &mut Transcript) -> Vec<Scalar> {
    let mut challenges: Vec<Scalar> = Vec::new();

    for (L, R) in proof.L_vec.iter().zip(proof.R_vec.iter()) {

        transcript.append_message(b"u_minus_one", L.as_bytes());
        transcript.append_message(b"u_plus_one", R.as_bytes());

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
    use sha3::Sha3_512;
    use std::iter;
    #[test]
    fn test_create_proof() {
        let mut rng = rand::thread_rng();

        let (a, b, t) = helper_dot_product(4);

        let G: Vec<RistrettoPoint> = (0..4).map(|_| RistrettoPoint::random(&mut rng)).collect();
        let H: Vec<RistrettoPoint> = (0..4).map(|_| RistrettoPoint::random(&mut rng)).collect();

        let mut transcript = Transcript::new(b"ip_no_zk");

        let mut Q = RistrettoPoint::hash_from_bytes::<Sha3_512>(b"test");
        let mut C = RistrettoPoint::vartime_multiscalar_mul(
            a.iter().chain(b.iter()).chain(iter::once(&t)),
            G.iter().chain(H.iter()).chain(iter::once(&Q)),
        );
        let proof = create(&mut transcript, &Q, C, G.clone(), H.clone(), a, b);

        transcript = Transcript::new(b"ip_no_zk");

        verify(&mut transcript, proof, C, &Q, t, &G, &H);
    }
}

fn helper_dot_product(n: usize) -> (Vec<Scalar>, Vec<Scalar>, Scalar) {
    let mut rng = rand::thread_rng();

    let a: Vec<Scalar> = (0..4).map(|_| Scalar::random(&mut rng)).collect();
    let b: Vec<Scalar> = (0..4).map(|_| Scalar::random(&mut rng)).collect();

    let t: Scalar = inner_product(&a, &b);
    (a, b, t)
}

fn inner_product(a: &[Scalar], b: &[Scalar]) -> Scalar {
    a.iter().zip(b.iter()).map(|(a, b)| a * b).sum()
}
