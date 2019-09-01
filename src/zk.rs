use crate::inner;
use crate::matrix::*;
use curve25519_dalek::{ristretto::RistrettoPoint, scalar::Scalar};
use merlin::Transcript;

pub struct Zk {
    inner: inner::Inner,
}

pub fn create(
    transcript: &mut Transcript,
    G_Vec: Vec<RistrettoPoint>,
    H_Vec: Vec<RistrettoPoint>,
    Q: RistrettoPoint,
    gamma_i: &block_matrix,
    w: Vec<Scalar>,
    r_prime: Vec<Scalar>,
) -> Zk {
    let proof = inner::create(transcript, G_Vec, H_Vec, &Q, gamma_i, w, r_prime);

    Zk { inner: proof }
}

impl Zk {
    pub fn verify(
        &self,
        transcript: &mut Transcript,
        G_Vec: Vec<RistrettoPoint>,
        H_Vec: Vec<RistrettoPoint>,
        Q: &RistrettoPoint,
        gamma_i: &block_matrix,
    ) -> bool {
        self.inner.verify(transcript, G_Vec, H_Vec, Q, gamma_i)
    }
}
