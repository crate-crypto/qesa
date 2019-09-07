#![allow(non_snake_case)]
use crate::inner_product_argument::qesa_inner;
use crate::matrix::*;
use curve25519_dalek::{ristretto::RistrettoPoint, scalar::Scalar};
use merlin::Transcript;

// QESA_ZK is a glorified wrpapper around Qesa_Inner
pub struct Zk {
    inner: qesa_inner::Inner,
}

pub fn create(
    transcript: &mut Transcript,
    G_Vec: Vec<RistrettoPoint>,
    H_Vec: Vec<RistrettoPoint>,
    Q: RistrettoPoint,
    gamma_i: &BlockMatrix,
    w: Vec<Scalar>,
) -> Zk {
    
    // 1. Compute r'
    //    
    let mut rng = rand::thread_rng();
    let r_prime = vec![Scalar::random(&mut rng),Scalar::random(&mut rng)]; 

    // The paper states that we should generate C'_w in Qesa_Zk
    // Then pass it to Qesa_Inner.
    // However, it makes more sense to generate it in Qesa_Inner 

    let qesa_inner_proof = qesa_inner::create(transcript, G_Vec, H_Vec, &Q, gamma_i, w, r_prime);
    Zk { inner: qesa_inner_proof }
}

impl Zk {
    pub fn verify(
        &self,
        transcript: &mut Transcript,
        G_Vec: Vec<RistrettoPoint>,
        H_Vec: Vec<RistrettoPoint>,
        Q: &RistrettoPoint,
        gamma_i: &BlockMatrix,
    ) -> bool {
        self.inner.verify(transcript, G_Vec, H_Vec, Q, gamma_i)
    }
}
