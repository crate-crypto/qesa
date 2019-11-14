use crate::lmpa::no_zk;
use crate::matrix::*;
use crate::transcript::TranscriptProtocol;
use curve25519_dalek::{ristretto::RistrettoPoint, scalar::Scalar};
use merlin::Transcript;

// Linear Map Pre-Image Argument without zero knowledge
struct SimpleZK {
    no_zk: no_zk::NoZK,
    a: Vec<RistrettoPoint>,
}

fn create(
    transcript: &mut Transcript,
    mut A: Vec<Vec<RistrettoPoint>>,
    G_Vec: Vec<RistrettoPoint>,
    H_Vec: Vec<RistrettoPoint>,
    w: Vec<Scalar>,
) -> SimpleZK {
    let mut n = G_Vec.len();

    let mut rng = rand::thread_rng();
    let r: Vec<Scalar> = (0..n).map(|_| Scalar::random(&mut rng)).collect();

    let a = matrixpoint_vector_mul(&A, &r);
    for element in a.iter() {
        transcript.append_message(b"a", element.compress().as_bytes());
    }

    let beta = transcript.challenge_scalar(b"beta");

    let w_prime = w
        .iter()
        .zip(r.iter())
        .map(|(witness, random)| beta * witness + random)
        .collect();

    let proof = no_zk::create(transcript, A, G_Vec, H_Vec, w_prime);

    SimpleZK { a: a, no_zk: proof }
}

impl SimpleZK {
    pub fn verify(
        &self,
        transcript: &mut Transcript,
        mut A: Vec<Vec<RistrettoPoint>>,
        G_Vec: Vec<RistrettoPoint>,
        H_Vec: Vec<RistrettoPoint>,
        mut n: usize,
        mut t: Vec<RistrettoPoint>,
    ) -> bool {
        assert_eq!(n, G_Vec.len());

        for element in self.a.iter() {
            transcript.append_message(b"a", element.compress().as_bytes());
        }
        let beta = transcript.challenge_scalar(b"beta");

        t = t
            .iter()
            .zip(self.a.iter())
            .map(|(t_i, a_i)| beta * t_i + a_i)
            .collect();

        self.no_zk.verify(transcript, A, G_Vec, H_Vec, n, t)
    }
}

use sha3::Sha3_512;

#[test]
fn test_lmpa_simple_zk_create_verify() {
    let n = 4;
    let mut rng = rand::thread_rng();

    let G: Vec<RistrettoPoint> = (0..n).map(|_| RistrettoPoint::random(&mut rng)).collect();
    let H: Vec<RistrettoPoint> = (0..n).map(|_| RistrettoPoint::random(&mut rng)).collect();

    let mut prover_transcript = Transcript::new(b"lmpa_simple_zk");

    // t = Aw
    let A = rand_matrix(n);
    let w: Vec<Scalar> = (0..n).map(|_| Scalar::random(&mut rng)).collect();
    let t = matrixpoint_vector_mul(&A, &w);

    let proof = create(&mut prover_transcript, A.clone(), G.clone(), H.clone(), w);

    let mut verifier_transcript = Transcript::new(b"lmpa_simple_zk");
    assert!(proof.verify(&mut verifier_transcript, A, G, H, n, t));
}

fn rand_matrix(n: usize) -> Vec<Vec<RistrettoPoint>> {
    let mut rng = rand::thread_rng();

    let mut matrix: Vec<Vec<RistrettoPoint>> = Vec::new();

    for _ in 0..n {
        let a: Vec<RistrettoPoint> = (0..n).map(|_| RistrettoPoint::random(&mut rng)).collect();

        matrix.push(a);
    }

    matrix
}
