use crate::matrix::*;
use crate::transcript::TranscriptProtocol;
use curve25519_dalek::{ristretto::RistrettoPoint, scalar::Scalar};
use merlin::Transcript;
// Linear Map Pre-Image Argument without zero knowledge
pub struct NoZK {
    L_Vec: Vec<Vec<RistrettoPoint>>,
    R_Vec: Vec<Vec<RistrettoPoint>>,
    w: Scalar,
}
pub fn create(
    transcript: &mut Transcript,
    mut A: Vec<Vec<RistrettoPoint>>,
    G_Vec: Vec<RistrettoPoint>,
    H_Vec: Vec<RistrettoPoint>,
    Q: &RistrettoPoint,
    mut w: Vec<Scalar>,
) -> NoZK {
    let mut n = G_Vec.len();

    let mut L_Vec: Vec<Vec<RistrettoPoint>> = Vec::new();
    let mut R_Vec: Vec<Vec<RistrettoPoint>> = Vec::new();

    while n != 1 {
        n = n / 2;
        let (A_0_1): Vec<_> = A
            .iter()
            .map(|row| {
                let (left, right) = row.split_at(row.len() / 2);

                (left.to_owned(), right.to_owned())
            })
            .collect();

        let (A_0): Vec<_> = A_0_1.iter().map(|(vec_left, _)| vec_left.clone()).collect();
        let (A_1): Vec<_> = A_0_1
            .iter()
            .map(|(_, vec_right)| vec_right.clone())
            .collect();

        let (w_0, w_1) = w.split_at(w.len() / 2);

        let L = matrixpoint_vector_mul(&A_0, w_1);
        for u_L in L.iter() {
            transcript.append_message(b"L", u_L.compress().as_bytes());
        }
        L_Vec.push(L);

        let R = matrixpoint_vector_mul(&A_1, w_0);
        for u_R in R.iter() {
            transcript.append_message(b"R", u_R.compress().as_bytes());
        }
        R_Vec.push(R);

        let challenge = transcript.challenge_scalar(b"challenge");

        let challenge_matrix = matrixpoint_scalar_mul(&A_1, &challenge);
        A = matrixpoint_matrix_add(&A_0, &challenge_matrix);

        w = w_0
            .iter()
            .zip(w_1.iter())
            .map(|(a, b)| a * challenge + b)
            .collect();
    }

    NoZK {
        L_Vec: L_Vec,
        R_Vec: R_Vec,
        w: w[0],
    }
}

impl NoZK {
    pub fn verify(
        &self,
        transcript: &mut Transcript,
        mut A: Vec<Vec<RistrettoPoint>>,
        G_Vec: Vec<RistrettoPoint>,
        H_Vec: Vec<RistrettoPoint>,
        Q: &RistrettoPoint,
        mut n: usize,
        mut t: Vec<RistrettoPoint>,
    ) -> bool {
        assert_eq!(n, G_Vec.len());

        let mut challenge_index = 0;

        while n != 1 {
            n = n / 2;

            let (A_0_1): Vec<_> = A
                .iter()
                .map(|row| {
                    let (left, right) = row.split_at(row.len() / 2);

                    (left.to_owned(), right.to_owned())
                })
                .collect();

            let (A_0): Vec<_> = A_0_1.iter().map(|(vec_left, _)| vec_left.clone()).collect();
            let (A_1): Vec<_> = A_0_1
                .iter()
                .map(|(_, vec_right)| vec_right.clone())
                .collect();

            // Add challenges to transcript
            for u_L in self.L_Vec[challenge_index].iter() {
                transcript.append_message(b"L", u_L.compress().as_bytes());
            }
            for u_R in self.R_Vec[challenge_index].iter() {
                transcript.append_message(b"R", u_R.compress().as_bytes());
            }

            let challenge = transcript.challenge_scalar(b"challenge");
            let challenge_matrix = matrixpoint_scalar_mul(&A_1, &challenge);
            A = matrixpoint_matrix_add(&A_0, &challenge_matrix);

            t = self.L_Vec[challenge_index]
                .iter()
                .zip(self.R_Vec[challenge_index].iter())
                .zip(t.iter())
                .map(|((l, r), t)| l + (challenge * t) + (challenge * challenge * r))
                .collect();

            challenge_index += 1;
        }

        let base_witnes = vec![self.w];

        let expected_t = matrixpoint_vector_mul(&A, &base_witnes);

        expected_t[0] == t[0]
    }
}
use sha3::Sha3_512;

#[test]
fn test_lmpa_no_zk_create_verify() {
    let n = 4;
    let mut rng = rand::thread_rng();

    let G: Vec<RistrettoPoint> = (0..n).map(|_| RistrettoPoint::random(&mut rng)).collect();
    let H: Vec<RistrettoPoint> = (0..n).map(|_| RistrettoPoint::random(&mut rng)).collect();
    let Q = RistrettoPoint::hash_from_bytes::<Sha3_512>(b"test point");

    let mut prover_transcript = Transcript::new(b"lmpa_no_zk");

    // t = Aw
    let A = rand_matrix(n);
    let w: Vec<Scalar> = (0..n).map(|_| Scalar::random(&mut rng)).collect();
    let t = matrixpoint_vector_mul(&A, &w);

    let proof = create(
        &mut prover_transcript,
        A.clone(),
        G.clone(),
        H.clone(),
        &Q,
        w,
    );

    let mut verifier_transcript = Transcript::new(b"lmpa_no_zk");
    assert!(proof.verify(&mut verifier_transcript, A, G, H, &Q, n, t));
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
