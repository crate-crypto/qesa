#![allow(non_snake_case)]
use crate::ipa::alm_zk;
use crate::math_utils::vandemonde_challenge;
use crate::matrix::*;
use crate::transcript::TranscriptProtocol;
use curve25519_dalek::{
    ristretto::{CompressedRistretto, RistrettoPoint},
    scalar::Scalar,
    traits::VartimeMultiscalarMul,
};
use merlin::Transcript;

#[derive(Clone)]
pub struct Inner {
    pub(crate) alm_zk: alm_zk::AlmZK,
    pub(crate) c_prime_w: CompressedRistretto,
    pub(crate) c_prime_prime_w: CompressedRistretto,
}

pub fn create(
    transcript: &mut Transcript,
    mut G_Vec: Vec<RistrettoPoint>,
    H_Vec: Vec<RistrettoPoint>,
    Q: &RistrettoPoint,
    gamma_i: &BlockMatrix,
    w: Vec<Scalar>,
    r_prime: Vec<Scalar>,
) -> Inner {
    let n = G_Vec.len();
    assert_eq!(w.len(), n - 2);
    assert_eq!(r_prime.len(), 2);

    // 1. Compute w' (extended witness)
    //
    let w_prime = [&w[..], &r_prime[..]].concat();

    //2. Compute commitment (C'_w) to w'
    //
    let c_prime_w = RistrettoPoint::vartime_multiscalar_mul(w_prime.iter(), G_Vec.iter());
    transcript.append_message(b"c_prime_w", c_prime_w.compress().as_bytes());

    //3. Compute challenges based on C'_w
    //
    let x_challenges: Vec<Scalar> = vandemonde_challenge(transcript.challenge_scalar(b"x"), n);
    assert_eq!(x_challenges.len(), n);

    // 4. batch each matrix into a single matrix using the challenges
    //
    // Given `k` matrices of size `n-2 x n-2`
    // We use the challenges to batch all of them into one `n-2 x n-2` matrix
    let gamma = gamma_i.block_matrix_batch(&x_challenges);

    //5. Fix the first element of the witness to be one (see paper)
    //
    let beta = x_challenges[1];
    G_Vec[0] = G_Vec[0] * beta.invert();

    //6. Compute r'' by rotating the r' by 90 degrees
    //
    let r_prime_prime = vec![-r_prime[1], r_prime[0]];

    // 7. Compute w'' (extended witness)
    //
    // Compute gamma * w
    let gamma_w = matrix_vector_mul(&gamma, &w);
    //
    let w_prime_prime = [&gamma_w[..], &r_prime_prime[..]].concat();

    //8. Compute commitment (C''_w) to w''
    //
    let c_prime_prime_w =
        RistrettoPoint::vartime_multiscalar_mul(w_prime_prime.iter(), H_Vec.iter());
    transcript.append_message(b"c_prime_prime_w", c_prime_prime_w.compress().as_bytes());

    //9. Compute `s` and `b` challenges
    //
    // `s` must be of size `n-2` and b must be of size `2`
    let s_challenges = vandemonde_challenge(transcript.challenge_scalar(b"s"), n - 2);
    let b_challenges = vandemonde_challenge(transcript.challenge_scalar(b"b"), 2);

    //10. compute s' by concatenating `s` and `b`
    //
    let s_prime = [&s_challenges[..], &b_challenges[..]].concat();
    assert_eq!(s_prime.len(), n);

    // 11. Compute `a` , ``b and `t`
    //
    // Recall that <a,b> = t
    //
    let gamma_prime = compute_gamma_prime(&gamma, n);
    let gamma_prime_t = matrix_transpose(&gamma_prime);
    let gamma_prime_t_s_prime = matrix_vector_mul(&gamma_prime_t, &s_prime);

    let a = row_row_sub(&w_prime, &s_prime);
    let b = row_row_add(&w_prime_prime, &gamma_prime_t_s_prime);
    let t = crate::math_utils::inner_product(&a, &b);

    // 12. Compute commitment (C_w) to witness
    //
    let C_w = RistrettoPoint::vartime_multiscalar_mul(
        a.iter().chain(b.iter()),
        G_Vec.iter().chain(H_Vec.iter()),
    );
    transcript.append_message(b"c_w", C_w.compress().as_bytes());

    //13. Call IPA_AlmZK as a sub-protocol
    //
    let alm_zk_proof = alm_zk::create(transcript, G_Vec, H_Vec, Q, C_w, a, b, t);
    Inner {
        alm_zk: alm_zk_proof,
        c_prime_w: c_prime_w.compress(),
        c_prime_prime_w: c_prime_prime_w.compress(),
    }
}

impl Inner {
    pub fn verify(
        &self,
        transcript: &mut Transcript,
        mut G_Vec: Vec<RistrettoPoint>,
        H_Vec: Vec<RistrettoPoint>,
        Q: &RistrettoPoint,
        gamma_i: &BlockMatrix,
    ) -> bool {
        let n = H_Vec.len();

        //1. Decompress the commitment (C'_w)
        //
        transcript.append_message(b"c_prime_w", self.c_prime_w.as_bytes());
        let c_prime_w = self.c_prime_w.decompress().unwrap();

        //2. Compute challenges based on C'_w
        //
        let x_challenges: Vec<Scalar> = vandemonde_challenge(transcript.challenge_scalar(b"x"), n);

        // 3. batch each matrix into a single matrix using the challenges
        //
        let gamma = gamma_i.block_matrix_batch(&x_challenges);

        //4. Fix the first element of the witness to be one (see paper)
        //
        let beta = x_challenges[1];
        G_Vec[0] = G_Vec[0] * beta.invert();

        //5. Decompress the commitment (C''_w)
        //
        let c_prime_prime_w = self.c_prime_prime_w.decompress().unwrap();
        let c_prime_w = c_prime_w - (beta - Scalar::one()) * G_Vec[0];
        transcript.append_message(b"c_prime_prime_w", c_prime_prime_w.compress().as_bytes());

        //6. Compute `s` and `b` challenges
        //
        let s_challenges = vandemonde_challenge(transcript.challenge_scalar(b"s"), n - 2);
        let b_challenges = vandemonde_challenge(transcript.challenge_scalar(b"b"), 2);

        //7. compute s' by concatenating `s` and `b`
        //
        let s_prime = [&s_challenges[..], &b_challenges[..]].concat();
        assert_eq!(s_prime.len(), n);

        // 8. Compute `t`
        //
        let gamma_prime = compute_gamma_prime(&gamma, n);

        let gamma_prime_t = matrix_transpose(&gamma_prime);
        let gamma_prime_t_s_prime = matrix_vector_mul(&gamma_prime_t, &s_prime);

        let gamma_t = matrix_transpose(&gamma);
        let gamma_t_s = matrix_vector_mul(&gamma_t, &s_challenges);

        let t = -crate::math_utils::inner_product(&s_challenges, &gamma_t_s);

        // 9. Compute commitment (C_w) to witness
        //
        let mut C_w = c_prime_w
            + c_prime_prime_w
            + RistrettoPoint::vartime_multiscalar_mul(gamma_prime_t_s_prime.iter(), H_Vec.iter());
        C_w = C_w - RistrettoPoint::vartime_multiscalar_mul(s_prime.iter(), G_Vec.iter());
        transcript.append_message(b"c_w", C_w.compress().as_bytes());

        //10. Call IPA_AlmZK as a sub-protocol
        //
        self.alm_zk
            .verify(transcript, &G_Vec, &H_Vec, &Q, n, C_w, t)
    }
}

fn compute_gamma_prime(gamma: &Vec<Vec<Scalar>>, n: usize) -> Vec<Vec<Scalar>> {
    // Pad the gamma rows with zeroes at the end (iss#1)
    let mut gamma_prime: Vec<Vec<Scalar>> = gamma
        .iter()
        .map(|row| {
            let mut padded_row = row.clone();
            padded_row.push(Scalar::zero());
            padded_row.push(Scalar::zero());
            padded_row
        })
        .collect();

    let mut row_n_minus_two = vec![Scalar::zero(); n - 2];
    row_n_minus_two.push(Scalar::zero());
    row_n_minus_two.push(-Scalar::from(1 as u8));
    gamma_prime.push(row_n_minus_two);

    let mut row_n_minus_one = vec![Scalar::zero(); n - 2];
    row_n_minus_one.push(Scalar::one());
    row_n_minus_one.push(Scalar::zero());

    gamma_prime.push(row_n_minus_one);

    gamma_prime
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math_utils::inner_product;
    use sha3::Sha3_512;
    #[test]
    fn test_create_qesa_inner_proof() {
        let mut rng = rand::thread_rng();

        let n = 4;

        let (witness, matrix) = helper_create_solutions(n - 2, 2);

        let G: Vec<RistrettoPoint> = (0..n).map(|_| RistrettoPoint::random(&mut rng)).collect();
        let H: Vec<RistrettoPoint> = (0..n).map(|_| RistrettoPoint::random(&mut rng)).collect();
        let Q = RistrettoPoint::hash_from_bytes::<Sha3_512>(b"test point");

        let r_prime: Vec<Scalar> = (0..2).map(|_| Scalar::random(&mut rng)).collect();

        let mut prover_transcript = Transcript::new(b"qesa_inner");

        let proof = create(
            &mut prover_transcript,
            G.clone(),
            H.clone(),
            &Q,
            &matrix,
            witness,
            r_prime,
        );

        let mut verifier_transcript = Transcript::new(b"qesa_inner");
        assert!(proof.verify(&mut verifier_transcript, G, H, &Q, &matrix))
    }
    // Creates a system of quadratic equations with solutions
    // and a witness
    fn helper_create_solutions(n: usize, num_of_matrices: usize) -> (Vec<Scalar>, BlockMatrix) {
        let mut rng = rand::thread_rng();
        let mut witness: Vec<Scalar> = (0..n).map(|_| Scalar::random(&mut rng)).collect();
        witness[0] = Scalar::one();

        let mut bm = BlockMatrix::new();

        for _ in 0..num_of_matrices {
            let mut gamma_i: Vec<Vec<Scalar>> = Vec::new();
            for _ in 0..n {
                // Use gram schmidt to create suitable solutions for each system of eqns
                let x: Vec<Scalar> = (0..n).map(|_| Scalar::random(&mut rng)).collect();
                let row_of_eqns = crate::ipa::gramschmidt::orth(&witness, &x);
                gamma_i.push(row_of_eqns)
            }
            bm.push(gamma_i);
        }

        // For now, we only use one set of system of equations
        (witness, bm)
    }

    #[test]
    fn test_helper_solutions() {
        let n = 4;

        let (witness, matrix) = helper_create_solutions(n - 2, 2);

        // Check that <w, gamma_i * w> = 0 for all i
        for gamma_i in matrix.block.iter() {
            let gamma_w = matrix_vector_mul(&gamma_i, &witness);
            assert_eq!(Scalar::zero(), inner_product(&witness, &gamma_w))
        }
    }
}
