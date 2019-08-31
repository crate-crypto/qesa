#![allow(non_snake_case)]
use crate::matrix::*;
use crate::transcript::TranscriptProtocol;
use curve25519_dalek::{
    ristretto::{CompressedRistretto, RistrettoPoint},
    scalar::Scalar,
    traits::VartimeMultiscalarMul,
};
use merlin::Transcript;
/*

Notes: This is a sub-section of qesa_zk.
- First six bulletpoints on protocol 4.5 before we Run protocol IPA_ALM_ZK
- Module structure will most likely change
*/
// XXX: We need to formalise the way data is added to the transcript
// XXX: The code currently does not make use of the efficiency of sparse matrices

struct Inner {
    // Set of equations we want the solution for
//    matrice : Vec<Scalar>
}

#[allow(dead_code)]
fn create(
    transcript: &mut Transcript,
    mut G_Vec: Vec<RistrettoPoint>,
    H_Vec: Vec<RistrettoPoint>,
    Q: &RistrettoPoint,
    gamma_i: block_matrix, // set of system of quadratic equations
    w: Vec<Scalar>,
    r_prime: Vec<Scalar>,
) {
    let n = G_Vec.len();

    assert_eq!(w.len(), n - 2);
    assert_eq!(r_prime.len(), 2);

    let w_prime = [&w[..], &r_prime[..]].concat();
    assert_eq!(w_prime.len(), n);

    let mut c_prime_w = RistrettoPoint::vartime_multiscalar_mul(w_prime.iter(), G_Vec.iter());
    transcript.append_message(b"c_prime_w", c_prime_w.compress().as_bytes());

    let x_challenges: Vec<Scalar> = vandemonde_challenge(transcript.challenge_scalar(b"x"), n);
    assert_eq!(x_challenges.len(), n);

    let gamma = gamma_i.block_matrix_batch(&x_challenges);

    let beta = x_challenges[0];

    // Change the first generator in g'
    G_Vec[0] = G_Vec[0] * beta.invert();

    c_prime_w = c_prime_w - (beta - Scalar::one()) * G_Vec[0];
    transcript.append_message(b"c_prime_w", c_prime_w.compress().as_bytes());

    // r_prime_prime is a rotation of r_prime by 90 degrees
    let r_prime_prime = vec![-r_prime[1], r_prime[0]];

    // -- DEBUG r_prime_prime
    let mut R: Vec<Vec<Scalar>> = Vec::new();
    R.push(vec![Scalar::zero(), -Scalar::one()]);
    R.push(vec![Scalar::one(), Scalar::zero()]);

    let expected_r_prime_prime = matrix_vector_mul(&R, &r_prime);
    assert_eq!(expected_r_prime_prime, r_prime_prime);
    // -- End DEBUG r_prime_prime

    let gamma_w = matrix_vector_mul(&gamma, &w);

    // DEBUG - gamma_w

    let should_be_zero = crate::math_utils::inner_product(&w, &gamma_w);
    assert_eq!(Scalar::zero(), should_be_zero);

    // End DEBUG - gamma_w

    let w_prime_prime = [&gamma_w[..], &r_prime_prime[..]].concat();

    let c_prime_prime_w =
        RistrettoPoint::vartime_multiscalar_mul(w_prime_prime.iter(), H_Vec.iter());

    transcript.append_message(b"c_prime_prime_w", c_prime_prime_w.compress().as_bytes());

    let s_challenges = vandemonde_challenge(transcript.challenge_scalar(b"s"), n - 2);
    let b_challenges = vandemonde_challenge(transcript.challenge_scalar(b"b"), 2);

    let s_prime = [&s_challenges[..], &b_challenges[..]].concat();
    assert_eq!(s_prime.len(), n);

    let gamma_prime = compute_gamma_prime(&gamma, n);

    let gamma_prime_t = matrix_transpose(&gamma_prime);
    let gamma_prime_t_s_prime = matrix_vector_mul(&gamma_prime_t, &s_prime);

    let a = row_row_sub(&w_prime, &s_prime);
    let b = row_row_add(&w_prime_prime, &gamma_prime_t_s_prime);

    let t = -crate::math_utils::inner_product(&a, &b);

    let gamma_t = matrix_transpose(&gamma);
    let gamma_t_s = matrix_vector_mul(&gamma_t, &s_challenges);

    let expected_t = crate::math_utils::inner_product(&s_challenges, &gamma_t_s);
    assert_eq!(t, expected_t);
}

// Creates a vector from the scalar `x`
// contents of vector = <x, x^2, x^3,.., x^n>

// XXX: double check that it is fine to use a vandermonde matrix to
// expand challenges instead of fetching each challenge from the distribution
// so we don't need `n` different challenges
fn vandemonde_challenge(x: Scalar, n: usize) -> Vec<Scalar> {
    let mut challenges: Vec<Scalar> = Vec::with_capacity(n);

    let mut x_n = x.clone();

    challenges.push(x_n);

    for i in 1..n {
        x_n = x_n * x_n;
        challenges.push(x_n)
    }

    assert_eq!(challenges.len(), n);

    challenges
}

fn compute_gamma_prime(gamma: &Vec<Vec<Scalar>>, n: usize) -> Vec<Vec<Scalar>> {
    let mut gamma_prime = gamma.clone();

    // Pad the gamma rows with zeroes at the end
    gamma_prime = gamma_prime
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

    let mut row_n_minus_one = vec![Scalar::zero(); n - 2];
    row_n_minus_one.push(Scalar::one());
    row_n_minus_one.push(Scalar::zero());

    gamma_prime.push(row_n_minus_two);
    gamma_prime.push(row_n_minus_one);

    gamma_prime
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math_utils::*;
    use sha3::Sha3_512;
    #[test]
    fn test_create_qesa_inner_proof() {
        let mut rng = rand::thread_rng();

        let n = 4;

        let (witness, matrix) = helper_create_solutions(n - 2, 2);

        let G: Vec<RistrettoPoint> = (0..n).map(|_| RistrettoPoint::random(&mut rng)).collect();
        let H: Vec<RistrettoPoint> = (0..n).map(|_| RistrettoPoint::random(&mut rng)).collect();

        let mut transcript = Transcript::new(b"qesa_inner");

        let Q = RistrettoPoint::hash_from_bytes::<Sha3_512>(b"test point");

        let r_prime = (0..2).map(|_| Scalar::random(&mut rng)).collect();

        create(&mut transcript, G, H, &Q, matrix, witness, r_prime);
    }
    // Creates a system of quadratic equations with solutions
    // and a witness
    fn helper_create_solutions(n: usize, num_of_matrices: usize) -> (Vec<Scalar>, block_matrix) {
        let mut rng = rand::thread_rng();
        let witness: Vec<Scalar> = (0..n).map(|_| Scalar::random(&mut rng)).collect();

        let mut bm = block_matrix::new();

        for i in 0..num_of_matrices {
            let mut gamma_i: Vec<Vec<Scalar>> = Vec::new();
            for _ in 0..n {
                // Use gram schmidt to create suitable solutions for each system of eqns
                let x: Vec<Scalar> = (0..n).map(|_| Scalar::random(&mut rng)).collect();
                let row_of_eqns = crate::inner_product_proof::gramschmidt::orth(&witness, &x);
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
