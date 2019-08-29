use crate::math_utils::inner_product;
use curve25519_dalek::scalar::Scalar;

fn sample_m_n_plus(n: usize) -> Vec<Scalar> {
    let mut r: Vec<Scalar> = vec![Scalar::zero(); n];

    let mut rng = rand::thread_rng();

    r[0] = Scalar::random(&mut rng);
    r[1] = Scalar::random(&mut rng);

    let mut i = 4;
    while i <= n {
        
        r[i / 2] = Scalar::random(&mut rng);
        r[i / 2 + 1] = Scalar::random(&mut rng);

        i = i * 2
    }

    r[n - 2] = Scalar::random(&mut rng);
    r[n - 1] = Scalar::random(&mut rng);

    assert_eq!(n, r.len());

    r
}

/// Returns an vector that is orthogonal to `a`
/// Formally, this function projects `b` orthogonally onto the
/// line spanned by `a`
pub fn orth(a: &[Scalar], b: &[Scalar]) -> Vec<Scalar> {
    let aa: Scalar = inner_product(&a, &a);
    let ab: Scalar = inner_product(&a, &b);

    assert_ne!(Scalar::zero(), ab);

    let x: Scalar = ab * aa.invert();

    let ax: Vec<Scalar> = a.iter().map(|k| k * x).collect();

    let b_minus_ax = b.iter().zip(ax.iter()).map(|(r, s)| r - s).collect();

    b_minus_ax
}

/// Samples a random matrix and returns a random vector
/// which is orthogonal to `a`
pub fn sample_gram_schmidt(a: &[Scalar]) -> Vec<Scalar> {
    let sampled_matrix = sample_m_n_plus(a.len());
    orth(a, &sampled_matrix)
}
/// Samples a random matrix and returns a random vector
/// which is orthogonal to `a` and to `b`
/// XXX: We can make a better API for this and naming convention
/// XXX: Naming it sample_gram_schmidt twice could imply that we sample twice
pub fn sample_gram_schmidt_twice(a: &[Scalar], b: &[Scalar]) -> Vec<Scalar> {
    let orth_a_b = orth(a, &b);
    let orth_a_sample = sample_gram_schmidt(a);

    orth(&orth_a_b, &orth_a_sample)
}

#[test]
fn test_orth() {
    let mut rng = rand::thread_rng();

    let a: Vec<Scalar> = (0..4).map(|_| Scalar::random(&mut rng)).collect();
    let b: Vec<Scalar> = (0..4).map(|_| Scalar::random(&mut rng)).collect();

    let orth_a = orth(&a, &b);

    assert_eq!(inner_product(&orth_a, &a), Scalar::zero());
}

#[test]
fn test_sample_gram_schmidt() {
    let mut rng = rand::thread_rng();

    let a: Vec<Scalar> = (0..4).map(|_| Scalar::random(&mut rng)).collect();

    let res = sample_gram_schmidt(&a);

    assert_eq!(inner_product(&res, &a), Scalar::zero());
}

#[test]
fn test_sample_gram_schmidt_twice() {
    let mut rng = rand::thread_rng();

    let a: Vec<Scalar> = (0..4).map(|_| Scalar::random(&mut rng)).collect();
    let b: Vec<Scalar> = (0..4).map(|_| Scalar::random(&mut rng)).collect();

    let res = sample_gram_schmidt_twice(&a, &b);

    assert_eq!(inner_product(&res, &a), Scalar::zero());
    assert_eq!(inner_product(&res, &b), Scalar::zero());
}
