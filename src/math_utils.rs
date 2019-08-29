use curve25519_dalek::scalar::Scalar;

/// Computes the inner product between two scalar vectors
pub fn inner_product(a: &[Scalar], b: &[Scalar]) -> Scalar {
    a.iter().zip(b.iter()).map(|(a, b)| a * b).sum()
}
/// Returns two random vectors `a` and `b` and a scalar `c`
/// such that <a,b> = c
pub fn helper_dot_product(n: usize) -> (Vec<Scalar>, Vec<Scalar>, Scalar) {
    let mut rng = rand::thread_rng();

    let a: Vec<Scalar> = (0..n).map(|_| Scalar::random(&mut rng)).collect();
    let b: Vec<Scalar> = (0..n).map(|_| Scalar::random(&mut rng)).collect();

    let t: Scalar = inner_product(&a, &b);
    (a, b, t)
}
