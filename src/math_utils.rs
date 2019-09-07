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

// Creates a vector from the scalar `x`
// contents of vector = <x, x^2, x^3,.., x^n>
// XXX: double check that it is fine to use a vandermonde matrix in regards to testing distributions to
// expand challenges instead of fetching each challenge from the distribution
// so we don't need `n` different challenges
pub fn vandemonde_challenge(x: Scalar, n: usize) -> Vec<Scalar> {
    let mut challenges: Vec<Scalar> = Vec::with_capacity(n);

    let mut x_n = x.clone();

    challenges.push(x_n);

    for _ in 1..n {
        x_n = x_n * x_n;
        challenges.push(x_n)
    }

    assert_eq!(challenges.len(), n);

    challenges
}
