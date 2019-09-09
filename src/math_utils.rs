use curve25519_dalek::scalar::Scalar;

/// Computes the inner product between two scalar vectors
pub fn inner_product(a: &[Scalar], b: &[Scalar]) -> Scalar {
    a.iter().zip(b.iter()).map(|(a, b)| a * b).sum()
}

// Creates a vector from the scalar `x`
// contents of vector = <x, x^2, x^3,.., x^n>
// XXX: double check that it is fine to use a vandermonde matrix in regards to testing distributions to
// expand challenges instead of fetching each challenge from the distribution
// so we don't need `n` different challenges
pub fn vandemonde_challenge(mut x: Scalar, n: usize) -> Vec<Scalar> {
    let mut challenges: Vec<Scalar> = Vec::with_capacity(n);
    for _ in 0..n {
        challenges.push(x);
        x = x * x;
    }
    challenges
}

#[test]
fn test_vandermonde() {
    let n = 32;

    let two = Scalar::from(2 as u8);

    let challenges = vandemonde_challenge(two, n);
    assert_eq!(challenges.len(), n);

    let mut two_n = two;

    for i in 0..challenges.len() {
        assert_eq!(two_n, challenges[i]);
        two_n = two_n * two_n;
    }
}
