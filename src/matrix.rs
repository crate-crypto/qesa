use crate::math_utils::inner_product;
use curve25519_dalek::scalar::Scalar;

// we will use dense matrix operations until sparse matrix functionality
// is implemented.

// multiplies every element in the matrix by the scalar
// Returns the new matrix
fn matrix_scalar_mul(matrix: &Vec<Vec<Scalar>>, x: &Scalar) -> Vec<Vec<Scalar>> {
    matrix
        .iter()
        .map(|row| row.iter().map(|elem| elem * x).collect())
        .collect()
}
// Takes a matrix and a vector
// Returns a new vector i.e. (Ax=b)
pub fn matrix_vector_mul(matrix: &Vec<Vec<Scalar>>, vec: &[Scalar]) -> Vec<Scalar> {
    matrix.iter().map(|row| inner_product(row, &vec)).collect()
}

// Takes two rows and adds them together
// The (i) entry of row a is added to the i entry of row b
pub fn row_row_add(row_a: &[Scalar], row_b: &[Scalar]) -> Vec<Scalar> {
    assert_eq!(row_a.len(), row_b.len());

    row_a.iter().zip(row_b.iter()).map(|(a, b)| a + b).collect()
}

// Takes two rows and subtracts them
// The (i) entry of row a is subtracted from the i entry of row b
pub fn row_row_sub(row_a: &[Scalar], row_b: &[Scalar]) -> Vec<Scalar> {
    assert_eq!(row_a.len(), row_b.len());

    row_a.iter().zip(row_b.iter()).map(|(a, b)| a - b).collect()
}

// Takes two matrices and adds them together
// The (i,j) entry of matrix A is added to the (i,j) entry of B
fn matrix_matrix_add(A: &Vec<Vec<Scalar>>, B: &Vec<Vec<Scalar>>) -> Vec<Vec<Scalar>> {
    assert_eq!(A.len(), B.len());

    A.iter()
        .zip(B.iter())
        .map(|(row_a, row_b)| row_row_add(row_a, row_b))
        .collect()
}

// Takes the transpose of a matrix
pub fn matrix_transpose(matrix: &Vec<Vec<Scalar>>) -> Vec<Vec<Scalar>> {
    let mut transpose: Vec<Vec<Scalar>> = vec![Vec::new(); matrix[0].len()];

    for (_, row) in matrix.iter().enumerate() {
        for (i, element) in row.iter().enumerate() {
            transpose[i].push(element.clone());
        }
    }

    transpose
}

pub struct block_matrix {
    pub block: Vec<Vec<Vec<Scalar>>>,
}

impl block_matrix {
    pub fn new() -> block_matrix {
        block_matrix { block: Vec::new() }
    }
    pub fn push(&mut self, matrix: Vec<Vec<Scalar>>) {
        self.block.push(matrix)
    }
    pub fn with(matrix: Vec<Vec<Scalar>>) -> block_matrix {
        let mut b = block_matrix { block: Vec::new() };

        b.block.push(matrix);
        b
    }
    // Takes a block matrix and uses challenges to batch all vectors together into one matrix
    // by summing them
    pub fn block_matrix_batch(&self, challenges: &[Scalar]) -> Vec<Vec<Scalar>> {
        // This returns the block of matrices with each matrix multiplied by challenge
        let block_challenge: Vec<Vec<Vec<Scalar>>> = self
            .block
            .iter()
            .zip(challenges.iter())
            .map(|(matrix, challenge)| matrix_scalar_mul(matrix, challenge))
            .collect();

        // sums up each matrix
        let mut sum_block = block_challenge[0].clone();
        for i in 1..block_challenge.len() {
            sum_block = matrix_matrix_add(&sum_block, &block_challenge[i])
        }

        sum_block
    }
}
