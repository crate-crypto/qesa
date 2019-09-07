use crate::copy::mapping::{Mapping, MessageMap};
use crate::inner;
use crate::math_utils::vandemonde_challenge;
use crate::matrix::*;
use crate::transcript::TranscriptProtocol;
use curve25519_dalek::{
    ristretto::{CompressedRistretto, RistrettoPoint},
    scalar::Scalar,
    traits::VartimeMultiscalarMul,
};
use merlin::Transcript;
// This module is in charge of allowing qesa to take commitments
// as input to the qesa_inner/zk protocol
#[derive(Clone)]
pub struct Copy {
    c_w: CompressedRistretto,
    inner: inner::Inner,
}
pub fn create(
    transcript: &mut Transcript,
    mut G_Vec: Vec<RistrettoPoint>,
    H_Vec: Vec<RistrettoPoint>,
    Q: &RistrettoPoint,
    map: Mapping,
    mut gamma_i: block_matrix,
    mut w: Vec<Scalar>,
    commitments: Vec<(Vec<Scalar>, Scalar)>, // (messages, rand)
) -> Copy {
    // 1. Setup Witness
    //
    // number of commitments stated in the map
    // should be equal to the number of commitments passed as an argument
    assert_eq!(map.num_commitments, commitments.len());

    // first element in the witness should be equal to one (see paper)
    assert_eq!(w[0], Scalar::one());

    let n = G_Vec.len();

    // Pad the witness with zeroes
    let pad_by = (n - 2) - (w.len());
    let zeroes: Vec<Scalar> = (0..pad_by).map(|_| Scalar::zero()).collect();
    w.extend(zeroes);

    // Create extended witness with two random values
    let mut rng = rand::thread_rng();
    let r_prime: Vec<Scalar> = (0..2).map(|_| Scalar::random(&mut rng)).collect();
    let mut w_prime = [&w[..], &r_prime[..]].concat();

    // Use the map to place the messages in the correct position
    // in the witness
    //
    // Collect all messages
    let all_messages :Vec<Scalar>= commitments.iter().map(|(msgs, _)|msgs).flatten().cloned().collect();
    for (msg, msg_map) in all_messages.iter().zip(map.msgs_map.iter()) {
        w_prime[msg_map.position as usize] = *msg;
    }

    // 2. Compute commitment to witness
    //
    let c_w = RistrettoPoint::vartime_multiscalar_mul(w_prime.iter(), G_Vec.iter());
    transcript.append_message(b"c_w", c_w.compress().as_bytes());

    // 3. Use challenge to adjust the witness value
    //
    let alpha = transcript.challenge_scalar(b"alpha");
    let alpha_challenges = vandemonde_challenge(alpha, commitments.len());

    // Compute aggregated random value
    //
    // Collect all of the random values from the openings
    let all_rand_openings: Vec<Scalar> = commitments.iter().map(|(_, rnd)| *rnd).collect();
    // use challenge to aggregate all values
    let aggregated_rand: Scalar = all_rand_openings
        .iter()
        .zip(alpha_challenges.iter())
        .map(|(rnd, alpha)| rnd * alpha)
        .sum();
    w_prime[n-1] += aggregated_rand; 

    // Compute aggregated witness values
    // XXX: This needs cleaning up and more commenting, so that readers know what is happening
    for msg_map in map.msgs_map.iter() {

        // For each message. We need to figure out what position 
        // it is in, in terms of it's original commitment 
        // This is the same as the key_id or generator indice (see paper - orange values)
        let key_id_position = msg_map.key_id;
        // This position is it's usable position in the witness vector. (see paper - green values)
        let msg_position = msg_map.position;
        // For each row of messages, the challenge value is the same
        // In order to determine what challenge is needed
        // We can fetch the crs_id 
        let crs_id = msg_map.crs_id;

        let w_msg_pos = w_prime[msg_position as usize];

        w_prime[key_id_position as usize] += w_msg_pos * alpha_challenges[crs_id as usize];
    }

    // 4. Modify gamma to match the new commited values
    // Make sure they are "correct" - expand on this
    //
    let zeroes: Vec<Scalar> = (0..n - 2).map(|_| Scalar::zero()).collect();

    for key_id in map.unique_key_indices {
        let mut matrix: Vec<Vec<Scalar>> = Vec::new();
        let mut new_eqn: Vec<Scalar> = (0..n - 2).map(|_| Scalar::zero()).collect();
        new_eqn[key_id as usize] = -Scalar::one();

        // Get all messages for this key_id
        let msgs_for_key_id: Vec<&MessageMap> = map
            .msgs_map
            .iter()
            .filter(|x| (x.key_id == key_id))
            .collect();

        for msg in msgs_for_key_id.iter() {
            new_eqn[msg.position as usize] = alpha_challenges[msg.crs_id as usize];
        }

        // Add to matrix
        matrix.push(new_eqn);
        // Add dummy lines, to make the matrix n-2 x n-2
        // XXX: Once sparse matrices are implemented, we can remove this.
        while matrix.len() != n - 2 {
            matrix.push(zeroes.clone());
        }

        // Add matrix to block matrix
        gamma_i.push(matrix);
    }

    // Decompose the extended witness by splitting vector back into r_prime and w_prime
    // Note that w_prime and r_prime are different than the original w and r_prime
    let r_prime = w_prime.split_off(n - 2);
    assert_eq!(r_prime.len(), 2);
    assert_eq!(w_prime.len(), n-2);
    //5. Run Qesa_inner as a sub protocol
    //
    let proof = inner::create(transcript, G_Vec, H_Vec, Q, &gamma_i, w_prime, r_prime);
    Copy {
        c_w: c_w.compress(),
        inner: proof,
    }
}
impl Copy {
    fn verify(
        &self,
        transcript: &mut Transcript,
        G_Vec: Vec<RistrettoPoint>,
        H_Vec: Vec<RistrettoPoint>,
        Q: &RistrettoPoint,
        mut gamma_i: block_matrix,
        map: Mapping,
        commitments: Vec<RistrettoPoint>,
    ) -> bool {
        let n = G_Vec.len();

        //1. Add commited witness to the transcript
        transcript.append_message(b"c_w", self.c_w.as_bytes());

        // 2. Use challenge to adjust the commitment values
        let alpha = transcript.challenge_scalar(b"alpha");
        let alpha_challenges = vandemonde_challenge(alpha, commitments.len());

        let mut commitment: RistrettoPoint = alpha_challenges
            .iter()
            .zip(commitments.iter())
            .map(|(alpha, commitment)| commitment * alpha)
            .sum();
        commitment = commitment + self.c_w.decompress().unwrap();

        // 4. Modify gamma to match the new commited values
        let zeroes: Vec<Scalar> = (0..n - 2).map(|_| Scalar::zero()).collect();

        for key_id in map.unique_key_indices {
            let mut matrix: Vec<Vec<Scalar>> = Vec::new();
            let mut new_eqn: Vec<Scalar> = (0..n - 2).map(|_| Scalar::zero()).collect();
            new_eqn[key_id as usize] = -Scalar::one();

            // Get all messages for this key_id
            let msgs_for_key_id: Vec<&MessageMap> = map
                .msgs_map
                .iter()
                .filter(|x| (x.key_id == key_id))
                .collect();

            for msg in msgs_for_key_id.iter() {
                new_eqn[msg.position as usize] = alpha_challenges[msg.crs_id as usize];
            }

            // Add to matrix
            matrix.push(new_eqn);
            // Add dummy lines, to make the matrix n-2 x n-2
            // Once sparse matrices are implemented, we can remove this.
            while matrix.len() != n - 2 {
                matrix.push(zeroes.clone());
            }
            // Add matrix to block matrix
            gamma_i.push(matrix);
        }
        let proof_inner = inner::Inner {
            alm_zk: self.inner.alm_zk.clone(),
            c_prime_w: commitment.compress(),
            c_prime_prime_w: self.inner.c_prime_prime_w,
        };

        proof_inner.verify(transcript, G_Vec, H_Vec, Q, &gamma_i)
        // self.inner.verify(transcript, G_Vec, H_Vec, Q, &gamma_i)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::copy::mapping;
    use sha3::Sha3_512;
    use std::collections::BTreeSet;
    #[test]
    fn test_copy_proof_creation() {
        // Create a witness with manadatory first element as one
        let mut witness = Vec::new();
        witness.push(Scalar::one());

        let n = 32;
        let mut rng = rand::thread_rng();

        // Create the CRS
        let G: Vec<RistrettoPoint> = (0..n).map(|_| RistrettoPoint::random(&mut rng)).collect();
        let H: Vec<RistrettoPoint> = (0..n).map(|_| RistrettoPoint::random(&mut rng)).collect();
        let Q = RistrettoPoint::hash_from_bytes::<Sha3_512>(b"test point");

        // Create a empty block matrix (no equations)
        let bm = block_matrix::new();

        // Generate the `messages` and random scalars for the commitments
        let num_commitments = 2;
        let (messages, rand_values) = helper_create_scalar_pair(num_commitments, [2, 3].to_vec());
        let total_messages: usize = messages.iter().map(|msg| msg.len()).sum();

        // The number of messages should be smaller than n-2-1
        // Because we need to map those onto the witness which will be at most `n-2-1` in length
        // XXX: even without this assert, the test will fail due to compute_mapping failing
        assert!((n - 2 - 1) > total_messages);

        // Generate the set of indices
        let sets_of_indices: Vec<_> = messages
            .iter()
            .map(|msg| {
                let num_indices = msg.len();

                let mut b: BTreeSet<u16> = BTreeSet::new();
                for i in 0..num_indices {
                    b.insert((i + 1) as u16);
                }
                b
            })
            .collect();

        // Generate actual commitments
        // XXX: Passing the actual EC points to the prover's creation is inconvenient because the prover also needs their
        // openings. For now, let's just pass the openings themselves.
        // What we can do is create an object that stores the `state`. We commit using this object, store the openings
        // and return the Pedersen Commitment
        // You can uncomment this code to convince yourself that the setup we did above, does indeed work without problems.
        let commitments: Vec<RistrettoPoint> = sets_of_indices
            .iter()
            .enumerate()
            .map(|(i, indices)| {
                let msg = &messages[i];
                let rnd = &rand_values[i];

                mapping::commit(&G, indices, msg, rnd)
            })
            .collect();

        let map = mapping::compute_mapping(&G, witness.len(), &sets_of_indices);
        let mut transcript = Transcript::new(b"copy");

        let openings = messages
            .into_iter()
            .zip(rand_values.into_iter())
            .map(|(message, rand)| (message, rand))
            .collect();

        let proof = create(
            &mut transcript,
            G.clone(),
            H.clone(),
            &Q,
            map.clone(),
            bm,
            witness,
            openings,
        );

        // Verify
        let mut transcript = Transcript::new(b"copy");
        let bm = block_matrix::new();
        assert!(proof.verify(&mut transcript, G, H, &Q, bm, map.clone(), commitments));
    }
    #[test]
    fn test_copy_proof_creation_single_commit() {
        // Create a witness with manadatory first element as one
        let mut witness = Vec::new();
        witness.push(Scalar::one());

        let n = 16;
        let mut rng = rand::thread_rng();

        // Create the CRS
        let G: Vec<RistrettoPoint> = (0..n).map(|_| RistrettoPoint::random(&mut rng)).collect();
        let H: Vec<RistrettoPoint> = (0..n).map(|_| RistrettoPoint::random(&mut rng)).collect();
        let Q = RistrettoPoint::hash_from_bytes::<Sha3_512>(b"test point");

        // Create a empty block matrix (no equations)
        let bm = block_matrix::new();

        // Generate the `message` and random scalar for the commitment
        let message : Vec<Scalar>= vec![Scalar::from(1 as u8),Scalar::from(1 as u8)];
        let rand = Scalar::from(1 as u8);

        // Generate indices for messages
        let mut set_of_indice : Vec<BTreeSet<u16>> = Vec::new();                
        let mut b: BTreeSet<u16> = BTreeSet::new();
                for i in 0..message.len() {
                    b.insert((i + 1) as u16);
                }
        set_of_indice.push(b);

        // Generate actual commitments
        let commitment = mapping::commit(&G, &set_of_indice[0], &message, &rand);
        let commitments: Vec<RistrettoPoint> = vec![commitment];
        
        let map = mapping::compute_mapping(&G, witness.len(), &set_of_indice);
        let mut transcript = Transcript::new(b"copy");

        let openings : Vec<(Vec<Scalar>, Scalar)>= vec![(message,rand)]; 

        let proof = create(
            &mut transcript,
            G.clone(),
            H.clone(),
            &Q,
            map.clone(),
            bm,
            witness,
            openings,
        );

        // Verify
        let mut transcript = Transcript::new(b"copy");
        let bm = block_matrix::new();
        assert!(proof.verify(&mut transcript, G, H, &Q, bm, map.clone(), commitments));
    }
    // Creates a pair of vectors, each of size `n` with random elements.
    // The first element represents the `messages` we want to commit to.
    // Each message size will be determined by the vector `k`. ie the size of the first message is determined by k[0]
    // The second element represents the randomness `r`
    fn helper_create_scalar_pair(n: usize, k: Vec<usize>) -> (Vec<Vec<Scalar>>, Vec<Scalar>) {
        assert!(n == k.len());

        // Generate all of the Commitment random values first.
        let randomness = generate_rand_scalars(n);

        let messages: Vec<Vec<Scalar>> = (0..n)
            .map(|i| {
                let num_messages = k[i];
                generate_rand_scalars(num_messages)
            })
            .collect();
        (messages, randomness)
    }

    fn generate_rand_scalars(n: usize) -> Vec<Scalar> {
        let mut rng = rand::thread_rng();
        (0..n).map(|_| Scalar::random(&mut rng)).collect()
    }
}
