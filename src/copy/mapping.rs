use curve25519_dalek::{ristretto::RistrettoPoint, scalar::Scalar, traits::VartimeMultiscalarMul};
use std::collections::BTreeSet;

//XXX: We can simplify this code, maybe we can compute the mapping more efficiently or in a more readable manner using iterators.

#[allow(dead_code)]
pub fn commit(
    G_Vec: &[RistrettoPoint],
    indices: BTreeSet<u16>,
    msg: Vec<Scalar>,
    rand: Scalar,
) -> RistrettoPoint {
    assert_eq!(indices.len(), msg.len());

    let points: Vec<RistrettoPoint> = indices.iter().map(|i| G_Vec[*i as usize]).collect();

    let mut commitment = RistrettoPoint::vartime_multiscalar_mul(msg.iter(), points.iter());

    // Add random scalar
    commitment = commitment + G_Vec.last().unwrap() * rand;

    commitment
}
#[derive(Debug)]
pub struct Mapping {
    num_commitments: usize,
    msg_map: Vec<MessageMap>,
    unique_key_indices: BTreeSet<u16>,
    last_index: u16,
}

#[derive(Debug)]
// Stores the relevant information about a specific message
// `key_id` refers to the indice of the point which committed the message
// If P = m * G[2] then key_id = 2
// crs_id; a commitment can have multiple messages. For example P = m_1 * G[1] + m_2 * G[2]
// The above commitment has two key_ids. The set of key_ids is called the commitmentCRS. We use crs for short.
// Each commitment has a crs. Both the verifier and the prover share a view of the set of crs's or the set of indices.
// The crs_id refers to the position of the specific `crs` that was used for a message in a commitment.
// position refers to the position of the message in the witness. This is computed in `compute_mapping`
pub struct MessageMap {
    msg_id: usize,
    key_id: u16,
    crs_id: u16,
    position: u16,
}

#[allow(dead_code)]
pub fn compute_mapping(
    G_Vec: &[RistrettoPoint],
    witness_size: usize,
    sets_of_indices: &[BTreeSet<u16>],
) -> Mapping {
    let n = G_Vec.len() as u16;

    let num_commitments = sets_of_indices.len();

    // Collect all of the indices into a set(automatically removing duplicates)
    let unique_indices: BTreeSet<u16> = sets_of_indices
        .iter()
        .map(|indices| indices.clone())
        .flatten()
        .collect();
    //XXX: assert unique_indices[i] cannot be equal to zero, as reserved for witness element (see paper)?

    // This is the amount of messages that are being committed to
    // For example, if a commitment = aG+bH+cQ ; the total for this commitment would be 3
    // Since (a,b,c) is the 3-tuple that is being committed to
    let total_messages: usize = sets_of_indices.iter().map(|indices| indices.len()).sum();

    let position = witness_size as u16;
    // These are the available positions in the witness vector for the messages
    // We filter out the generator indices as they are `reserved` positions in the witness (see paper)
    let available_positions: Vec<u16> = (position..n - 2)
        .filter(|x| !unique_indices.contains(&x))
        .take(total_messages)
        .map(|x| x)
        .collect();

    let last_index = *available_positions.last().unwrap();

    // This checks if there are more messages than available positions
    // The take() method will `take` as many items as it can
    // If there is not enough, it will not panic
    // Hence, the following assert method.
    assert_eq!(total_messages, available_positions.len());

    // In order to perform Qesa_Copy, we need various information
    // This information, will be stored in a collection of maps.
    //
    // The first map, we will need is a map from each message to an available position in the witness - map(a)
    // To visualise simple example:
    // Imagine we had two commitments: P_1 = a * G[1] + b * G[2] and P_2 = c * G[3] + d * G[4]
    // We assign each committed item a id
    // a = 1 , b = 2, c = 3, d = 4
    // We then assign each id to an available position
    // `a` would go into the first available position
    // `b` would go in the second available position
    // `c` would go in the third available position
    // `d` would go in the fourth available position
    // Note that the total_messages = 4 in this case
    //
    // Another mapping we need is the mapping from a message to it's assosciated point - map(b)
    // For example:
    // Looking solely at P_1 = a * G[1] + b * G[2]
    // the message `a` would map to G[1]
    // the message `b` would map to G[2]
    // We may refer to this assosciated point as the `key_id`
    // the `key_id` refers to the position that the point is in the vector of points; G[1] has a key_id of 1.
    //
    //
    // Another mapping we need is the mapping from a message to the assosciated set_of_indices
    // For example:
    // When we `commit()` we need to supply the indices that we want to use to commit the messages
    // Looking at P_1 = a * G[1] + b * G[2]
    // The commit method would look something like commit(msgs:[a,b], indices: [1,2], point_vector: G_Vec)
    // As you can see, each commitment, has an assosciated set of indices. We want to create a map for this.
    // We assume that the prover and verifier share a common ordered view of the set of indices, or this mapping will be off.
    // For a protocol like rangeproof; The set of indices will always be the same, so we do not need to send this.
    // I believe the general use-case will be similar to rangeproofs and so sending the set of indices, will be special.

    // Create all message maps and assign each message map an ID
    // This ID corresponds to a specific message and depends on what order the commitments are in
    let mut msgs_map: Vec<MessageMap> = (0..total_messages)
        .map(|i| MessageMap {
            msg_id: i,
            key_id: 0,
            crs_id: 0,
            position: 0,
        })
        .collect();

    // Give each message an available position in the witness vector - map(a)
    for (msg_map, av_pos) in msgs_map.iter_mut().zip(available_positions.iter()) {
        msg_map.position = *av_pos;
    }

    // Now we map each msg to it's assosciated point - map(b)
    let mut msgs_map_iter = msgs_map.iter_mut();
    for set_of_indices in sets_of_indices.iter() {
        for indice in set_of_indices {
            msgs_map_iter.next().unwrap().key_id = *indice;
        }
    }

    // Now we map the messages to their assosciated set of indices - map(c)
    // XXX: We can combine this into the above for loop, for commit history we commit like this first as it's a bit clearer
    let mut msgs_map_iter = msgs_map.iter_mut();
    for (crs_id, set_of_indices) in sets_of_indices.iter().enumerate() {
        for _ in set_of_indices {
            msgs_map_iter.next().unwrap().crs_id = crs_id as u16;
        }
    }

    Mapping {
        num_commitments: num_commitments,
        msg_map: msgs_map,
        unique_key_indices: unique_indices,
        last_index: last_index,
    }
}

#[test]
fn test_random_compute_mapping() {
    use super::*;
    use rand::Rng;
    let mut rng = rand::thread_rng();

    let n = 18;
    let num_of_commitments = 4;

    let G: Vec<RistrettoPoint> = (0..n).map(|_| RistrettoPoint::random(&mut rng)).collect();

    // We will hae four commitments and so four sets of indices
    let sets_of_indices: Vec<BTreeSet<u16>> = (0..num_of_commitments)
        .map(|_| {
            let mut set = BTreeSet::new();

            // Each commitment will only commit to two elements at most(excluding the random value)
            let num_of_indices = 2;

            for i in 0..num_of_indices {
                // XXX: Need to confirm as there was a disagreement; we cannot use the first element as that is reserved for the fixed witness element
                let num: u16 = rand::thread_rng().gen_range(2, n - 3); // [2,n-3]
                set.insert(num);
            }
            set
        })
        .collect();
    // The above code tells us that we are committing to (num_of_commitments * num_of_indices) elements
    // In this case we are committing to max 8 elements,if each element is commmited to with a different indice, then
    // we have 8 x 2 = 16 reserved spots in the witness for the generator points, and 2 spots reserved for the randomness.
    // `n` therefore needs to be _atleast_ 18
    // Note that since, randomness is used to generate the indices, the tests may sporadically pass with `n` <  18
    // if the same indice is added to the set

    let map = compute_mapping(&G, 1, &sets_of_indices);
    assert_eq!(num_of_commitments, map.num_commitments);

    // Check that each message has been assigned a position
    let total_messages: usize = sets_of_indices.iter().map(|indices| indices.len()).sum();
    for msg_id in 0..total_messages {
        let msg_map = map.msg_map.iter().find(|&x| x.msg_id == msg_id);
        assert!(msg_map.is_some());
    }

    // Check that the last index corresponds to the position of the last message
    let last_msg_id = total_messages - 1;
    assert_eq!(map.msg_map[last_msg_id].position, map.last_index);

    // There should only be 4 crs_ids from 0..3 since we only have four commitments
    // They should also be ordered by msg_id and ascending
    // This code essentially `splits` the vector whereever the crs_id should be different
    // We will then have a vector of vectors of messagemaps
    // XXX: Maybe this is too complicated and we should refactor to simplify?
    let mut map_msg : &[MessageMap] = &map.msg_map;
    let chunked_maps : Vec<_>= sets_of_indices.iter().enumerate().map(|(i,set_of_indice)| {
        let num_of_indice = set_of_indice.len();
        let (head, tail) = map_msg.split_at(num_of_indice);
        map_msg = tail;
        head
    }).collect();
    for (crs_id,maps) in chunked_maps.iter().enumerate(){
        for map in maps.iter() {
            assert!(map.crs_id == crs_id as u16)
        }
    }
}
