#![allow(non_snake_case)]
use curve25519_dalek::{ristretto::RistrettoPoint, scalar::Scalar, traits::VartimeMultiscalarMul};
use std::collections::BTreeSet;

/// The purpose of the mapping code is to map each message into the witness vector
///
/// What is a message?
/// A message is an element in a Pedersen commitment
/// For example, if P = aG + bH , then we have 2 messages `a` and `b`
///
/// In particular, we want to store three pieces of meta-data about a message; the commmitment_index, the position_in_witness and the commitment_crs_id
///
/// The `commitment_index` is the indice of the point that was used to commit that message
/// For example, Imagine we have a vector of points that we use to perform a Pedersen commitment `G`
/// We can commit to a commitment using P = a * G[0] + b * G[1] + c * G[2]
/// In this example, we have three message `a` , `b` and `c`.
/// The commitment_index for `a` is `0` because we used the point at the 0th index from the vector of points `G`
///
/// The `position_in_witness` is a value that is computed `compute_mapping`
/// It refers to an available position in the witness that this message can be placed into
/// This file only computes the mapping, Qesa_Copy places the messages in the relevant positions
///
/// The `commitment_crs_id` is the id of the set that was used to commit the message
/// In this proving system, one can commit using a continuos sequence of points: P = a * G[0] + b * G[1] + c * G[2]
/// Or you can commit to a combination of points: P = a * G[3] + b * G[5] + c * G[2] + d * G[4]
/// Therefore, when we commit messages we also need to specify what set of points are used to commit the message
/// This set is what we call a `commitment_crs` . In the above example, the commitment_crs would be {3,5,2,4}
/// The `commitment_crs_id` is number that refers to a specific `commitment_crs`. In particular, the `commitment_crs_id` for a specific message
/// is the id of the commitment_crs that was used to commit that message.

pub fn commit(
    G_Vec: &[RistrettoPoint],
    indices: &BTreeSet<u16>,
    msg: &[Scalar],
    rand: &Scalar,
) -> RistrettoPoint {
    assert_eq!(indices.len(), msg.len());

    let points: Vec<RistrettoPoint> = indices.iter().map(|i| G_Vec[*i as usize]).collect();

    let mut commitment = RistrettoPoint::vartime_multiscalar_mul(msg.iter(), points.iter());

    // Add random scalar
    commitment = commitment + G_Vec.last().unwrap() * rand;

    commitment
}
#[derive(Debug, Clone)]
pub struct Mapping {
    pub(crate) num_commitments: usize,
    pub(crate) msgs_map: Vec<MessageMap>,
    pub(crate) unique_commitment_indices: BTreeSet<u16>,
    pub(crate) last_index: u16,
}

#[derive(Debug, Clone)]
// Stores the relevant information about a specific message
pub struct MessageMap {
    pub(crate) msg_id: usize,
    pub(crate) commitment_index: u16,
    pub(crate) commitment_crs_id: u16,
    pub(crate) position_in_witness: u16,
}
pub fn compute_mapping(
    n: usize,
    witness_size: usize,
    sets_of_indices: &[BTreeSet<u16>],
) -> Mapping {
    let n = n as u16;

    let num_commitments = sets_of_indices.len();

    // Collect all of the indices into a set(automatically removing duplicates)
    let unique_indices: BTreeSet<u16> = sets_of_indices
        .iter()
        .map(|indices| indices.clone())
        .flatten()
        .collect();

    // The first element in the witness is reserved and so users cannot commit to it.
    unique_indices.contains(&0);

    // This is the total amount of messages that are being committed to
    // Which is equal to the total amount of indices
    let total_messages: usize = sets_of_indices.iter().map(|indices| indices.len()).sum();

    let position = witness_size as u16;
    // These are the available positions in the witness vector for the messages
    // We filter out the indices of points that are being used to commit the messages as they are `reserved` positions in the witness (see paper)
    let available_positions: Vec<u16> = (position..n - 2)
        .filter(|x| !unique_indices.contains(&x))
        .take(total_messages)
        .map(|x| x)
        .collect();
    let last_index = *available_positions.last().unwrap();

    // This checks if there are more messages than available positions in the witness
    // The take() method will `take` as many items as it can
    // If there is not enough, it will not panic
    // Hence, the following assert method.
    assert_eq!(total_messages, available_positions.len());

    // Create all message maps and assign each message map an ID
    // ID = 0 corresponds to the first message, ID = 1 corresponds to the second message and so on.
    // The first message will the the first element from the commitment, which depends on the commitment
    // Therefore, permuting the commitment after a mapping will fail verification
    let mut msgs_map: Vec<MessageMap> = (0..total_messages)
        .map(|i| MessageMap {
            msg_id: i,
            commitment_index: 0,
            commitment_crs_id: 0,
            position_in_witness: 0,
        })
        .collect();

    // Map each message to an available position in the witness vector
    for (msg_map, av_pos) in msgs_map.iter_mut().zip(available_positions.iter()) {
        msg_map.position_in_witness = *av_pos;
    }

    // Map each msg to it's assosciated commitment index.
    // Recall this is the index of the point which is used to commit the message in the commmitment
    let mut msgs_map_iter = msgs_map.iter_mut();
    for set_of_indices in sets_of_indices.iter() {
        for index in set_of_indices {
            msgs_map_iter.next().unwrap().commitment_index = *index;
        }
    }

    // Map the messages to it's assosciated commitment_crs_id
    // XXX: We can combine this into the above for loop, for commit history we commit like this first as it's a bit clearer
    let mut msgs_map_iter = msgs_map.iter_mut();
    for (crs_id, set_of_indices) in sets_of_indices.iter().enumerate() {
        for _ in set_of_indices {
            msgs_map_iter.next().unwrap().commitment_crs_id = crs_id as u16;
        }
    }

    Mapping {
        num_commitments: num_commitments,
        msgs_map: msgs_map,
        unique_commitment_indices: unique_indices,
        last_index: last_index,
    }
}

#[test]
fn test_random_compute_mapping() {
    use rand::Rng;
    let mut rng = rand::thread_rng();

    let n = 18;
    let num_of_commitments = 4;

    let G: Vec<RistrettoPoint> = (0..n).map(|_| RistrettoPoint::random(&mut rng)).collect();

    // We will have four commitments and so four sets of indices
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

    let map = compute_mapping(G.len(), 1, &sets_of_indices);
    assert_eq!(num_of_commitments, map.num_commitments);

    // Check that each message has been assigned a position
    let total_messages: usize = sets_of_indices.iter().map(|indices| indices.len()).sum();
    for msg_id in 0..total_messages {
        let msg_map = map.msgs_map.iter().find(|&x| x.msg_id == msg_id);
        assert!(msg_map.is_some());
    }

    // Check that the last index corresponds to the position of the last message
    let last_msg_id = total_messages - 1;
    assert_eq!(
        map.msgs_map[last_msg_id].position_in_witness,
        map.last_index
    );

    // There should only be 4 crs_ids from 0..3 since we only have four commitments
    // They should also be ordered by msg_id and ascending
    // This code essentially `splits` the vector whereever the crs_id should be different
    // We will then have a vector of vectors of messagemaps
    // XXX: Maybe this is too complicated and we should refactor to simplify?
    let mut map_msg: &[MessageMap] = &map.msgs_map;
    let chunked_maps: Vec<_> = sets_of_indices
        .iter()
        .enumerate()
        .map(|(_, set_of_indice)| {
            let num_of_indice = set_of_indice.len();
            let (head, tail) = map_msg.split_at(num_of_indice);
            map_msg = tail;
            head
        })
        .collect();

    for (crs_id, maps) in chunked_maps.iter().enumerate() {
        for map in maps.iter() {
            assert!(map.commitment_crs_id == crs_id as u16)
        }
    }
}
