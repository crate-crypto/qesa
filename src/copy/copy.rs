// This module is in charge of allowing qesa to take commitments
// as input to the qesa_inner/zk protocol
//XXX: This file in particular should be cleaned up, as first pass is very rough.
// XXX: Is a u16 good enough as the type for the mapping? For: sha256 was 26k, can we expect circuits to be bigger than sha256? Probably not.
pub struct Copy {}

