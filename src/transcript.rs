//! Defines a `TranscriptProtocol` trait for using a Merlin transcript.

use curve25519_dalek::scalar::Scalar;
use merlin::Transcript;

pub trait TranscriptProtocol {
    /// Compute a `label`ed challenge variable.
    fn challenge_scalar(&mut self, label: &'static [u8]) -> Scalar;
}

impl TranscriptProtocol for Transcript {
    fn challenge_scalar(&mut self, label: &'static [u8]) -> Scalar {
        let mut buf = [0u8; 64];
        self.challenge_bytes(label, &mut buf);

        Scalar::from_bytes_mod_order_wide(&buf)
    }
}
