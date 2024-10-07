use super::public::PublicKey;
use crate::bn254::{Fr, G1Affine};

#[derive(Clone, Copy, Debug)]
pub(crate) struct PrivateKey(pub(crate) Fr);

impl PrivateKey {
    pub(crate) fn new(value: Fr) -> Self {
        Self(value)
    }

    pub(crate) fn sign(self) {}

    pub(crate) fn public_key(self) -> PublicKey {
        let projective = G1Affine::generator() * self.0;
        let value = G1Affine::from(projective);

        PublicKey(value)
    }
}
