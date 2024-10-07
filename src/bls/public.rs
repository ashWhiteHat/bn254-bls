use crate::bn254::G1Affine;

#[derive(Clone, Copy, Debug)]
pub(crate) struct PublicKey(pub(crate) G1Affine);
