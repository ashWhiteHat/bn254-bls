//! scalar field

use core::fmt::{Debug, Formatter, Result};
use core::ops::Mul;
use rand_core::RngCore;

use super::limbs::{mont, mul, random_limbs, to_nafs, Nafs};

pub(crate) const MODULUS: [u64; 4] = [
    0x43e1f593f0000001,
    0x2833e84879b97091,
    0xb85045b68181585d,
    0x30644e72e131a029,
];

/// `R = 2^256 mod r`
pub(crate) const R: [u64; 4] = [
    0xac96341c4ffffffb,
    0x36fc76959f60cd29,
    0x666ea36f7879462e,
    0x0e0a77c19a07df2f,
];

/// `R^2 = 2^512 mod r`
pub(crate) const R2: [u64; 4] = [
    0x1bb8e645ae216da7,
    0x53fe3ab1e35c59e3,
    0x8c49833d53bb8085,
    0x0216d0b17f4e44a5,
];

/// `R^3 = 2^768 mod r`
pub(crate) const R3: [u64; 4] = [
    0x5e94d8e1b4bf0040,
    0x2a489cbe1cfbb6b8,
    0x893cc664a19fcfed,
    0x0cf8594b7fcc657c,
];

/// INV = -(r^{-1} mod 2^64) mod 2^64
pub const INV: u64 = 0xc2e1f593efffffff;

#[derive(Clone, Copy)]
pub struct Fr(pub [u64; 4]);

impl Fr {
    pub fn random<R: RngCore>(rand: &mut R) -> Self {
        Self(random_limbs(rand, R2, R3, MODULUS, INV))
    }

    pub(crate) const fn montgomery_reduce(self) -> [u64; 4] {
        mont(
            [self.0[0], self.0[1], self.0[2], self.0[3], 0, 0, 0, 0],
            MODULUS,
            INV,
        )
    }

    pub(crate) fn to_nafs(self) -> Nafs {
        to_nafs(self.montgomery_reduce())
    }
}

impl Mul<Fr> for Fr {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self(mul(self.0, rhs.0, MODULUS, INV))
    }
}

impl Debug for Fr {
    fn fmt(&self, f: &mut Formatter) -> Result {
        write!(f, "0x")?;
        for limb in self.montgomery_reduce().iter().rev() {
            for byte in limb.to_be_bytes() {
                write!(f, "{:02x}", byte)?;
            }
        }
        Ok(())
    }
}
