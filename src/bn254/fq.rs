//! base field

use super::limbs::{
    add, double, from_u64, invert, little_fermat, mont, mul, neg, random_limbs, square, sub,
};
use super::math::sbb;
use core::fmt::{Debug, Formatter, Result};
use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub};
use rand_core::RngCore;

pub(crate) const MODULUS: [u64; 4] = [
    0x3c208c16d87cfd47,
    0x97816a916871ca8d,
    0xb85045b68181585d,
    0x30644e72e131a029,
];

/// R = 2^256 mod q
pub(crate) const R: [u64; 4] = [
    0xd35d438dc58f0d9d,
    0x0a78eb28f5c70b3d,
    0x666ea36f7879462c,
    0x0e0a77c19a07df2f,
];

/// R^2 = 2^512 mod q
pub(crate) const R2: [u64; 4] = [
    0xf32cfc5b538afa89,
    0xb5e71911d44501fb,
    0x47ab1eff0a417ff6,
    0x06d89f71cab8351f,
];

/// R^3 = 2^768 mod q
pub(crate) const R3: [u64; 4] = [
    0xb1cd6dafda1530df,
    0x62f210e6a7283db6,
    0xef7f0b0c0ada0afb,
    0x20fd6e902d592544,
];

/// INV = -(q^{-1} mod 2^64) mod 2^64
pub(crate) const INV: u64 = 0x87d20782e4866389;

#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) struct Fq(pub(crate) [u64; 4]);

impl Fq {
    pub(crate) const fn zero() -> Self {
        Self([0; 4])
    }

    pub(crate) const fn is_zero(self) -> bool {
        self.0[0] == 0 && self.0[1] == 0 && self.0[2] == 0 && self.0[3] == 0
    }

    pub(crate) const fn one() -> Self {
        Self(R)
    }

    pub(crate) const fn double(self) -> Self {
        Self(double(self.0, MODULUS))
    }

    pub(crate) const fn square(self) -> Self {
        Self(square(self.0, MODULUS, INV))
    }

    pub(crate) fn invert(self) -> Option<Self> {
        match invert(self.0, little_fermat(MODULUS), R, MODULUS, INV) {
            Some(x) => Some(Self(x)),
            None => None,
        }
    }

    pub(crate) const fn to_mont_form(val: [u64; 4]) -> Self {
        Self(mul(val, R2, MODULUS, INV))
    }

    pub(crate) const fn montgomery_reduce(self) -> [u64; 4] {
        mont(
            [self.0[0], self.0[1], self.0[2], self.0[3], 0, 0, 0, 0],
            MODULUS,
            INV,
        )
    }

    pub fn random<R: RngCore>(rand: &mut R) -> Self {
        Self(random_limbs(rand, R2, R3, MODULUS, INV))
    }
}

impl Add for Fq {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(add(self.0, rhs.0, MODULUS))
    }
}

impl AddAssign for Fq {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl Sub for Fq {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self(sub(self.0, rhs.0, MODULUS))
    }
}

impl Neg for Fq {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(neg(self.0, MODULUS))
    }
}

impl Mul<Fq> for Fq {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self(mul(self.0, rhs.0, MODULUS, INV))
    }
}

impl MulAssign for Fq {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl From<u64> for Fq {
    fn from(val: u64) -> Fq {
        Fq(from_u64(val, R2, MODULUS, INV))
    }
}

impl Debug for Fq {
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

impl Fq {
    pub fn to_bytes(self) -> [u8; 32] {
        let mut res = [0; 32];
        let tmp = self.montgomery_reduce();
        res[0..8].copy_from_slice(&tmp[0].to_le_bytes());
        res[8..16].copy_from_slice(&tmp[1].to_le_bytes());
        res[16..24].copy_from_slice(&tmp[2].to_le_bytes());
        res[24..32].copy_from_slice(&tmp[3].to_le_bytes());

        res
    }

    pub fn from_bytes(bytes: [u8; 32]) -> Option<Self> {
        let l0 = u64::from_le_bytes(bytes[0..8].try_into().unwrap());
        let l1 = u64::from_le_bytes(bytes[8..16].try_into().unwrap());
        let l2 = u64::from_le_bytes(bytes[16..24].try_into().unwrap());
        let l3 = u64::from_le_bytes(bytes[24..32].try_into().unwrap());

        let (_, borrow) = sbb(l0, MODULUS[0], 0);
        let (_, borrow) = sbb(l1, MODULUS[1], borrow);
        let (_, borrow) = sbb(l2, MODULUS[2], borrow);
        let (_, borrow) = sbb(l3, MODULUS[3], borrow);

        if borrow & 1 == 1 {
            Some(Self([l0, l1, l2, l3]) * Self(R2))
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_core::OsRng;

    #[test]
    fn test_montgomery() {
        let mut rng = OsRng;
        for _ in 0..100000 {
            let s = Fq::random(&mut rng);
            let r = s.montgomery_reduce();
            let s_prime = Fq::to_mont_form(r);
            assert_eq!(s, s_prime);
        }
    }
}
