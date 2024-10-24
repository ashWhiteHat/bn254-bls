use super::math::{adb, adc, mac, sba, sbb};
use rand_core::RngCore;

#[inline(always)]
pub(crate) const fn add(a: [u64; 4], b: [u64; 4], p: [u64; 4]) -> [u64; 4] {
    let (l0, c) = adb(a[0], b[0]);
    let (l1, c) = adc(a[1], b[1], c);
    let (l2, c) = adc(a[2], b[2], c);
    let l3 = a[3].wrapping_add(b[3]).wrapping_add(c);

    let (l0, brw) = sba(l0, p[0]);
    let (l1, brw) = sbb(l1, p[1], brw);
    let (l2, brw) = sbb(l2, p[2], brw);
    let (l3, brw) = sbb(l3, p[3], brw);

    let (l0, c) = adb(l0, p[0] & brw);
    let (l1, c) = adc(l1, p[1] & brw, c);
    let (l2, c) = adc(l2, p[2] & brw, c);
    let l3 = l3.wrapping_add(p[3] & brw).wrapping_add(c);

    [l0, l1, l2, l3]
}

#[inline(always)]
pub(crate) const fn sub(a: [u64; 4], b: [u64; 4], p: [u64; 4]) -> [u64; 4] {
    let (l0, brw) = sba(a[0], b[0]);
    let (l1, brw) = sbb(a[1], b[1], brw);
    let (l2, brw) = sbb(a[2], b[2], brw);
    let (l3, brw) = sbb(a[3], b[3], brw);

    let (l0, c) = adb(l0, p[0] & brw);
    let (l1, c) = adc(l1, p[1] & brw, c);
    let (l2, c) = adc(l2, p[2] & brw, c);
    let l3 = l3.wrapping_add(p[3] & brw).wrapping_add(c);

    [l0, l1, l2, l3]
}

#[inline(always)]
pub(crate) const fn double(a: [u64; 4], p: [u64; 4]) -> [u64; 4] {
    let l0 = a[0] << 1;
    let l1 = a[1] << 1 | a[0] >> 63;
    let l2 = a[2] << 1 | a[1] >> 63;
    let l3 = a[3] << 1 | a[2] >> 63;

    let (l0, brw) = sba(l0, p[0]);
    let (l1, brw) = sbb(l1, p[1], brw);
    let (l2, brw) = sbb(l2, p[2], brw);
    let (l3, brw) = sbb(l3, p[3], brw);

    let (l0, c) = adb(l0, p[0] & brw);
    let (l1, c) = adc(l1, p[1] & brw, c);
    let (l2, c) = adc(l2, p[2] & brw, c);
    let l3 = l3.wrapping_add(p[3] & brw).wrapping_add(c);

    [l0, l1, l2, l3]
}

#[inline(always)]
pub(crate) const fn mul(a: [u64; 4], b: [u64; 4], p: [u64; 4], inv: u64) -> [u64; 4] {
    let s = a[0] as u128 * b[0] as u128;
    let (l0, c) = (s as u64, (s >> 64) as u64);
    let s = a[0] as u128 * b[1] as u128 + c as u128;
    let (l1, c) = (s as u64, (s >> 64) as u64);
    let s = a[0] as u128 * b[2] as u128 + c as u128;
    let (l2, c) = (s as u64, (s >> 64) as u64);
    let s = a[0] as u128 * b[3] as u128 + c as u128;
    let (l3, l4) = (s as u64, (s >> 64) as u64);

    let s = a[1] as u128 * b[0] as u128 + l1 as u128;
    let (l1, c) = (s as u64, (s >> 64) as u64);
    let (l2, c) = mac(l2, a[1], b[1], c);
    let (l3, c) = mac(l3, a[1], b[2], c);
    let (l4, l5) = mac(l4, a[1], b[3], c);

    let s = a[2] as u128 * b[0] as u128 + l2 as u128;
    let (l2, c) = (s as u64, (s >> 64) as u64);
    let (l3, c) = mac(l3, a[2], b[1], c);
    let (l4, c) = mac(l4, a[2], b[2], c);
    let (l5, l6) = mac(l5, a[2], b[3], c);

    let s = a[3] as u128 * b[0] as u128 + l3 as u128;
    let (l3, c) = (s as u64, (s >> 64) as u64);
    let (l4, c) = mac(l4, a[3], b[1], c);
    let (l5, c) = mac(l5, a[3], b[2], c);
    let (l6, l7) = mac(l6, a[3], b[3], c);

    mont([l0, l1, l2, l3, l4, l5, l6, l7], p, inv)
}

#[inline(always)]
pub(crate) const fn square(a: [u64; 4], p: [u64; 4], inv: u64) -> [u64; 4] {
    let s = a[1] as u128 * a[0] as u128;
    let (l1, c) = (s as u64, (s >> 64) as u64);
    let s = a[2] as u128 * a[0] as u128 + c as u128;
    let (l2, c) = (s as u64, (s >> 64) as u64);
    let s = a[3] as u128 * a[0] as u128 + c as u128;
    let (l3, c) = (s as u64, (s >> 64) as u64);
    let s = a[1] as u128 * a[3] as u128 + c as u128;
    let (l4, c) = (s as u64, (s >> 64) as u64);
    let s = a[2] as u128 * a[3] as u128 + c as u128;
    let (l5, l6) = (s as u64, (s >> 64) as u64);
    let s = a[1] as u128 * a[2] as u128 + l3 as u128;
    let (l3, c) = (s as u64, (s >> 64) as u64);
    let (l4, c) = adb(l4, c);
    let l5 = l5.wrapping_add(c);

    let l7 = l6 >> 63;
    let l6 = (l6 << 1) | (l5 >> 63);
    let l5 = (l5 << 1) | (l4 >> 63);
    let l4 = (l4 << 1) | (l3 >> 63);
    let l3 = (l3 << 1) | (l2 >> 63);
    let l2 = (l2 << 1) | (l1 >> 63);
    let l1 = l1 << 1;

    let s = a[0] as u128 * a[0] as u128;
    let (l0, c) = (s as u64, (s >> 64) as u64);
    let (l1, c) = adb(l1, c);
    let (l2, c) = mac(l2, a[1], a[1], c);
    let (l3, c) = adb(l3, c);
    let (l4, c) = mac(l4, a[2], a[2], c);
    let (l5, c) = adb(l5, c);
    let (l6, c) = mac(l6, a[3], a[3], c);
    let l7 = l7.wrapping_add(c);

    mont([l0, l1, l2, l3, l4, l5, l6, l7], p, inv)
}

#[inline(always)]
pub(crate) const fn neg(a: [u64; 4], p: [u64; 4]) -> [u64; 4] {
    if (a[0] | a[1] | a[2] | a[3]) == 0 {
        a
    } else {
        let (l0, b) = sba(p[0], a[0]);
        let (l1, b) = sbb(p[1], a[1], b);
        let (l2, b) = sbb(p[2], a[2], b);
        let l3 = (p[3]).wrapping_sub(a[3]).wrapping_sub(b >> 63);

        [l0, l1, l2, l3]
    }
}

#[inline(always)]
pub(crate) const fn mont(a: [u64; 8], p: [u64; 4], inv: u64) -> [u64; 4] {
    let rhs = a[0].wrapping_mul(inv);
    let s = rhs as u128 * p[0] as u128 + a[0] as u128;
    let d = (s >> 64) as u64;
    let (l1, d) = mac(a[1], rhs, p[1], d);
    let (l2, d) = mac(a[2], rhs, p[2], d);
    let (l3, d) = mac(a[3], rhs, p[3], d);
    let (l4, e) = adb(a[4], d);

    let rhs = l1.wrapping_mul(inv);
    let s = rhs as u128 * p[0] as u128 + l1 as u128;
    let d = (s >> 64) as u64;
    let (l2, d) = mac(l2, rhs, p[1], d);
    let (l3, d) = mac(l3, rhs, p[2], d);
    let (l4, d) = mac(l4, rhs, p[3], d);
    let (l5, e) = adc(a[5], e, d);

    let rhs = l2.wrapping_mul(inv);
    let s = rhs as u128 * p[0] as u128 + l2 as u128;
    let d = (s >> 64) as u64;
    let (l3, d) = mac(l3, rhs, p[1], d);
    let (l4, d) = mac(l4, rhs, p[2], d);
    let (l5, d) = mac(l5, rhs, p[3], d);
    let (l6, e) = adc(a[6], e, d);

    let rhs = l3.wrapping_mul(inv);
    let s = rhs as u128 * p[0] as u128 + l3 as u128;
    let d = (s >> 64) as u64;
    let (l4, d) = mac(l4, rhs, p[1], d);
    let (l5, d) = mac(l5, rhs, p[2], d);
    let (l6, d) = mac(l6, rhs, p[3], d);
    let l7 = a[7].wrapping_add(e).wrapping_add(d);

    let (l0, brw) = sba(l4, p[0]);
    let (l1, brw) = sbb(l5, p[1], brw);
    let (l2, brw) = sbb(l6, p[2], brw);
    let (l3, brw) = sbb(l7, p[3], brw);

    let (l0, c) = adb(l0, p[0] & brw);
    let (l1, c) = adc(l1, p[1] & brw, c);
    let (l2, c) = adc(l2, p[2] & brw, c);
    let l3 = l3.wrapping_add(p[3] & brw).wrapping_add(c);

    [l0, l1, l2, l3]
}

#[inline(always)]
pub fn invert(
    a: [u64; 4],
    little_fermat: [u64; 4],
    identity: [u64; 4],
    p: [u64; 4],
    inv: u64,
) -> Option<[u64; 4]> {
    let zero: [u64; 4] = [0; 4];
    if a == zero {
        None
    } else {
        Some(pow(a, little_fermat, identity, p, inv))
    }
}

pub type Bits = Vec<u8>;

#[inline(always)]
pub fn to_bits(val: [u64; 4]) -> Bits {
    let mut index = 256;
    let mut bits: [u8; 256] = [0; 256];
    for limb in val {
        for byte in limb.to_le_bytes().iter() {
            for i in 0..8 {
                index -= 1;
                bits[index] = byte >> i & 1;
            }
        }
    }
    bits.to_vec()
}

#[inline(always)]
pub(crate) fn pow(
    a: [u64; 4],
    b: [u64; 4],
    mut identity: [u64; 4],
    p: [u64; 4],
    inv: u64,
) -> [u64; 4] {
    let zero: [u64; 4] = [0; 4];
    if b == zero {
        return identity;
    } else if a == zero {
        return zero;
    }
    let bits = to_bits(b);
    for &bit in bits.iter() {
        identity = square(identity, p, inv);
        if bit == 1 {
            identity = mul(identity, a, p, inv);
        }
    }
    identity
}

#[inline(always)]
pub(crate) const fn little_fermat(p: [u64; 4]) -> [u64; 4] {
    sub([0; 4], [2, 0, 0, 0], p)
}

#[inline(always)]
const fn from_u512(limbs: [u64; 8], r2: [u64; 4], r3: [u64; 4], p: [u64; 4], inv: u64) -> [u64; 4] {
    let a = mul([limbs[0], limbs[1], limbs[2], limbs[3]], r2, p, inv);
    let b = mul([limbs[4], limbs[5], limbs[6], limbs[7]], r3, p, inv);
    add(a, b, p)
}

#[inline(always)]
pub fn random_limbs(
    mut rand: impl RngCore,
    r2: [u64; 4],
    r3: [u64; 4],
    p: [u64; 4],
    inv: u64,
) -> [u64; 4] {
    from_u512(
        [
            rand.next_u64(),
            rand.next_u64(),
            rand.next_u64(),
            rand.next_u64(),
            rand.next_u64(),
            rand.next_u64(),
            rand.next_u64(),
            rand.next_u64(),
        ],
        r2,
        r3,
        p,
        inv,
    )
}

pub const fn from_u64(val: u64, r2: [u64; 4], p: [u64; 4], inv: u64) -> [u64; 4] {
    mul([val, 0, 0, 0], r2, p, inv)
}

#[derive(Clone, Copy, Eq, PartialEq, Debug)]
pub(crate) enum Naf {
    Zero = 0,
    Plus = 1,
    Minus = 2,
}

impl From<i8> for Naf {
    fn from(value: i8) -> Self {
        match value {
            0 => Naf::Zero,
            1 => Naf::Plus,
            -1 => Naf::Minus,
            _ => unimplemented!(),
        }
    }
}

pub(crate) type Nafs = Vec<Naf>;

#[inline(always)]
pub(crate) fn to_nafs(val: [u64; 4]) -> Nafs {
    let mut index = 0;
    let mut bits: [u8; 258] = [0; 258];
    for limb in val {
        for byte in limb.to_le_bytes().iter() {
            for i in 0..8 {
                bits[index] = byte >> i & 1;
                index += 1;
            }
        }
    }
    let mut carry = 0;
    let mut nafs: Vec<Naf> = bits
        .iter()
        .map(|bit| {
            let triple = bit * 3;
            let bit_3 = (triple + carry) % 2;
            carry = (triple + carry) / 2;
            (bit_3 as i8 - *bit as i8).into()
        })
        .collect::<Vec<_>>()
        .into_iter()
        .rev()
        .skip_while(|x| x == &Naf::Zero)
        .collect::<Vec<_>>();
    nafs.pop();
    nafs
}
