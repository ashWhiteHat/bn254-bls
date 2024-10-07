use crate::bn254::fq::Fq;
use crate::bn254::fq2::Fq2;
use crate::bn254::fr::Fr;
use crate::bn254::limbs::Naf;
use crate::bn254::params::{FROBENIUS_COEFF_FQ6_C1, SIX_U_PLUS_2_NAF, XI_TO_Q_MINUS_1_OVER_2};
use core::ops::{Add, AddAssign, Mul, Neg, SubAssign};

#[derive(Clone, Copy, Debug)]
pub struct G2Affine {
    x: Fq2,
    y: Fq2,
    is_infinity: bool,
}

pub(crate) const G2_GENERATOR_X: Fq2 = Fq2([
    Fq::to_mont_form([
        0x46debd5cd992f6ed,
        0x674322d4f75edadd,
        0x426a00665e5c4479,
        0x1800deef121f1e76,
    ]),
    Fq::to_mont_form([
        0x97e485b7aef312c2,
        0xf1aa493335a9e712,
        0x7260bfb731fb5d25,
        0x198e9393920d483a,
    ]),
]);

pub(crate) const G2_GENERATOR_Y: Fq2 = Fq2([
    Fq::to_mont_form([
        0x4ce6cc0166fa7daa,
        0xe3d1e7690c43d37b,
        0x4aab71808dcb408f,
        0x12c85ea5db8c6deb,
    ]),
    Fq::to_mont_form([
        0x55acdadcd122975b,
        0xbc4b313370b38ef3,
        0xec9e99ad690c3395,
        0x090689d0585ff075,
    ]),
]);

pub(crate) const G2_PARAM_B: Fq2 = Fq2([
    Fq::to_mont_form([
        0x3267e6dc24a138e5,
        0xb5b4c5e559dbefa3,
        0x81be18991be06ac3,
        0x2b149d40ceb8aaae,
    ]),
    Fq::to_mont_form([
        0xe4a2bd0685c315d2,
        0xa74fa084e52d1852,
        0xcd2cafadeed8fdf4,
        0x009713b03af0fed4,
    ]),
]);

impl G2Affine {
    fn is_identity(self) -> bool {
        self.is_infinity
    }

    pub const fn generator() -> Self {
        Self {
            x: G2_GENERATOR_X,
            y: G2_GENERATOR_Y,
            is_infinity: false,
        }
    }

    pub const fn identity() -> Self {
        Self {
            x: Fq2::zero(),
            y: Fq2::one(),
            is_infinity: true,
        }
    }

    pub(crate) fn param_3b() -> Fq2 {
        G2_PARAM_B.double() + G2_PARAM_B
    }

    pub const fn to_projective(self) -> G2Projective {
        G2Projective {
            x: self.x,
            y: self.y,
            z: Fq2::one(),
        }
    }
}

impl Neg for G2Affine {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            x: self.x,
            y: -self.y,
            is_infinity: self.is_infinity,
        }
    }
}

impl Mul<Fr> for G2Affine {
    type Output = G2Projective;

    fn mul(self, rhs: Fr) -> Self::Output {
        scalar_point(self, rhs)
    }
}

impl From<G2Projective> for G2Affine {
    fn from(value: G2Projective) -> Self {
        match value.z.invert() {
            Some(z_inv) => Self {
                x: value.x * z_inv,
                y: value.y * z_inv,
                is_infinity: false,
            },
            None => Self::identity(),
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct G2Projective {
    x: Fq2,
    y: Fq2,
    z: Fq2,
}

impl G2Projective {
    pub const fn identity() -> Self {
        Self {
            x: Fq2::zero(),
            y: Fq2::one(),
            z: Fq2::zero(),
        }
    }

    pub fn is_identity(self) -> bool {
        self.z == Fq2::zero()
    }

    pub(crate) fn double_eval(&mut self) -> PairingCoeff {
        // Adaptation of Algorithm 26, https://eprint.iacr.org/2010/354.pdf
        let tmp0 = self.x.square();
        let tmp1 = self.y.square();
        let tmp2 = tmp1.square();
        let tmp3 = (tmp1 + self.x).square() - tmp0 - tmp2;
        let tmp3 = tmp3.double();
        let tmp4 = tmp0.double() + tmp0;
        let tmp6 = self.x + tmp4;
        let tmp5 = tmp4.square();
        let zsquared = self.z.square();
        self.x = tmp5 - tmp3.double();
        self.z = (self.z + self.y).square() - tmp1 - zsquared;
        self.y = (tmp3 - self.x) * tmp4 - tmp2.double().double().double();
        let tmp3 = -(tmp4 * zsquared).double();
        let tmp6 = tmp6.square() - tmp0 - tmp5;
        let tmp1 = tmp1.double().double();
        let tmp6 = tmp6 - tmp1;
        let tmp0 = self.z * zsquared;
        let tmp0 = tmp0.double();

        PairingCoeff(tmp0, tmp3, tmp6)
    }

    pub(crate) fn add_eval(&mut self, rhs: G2Affine) -> PairingCoeff {
        // Adaptation of Algorithm 27, https://eprint.iacr.org/2010/354.pdf
        let zsquared = self.z.square();
        let ysquared = rhs.y.square();
        let t0 = zsquared * rhs.x;
        let t1 = ((rhs.y + self.z).square() - ysquared - zsquared) * zsquared;
        let t2 = t0 - self.x;
        let t3 = t2.square();
        let t4 = t3.double().double();
        let t5 = t4 * t2;
        let t6 = t1 - self.y.double();
        let t9 = t6 * rhs.x;
        let t7 = t4 * self.x;
        self.x = t6.square() - t5 - t7.double();
        self.z = (self.z + t2).square() - zsquared - t3;
        let t10 = rhs.y + self.z;
        let t8 = (t7 - self.x) * t6;
        let t0 = self.y * t5;
        self.y = t8 - t0.double();
        let t10 = t10.square() - ysquared;
        let ztsquared = self.z.square();
        let t10 = t10 - ztsquared;
        let t9 = t9.double() - t10;
        let t10 = self.z.double();
        let t1 = -t6.double();

        PairingCoeff(t10, t1, t9)
    }
}

impl Add for G2Projective {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        add_projective(self, rhs)
    }
}

impl AddAssign<G2Affine> for G2Projective {
    fn add_assign(&mut self, rhs: G2Affine) {
        *self = add_mixed(rhs, *self)
    }
}

impl SubAssign<G2Affine> for G2Projective {
    fn sub_assign(&mut self, rhs: G2Affine) {
        *self = add_mixed(-rhs, *self)
    }
}

impl From<G2Affine> for G2Projective {
    fn from(affine: G2Affine) -> Self {
        if affine.is_identity() {
            Self {
                x: Fq2::zero(),
                y: Fq2::one(),
                z: Fq2::zero(),
            }
        } else {
            Self {
                x: affine.x,
                y: affine.y,
                z: Fq2::one(),
            }
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct PairingCoeff(pub(crate) Fq2, pub(crate) Fq2, pub(crate) Fq2);

pub struct G2PairingAffine {
    pub(crate) coeffs: Vec<PairingCoeff>,
    is_infinity: bool,
}

impl G2PairingAffine {
    pub fn is_identity(&self) -> bool {
        self.is_infinity
    }
}

impl From<G2Affine> for G2PairingAffine {
    fn from(g2: G2Affine) -> Self {
        if g2.is_identity() {
            Self {
                coeffs: vec![],
                is_infinity: true,
            }
        } else {
            let mut coeffs = vec![];
            let mut g2_projective = G2Projective::from(g2);
            let neg = -g2;

            for i in (1..SIX_U_PLUS_2_NAF.len()).rev() {
                coeffs.push(g2_projective.double_eval());
                let x = SIX_U_PLUS_2_NAF[i - 1];
                match x {
                    1 => {
                        coeffs.push(g2_projective.add_eval(g2));
                    }
                    -1 => {
                        coeffs.push(g2_projective.add_eval(neg));
                    }
                    _ => continue,
                }
            }

            let mut q = g2;

            q.x.0[1] = -q.x.0[1];
            q.x *= FROBENIUS_COEFF_FQ6_C1[1];

            q.y.0[1] = -q.y.0[1];
            q.y *= XI_TO_Q_MINUS_1_OVER_2;

            coeffs.push(g2_projective.add_eval(q));

            let mut minusq2 = g2;
            minusq2.x *= FROBENIUS_COEFF_FQ6_C1[2];

            coeffs.push(g2_projective.add_eval(minusq2));

            Self {
                coeffs,
                is_infinity: false,
            }
        }
    }
}

/// add affine
#[inline(always)]
pub fn add_affine(lhs: G2Affine, rhs: G2Affine) -> G2Projective {
    if lhs.is_identity() {
        return rhs.to_projective();
    } else if rhs.is_identity() {
        return lhs.to_projective();
    }

    let (x0, y0) = (lhs.x, lhs.y);
    let (x1, y1) = (rhs.x, rhs.y);

    if x0 == x1 {
        if y0 == y1 {
            return double_affine(lhs);
        } else {
            return G2Projective::identity();
        }
    }

    let s = y0 - y1;
    let u = x0 - x1;
    let uu = u.square();
    let w = s.square() - uu * (x0 + x1);
    let uuu = uu * u;

    let x = u * w;
    let y = s * (x0 * uu - w) - y0 * uuu;
    let z = uuu;

    G2Projective { x, y, z }
}

/// double affine
#[inline(always)]
pub fn double_affine(point: G2Affine) -> G2Projective {
    // Algorithm 9, https://eprint.iacr.org/2015/1060.pdf
    let b3 = G2Affine::param_3b();
    let (x, y) = (point.x, point.y);

    let t0 = y.square();
    let z3 = t0.double().double().double();
    let x3 = b3 * z3;
    let y3 = t0 + b3;
    let z3 = y * z3;
    let t1 = b3.double();
    let t2 = t1 + b3;
    let t0 = t0 - t2;
    let y3 = t0 * y3;
    let y3 = x3 + y3;
    let t1 = x * y;
    let x3 = t0 * t1;
    let x3 = x3.double();

    G2Projective {
        x: x3,
        y: y3,
        z: z3,
    }
}

/// add affine to projective
#[inline(always)]
pub fn add_mixed(lhs: G2Affine, rhs: G2Projective) -> G2Projective {
    if lhs.is_identity() {
        return rhs;
    } else if rhs.is_identity() {
        return lhs.to_projective();
    }

    let (x1, y1, z1) = (rhs.x, rhs.y, rhs.z);
    let (x2, y2) = (lhs.x, lhs.y);

    let s1 = y2 * z1;
    let u1 = x2 * z1;

    if u1 == x1 {
        if s1 == y1 {
            return double_affine(lhs);
        } else {
            return G2Projective::identity();
        }
    }

    let u = s1 - y1;
    let uu = u.square();
    let v = u1 - x1;
    let vv = v.square();
    let vvv = vv * v;
    let r = vv * x1;
    let a = uu * z1 - vvv - r.double();

    let x = v * a;
    let y = u * (r - a) - vvv * y1;
    let z = vvv * z1;

    G2Projective { x, y, z }
}

/// add projective
#[inline(always)]
pub fn add_projective(lhs: G2Projective, rhs: G2Projective) -> G2Projective {
    if lhs.is_identity() {
        return rhs;
    } else if rhs.is_identity() {
        return lhs;
    }

    let (x0, y0, z0) = (lhs.x, lhs.y, lhs.z);
    let (x1, y1, z1) = (rhs.x, rhs.y, rhs.z);

    let s1 = y0 * z1;
    let s2 = y1 * z0;
    let u1 = x0 * z1;
    let u2 = x1 * z0;

    if u1 == u2 {
        if s1 == s2 {
            return double_projective(lhs);
        } else {
            return G2Projective::identity();
        }
    }

    let s = s1 - s2;
    let u = u1 - u2;
    let uu = u.square();
    let v = z0 * z1;
    let w = s.square() * v - uu * (u1 + u2);
    let uuu = uu * u;

    let x = u * w;
    let y = s * (u1 * uu - w) - s1 * uuu;
    let z = uuu * v;

    G2Projective { x, y, z }
}

/// double projective
#[inline(always)]
pub fn double_projective(lhs: G2Projective) -> G2Projective {
    // Algorithm 9, https://eprint.iacr.org/2015/1060.pdf
    let (x, y, z) = (lhs.x, lhs.y, lhs.z);
    let b3 = G2Affine::param_3b();

    let t0 = y.square();
    let z3 = t0.double().double().double();
    let t1 = y * z;
    let t2 = z.square();
    let t2 = t2 * b3;
    let x3 = t2 * z3;
    let y3 = t0 + t2;
    let z3 = t1 * z3;
    let t1 = t2.double();
    let t2 = t1 + t2;
    let t0 = t0 - t2;
    let y3 = t0 * y3;
    let y3 = x3 + y3;
    let t1 = x * y;
    let x3 = t0 * t1;
    let x3 = x3.double();

    G2Projective {
        x: x3,
        y: y3,
        z: z3,
    }
}

/// scalar
#[inline(always)]
pub fn scalar_point(point: G2Affine, scalar: Fr) -> G2Projective {
    let mut res = G2Projective::identity();
    for &naf in scalar.to_nafs().iter() {
        res = double_projective(res);
        if naf == Naf::Plus {
            res += point;
        } else if naf == Naf::Minus {
            res -= point;
        }
    }
    res
}
