use crate::fq::Fq;
use crate::fr::Fr;
use crate::limbs::Naf;
use core::ops::{Add, AddAssign, Mul, Neg, SubAssign};

pub(crate) const G1_GENERATOR_X: Fq = Fq::one();
pub(crate) const G1_GENERATOR_Y: Fq = Fq::to_mont_form([2, 0, 0, 0]);
pub(crate) const G1_PARAM_B: Fq = Fq::to_mont_form([3, 0, 0, 0]);

#[derive(Clone, Copy, Debug)]
pub struct G1Affine {
    pub(crate) x: Fq,
    pub(crate) y: Fq,
    is_infinity: bool,
}

impl G1Affine {
    pub(crate) fn is_identity(self) -> bool {
        self.is_infinity
    }

    pub const fn generator() -> Self {
        Self {
            x: G1_GENERATOR_X,
            y: G1_GENERATOR_Y,
            is_infinity: false,
        }
    }

    pub const fn identity() -> Self {
        Self {
            x: Fq::zero(),
            y: Fq::one(),
            is_infinity: true,
        }
    }

    pub(crate) fn param_3b() -> Fq {
        G1_PARAM_B.double() + G1_PARAM_B
    }

    pub const fn to_projective(self) -> G1Projective {
        G1Projective {
            x: self.x,
            y: self.y,
            z: Fq::one(),
        }
    }
}

impl Neg for G1Affine {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            x: self.x,
            y: -self.y,
            is_infinity: self.is_infinity,
        }
    }
}

impl Add for G1Affine {
    type Output = G1Projective;

    fn add(self, rhs: G1Affine) -> Self::Output {
        add_affine(self, rhs)
    }
}

impl Add<G1Projective> for G1Affine {
    type Output = G1Projective;

    fn add(self, rhs: G1Projective) -> Self::Output {
        add_mixed(self, rhs)
    }
}

impl Mul<Fr> for G1Affine {
    type Output = G1Projective;

    fn mul(self, rhs: Fr) -> Self::Output {
        scalar_point(self, rhs)
    }
}

impl From<G1Projective> for G1Affine {
    fn from(value: G1Projective) -> Self {
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

#[derive(Debug, Clone, Copy)]
pub struct G1Projective {
    pub(crate) x: Fq,
    pub(crate) y: Fq,
    pub(crate) z: Fq,
}

impl G1Projective {
    pub const fn identity() -> Self {
        Self {
            x: Fq::zero(),
            y: Fq::one(),
            z: Fq::zero(),
        }
    }

    pub fn is_identity(self) -> bool {
        self.z == Fq::zero()
    }
}

impl Add for G1Projective {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        add_projective(self, rhs)
    }
}

impl AddAssign<G1Affine> for G1Projective {
    fn add_assign(&mut self, rhs: G1Affine) {
        *self = add_mixed(rhs, *self)
    }
}

impl SubAssign<G1Affine> for G1Projective {
    fn sub_assign(&mut self, rhs: G1Affine) {
        *self = add_mixed(-rhs, *self)
    }
}

/// add affine
#[inline(always)]
pub fn add_affine(lhs: G1Affine, rhs: G1Affine) -> G1Projective {
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
            return G1Projective::identity();
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

    G1Projective { x, y, z }
}

/// double affine
#[inline(always)]
pub fn double_affine(point: G1Affine) -> G1Projective {
    // Algorithm 9, https://eprint.iacr.org/2015/1060.pdf
    let b3 = G1Affine::param_3b();
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

    G1Projective {
        x: x3,
        y: y3,
        z: z3,
    }
}

/// add affine to projective
#[inline(always)]
pub fn add_mixed(lhs: G1Affine, rhs: G1Projective) -> G1Projective {
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
            return G1Projective::identity();
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

    G1Projective { x, y, z }
}

/// add projective
#[inline(always)]
pub fn add_projective(lhs: G1Projective, rhs: G1Projective) -> G1Projective {
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
            return G1Projective::identity();
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

    G1Projective { x, y, z }
}

/// double projective
#[inline(always)]
pub fn double_projective(lhs: G1Projective) -> G1Projective {
    // Algorithm 9, https://eprint.iacr.org/2015/1060.pdf
    let (x, y, z) = (lhs.x, lhs.y, lhs.z);
    let b3 = G1Affine::param_3b();

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

    G1Projective {
        x: x3,
        y: y3,
        z: z3,
    }
}

/// scalar
#[inline(always)]
pub fn scalar_point(point: G1Affine, scalar: Fr) -> G1Projective {
    let mut res = G1Projective::identity();
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
