mod fq;
mod fq12;
mod fq2;
mod fq6;
mod fr;
mod g1;
mod g2;
mod gt;
mod limbs;
mod math;
mod pairing;
mod params;

pub use fq12::Fq12;
pub use fr::Fr;
pub use g1::G1Affine;
pub use g2::G2Affine;
pub use gt::Gt;
pub use pairing::AteParing;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_final_exp() {
        assert_eq!(Fq12::one().final_exp(), Gt::identity());
    }

    #[test]
    fn test_generator_pairing() {
        let g1 = G1Affine::generator();
        let g2 = G2Affine::generator();
        let gt = Gt::generator();

        assert_eq!(gt, AteParing::pairing(g1, g2));
    }

    #[test]
    fn test_signed_pairing() {
        let g = G1Affine::generator();
        let h = G2Affine::generator();

        let p = -AteParing::pairing(g, h);
        let q = AteParing::pairing(g, -h);
        let r = AteParing::pairing(-g, h);

        assert_eq!(p, q);
        assert_eq!(q, r);
    }
}
