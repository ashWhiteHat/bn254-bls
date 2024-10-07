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
    use rand_core::OsRng;

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

    #[test]
    fn test_normal_pairing() {
        let g1 = G1Affine::generator();
        let g2 = G2Affine::generator();
        let mut rng = OsRng;

        for _ in 0..10 {
            let a = Fr::random(&mut rng);
            let b = Fr::random(&mut rng);
            let c = a * b;

            let g = G1Affine::from(g1 * a);
            let h = G2Affine::from(g2 * b);
            let p = AteParing::pairing(g, h);

            let expected = G1Affine::from(g1 * c);
            let test = G2Affine::from(g2 * c);

            assert_eq!(p, AteParing::pairing(expected, g2));
            assert_eq!(p, AteParing::pairing(g1, test));
        }
    }
}
