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

pub use g1::G1Affine;
pub use g2::G2Affine;
pub use gt::Gt;
pub use pairing::AteParing;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generator_pairing() {
        let g1 = G1Affine::generator();
        let g2 = G2Affine::generator();
        let gt = Gt::generator();

        assert_eq!(gt, AteParing::pairing(g1, g2));
    }
}
