use crate::fq::Fq;
use crate::fq12::Fq12;
use crate::fq2::Fq2;
use crate::g1::G1Affine;
use crate::g2::{G2Affine, G2PairingAffine};
use crate::gt::Gt;
use crate::params::SIX_U_PLUS_2_NAF;

/// Ate pairing struct holds necessary components for pairing.
/// `pairing` function takes G1 and G2 group elements and output
/// GT target group element.
pub struct AteParing;

impl AteParing {
    pub fn pairing(g1: G1Affine, g2: G2Affine) -> Gt {
        let g2 = G2PairingAffine::from(g2);
        Self::multi_miller_loop(&[(g1, g2)]).final_exp()
    }

    pub fn multi_miller_loop(pairs: &[(G1Affine, G2PairingAffine)]) -> Fq12 {
        let mut pairs = pairs
            .iter()
            .filter(|(a, b)| !(a.is_identity()) && !b.is_identity())
            .map(|(g1, g2)| (g1, g2.coeffs.iter()))
            .collect::<Vec<_>>();

        let mut acc = Fq12::one();

        for i in (1..SIX_U_PLUS_2_NAF.len()).rev() {
            if i != SIX_U_PLUS_2_NAF.len() - 1 {
                acc.square_assign();
            }
            for &mut (p, ref mut coeffs) in &mut pairs {
                acc = acc.untwist(*coeffs.next().unwrap(), *p);
            }
            let x = SIX_U_PLUS_2_NAF[i - 1];
            match x {
                1 => {
                    for &mut (p, ref mut coeffs) in &mut pairs {
                        acc = acc.untwist(*coeffs.next().unwrap(), *p);
                    }
                }
                -1 => {
                    for &mut (p, ref mut coeffs) in &mut pairs {
                        acc = acc.untwist(*coeffs.next().unwrap(), *p);
                    }
                }
                _ => continue,
            }
        }

        for &mut (p, ref mut coeffs) in &mut pairs {
            acc = acc.untwist(*coeffs.next().unwrap(), *p);
        }

        for &mut (p, ref mut coeffs) in &mut pairs {
            acc = acc.untwist(*coeffs.next().unwrap(), *p);
        }

        for &mut (_p, ref mut coeffs) in &mut pairs {
            assert_eq!(coeffs.next(), None);
        }

        acc
    }
}
