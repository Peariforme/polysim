//! Solubility parameter calculations.
//!
//! All solubility parameters are in **(MPa)^0.5**.

use crate::error::PolySimError;
use crate::polymer::PolymerChain;
use crate::properties::group_contribution::{total_ecoh, total_vw, GroupDatabase};

/// Hildebrand solubility parameter (MPa)^0.5.
///
/// Computed as:
///
/// **delta = sqrt(Ecoh / Vw)**
///
/// where `Ecoh` is the total cohesive energy (J/mol) and `Vw` is the total
/// Van der Waals volume (cm^3/mol) from group contributions. The result is
/// converted from (J/cm^3)^0.5 to (MPa)^0.5 (numerically equivalent since
/// 1 J/cm^3 = 1 MPa).
///
/// # Errors
///
/// Returns `PolySimError::GroupDecomposition` if the SMILES cannot be decomposed.
///
/// # Reference
///
/// Hildebrand, J. H. & Scott, R. L. (1950). *The Solubility of
/// Non-Electrolytes*, Reinhold.
///
/// Van Krevelen, D. W. (1990). *Properties of Polymers*, 3rd ed.,
/// Elsevier. Chapter 7.
pub fn hildebrand_solubility_parameter(chain: &PolymerChain) -> Result<f64, PolySimError> {
    let groups = GroupDatabase::decompose(chain)?;

    let ecoh = total_ecoh(&groups);
    let vw = total_vw(&groups);

    if vw < f64::EPSILON {
        return Err(PolySimError::GroupDecomposition(
            "Van der Waals volume is zero".into(),
        ));
    }

    // Ecoh is in J/mol, Vw is in cm^3/mol.
    // delta = sqrt(Ecoh / Vw) in (J/cm^3)^0.5 = (MPa)^0.5
    Ok((ecoh / vw).sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::{linear::LinearBuilder, BuildStrategy};

    fn build_chain(bigsmiles: &str, n: usize) -> PolymerChain {
        let bs = crate::parse(bigsmiles).expect("valid BigSMILES");
        LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(n))
            .homopolymer()
            .expect("build should succeed")
    }

    #[test]
    fn hildebrand_pe() {
        // PE: delta exp ~ 16.2 (MPa)^0.5
        let chain = build_chain("{[]CC[]}", 50);
        let delta = hildebrand_solubility_parameter(&chain).unwrap();
        // PE: all CH2 groups, Ecoh=4100, Vw=10.23 per CH2
        // delta = sqrt(4100/10.23) = sqrt(400.8) = 20.0
        // With terminal CH3: slightly different
        assert!(
            delta > 15.0 && delta < 25.0,
            "PE delta = {delta:.1}, expected ~16-22 range"
        );
    }

    #[test]
    fn hildebrand_pvc() {
        // PVC: delta exp ~ 19.5 (MPa)^0.5
        let chain = build_chain("{[]C(Cl)C[]}", 50);
        let delta = hildebrand_solubility_parameter(&chain).unwrap();
        assert!(delta > 15.0 && delta < 30.0, "PVC delta = {delta:.1}");
    }

    #[test]
    fn hildebrand_positive() {
        let chain = build_chain("{[]CC[]}", 10);
        let delta = hildebrand_solubility_parameter(&chain).unwrap();
        assert!(delta > 0.0, "delta should be positive, got {delta}");
    }

    #[test]
    fn hildebrand_independent_of_chain_length() {
        // Solubility parameter is intensive -- should be similar for
        // different chain lengths (up to end-group effects).
        let chain_short = build_chain("{[]CC[]}", 10);
        let chain_long = build_chain("{[]CC[]}", 100);
        let d_short = hildebrand_solubility_parameter(&chain_short).unwrap();
        let d_long = hildebrand_solubility_parameter(&chain_long).unwrap();
        let diff = (d_short - d_long).abs();
        assert!(
            diff < 2.0,
            "delta should be similar: short={d_short:.2}, long={d_long:.2}"
        );
    }
}
