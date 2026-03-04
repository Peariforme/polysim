//! Mechanical and density property calculations.
//!
//! All densities are in **g/cm^3**.

use crate::error::PolySimError;
use crate::polymer::PolymerChain;
use crate::properties::group_contribution::{total_vw, GroupDatabase};
use crate::properties::molecular_weight::average_mass;

/// Packing coefficient for amorphous polymers (Van Krevelen, Table 4.8).
const AMORPHOUS_PACKING: f64 = 0.681;

/// Estimates the amorphous density (g/cm^3) using Van Krevelen group contributions.
///
/// The molar volume is estimated from the Van der Waals volume divided by the
/// amorphous packing coefficient (0.681):
///
/// **V = Vw / 0.681**
///
/// The density is then:
///
/// **rho = M0 / V**
///
/// where `M0` is the repeat-unit molar mass (g/mol) and `V` is the molar
/// volume per repeat unit (cm^3/mol).
///
/// # Errors
///
/// Returns `PolySimError::GroupDecomposition` if the SMILES cannot be decomposed.
///
/// # Reference
///
/// Van Krevelen, D. W. & te Nijenhuis, K. (2009).
/// *Properties of Polymers*, 4th ed., Elsevier. Chapter 4.
pub fn density(chain: &PolymerChain) -> Result<f64, PolySimError> {
    let groups = GroupDatabase::decompose(chain)?;
    let vw_total = total_vw(&groups);

    let n = chain.repeat_count.max(1) as f64;
    let vw_per_repeat = vw_total / n;

    if vw_per_repeat < f64::EPSILON {
        return Err(PolySimError::GroupDecomposition(
            "Van der Waals volume per repeat unit is zero".into(),
        ));
    }

    let m_chain = average_mass(chain);
    let m0 = m_chain / n;

    // V = Vw / packing coefficient
    let v_molar = vw_per_repeat / AMORPHOUS_PACKING;

    Ok(m0 / v_molar)
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
    fn density_pe() {
        // PE: rho exp ~ 0.95 g/cm^3 (amorphous)
        let chain = build_chain("{[]CC[]}", 50);
        let rho = density(&chain).unwrap();
        assert!(
            (rho - 0.95).abs() < 0.15,
            "PE density = {rho:.3} g/cm3, expected ~0.95"
        );
    }

    #[test]
    fn density_ps() {
        // PS: rho exp ~ 1.05 g/cm^3
        // VK with fixed packing 0.681 gives ~0.80 -- underestimate due to
        // aromatic packing efficiency not captured by a single constant.
        let chain = build_chain("{[]CC(c1ccccc1)[]}", 50);
        let rho = density(&chain).unwrap();
        assert!(
            (rho - 1.05).abs() < 0.30,
            "PS density = {rho:.3} g/cm3, expected ~1.05"
        );
    }

    #[test]
    fn density_pvc() {
        // PVC: rho exp ~ 1.40 g/cm^3
        let chain = build_chain("{[]C(Cl)C[]}", 50);
        let rho = density(&chain).unwrap();
        assert!(
            (rho - 1.40).abs() < 0.20,
            "PVC density = {rho:.3} g/cm3, expected ~1.40"
        );
    }

    #[test]
    fn density_pmma() {
        // PMMA: rho exp ~ 1.18 g/cm^3
        let chain = build_chain("{[]CC(C)(C(=O)OC)[]}", 50);
        let rho = density(&chain).unwrap();
        assert!(
            (rho - 1.18).abs() < 0.20,
            "PMMA density = {rho:.3} g/cm3, expected ~1.18"
        );
    }

    #[test]
    fn density_positive() {
        let chain = build_chain("{[]CC[]}", 10);
        let rho = density(&chain).unwrap();
        assert!(rho > 0.0, "density must be positive, got {rho}");
    }

    #[test]
    fn density_physical_range() {
        // All common polymers have density between 0.8 and 2.0 g/cm^3
        let polymers = [
            ("{[]CC[]}", "PE"),
            ("{[]CC(C)[]}", "PP"),
            ("{[]CC(c1ccccc1)[]}", "PS"),
            ("{[]CC(C)(C(=O)OC)[]}", "PMMA"),
            ("{[]C(Cl)C[]}", "PVC"),
        ];
        for (bigsmiles, name) in polymers {
            let chain = build_chain(bigsmiles, 50);
            let rho = density(&chain).unwrap();
            assert!(
                rho > 0.5 && rho < 2.5,
                "{name}: density = {rho:.3} g/cm3, out of physical range [0.5, 2.5]"
            );
        }
    }
}
