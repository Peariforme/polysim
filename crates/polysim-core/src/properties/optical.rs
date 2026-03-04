//! Optical property calculations.

use crate::error::PolySimError;
use crate::polymer::PolymerChain;
use crate::properties::group_contribution::{total_ri, GroupDatabase};
use crate::properties::mechanical::density;
use crate::properties::molecular_weight::average_mass;

/// Estimates the refractive index using the Lorentz-Lorenz equation.
///
/// The molar refraction `Rm` is summed from group contributions, and the
/// molar volume `V = M0 / rho` is obtained from the predicted density.
///
/// **n = sqrt((1 + 2*Rm/V) / (1 - Rm/V))**
///
/// # Errors
///
/// Returns `PolySimError::GroupDecomposition` if the SMILES cannot be decomposed.
///
/// # Reference
///
/// Van Krevelen, D. W. & te Nijenhuis, K. (2009).
/// *Properties of Polymers*, 4th ed., Elsevier. Chapter 11.
pub fn refractive_index(chain: &PolymerChain) -> Result<f64, PolySimError> {
    let groups = GroupDatabase::decompose(chain)?;

    let n_repeat = chain.repeat_count.max(1) as f64;

    // Molar refraction per repeat unit (cm^3/mol)
    let rm_total = total_ri(&groups);
    let rm = rm_total / n_repeat;

    // Molar volume per repeat unit: V = M0 / rho
    let rho = density(chain)?;
    let m0 = average_mass(chain) / n_repeat;
    let v = m0 / rho; // cm^3/mol (since rho is g/cm^3 and M0 is g/mol)

    if v < f64::EPSILON {
        return Err(PolySimError::GroupDecomposition(
            "molar volume is zero".into(),
        ));
    }

    let ratio = rm / v;

    // Guard against unphysical Rm/V >= 1 (would give negative under sqrt)
    if ratio >= 1.0 {
        return Err(PolySimError::GroupDecomposition(
            "Rm/V ratio >= 1, unphysical".into(),
        ));
    }

    // Lorentz-Lorenz: n = sqrt((1 + 2*r) / (1 - r))
    let n = ((1.0 + 2.0 * ratio) / (1.0 - ratio)).sqrt();

    Ok(n)
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
    fn ri_pe() {
        // PE: n exp ~ 1.49
        let chain = build_chain("{[]CC[]}", 50);
        let n = refractive_index(&chain).unwrap();
        assert!(
            (n - 1.49).abs() < 0.10,
            "PE refractive index = {n:.3}, expected ~1.49"
        );
    }

    #[test]
    fn ri_pvc() {
        let chain = build_chain("{[]C(Cl)C[]}", 50);
        let n = refractive_index(&chain).unwrap();
        assert!(
            n > 1.3 && n < 1.8,
            "PVC refractive index = {n:.3}, expected 1.3-1.8"
        );
    }

    #[test]
    fn ri_positive_and_greater_than_one() {
        let chain = build_chain("{[]CC[]}", 10);
        let n = refractive_index(&chain).unwrap();
        assert!(n > 1.0, "refractive index must be > 1.0, got {n:.3}");
    }

    #[test]
    fn ri_physical_range() {
        let polymers = [
            ("{[]CC[]}", "PE"),
            ("{[]CC(C)[]}", "PP"),
            ("{[]CC(c1ccccc1)[]}", "PS"),
            ("{[]CC(C)(C(=O)OC)[]}", "PMMA"),
            ("{[]C(Cl)C[]}", "PVC"),
        ];
        for (bigsmiles, name) in polymers {
            let chain = build_chain(bigsmiles, 50);
            let n = refractive_index(&chain).unwrap();
            assert!(
                n > 1.3 && n < 1.8,
                "{name}: refractive index = {n:.3}, out of [1.3, 1.8]"
            );
        }
    }
}
