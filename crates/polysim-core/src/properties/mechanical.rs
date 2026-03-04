//! Mechanical and density property calculations.
//!
//! All densities are in **g/cm^3**.

use crate::error::PolySimError;
use crate::polymer::PolymerChain;
use crate::properties::group_contribution::{total_rao, total_vw, GroupDatabase, GroupMatch};
use crate::properties::molecular_weight::average_mass;

/// Packing coefficient for aliphatic amorphous polymers (Van Krevelen, Table 4.3).
const K_ALIPHATIC: f64 = 0.681;

/// Packing coefficient for aromatic polymers like PS, PC (Van Krevelen, Table 4.3).
const K_AROMATIC: f64 = 0.895;

/// Computes a weighted packing coefficient based on the aromatic volume fraction.
///
/// VK Table 4.3 gives k=0.681 for aliphatic and k=0.895 for aromatic polymers
/// (pi-stacking of phenyl rings). The coefficient is linearly interpolated
/// based on the fraction of Van der Waals volume from aromatic groups.
fn packing_coefficient(groups: &[GroupMatch]) -> f64 {
    let vw_total: f64 = groups.iter().map(|gm| gm.group.vw * gm.count as f64).sum();
    let vw_aromatic: f64 = groups
        .iter()
        .filter(|gm| gm.group.name == "-C6H5" || gm.group.name == "-C6H4-")
        .map(|gm| gm.group.vw * gm.count as f64)
        .sum();
    let f_arom = if vw_total > f64::EPSILON {
        vw_aromatic / vw_total
    } else {
        0.0
    };
    K_ALIPHATIC + f_arom * (K_AROMATIC - K_ALIPHATIC)
}

/// Estimates the amorphous density (g/cm^3) using Van Krevelen group contributions.
///
/// The molar volume is estimated from the Van der Waals volume divided by
/// a packing coefficient that depends on the aromatic content:
///
/// **V = Vw / k**
///
/// where k is interpolated between 0.681 (aliphatic) and 0.895 (aromatic)
/// based on the aromatic volume fraction. The density is then:
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
/// *Properties of Polymers*, 4th ed., Elsevier. Chapter 4, Table 4.3.
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

    let k = packing_coefficient(&groups);
    let v_molar = vw_per_repeat / k;

    Ok(m0 / v_molar)
}

/// Estimates the Young's modulus (GPa) using the Rao function.
///
/// The Rao function (sound velocity increment) is summed from group
/// contributions and combined with the molar volume to estimate
/// the longitudinal sound velocity. The modulus is then:
///
/// **E = rho * U^2 * 10^-9**
///
/// where `U = Rao_total / V_molar` is the sound velocity (m/s), `rho` is the
/// density (kg/m^3), and the factor 10^-9 converts Pa to GPa.
///
/// # Errors
///
/// Returns `PolySimError::GroupDecomposition` if the SMILES cannot be decomposed.
///
/// # Reference
///
/// Van Krevelen, D. W. & te Nijenhuis, K. (2009).
/// *Properties of Polymers*, 4th ed., Elsevier. Chapter 13.
pub fn youngs_modulus(chain: &PolymerChain) -> Result<f64, PolySimError> {
    let groups = GroupDatabase::decompose(chain)?;
    let vw_total = total_vw(&groups);
    let rao_total = total_rao(&groups);

    let n = chain.repeat_count.max(1) as f64;
    let vw_per_repeat = vw_total / n;

    if vw_per_repeat < f64::EPSILON {
        return Err(PolySimError::GroupDecomposition(
            "Van der Waals volume per repeat unit is zero".into(),
        ));
    }

    let m_chain = average_mass(chain);
    let m0 = m_chain / n;
    let rao_per_repeat = rao_total / n;

    let k = packing_coefficient(&groups);
    // V_molar in cm^3/mol
    let v_molar = vw_per_repeat / k;

    // rho in g/cm^3
    let rho = m0 / v_molar;
    // Convert to kg/m^3: rho_si = rho * 1000
    let rho_si = rho * 1000.0;

    // Sound velocity U = Rao / V_molar (in cm/s units from VK)
    // Rao is in (cm^3/mol)(cm/s)^(1/3), V_molar in cm^3/mol
    // U = (Rao / V_molar)^3 in (cm/s)^1 -- VK convention
    // E = rho * U^2
    // U = (Rao / V)^3 in cm/s (VK convention: Rao in (cm^3/mol)(cm/s)^(1/3),
    // V in cm^3/mol). Convert to m/s by dividing by 100.
    let u_ratio = rao_per_repeat / v_molar;
    let u_cm_s = u_ratio.powi(3);
    let u_m_s = u_cm_s / 100.0;
    // E = rho (kg/m^3) * U^2 (m/s)^2 = Pa, convert to GPa
    let e_pa = rho_si * u_m_s * u_m_s;
    let e_gpa = e_pa / 1.0e9;

    Ok(e_gpa)
}

/// Estimates the tensile strength (MPa) from the Young's modulus.
///
/// Uses the empirical correlation:
///
/// **sigma = E / 35**
///
/// This gives a rough estimate; typical polymers have sigma_y ~ 30-100 MPa
/// and E ~ 1-4 GPa. The factor 35 is a typical strain-to-yield ratio
/// for glassy amorphous polymers.
///
/// # Errors
///
/// Returns `PolySimError::GroupDecomposition` if the SMILES cannot be decomposed.
///
/// # Reference
///
/// Van Krevelen, D. W. & te Nijenhuis, K. (2009).
/// *Properties of Polymers*, 4th ed., Elsevier. Chapter 13.
pub fn tensile_strength(chain: &PolymerChain) -> Result<f64, PolySimError> {
    let e_gpa = youngs_modulus(chain)?;
    // sigma in MPa: E(GPa) * 1000 / 35
    Ok(e_gpa * 1000.0 / 35.0)
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
        // VK Table 4.3 gives k=0.895 for aromatic polymers (pi-stacking).
        let chain = build_chain("{[]CC(c1ccccc1)[]}", 50);
        let rho = density(&chain).unwrap();
        assert!(
            (rho - 1.05).abs() < 0.15,
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

    #[test]
    fn youngs_modulus_pe() {
        // PE: E exp ~ 1.0 GPa (amorphous)
        let chain = build_chain("{[]CC[]}", 50);
        let e = youngs_modulus(&chain).unwrap();
        assert!(
            e > 0.1 && e < 10.0,
            "PE Young's modulus = {e:.2} GPa, expected 0.1-10 range"
        );
    }

    #[test]
    fn youngs_modulus_ps() {
        // PS: E exp ~ 3.0-3.5 GPa
        let chain = build_chain("{[]CC(c1ccccc1)[]}", 50);
        let e = youngs_modulus(&chain).unwrap();
        assert!(
            e > 0.5 && e < 15.0,
            "PS Young's modulus = {e:.2} GPa, expected 0.5-15 range"
        );
    }

    #[test]
    fn youngs_modulus_positive() {
        let chain = build_chain("{[]CC[]}", 10);
        let e = youngs_modulus(&chain).unwrap();
        assert!(e > 0.0, "Young's modulus must be positive, got {e}");
    }

    #[test]
    fn tensile_strength_positive() {
        let chain = build_chain("{[]CC[]}", 50);
        let sigma = tensile_strength(&chain).unwrap();
        assert!(
            sigma > 0.0,
            "tensile strength must be positive, got {sigma}"
        );
    }

    #[test]
    fn tensile_strength_physical_range() {
        // Typical polymer tensile strength: 10-200 MPa
        let polymers = [
            ("{[]CC[]}", "PE"),
            ("{[]CC(c1ccccc1)[]}", "PS"),
            ("{[]C(Cl)C[]}", "PVC"),
        ];
        for (bigsmiles, name) in polymers {
            let chain = build_chain(bigsmiles, 50);
            let sigma = tensile_strength(&chain).unwrap();
            assert!(
                sigma > 1.0 && sigma < 500.0,
                "{name}: tensile strength = {sigma:.1} MPa, out of physical range [1, 500]"
            );
        }
    }
}
