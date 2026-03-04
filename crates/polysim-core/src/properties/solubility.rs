//! Solubility parameter calculations.
//!
//! All solubility parameters are in **(MPa)^0.5**.

use crate::error::PolySimError;
use crate::polymer::PolymerChain;
use crate::properties::group_contribution::{
    total_ecoh, total_ed, total_eh, total_ep, total_vw, GroupDatabase,
};

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

    // Vmol = Vw / 0.667 for flexible-chain polymers (VK Table 4.3).
    // delta = sqrt(Ecoh / Vmol) in (J/cm^3)^0.5 = (MPa)^0.5
    let vmol = vw / 0.667;
    Ok((ecoh / vmol).sqrt())
}

/// Hansen solubility parameters in (MPa)^0.5.
#[derive(Debug, Clone, PartialEq)]
pub struct HansenParams {
    /// Dispersive component (MPa)^0.5.
    pub delta_d: f64,
    /// Polar component (MPa)^0.5.
    pub delta_p: f64,
    /// Hydrogen-bonding component (MPa)^0.5.
    pub delta_h: f64,
    /// Total solubility parameter (MPa)^0.5 = sqrt(delta_d^2 + delta_p^2 + delta_h^2).
    pub delta_t: f64,
}

/// Estimates the Hansen solubility parameters from group contributions.
///
/// Each component is computed as:
///
/// - **delta_d = sqrt(Ed / V)**
/// - **delta_p = sqrt(Ep / V)**
/// - **delta_h = sqrt(Eh / V)**
/// - **delta_t = sqrt(delta_d^2 + delta_p^2 + delta_h^2)**
///
/// where `Ed`, `Ep`, `Eh` are the dispersive, polar, and hydrogen-bonding
/// cohesive energy contributions (J/mol), and `V` is the Van der Waals
/// volume (cm^3/mol).
///
/// # Errors
///
/// Returns `PolySimError::GroupDecomposition` if the SMILES cannot be decomposed.
///
/// # Reference
///
/// Hansen, C. M. (2007). *Hansen Solubility Parameters: A User's Handbook*,
/// 2nd ed., CRC Press.
///
/// Van Krevelen, D. W. & te Nijenhuis, K. (2009).
/// *Properties of Polymers*, 4th ed., Elsevier. Chapter 7.
pub fn hansen_solubility_parameters(chain: &PolymerChain) -> Result<HansenParams, PolySimError> {
    let groups = GroupDatabase::decompose(chain)?;
    let vw = total_vw(&groups);

    if vw < f64::EPSILON {
        return Err(PolySimError::GroupDecomposition(
            "Van der Waals volume is zero".into(),
        ));
    }

    let ed = total_ed(&groups);
    let ep = total_ep(&groups);
    let eh = total_eh(&groups);

    // Vmol = Vw / 0.667 for flexible-chain polymers (VK Table 4.3).
    let vmol = vw / 0.667;
    let delta_d = (ed / vmol).sqrt();
    let delta_p = (ep / vmol).sqrt();
    let delta_h = (eh / vmol).sqrt();
    let delta_t = (delta_d.powi(2) + delta_p.powi(2) + delta_h.powi(2)).sqrt();

    Ok(HansenParams {
        delta_d,
        delta_p,
        delta_h,
        delta_t,
    })
}

/// Relative Energy Difference (RED) between a polymer and a solvent.
///
/// **RED = Ra / R0**
///
/// where `Ra = sqrt(4*(delta_d1 - delta_d2)^2 + (delta_p1 - delta_p2)^2 + (delta_h1 - delta_h2)^2)`
/// is the Hansen distance and `R0` is the interaction radius of the solvent sphere.
///
/// - RED < 1: polymer is likely soluble in the solvent.
/// - RED > 1: polymer is likely insoluble.
pub fn red_distance(polymer: &HansenParams, solvent: &HansenParams, r0: f64) -> f64 {
    let ra = ((4.0 * (polymer.delta_d - solvent.delta_d).powi(2))
        + (polymer.delta_p - solvent.delta_p).powi(2)
        + (polymer.delta_h - solvent.delta_h).powi(2))
    .sqrt();
    ra / r0
}

/// Gas constant R in J/(mol*K).
const R_GAS: f64 = 8.314;

/// Reference molar volume for the Flory-Huggins interaction parameter (cm^3/mol).
///
/// By convention, Vr = 100 cm^3/mol is used as the lattice-site volume in the
/// Flory-Huggins theory for polymer blends.
const VR: f64 = 100.0;

/// Estimates the Flory-Huggins interaction parameter chi between two polymers.
///
/// The Hildebrand-based approximation is:
///
/// **chi = Vr * (delta_A - delta_B)^2 / (R * T)**
///
/// where `delta_A` and `delta_B` are the Hildebrand solubility parameters
/// in (MPa)^0.5 = (J/cm^3)^0.5, `Vr` = 100 cm^3/mol (reference volume),
/// `R` = 8.314 J/(mol*K), and `T` is the temperature in Kelvin.
///
/// - chi < 0.5: polymers are likely miscible.
/// - chi > 0.5: polymers are likely immiscible.
///
/// # Errors
///
/// Returns `PolySimError::GroupDecomposition` if either SMILES cannot be decomposed.
///
/// # Reference
///
/// Flory, P. J. (1953). *Principles of Polymer Chemistry*, Cornell University Press.
///
/// Van Krevelen, D. W. & te Nijenhuis, K. (2009).
/// *Properties of Polymers*, 4th ed., Elsevier. Chapter 7.
pub fn flory_huggins_chi(
    polymer_a: &PolymerChain,
    polymer_b: &PolymerChain,
    temperature: f64,
) -> Result<f64, PolySimError> {
    let delta_a = hildebrand_solubility_parameter(polymer_a)?;
    let delta_b = hildebrand_solubility_parameter(polymer_b)?;
    // delta in (MPa)^0.5 = (J/cm^3)^0.5, so (delta_a - delta_b)^2 is J/cm^3.
    // Vr in cm^3/mol, R in J/(mol*K), T in K.
    // chi = Vr * (delta_a - delta_b)^2 / (R * T) is dimensionless.
    let chi = VR * (delta_a - delta_b).powi(2) / (R_GAS * temperature);
    Ok(chi)
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

    #[test]
    fn hansen_pe_mostly_dispersive() {
        let chain = build_chain("{[]CC[]}", 50);
        let h = hansen_solubility_parameters(&chain).unwrap();
        assert!(h.delta_d > 15.0, "PE delta_d = {:.1}", h.delta_d);
        assert!(
            h.delta_p < 1.0,
            "PE delta_p should be ~0, got {:.1}",
            h.delta_p
        );
        assert!(
            h.delta_h < 1.0,
            "PE delta_h should be ~0, got {:.1}",
            h.delta_h
        );
        let hild = hildebrand_solubility_parameter(&chain).unwrap();
        assert!(
            (h.delta_t - hild).abs() < 1.0,
            "PE Hansen delta_t={:.1} should match Hildebrand={:.1}",
            h.delta_t,
            hild
        );
    }

    #[test]
    fn hansen_pvc_has_polar_component() {
        let chain = build_chain("{[]C(Cl)C[]}", 50);
        let h = hansen_solubility_parameters(&chain).unwrap();
        assert!(h.delta_d > 10.0, "PVC delta_d = {:.1}", h.delta_d);
        assert!(
            h.delta_p > 1.0,
            "PVC should have polar component, got {:.1}",
            h.delta_p
        );
    }

    #[test]
    fn hansen_components_sum_to_total() {
        let chain = build_chain("{[]C(Cl)C[]}", 50);
        let h = hansen_solubility_parameters(&chain).unwrap();
        let sum = (h.delta_d.powi(2) + h.delta_p.powi(2) + h.delta_h.powi(2)).sqrt();
        assert!(
            (h.delta_t - sum).abs() < 0.01,
            "delta_t={:.2} should equal sqrt(d^2+p^2+h^2)={sum:.2}",
            h.delta_t
        );
    }

    #[test]
    fn red_distance_same_polymer_is_zero() {
        let chain = build_chain("{[]CC[]}", 50);
        let h = hansen_solubility_parameters(&chain).unwrap();
        let red = red_distance(&h, &h, 5.0);
        assert!(
            red < 0.01,
            "RED of polymer with itself should be ~0, got {red:.3}"
        );
    }

    #[test]
    fn flory_huggins_same_polymer_is_zero() {
        let pe = build_chain("{[]CC[]}", 50);
        let chi = flory_huggins_chi(&pe, &pe, 298.0).unwrap();
        assert!(
            chi < 0.001,
            "chi of polymer with itself should be ~0, got {chi:.4}"
        );
    }

    #[test]
    fn flory_huggins_pe_pvc_larger_than_similar_pair() {
        // PE/PVC should have larger chi than PE/PP (more different delta)
        let pe = build_chain("{[]CC[]}", 50);
        let pp = build_chain("{[]CC(C)[]}", 50);
        let pvc = build_chain("{[]C(Cl)C[]}", 50);
        let chi_pe_pp = flory_huggins_chi(&pe, &pp, 298.0).unwrap();
        let chi_pe_pvc = flory_huggins_chi(&pe, &pvc, 298.0).unwrap();
        assert!(
            chi_pe_pvc > chi_pe_pp,
            "PE/PVC chi={chi_pe_pvc:.3} should be > PE/PP chi={chi_pe_pp:.3}"
        );
    }

    #[test]
    fn flory_huggins_positive() {
        let pe = build_chain("{[]CC[]}", 50);
        let pvc = build_chain("{[]C(Cl)C[]}", 50);
        let chi = flory_huggins_chi(&pe, &pvc, 298.0).unwrap();
        assert!(chi >= 0.0, "chi must be non-negative, got {chi}");
    }

    #[test]
    fn flory_huggins_decreases_with_temperature() {
        // chi ~ 1/T, so higher T -> lower chi
        let pe = build_chain("{[]CC[]}", 50);
        let ps = build_chain("{[]CC(c1ccccc1)[]}", 50);
        let chi_low = flory_huggins_chi(&pe, &ps, 300.0).unwrap();
        let chi_high = flory_huggins_chi(&pe, &ps, 500.0).unwrap();
        assert!(
            chi_high < chi_low,
            "chi should decrease with T: chi(300K)={chi_low:.3}, chi(500K)={chi_high:.3}"
        );
    }
}
