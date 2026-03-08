//! Gas permeability prediction via the Permachor method (Van Krevelen, Ch. 22).
//!
//! The Permachor method estimates gas permeability by summing structural
//! increments τᵢ for each functional group:
//!
//! **log₁₀(P) = Σ τᵢ · countᵢ**
//!
//! where P is in Barrer (1 Barrer = 10⁻¹⁰ cm³(STP)·cm / cm²·s·cmHg).
//!
//! # Reference
//!
//! Van Krevelen, D. W. & te Nijenhuis, K. (2009).
//! *Properties of Polymers*, 4th ed., Elsevier. Chapter 22.

use crate::error::PolySimError;
use crate::polymer::PolymerChain;
use crate::properties::group_contribution::GroupDatabase;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Gas species for permeability prediction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Gas {
    /// Oxygen (O₂).
    O2,
    /// Carbon dioxide (CO₂).
    CO2,
    /// Nitrogen (N₂).
    N2,
    /// Water vapour (H₂O).
    H2O,
    /// Helium (He).
    He,
    /// Hydrogen (H₂).
    H2,
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Permeability in Barrer (1 Barrer = 10⁻¹⁰ cm³(STP)·cm / cm²·s·cmHg).
///
/// Uses the Permachor method: log₁₀(P) = Σ τᵢ · countᵢ
///
/// # Errors
///
/// Returns `PolySimError::GroupDecomposition` if the SMILES cannot be parsed.
pub fn gas_permeability(chain: &PolymerChain, gas: Gas) -> Result<f64, PolySimError> {
    let groups = GroupDatabase::decompose(chain)?;
    let log_p: f64 = groups
        .iter()
        .map(|gm| permachor_contribution(gm.group.name, gas) * gm.count as f64)
        .sum();
    Ok(10_f64.powf(log_p))
}

// ---------------------------------------------------------------------------
// Internal: Permachor contributions
// ---------------------------------------------------------------------------

/// Returns the Permachor increment τ for a given group and gas.
///
/// Values from Van Krevelen Table 22.2 (per structural group).
fn permachor_contribution(name: &str, gas: Gas) -> f64 {
    match (name, gas) {
        ("-CH2-", Gas::O2) => 0.23,
        ("-CH2-", Gas::CO2) => 0.47,
        ("-CH2-", Gas::N2) => 0.17,
        ("-CH2-", Gas::H2O) => -0.20,
        ("-CH2-", Gas::He) => 0.40,
        ("-CH2-", Gas::H2) => 0.35,
        ("-CH3", Gas::O2) => 0.31,
        ("-CH3", Gas::CO2) => 0.62,
        ("-CH3", Gas::N2) => 0.23,
        ("-CH3", Gas::H2O) => -0.10,
        ("-CH3", Gas::He) => 0.50,
        ("-CH3", Gas::H2) => 0.45,
        ("-COO-", Gas::O2) => -0.65,
        ("-COO-", Gas::CO2) => -0.50,
        ("-COO-", Gas::N2) => -0.70,
        ("-COO-", Gas::H2O) => 0.80,
        ("-COO-", Gas::He) => -0.30,
        ("-COO-", Gas::H2) => -0.40,
        ("-O-", Gas::O2) => -0.30,
        ("-O-", Gas::CO2) => -0.20,
        ("-O-", Gas::N2) => -0.35,
        ("-O-", Gas::H2O) => 0.50,
        ("-O-", Gas::He) => -0.15,
        ("-O-", Gas::H2) => -0.20,
        ("-Cl", Gas::O2) => -0.80,
        ("-Cl", Gas::CO2) => -0.60,
        ("-Cl", Gas::N2) => -0.85,
        ("-Cl", Gas::H2O) => -0.10,
        ("-Cl", Gas::He) => -0.40,
        ("-Cl", Gas::H2) => -0.50,
        ("-C6H5", Gas::O2) => 0.05,
        ("-C6H5", Gas::CO2) => 0.10,
        ("-C6H5", Gas::N2) => 0.03,
        ("-C6H5", Gas::H2O) => -0.30,
        ("-C6H5", Gas::He) => 0.15,
        ("-C6H5", Gas::H2) => 0.12,
        ("-C6H4-", Gas::O2) => 0.04,
        ("-C6H4-", Gas::CO2) => 0.08,
        ("-C6H4-", Gas::N2) => 0.02,
        ("-C6H4-", Gas::H2O) => -0.25,
        ("-C6H4-", Gas::He) => 0.12,
        ("-C6H4-", Gas::H2) => 0.10,
        // Default neutral contribution for unrecognised groups.
        (_, Gas::O2) => 0.10,
        (_, Gas::CO2) => 0.20,
        (_, Gas::N2) => 0.07,
        (_, Gas::H2O) => 0.00,
        (_, Gas::He) => 0.20,
        (_, Gas::H2) => 0.18,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

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
    fn pe_o2_permeability_positive() {
        let chain = build_chain("{[]CC[]}", 20);
        let p = gas_permeability(&chain, Gas::O2).unwrap();
        assert!(p > 0.0, "O2 permeability should be positive, got {p}");
    }

    #[test]
    fn pvc_lower_o2_than_pe() {
        // PVC has Cl groups which strongly reduce O2 permeability.
        let pe = build_chain("{[]CC[]}", 20);
        let pvc = build_chain("{[]C(Cl)C[]}", 20);
        let p_pe = gas_permeability(&pe, Gas::O2).unwrap();
        let p_pvc = gas_permeability(&pvc, Gas::O2).unwrap();
        assert!(
            p_pvc < p_pe,
            "PVC O2 permeability ({p_pvc:.4}) should be < PE ({p_pe:.4})"
        );
    }

    #[test]
    fn all_gases_return_ok() {
        let chain = build_chain("{[]CC[]}", 10);
        for gas in [Gas::O2, Gas::CO2, Gas::N2, Gas::H2O, Gas::He, Gas::H2] {
            assert!(gas_permeability(&chain, gas).is_ok());
        }
    }
}
