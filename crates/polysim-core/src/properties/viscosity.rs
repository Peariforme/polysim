//! Intrinsic viscosity prediction using the Mark-Houwink equation.
//!
//! All intrinsic viscosities are in **dL/g**.

use crate::error::PolySimError;
use crate::polymer::PolymerChain;
use crate::properties::molecular_weight::average_mass;

// ---------------------------------------------------------------------------
// Mark-Houwink parameters
// ---------------------------------------------------------------------------

/// Mark-Houwink constants for a given polymer/solvent/temperature system.
///
/// The constants relate the intrinsic viscosity to the number-average molar
/// mass via the Mark-Houwink equation:
///
/// **[η] = K · Mn^a**
///
/// # Units
///
/// `k` is given in **mL/g** (cm³/g), consistent with the *Polymer Handbook*
/// (Brandrup & Immergut, 4th ed., Section VII). The function
/// [`intrinsic_viscosity`] converts the result to **dL/g** (1 dL = 100 mL).
#[derive(Debug, Clone, PartialEq)]
pub struct MarkHouwinkParams {
    /// Pre-exponential constant K (mL/g).
    pub k: f64,
    /// Mark-Houwink exponent a (dimensionless, typical range 0.5–0.8).
    pub a: f64,
    /// Solvent name (informational).
    pub solvent: String,
    /// Temperature (K) at which the constants were determined.
    pub temperature: f64,
}

impl MarkHouwinkParams {
    /// Mark-Houwink parameters for polystyrene in toluene at 25 °C.
    ///
    /// Source: Brandrup, J. & Immergut, E. H. (1999). *Polymer Handbook*, 4th ed.,
    /// Wiley, Section VII, p. VII-1. Also confirmed by Berry (1967),
    /// *J. Chem. Phys.* **46**, 1338.
    pub fn ps_toluene_25c() -> Self {
        Self {
            k: 1.16e-2,
            a: 0.73,
            solvent: "toluene".into(),
            temperature: 298.15,
        }
    }

    /// Mark-Houwink parameters for poly(methyl methacrylate) in acetone at 25 °C.
    ///
    /// Source: Brandrup, J. & Immergut, E. H. (1999). *Polymer Handbook*, 4th ed.,
    /// Wiley, Section VII, p. VII-32. Constants apply to atactic PMMA.
    pub fn pmma_acetone_25c() -> Self {
        Self {
            k: 7.5e-3,
            a: 0.70,
            solvent: "acetone".into(),
            temperature: 298.15,
        }
    }

    /// Mark-Houwink parameters for polyethylene in decalin at 135 °C.
    ///
    /// Source: Brandrup, J. & Immergut, E. H. (1999). *Polymer Handbook*, 4th ed.,
    /// Wiley, Section VII, p. VII-14. Measured at 135 °C to ensure full dissolution
    /// of semicrystalline PE; decalin (decahydronaphthalene) is a good solvent for
    /// PE at elevated temperature.
    pub fn pe_decalin_135c() -> Self {
        Self {
            k: 6.2e-2,
            a: 0.70,
            solvent: "decalin".into(),
            temperature: 408.15,
        }
    }
}

// ---------------------------------------------------------------------------
// Core calculation
// ---------------------------------------------------------------------------

/// Estimates the intrinsic viscosity [η] (dL/g) using the Mark-Houwink equation.
///
/// The intrinsic viscosity is computed as:
///
/// **[η] = K · Mn^a / 100**
///
/// where `K` is in mL/g (Brandrup convention), `Mn` is the number-average
/// molar mass in g/mol computed from [`average_mass`], `a` is the dimensionless
/// Mark-Houwink exponent, and the factor 1/100 converts from mL/g to dL/g
/// (1 dL = 100 mL = 100 cm³).
///
/// # Errors
///
/// Returns [`PolySimError::GroupDecomposition`] if:
/// - `params.k` is non-positive.
/// - `params.a` is outside the range [0.0, 2.0].
///
/// # Reference
///
/// Mark, H. (1938). In *Der feste Körper*, Hirzel, Leipzig, p. 103.
/// Houwink, R. (1940). *J. Prakt. Chem.* **157**, 15.
/// Brandrup, J. & Immergut, E. H. (1999). *Polymer Handbook*, 4th ed., Wiley,
/// Section VII.
pub fn intrinsic_viscosity(
    chain: &PolymerChain,
    params: &MarkHouwinkParams,
) -> Result<f64, PolySimError> {
    if params.k <= 0.0 {
        return Err(PolySimError::GroupDecomposition(
            "Mark-Houwink K must be positive".into(),
        ));
    }
    if params.a < 0.0 || params.a > 2.0 {
        return Err(PolySimError::GroupDecomposition(format!(
            "Mark-Houwink exponent a = {} is outside the range [0, 2]",
            params.a
        )));
    }

    let mn = average_mass(chain);

    // [η] in mL/g = K (mL/g) · Mn^a
    // Convert to dL/g: divide by 100 (1 dL = 100 mL)
    let eta_ml_g = params.k * mn.powf(params.a);
    Ok(eta_ml_g / 100.0)
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

    // PS repeat unit: ~104.15 g/mol. n=1000 → Mn ≈ 104 150 g/mol
    // [η] = 1.16e-2 × 104150^0.73 / 100
    // 104150^0.73 = exp(0.73 × ln(104150)) ≈ exp(0.73 × 11.553) ≈ exp(8.434) ≈ 4585
    // [η] ≈ 1.16e-2 × 4585 / 100 ≈ 0.532 dL/g
    #[test]
    fn ps_toluene_physical_range() {
        let chain = build_chain("{[]CC(c1ccccc1)[]}", 1000);
        let eta = intrinsic_viscosity(&chain, &MarkHouwinkParams::ps_toluene_25c()).unwrap();
        assert!(
            eta > 0.3 && eta < 2.0,
            "PS [η] = {eta:.3} dL/g, expected 0.3–2.0"
        );
    }

    // PE repeat unit: ~28.05 g/mol. n=1000 → Mn ≈ 28 050 g/mol
    // [η] = 6.2e-2 × 28050^0.70 / 100
    // 28050^0.70 = exp(0.70 × ln(28050)) ≈ exp(0.70 × 10.242) ≈ exp(7.169) ≈ 1293
    // [η] ≈ 6.2e-2 × 1293 / 100 ≈ 0.802 dL/g
    #[test]
    fn pe_decalin_physical_range() {
        let chain = build_chain("{[]CC[]}", 1000);
        let eta = intrinsic_viscosity(&chain, &MarkHouwinkParams::pe_decalin_135c()).unwrap();
        assert!(
            eta > 0.3 && eta < 3.0,
            "PE [η] = {eta:.3} dL/g, expected 0.3–3.0"
        );
    }

    // PMMA repeat unit: ~100.12 g/mol. n=1000 → Mn ≈ 100 120 g/mol
    // [η] = 7.5e-3 × 100120^0.70 / 100
    // 100120^0.70 = exp(0.70 × ln(100120)) ≈ exp(0.70 × 11.514) ≈ exp(8.060) ≈ 3163
    // [η] ≈ 7.5e-3 × 3163 / 100 ≈ 0.237 dL/g
    #[test]
    fn pmma_acetone_physical_range() {
        let chain = build_chain("{[]CC(C)(C(=O)OC)[]}", 1000);
        let eta = intrinsic_viscosity(&chain, &MarkHouwinkParams::pmma_acetone_25c()).unwrap();
        assert!(
            eta > 0.1 && eta < 1.5,
            "PMMA [η] = {eta:.3} dL/g, expected 0.1–1.5"
        );
    }

    #[test]
    fn positive_result() {
        let chain = build_chain("{[]CC[]}", 100);
        let eta = intrinsic_viscosity(&chain, &MarkHouwinkParams::pe_decalin_135c()).unwrap();
        assert!(eta > 0.0, "[η] must be positive, got {eta}");
    }

    #[test]
    fn increases_with_chain_length() {
        // [η] ∝ Mn^a with a > 0, so longer chains → higher [η]
        let chain_short = build_chain("{[]CC(c1ccccc1)[]}", 100);
        let chain_long = build_chain("{[]CC(c1ccccc1)[]}", 1000);
        let eta_short =
            intrinsic_viscosity(&chain_short, &MarkHouwinkParams::ps_toluene_25c()).unwrap();
        let eta_long =
            intrinsic_viscosity(&chain_long, &MarkHouwinkParams::ps_toluene_25c()).unwrap();
        assert!(
            eta_long > eta_short,
            "[η] must increase with chain length: short={eta_short:.3}, long={eta_long:.3}"
        );
    }

    #[test]
    fn invalid_k_returns_error() {
        let chain = build_chain("{[]CC[]}", 50);
        let bad_params = MarkHouwinkParams {
            k: -1.0,
            a: 0.70,
            solvent: String::new(),
            temperature: 298.15,
        };
        assert!(intrinsic_viscosity(&chain, &bad_params).is_err());
    }

    #[test]
    fn invalid_a_returns_error() {
        let chain = build_chain("{[]CC[]}", 50);
        let bad_params = MarkHouwinkParams {
            k: 1.0e-2,
            a: 3.0,
            solvent: String::new(),
            temperature: 298.15,
        };
        assert!(intrinsic_viscosity(&chain, &bad_params).is_err());
    }
}
