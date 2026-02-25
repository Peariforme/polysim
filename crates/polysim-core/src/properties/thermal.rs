use crate::polymer::PolymerChain;

/// Estimates the glass transition temperature (K) using the Fox equation.
///
/// # Arguments
///
/// `components` — slice of `(weight_fraction, Tg_homopolymer_K)` pairs,
/// one per distinct repeat unit. Weight fractions must sum to 1.0.
///
/// # Reference
///
/// Fox, T. G. (1956). *Bull. Am. Phys. Soc.* **1**, 123.
///
/// # Example
///
/// ```rust
/// use polysim_core::properties::thermal::tg_fox;
///
/// // 50/50 blend of PS (Tg ≈ 373 K) and PMMA (Tg ≈ 378 K)
/// let tg = tg_fox(&[(0.5, 373.0), (0.5, 378.0)]);
/// assert!((tg - 375.4).abs() < 0.2);
/// ```
pub fn tg_fox(components: &[(f64, f64)]) -> f64 {
    let inv_tg: f64 = components.iter().map(|(wi, tgi)| wi / tgi).sum();
    1.0 / inv_tg
}

/// Estimates Tg (K) using the Van Krevelen group-contribution method.
///
/// # Reference
///
/// Van Krevelen, D. W. & te Nijenhuis, K. (2009).
/// *Properties of Polymers*, 4th ed., Elsevier. Chapter 6.
pub fn tg_van_krevelen(_chain: &PolymerChain) -> f64 {
    todo!("Van Krevelen group-contribution Tg")
}

/// Qualitative tendency of a polymer chain to crystallise.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CrystallizationTendency {
    /// Highly regular chain — expected to crystallise readily (e.g. PE, isotactic PP).
    High,
    /// Moderate regularity — partial crystallisation possible (e.g. syndiotactic PS).
    Medium,
    /// Low regularity — unlikely to crystallise significantly.
    Low,
    /// Fully amorphous — no crystallisation expected (e.g. atactic PS, PMMA).
    Amorphous,
}

/// Estimates the crystallisation tendency of a polymer chain based on its
/// structural regularity and symmetry.
pub fn crystallization_tendency(_chain: &PolymerChain) -> CrystallizationTendency {
    todo!("estimate crystallisation tendency from SMILES regularity/symmetry")
}
