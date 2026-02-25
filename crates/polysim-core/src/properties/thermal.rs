use crate::polymer::PolymerChain;

/// Estimate the glass transition temperature (K) using the Fox equation.
///
/// `components` — `(weight_fraction, Tg_homopolymer_K)` for each repeat unit.
///
/// # Reference
/// Fox, T. G. (1956). *Bull. Am. Phys. Soc.* **1**, 123.
pub fn tg_fox(components: &[(f64, f64)]) -> f64 {
    let inv_tg: f64 = components.iter().map(|(wi, tgi)| wi / tgi).sum();
    1.0 / inv_tg
}

/// Estimate Tg (K) using Van Krevelen group-contribution method.
///
/// # Reference
/// Van Krevelen, D. W. (1990). *Properties of Polymers*, 3rd ed., Elsevier.
pub fn tg_van_krevelen(_chain: &PolymerChain) -> f64 {
    todo!("Van Krevelen group-contribution Tg")
}

/// Tendency of a polymer chain to crystallise, based on structural regularity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CrystallizationTendency {
    High,
    Medium,
    Low,
    Amorphous,
}

/// Estimate the crystallisation tendency of a polymer chain.
pub fn crystallization_tendency(_chain: &PolymerChain) -> CrystallizationTendency {
    todo!("estimate crystallisation tendency from SMILES regularity/symmetry")
}
