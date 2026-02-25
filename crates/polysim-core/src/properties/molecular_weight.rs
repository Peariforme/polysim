use crate::polymer::PolymerChain;

/// Computes the monoisotopic mass of a polymer chain from its SMILES (g/mol).
///
/// Uses the most abundant isotope for each element (e.g. ¹²C = 12.000,
/// ¹H = 1.00783, ¹⁶O = 15.9949, …).
pub fn monoisotopic_mass(_chain: &PolymerChain) -> f64 {
    todo!("parse SMILES atoms and sum monoisotopic masses")
}

/// Computes the average molecular mass of a polymer chain from its SMILES (g/mol).
///
/// Uses IUPAC standard atomic weights (isotopically averaged).
pub fn average_mass(_chain: &PolymerChain) -> f64 {
    todo!("parse SMILES atoms and sum average atomic weights")
}
