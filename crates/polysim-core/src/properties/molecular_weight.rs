use crate::polymer::PolymerChain;

/// Compute the monoisotopic mass of a polymer chain from its SMILES (g/mol).
pub fn monoisotopic_mass(_chain: &PolymerChain) -> f64 {
    todo!("parse SMILES atoms and sum monoisotopic masses")
}

/// Compute the average molecular mass of a polymer chain from its SMILES (g/mol).
pub fn average_mass(_chain: &PolymerChain) -> f64 {
    todo!("parse SMILES atoms and sum average atomic weights")
}
