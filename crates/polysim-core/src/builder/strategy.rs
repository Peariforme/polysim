/// Determines how many repeat units are incorporated into a generated chain.
///
/// All mass-based variants use SI/chemistry conventions:
/// - molecular weights in **g/mol**
/// - monoisotopic masses in **g/mol**
#[derive(Debug, Clone)]
pub enum BuildStrategy {
    /// Generate exactly `n` repeat units.
    ByRepeatCount(usize),

    /// Target number-average molecular weight (Mn) in g/mol.
    ///
    /// The repeat count is chosen so that the chain Mn is as close as possible
    /// to the given target. Requires molecular weight calculation to be
    /// implemented (see `properties::molecular_weight`).
    ByTargetMn(f64),

    /// Target an exact (monoisotopic) chain mass in g/mol.
    ///
    /// The repeat count is chosen so that the monoisotopic mass is as close as
    /// possible to the given target. Requires molecular weight calculation to be
    /// implemented (see `properties::molecular_weight`).
    ByExactMass(f64),
}
