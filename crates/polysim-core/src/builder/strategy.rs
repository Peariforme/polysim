/// Determines the length of a generated polymer chain.
#[derive(Debug, Clone)]
pub enum BuildStrategy {
    /// Generate exactly `n` repeat units.
    ByRepeatCount(usize),

    /// Target number-average molecular weight in g/mol.
    /// Chain length is chosen to approximate the target as closely as possible.
    ByTargetMn(f64),

    /// Generate a chain whose molecular weight is as close as possible to the
    /// given exact mass in g/mol (monoisotopic).
    ByExactMass(f64),
}
