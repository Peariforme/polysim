use thiserror::Error;

/// All errors that can be produced by polysim operations.
#[derive(Debug, Error)]
pub enum PolySimError {
    /// A BigSMILES string could not be parsed.
    #[error("BigSMILES parse error: {0}")]
    Parse(#[from] bigsmiles::ParseError),

    /// The [`BuildStrategy`](crate::BuildStrategy) is invalid or not yet supported.
    #[error("Invalid build strategy: {0}")]
    BuildStrategy(String),

    /// The BigSMILES contains no stochastic object (`{...}`), so no repeat units
    /// are available for chain generation.
    #[error("No stochastic object (repeat units) found in BigSMILES")]
    NoStochasticObject,

    /// The stochastic object contains the wrong number of repeat units for the
    /// requested architecture.
    #[error("Incompatible repeat unit count for {architecture}: got {got}, need {need}")]
    RepeatUnitCount {
        architecture: &'static str,
        got: usize,
        need: usize,
    },

    /// The weight fractions supplied to a copolymer builder do not sum to 1.0.
    #[error("Weight fractions must sum to 1.0 (got {sum:.4})")]
    InvalidFractions { sum: f64 },

    /// A single repeat unit already uses more than 99 distinct ring-closure numbers,
    /// which exceeds the SMILES specification.
    #[error(
        "Ring number overflow: the repeat unit uses {max_ring} ring closure(s), \
         SMILES maximum is {max_supported}"
    )]
    RingNumberOverflow { max_ring: u32, max_supported: u32 },
}
