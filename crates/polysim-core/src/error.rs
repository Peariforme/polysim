use thiserror::Error;

#[derive(Debug, Error)]
pub enum PolySimError {
    #[error("BigSMILES parse error: {0}")]
    Parse(#[from] bigsmiles::ParseError),

    #[error("Invalid build strategy: {0}")]
    BuildStrategy(String),

    #[error("No stochastic object (repeat units) found in BigSMILES")]
    NoStochasticObject,

    #[error("Incompatible repeat unit count for {architecture}: got {got}, need {need}")]
    RepeatUnitCount {
        architecture: &'static str,
        got: usize,
        need: usize,
    },

    #[error("Weight fractions must sum to 1.0 (got {sum:.4})")]
    InvalidFractions { sum: f64 },

    #[error(
        "ring number overflow: l'unité de répétition utilise {max_ring} ring(s), \
         maximum supporté par SMILES = {max_supported}"
    )]
    RingNumberOverflow { max_ring: u32, max_supported: u32 },
}
