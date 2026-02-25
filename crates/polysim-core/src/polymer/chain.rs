/// A single, fully resolved polymer chain instance.
///
/// A `PolymerChain` is the output of a builder: it holds the concrete SMILES
/// string for the generated chain together with metadata computed at build time.
#[derive(Debug, Clone)]
pub struct PolymerChain {
    /// SMILES string representing this specific chain.
    pub smiles: String,
    /// Number of repeat units incorporated into the chain.
    pub repeat_count: usize,
    /// Number-average molecular weight in g/mol.
    ///
    /// Currently `0.0` until `properties::molecular_weight` is implemented.
    pub mn: f64,
}

impl PolymerChain {
    /// Creates a new `PolymerChain` with the given SMILES, repeat count, and Mn.
    pub fn new(smiles: String, repeat_count: usize, mn: f64) -> Self {
        Self {
            smiles,
            repeat_count,
            mn,
        }
    }
}

impl std::fmt::Display for PolymerChain {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.smiles)
    }
}
