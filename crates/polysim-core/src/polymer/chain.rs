/// A single, fully resolved polymer chain instance.
#[derive(Debug, Clone)]
pub struct PolymerChain {
    /// SMILES string representing this specific chain.
    pub smiles: String,
    /// Number of repeat units incorporated.
    pub repeat_count: usize,
    /// Computed number-average molecular weight (g/mol).
    pub mn: f64,
}

impl PolymerChain {
    pub fn new(smiles: String, repeat_count: usize, mn: f64) -> Self {
        Self { smiles, repeat_count, mn }
    }
}

impl std::fmt::Display for PolymerChain {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.smiles)
    }
}
