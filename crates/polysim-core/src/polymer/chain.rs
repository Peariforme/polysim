/// Composition unit for copolymer chains.
///
/// Stores a single repeat unit type with its molar fraction in the chain.
#[derive(Debug, Clone, PartialEq)]
pub struct MonomerUnit {
    /// SMILES string of the repeat unit (e.g. "CC" for ethylene).
    pub smiles: String,
    /// Molar fraction of this unit in the chain (0.0–1.0).
    pub fraction: f64,
}

impl MonomerUnit {
    /// Creates a new `MonomerUnit`.
    pub fn new(smiles: impl Into<String>, fraction: f64) -> Self {
        Self {
            smiles: smiles.into(),
            fraction,
        }
    }
}

/// Polymer chain architecture classification.
#[derive(Debug, Clone, PartialEq, Default)]
pub enum Architecture {
    /// Simple linear chain (default).
    #[default]
    Linear,
    /// Star polymer with `arms` number of arms radiating from a central core.
    Star { arms: usize },
    /// Comb polymer with branches every `branch_spacing` backbone units.
    Comb { branch_spacing: usize },
    /// Dendrimer of the given `generation`.
    Dendrimer { generation: usize },
    /// Cyclic polymer (no chain ends).
    Cyclic,
    /// Gradient copolymer with composition varying along the chain.
    Gradient,
    /// Graft copolymer with randomly placed branches at `graft_fraction`.
    Graft { graft_fraction: f64 },
}

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
    pub mn: f64,
    /// Monomer composition: each unit type with its molar fraction.
    ///
    /// Homopolymers have a single element with fraction 1.0.
    /// Empty when composition was not tracked by the builder.
    pub composition: Vec<MonomerUnit>,
    /// Polymer architecture (linear by default).
    pub architecture: Architecture,
}

impl PolymerChain {
    /// Creates a new `PolymerChain` with the given SMILES, repeat count, and Mn.
    ///
    /// `composition` defaults to empty and `architecture` to `Linear`.
    /// Use the builder methods [`Self::with_composition`] and
    /// [`Self::with_architecture`] to populate these fields.
    pub fn new(smiles: String, repeat_count: usize, mn: f64) -> Self {
        Self {
            smiles,
            repeat_count,
            mn,
            composition: Vec::new(),
            architecture: Architecture::default(),
        }
    }

    /// Attaches monomer composition metadata to this chain.
    pub fn with_composition(mut self, composition: Vec<MonomerUnit>) -> Self {
        self.composition = composition;
        self
    }

    /// Attaches architecture metadata to this chain.
    pub fn with_architecture(mut self, architecture: Architecture) -> Self {
        self.architecture = architecture;
        self
    }
}

impl std::fmt::Display for PolymerChain {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.smiles)
    }
}
