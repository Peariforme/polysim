use crate::error::PolySimError;

use super::PolymerChain;

/// A collection of polymer chains representing a polydisperse sample.
#[derive(Debug, Clone)]
pub struct PolymerEnsemble {
    chains: Vec<PolymerChain>,
}

impl PolymerEnsemble {
    /// Creates a new ensemble from a vector of chains.
    ///
    /// # Errors
    ///
    /// Returns [`PolySimError::EmptyEnsemble`] if `chains` is empty.
    pub fn new(chains: Vec<PolymerChain>) -> Result<Self, PolySimError> {
        if chains.is_empty() {
            return Err(PolySimError::EmptyEnsemble);
        }
        Ok(Self { chains })
    }

    /// Returns a reference to the individual chains.
    pub fn chains(&self) -> &[PolymerChain] {
        &self.chains
    }

    /// Number of chains in the ensemble (always ≥ 1).
    pub fn len(&self) -> usize {
        self.chains.len()
    }

    /// Always returns `false` — an ensemble is guaranteed non-empty by construction.
    pub fn is_empty(&self) -> bool {
        false
    }

    /// Number-average molecular weight: Mn = Σ Mi / N
    pub fn mn(&self) -> f64 {
        let sum: f64 = self.chains.iter().map(|c| c.mn).sum();
        sum / self.chains.len() as f64
    }

    /// Weight-average molecular weight: Mw = Σ Mi² / Σ Mi
    pub fn mw(&self) -> f64 {
        let sum_mi: f64 = self.chains.iter().map(|c| c.mn).sum();
        let sum_mi2: f64 = self.chains.iter().map(|c| c.mn * c.mn).sum();
        sum_mi2 / sum_mi
    }

    /// Polydispersity index: PDI = Mw / Mn
    pub fn pdi(&self) -> f64 {
        self.mw() / self.mn()
    }
}
