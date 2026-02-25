use crate::{error::PolySimError, polymer::PolymerChain};
use super::strategy::BuildStrategy;
use bigsmiles::BigSmiles;

/// Builder for linear polymer architectures derived from a single BigSMILES.
pub struct LinearBuilder {
    bigsmiles: BigSmiles,
    strategy: BuildStrategy,
}

impl LinearBuilder {
    pub fn new(bigsmiles: BigSmiles, strategy: BuildStrategy) -> Self {
        Self { bigsmiles, strategy }
    }

    /// Build a linear homopolymer (single repeat unit, repeated n times).
    pub fn homopolymer(&self) -> Result<PolymerChain, PolySimError> {
        todo!("implement homopolymer generation")
    }

    /// Build a random (statistical) copolymer.
    ///
    /// `fractions` — weight fraction of each repeat unit (must sum to 1.0).
    /// The BigSMILES must contain exactly `fractions.len()` repeat units.
    pub fn random_copolymer(&self, fractions: &[f64]) -> Result<PolymerChain, PolySimError> {
        let sum: f64 = fractions.iter().sum();
        if (sum - 1.0).abs() > 1e-6 {
            return Err(PolySimError::InvalidFractions { sum });
        }
        todo!("implement random copolymer generation")
    }

    /// Build an alternating copolymer (–A–B–A–B–).
    ///
    /// The BigSMILES must contain exactly 2 repeat units.
    pub fn alternating_copolymer(&self) -> Result<PolymerChain, PolySimError> {
        todo!("implement alternating copolymer generation")
    }

    /// Build a block copolymer (–AAAA–BBBB–).
    ///
    /// `block_lengths` — number of repeat units per block, in order.
    /// The BigSMILES must contain exactly `block_lengths.len()` repeat units.
    pub fn block_copolymer(&self, block_lengths: &[usize]) -> Result<PolymerChain, PolySimError> {
        todo!("implement block copolymer generation")
    }
}
