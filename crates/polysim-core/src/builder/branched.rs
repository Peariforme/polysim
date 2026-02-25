use crate::{error::PolySimError, polymer::PolymerChain};
use super::strategy::BuildStrategy;
use bigsmiles_core::BigSmiles;

/// Builder for non-linear polymer architectures (branched, graft, macromonomer).
pub struct BranchedBuilder {
    /// BigSMILES of the backbone.
    backbone: BigSmiles,
    /// BigSMILES of the branch / side chain.
    branch: BigSmiles,
    strategy: BuildStrategy,
}

impl BranchedBuilder {
    pub fn new(backbone: BigSmiles, branch: BigSmiles, strategy: BuildStrategy) -> Self {
        Self { backbone, branch, strategy }
    }

    /// Build a comb (regularly branched) polymer.
    ///
    /// `branch_every` — attach one branch every N backbone repeat units.
    pub fn comb_polymer(&self, branch_every: usize) -> Result<PolymerChain, PolySimError> {
        todo!("implement comb/branched polymer generation")
    }

    /// Build a graft copolymer (random branch-point placement).
    ///
    /// `graft_fraction` — fraction of backbone repeat units that carry a branch.
    pub fn graft_copolymer(&self, graft_fraction: f64) -> Result<PolymerChain, PolySimError> {
        todo!("implement graft copolymer generation")
    }

    /// Build a macromonomer (a single branch/side chain with a polymerisable end group).
    pub fn macromonomer(&self) -> Result<PolymerChain, PolySimError> {
        todo!("implement macromonomer generation")
    }
}
