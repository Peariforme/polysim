use bigsmiles::BigSmiles;

use crate::{error::PolySimError, polymer::PolymerChain};

use super::strategy::BuildStrategy;

/// Builder for non-linear polymer architectures (branched, graft, macromonomer).
///
/// Unlike [`LinearBuilder`](super::linear::LinearBuilder), this builder takes
/// two BigSMILES strings: one for the **backbone** and one for the **branch**.
// Fields are stored for future use once the builder methods are implemented.
#[allow(dead_code)]
pub struct BranchedBuilder {
    /// BigSMILES of the backbone chain.
    backbone: BigSmiles,
    /// BigSMILES of the branch / side chain.
    branch: BigSmiles,
    /// Strategy that controls backbone chain length.
    strategy: BuildStrategy,
}

impl BranchedBuilder {
    /// Creates a new builder from backbone and branch BigSMILES strings plus a
    /// build strategy that governs the backbone length.
    pub fn new(backbone: BigSmiles, branch: BigSmiles, strategy: BuildStrategy) -> Self {
        Self {
            backbone,
            branch,
            strategy,
        }
    }

    /// Generates a comb (regularly branched) polymer.
    ///
    /// `branch_every` — attach one branch every N backbone repeat units.
    pub fn comb_polymer(&self, _branch_every: usize) -> Result<PolymerChain, PolySimError> {
        todo!("implement comb/branched polymer generation")
    }

    /// Generates a graft copolymer (random branch-point placement).
    ///
    /// `graft_fraction` — fraction of backbone repeat units that carry a branch
    /// (0.0 = no grafting, 1.0 = every backbone unit is grafted).
    pub fn graft_copolymer(&self, _graft_fraction: f64) -> Result<PolymerChain, PolySimError> {
        todo!("implement graft copolymer generation")
    }

    /// Generates a macromonomer: a single branch/side chain with a
    /// polymerisable end group.
    pub fn macromonomer(&self) -> Result<PolymerChain, PolySimError> {
        todo!("implement macromonomer generation")
    }
}
