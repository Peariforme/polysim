//! Chain length distribution models for polydisperse polymer ensembles.
//!
//! Each distribution implements [`ChainLengthDistribution`] and can sample
//! repeat-unit counts given target Mn, PDI, and repeat-unit molar mass.

pub mod flory;
pub mod log_normal;
pub mod schulz_zimm;

pub use flory::Flory;
pub use log_normal::LogNormal;
pub use schulz_zimm::SchulzZimm;

use rand::RngCore;

/// A distribution that can sample chain lengths (as repeat-unit counts).
pub trait ChainLengthDistribution {
    /// Sample `num_chains` chain lengths (each ≥ 1).
    ///
    /// - `mn` — target number-average molecular weight (g/mol).
    /// - `pdi` — target polydispersity index (Mw/Mn, ≥ 1.0).
    /// - `m0` — molar mass of a single repeat unit (g/mol).
    /// - `rng` — random number generator (pass a seeded RNG for reproducibility).
    fn sample(
        &self,
        mn: f64,
        pdi: f64,
        m0: f64,
        num_chains: usize,
        rng: &mut dyn RngCore,
    ) -> Vec<usize>;

    /// Human-readable name for display purposes.
    fn name(&self) -> &'static str;
}
