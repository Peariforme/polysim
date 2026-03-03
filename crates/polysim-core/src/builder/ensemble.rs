use bigsmiles::BigSmiles;
use rand::{rngs::StdRng, SeedableRng};

use crate::{
    distribution::ChainLengthDistribution,
    error::PolySimError,
    polymer::{PolymerChain, PolymerEnsemble},
    properties::molecular_weight::average_mass,
};

use super::linear::build_linear_smiles;

/// Default number of chains in an ensemble.
pub const DEFAULT_NUM_CHAINS: usize = 100;

/// Builder that generates a polydisperse ensemble of polymer chains.
pub struct EnsembleBuilder<D: ChainLengthDistribution> {
    bigsmiles: BigSmiles,
    distribution: D,
    mn: f64,
    pdi: f64,
    num_chains: usize,
    seed: Option<u64>,
}

impl<D: ChainLengthDistribution> EnsembleBuilder<D> {
    /// Creates a new ensemble builder.
    ///
    /// - `bigsmiles` — parsed BigSMILES describing the polymer.
    /// - `distribution` — chain length distribution model.
    /// - `mn` — target number-average molecular weight (g/mol).
    /// - `pdi` — target polydispersity index (Mw/Mn).
    pub fn new(bigsmiles: BigSmiles, distribution: D, mn: f64, pdi: f64) -> Self {
        Self {
            bigsmiles,
            distribution,
            mn,
            pdi,
            num_chains: DEFAULT_NUM_CHAINS,
            seed: None,
        }
    }

    /// Override the number of chains to generate (default: 100).
    pub fn num_chains(mut self, n: usize) -> Self {
        self.num_chains = n;
        self
    }

    /// Set a seed for reproducible sampling.
    pub fn seed(mut self, seed: u64) -> Self {
        self.seed = Some(seed);
        self
    }

    /// Build a polydisperse ensemble of homopolymer chains.
    ///
    /// # Errors
    ///
    /// - [`PolySimError::NoStochasticObject`] if no `{...}` block found.
    /// - [`PolySimError::RepeatUnitCount`] if ≠ 1 repeat unit.
    /// - [`PolySimError::EmptyEnsemble`] if `num_chains` is 0.
    pub fn homopolymer_ensemble(&self) -> Result<PolymerEnsemble, PolySimError> {
        let stoch = self
            .bigsmiles
            .first_stochastic()
            .ok_or(PolySimError::NoStochasticObject)?;

        if stoch.repeat_units.len() != 1 {
            return Err(PolySimError::RepeatUnitCount {
                architecture: "homopolymer",
                got: stoch.repeat_units.len(),
                need: 1,
            });
        }

        let smiles_raw = &stoch.repeat_units[0].smiles_raw;

        // Two-point calibration to separate repeat-unit mass (m0) from
        // end-group mass (m_end): MW(n) = n × m0 + m_end.
        let mw1 = average_mass(&PolymerChain::new(
            build_linear_smiles(smiles_raw, 1)?,
            1,
            0.0,
        ));
        let mw2 = average_mass(&PolymerChain::new(
            build_linear_smiles(smiles_raw, 2)?,
            2,
            0.0,
        ));
        let m0 = mw2 - mw1;
        let m_end = mw1 - m0;

        // Adjust target Mn to account for end groups: Mn_target = Xn × m0 + m_end
        let target_mn_corrected = self.mn - m_end;

        // Sample chain lengths from distribution.
        let mut rng: Box<dyn rand::RngCore> = match self.seed {
            Some(s) => Box::new(StdRng::seed_from_u64(s)),
            None => Box::new(rand::rng()),
        };
        let lengths = self.distribution.sample(
            target_mn_corrected,
            self.pdi,
            m0,
            self.num_chains,
            &mut *rng,
        );

        // Build each chain.
        let chains: Result<Vec<PolymerChain>, PolySimError> = lengths
            .into_iter()
            .map(|n| {
                let smiles = build_linear_smiles(smiles_raw, n)?;
                let chain = PolymerChain::new(smiles, n, 0.0);
                let mn = average_mass(&chain);
                Ok(PolymerChain::new(chain.smiles, n, mn))
            })
            .collect();

        PolymerEnsemble::new(chains?)
    }
}
