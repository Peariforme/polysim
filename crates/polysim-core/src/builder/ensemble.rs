use bigsmiles::BigSmiles;
use rand::distr::weighted::WeightedIndex;
use rand::prelude::*;
use rand::rngs::StdRng;

use crate::{
    distribution::ChainLengthDistribution,
    error::PolySimError,
    polymer::{PolymerChain, PolymerEnsemble},
    properties::molecular_weight::average_mass,
};

use super::linear::{
    build_copolymer_smiles, build_linear_smiles, gradient_fraction, GradientProfile,
};

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
                need_min: 1,
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
        let mut rng = self.make_rng();
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

    /// Build a polydisperse ensemble of random copolymer chains.
    ///
    /// `fractions` — weight fraction of each repeat unit (must sum to 1.0).
    pub fn random_copolymer_ensemble(
        &self,
        fractions: &[f64],
    ) -> Result<PolymerEnsemble, PolySimError> {
        let sum: f64 = fractions.iter().sum();
        if (sum - 1.0).abs() > 1e-6 {
            return Err(PolySimError::InvalidFractions { sum });
        }

        let (units, m0_avg, m_end) = self.copolymer_calibration(fractions)?;

        let target_mn_corrected = self.mn - m_end;
        let mut rng = self.make_rng();
        let lengths = self.distribution.sample(
            target_mn_corrected,
            self.pdi,
            m0_avg,
            self.num_chains,
            &mut *rng,
        );

        let dist = WeightedIndex::new(fractions)
            .map_err(|e| PolySimError::BuildStrategy(format!("invalid weight fractions: {e}")))?;

        let chains: Result<Vec<PolymerChain>, PolySimError> = lengths
            .into_iter()
            .map(|n| {
                let sequence: Vec<&str> = (0..n).map(|_| units[dist.sample(&mut *rng)]).collect();
                let smiles = build_copolymer_smiles(&sequence)?;
                let chain = PolymerChain::new(smiles, n, 0.0);
                let mn = average_mass(&chain);
                Ok(PolymerChain::new(chain.smiles, n, mn))
            })
            .collect();

        PolymerEnsemble::new(chains?)
    }

    /// Build a polydisperse ensemble of alternating copolymer chains.
    pub fn alternating_copolymer_ensemble(&self) -> Result<PolymerEnsemble, PolySimError> {
        let stoch = self
            .bigsmiles
            .first_stochastic()
            .ok_or(PolySimError::NoStochasticObject)?;

        if stoch.repeat_units.len() < 2 {
            return Err(PolySimError::RepeatUnitCount {
                architecture: "alternating copolymer",
                got: stoch.repeat_units.len(),
                need_min: 2,
            });
        }

        let units: Vec<&str> = stoch
            .repeat_units
            .iter()
            .map(|f| f.smiles_raw.as_str())
            .collect();
        let k = units.len();

        // Calibrate using one full cycle as the "composite repeat unit".
        let cycle_smiles: Vec<&str> = units.to_vec();
        let cycle1 = build_copolymer_smiles(&cycle_smiles)?;
        let cycle2_seq: Vec<&str> = units.iter().chain(units.iter()).copied().collect();
        let cycle2 = build_copolymer_smiles(&cycle2_seq)?;

        let mw1 = average_mass(&PolymerChain::new(cycle1, k, 0.0));
        let mw2 = average_mass(&PolymerChain::new(cycle2, k * 2, 0.0));
        let m0_cycle = mw2 - mw1; // mass per full cycle
        let m_end = mw1 - m0_cycle;

        let target_mn_corrected = self.mn - m_end;
        // m0 for distribution = mass per individual unit (average)
        let m0_per_unit = m0_cycle / k as f64;

        let mut rng = self.make_rng();
        let lengths = self.distribution.sample(
            target_mn_corrected,
            self.pdi,
            m0_per_unit,
            self.num_chains,
            &mut *rng,
        );

        let chains: Result<Vec<PolymerChain>, PolySimError> = lengths
            .into_iter()
            .map(|n| {
                let sequence: Vec<&str> = (0..n).map(|i| units[i % k]).collect();
                let smiles = build_copolymer_smiles(&sequence)?;
                let chain = PolymerChain::new(smiles, n, 0.0);
                let mn = average_mass(&chain);
                Ok(PolymerChain::new(chain.smiles, n, mn))
            })
            .collect();

        PolymerEnsemble::new(chains?)
    }

    /// Build a polydisperse ensemble of block copolymer chains.
    ///
    /// `block_ratios` — relative proportion of each block (must sum to 1.0).
    /// Each chain has a different total n, distributed according to the ratios.
    pub fn block_copolymer_ensemble(
        &self,
        block_ratios: &[f64],
    ) -> Result<PolymerEnsemble, PolySimError> {
        let sum: f64 = block_ratios.iter().sum();
        if (sum - 1.0).abs() > 1e-6 {
            return Err(PolySimError::InvalidBlockRatios { sum });
        }

        let stoch = self
            .bigsmiles
            .first_stochastic()
            .ok_or(PolySimError::NoStochasticObject)?;

        if stoch.repeat_units.len() < 2 {
            return Err(PolySimError::RepeatUnitCount {
                architecture: "block copolymer",
                got: stoch.repeat_units.len(),
                need_min: 2,
            });
        }

        if block_ratios.len() != stoch.repeat_units.len() {
            return Err(PolySimError::RepeatUnitCount {
                architecture: "block copolymer (ratios count mismatch)",
                got: block_ratios.len(),
                need_min: stoch.repeat_units.len(),
            });
        }

        let units: Vec<&str> = stoch
            .repeat_units
            .iter()
            .map(|f| f.smiles_raw.as_str())
            .collect();

        // Calibrate: weighted average mass per unit
        let mut unit_masses = Vec::with_capacity(units.len());
        let mut m_end_sum = 0.0;
        for &unit in &units {
            let mw1 = average_mass(&PolymerChain::new(build_linear_smiles(unit, 1)?, 1, 0.0));
            let mw2 = average_mass(&PolymerChain::new(build_linear_smiles(unit, 2)?, 2, 0.0));
            let m0 = mw2 - mw1;
            unit_masses.push(m0);
            m_end_sum += mw1 - m0;
        }
        let m_end = m_end_sum / units.len() as f64;
        let m0_avg: f64 = block_ratios
            .iter()
            .zip(unit_masses.iter())
            .map(|(r, m)| r * m)
            .sum();

        let target_mn_corrected = self.mn - m_end;
        let mut rng = self.make_rng();
        let lengths = self.distribution.sample(
            target_mn_corrected,
            self.pdi,
            m0_avg,
            self.num_chains,
            &mut *rng,
        );

        let chains: Result<Vec<PolymerChain>, PolySimError> = lengths
            .into_iter()
            .map(|n| {
                // Distribute n across blocks proportionally to ratios
                let block_lengths: Vec<usize> = distribute_n_by_ratios(n, block_ratios);
                let sequence: Vec<&str> = block_lengths
                    .iter()
                    .zip(units.iter())
                    .flat_map(|(&len, &unit)| std::iter::repeat_n(unit, len))
                    .collect();
                let smiles = build_copolymer_smiles(&sequence)?;
                let total = sequence.len();
                let chain = PolymerChain::new(smiles, total, 0.0);
                let mn = average_mass(&chain);
                Ok(PolymerChain::new(chain.smiles, total, mn))
            })
            .collect();

        PolymerEnsemble::new(chains?)
    }

    /// Build a polydisperse ensemble of gradient copolymer chains.
    ///
    /// Each chain has composition that varies from `f_start` to `f_end`
    /// according to `profile`. The BigSMILES must contain exactly 2 repeat units.
    pub fn gradient_copolymer_ensemble(
        &self,
        profile: &GradientProfile,
    ) -> Result<PolymerEnsemble, PolySimError> {
        let stoch = self
            .bigsmiles
            .first_stochastic()
            .ok_or(PolySimError::NoStochasticObject)?;

        if stoch.repeat_units.len() != 2 {
            return Err(PolySimError::RepeatUnitCount {
                architecture: "gradient copolymer",
                got: stoch.repeat_units.len(),
                need_min: 2,
            });
        }

        let units: Vec<&str> = stoch
            .repeat_units
            .iter()
            .map(|f| f.smiles_raw.as_str())
            .collect();

        // Calibrate mass per unit for each type
        let mut unit_masses = Vec::with_capacity(2);
        let mut m_end_sum = 0.0;
        for &unit in &units {
            let mw1 = average_mass(&PolymerChain::new(build_linear_smiles(unit, 1)?, 1, 0.0));
            let mw2 = average_mass(&PolymerChain::new(build_linear_smiles(unit, 2)?, 2, 0.0));
            let m0 = mw2 - mw1;
            unit_masses.push(m0);
            m_end_sum += mw1 - m0;
        }
        let m_end = m_end_sum / 2.0;

        // Average mass per unit using the mean gradient fraction across the chain
        // For a chain of variable length, approximate with 100-point average.
        let n_sample = 100usize;
        let avg_f_a: f64 = (0..n_sample)
            .map(|i| gradient_fraction(profile, i, n_sample))
            .sum::<f64>()
            / n_sample as f64;
        let m0_avg = avg_f_a * unit_masses[0] + (1.0 - avg_f_a) * unit_masses[1];

        let target_mn_corrected = self.mn - m_end;
        let mut rng = self.make_rng();
        let lengths = self.distribution.sample(
            target_mn_corrected,
            self.pdi,
            m0_avg,
            self.num_chains,
            &mut *rng,
        );

        let chains: Result<Vec<PolymerChain>, PolySimError> = lengths
            .into_iter()
            .map(|n| {
                let sequence: Vec<&str> = (0..n)
                    .map(|i| {
                        let f_a = gradient_fraction(profile, i, n);
                        let pick: f64 = rng.random();
                        if pick < f_a {
                            units[0]
                        } else {
                            units[1]
                        }
                    })
                    .collect();
                let smiles = build_copolymer_smiles(&sequence)?;
                let chain = PolymerChain::new(smiles, n, 0.0);
                let mn = average_mass(&chain);
                Ok(PolymerChain::new(chain.smiles, n, mn))
            })
            .collect();

        PolymerEnsemble::new(chains?)
    }

    // --- private helpers ---

    fn make_rng(&self) -> Box<dyn rand::RngCore> {
        match self.seed {
            Some(s) => Box::new(StdRng::seed_from_u64(s)),
            None => Box::new(rand::rng()),
        }
    }

    /// Extracts units and computes weighted average mass for copolymer calibration.
    fn copolymer_calibration(
        &self,
        fractions: &[f64],
    ) -> Result<(Vec<&str>, f64, f64), PolySimError> {
        let stoch = self
            .bigsmiles
            .first_stochastic()
            .ok_or(PolySimError::NoStochasticObject)?;

        if stoch.repeat_units.len() < 2 {
            return Err(PolySimError::RepeatUnitCount {
                architecture: "random copolymer",
                got: stoch.repeat_units.len(),
                need_min: 2,
            });
        }

        if fractions.len() != stoch.repeat_units.len() {
            return Err(PolySimError::RepeatUnitCount {
                architecture: "random copolymer (fractions count mismatch)",
                got: fractions.len(),
                need_min: stoch.repeat_units.len(),
            });
        }

        let units: Vec<&str> = stoch
            .repeat_units
            .iter()
            .map(|f| f.smiles_raw.as_str())
            .collect();

        let mut unit_masses = Vec::with_capacity(units.len());
        let mut m_end_sum = 0.0;
        for &unit in &units {
            let mw1 = average_mass(&PolymerChain::new(build_linear_smiles(unit, 1)?, 1, 0.0));
            let mw2 = average_mass(&PolymerChain::new(build_linear_smiles(unit, 2)?, 2, 0.0));
            let m0 = mw2 - mw1;
            unit_masses.push(m0);
            m_end_sum += mw1 - m0;
        }

        let m_end = m_end_sum / units.len() as f64;
        let m0_avg: f64 = fractions
            .iter()
            .zip(unit_masses.iter())
            .map(|(f, m)| f * m)
            .sum();

        Ok((units, m0_avg, m_end))
    }
}

/// Distributes total n across blocks proportionally to ratios.
/// Ensures sum of block lengths equals n (uses largest-remainder method).
fn distribute_n_by_ratios(n: usize, ratios: &[f64]) -> Vec<usize> {
    let n_f = n as f64;
    let raw: Vec<f64> = ratios.iter().map(|r| r * n_f).collect();
    let mut floors: Vec<usize> = raw.iter().map(|&x| x.floor() as usize).collect();
    let mut remainder: Vec<(usize, f64)> = raw
        .iter()
        .enumerate()
        .map(|(i, &x)| (i, x - x.floor()))
        .collect();

    let allocated: usize = floors.iter().sum();
    let mut deficit = n.saturating_sub(allocated);

    // Sort by fractional part descending, allocate remaining units
    remainder.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    for (i, _) in &remainder {
        if deficit == 0 {
            break;
        }
        floors[*i] += 1;
        deficit -= 1;
    }

    floors
}
