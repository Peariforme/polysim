use bigsmiles::{BigSmiles, BigSmilesSegment};
use rand::distr::weighted::WeightedIndex;
use rand::prelude::*;
use rand::rngs::StdRng;

use crate::{
    error::PolySimError,
    polymer::{Architecture, MonomerUnit, PolymerChain},
    properties::molecular_weight::{average_mass, monoisotopic_mass},
};

use super::strategy::BuildStrategy;

/// Gradient composition profile for gradient copolymers.
#[derive(Debug, Clone)]
pub enum GradientProfile {
    /// Linear gradient: f_a(i) = f_start + (f_end - f_start) * i / (n-1)
    Linear { f_start: f64, f_end: f64 },
    /// Sigmoid gradient: f_a(i) = f_start + (f_end - f_start) * sigmoid((i/n - 0.5) * 10)
    Sigmoid { f_start: f64, f_end: f64 },
}

/// Builder for linear polymer architectures.
///
/// Supports homopolymers, random/alternating/block copolymers — all derived
/// from a single BigSMILES string.
pub struct LinearBuilder {
    bigsmiles: BigSmiles,
    strategy: BuildStrategy,
    seed: Option<u64>,
}

impl LinearBuilder {
    /// Creates a new builder from a parsed BigSMILES and a build strategy.
    pub fn new(bigsmiles: BigSmiles, strategy: BuildStrategy) -> Self {
        Self {
            bigsmiles,
            strategy,
            seed: None,
        }
    }

    /// Set a random seed for reproducible copolymer generation.
    pub fn seed(mut self, seed: u64) -> Self {
        self.seed = Some(seed);
        self
    }

    /// Generates a linear homopolymer (single repeat unit, repeated *n* times).
    ///
    /// # Errors
    ///
    /// - [`PolySimError::NoStochasticObject`] if the BigSMILES contains no
    ///   stochastic object (`{...}`).
    /// - [`PolySimError::RepeatUnitCount`] if the stochastic object contains ≠ 1
    ///   repeat unit.
    /// - [`PolySimError::BuildStrategy`] if the strategy yields *n* = 0.
    ///
    /// # Example
    ///
    /// ```rust
    /// use polysim_core::{parse, builder::{linear::LinearBuilder, BuildStrategy}};
    ///
    /// let bs = parse("{[]CC(C)[]}").unwrap(); // polypropylene
    /// let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(3))
    ///     .homopolymer()
    ///     .unwrap();
    ///
    /// assert_eq!(chain.smiles, "CC(C)CC(C)CC(C)");
    /// assert_eq!(chain.repeat_count, 3);
    /// ```
    pub fn homopolymer(&self) -> Result<PolymerChain, PolySimError> {
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

        let fragment = &stoch.repeat_units[0];
        let n = self.resolve_n(&fragment.smiles_raw)?;

        if n == 0 {
            return Err(PolySimError::BuildStrategy(
                "repeat count must be ≥ 1".to_string(),
            ));
        }

        let body = build_linear_smiles(&fragment.smiles_raw, n)?;
        let smiles = self.with_end_groups(&body);
        let chain = PolymerChain::new(smiles, n, 0.0);
        let mn = average_mass(&chain);
        Ok(PolymerChain::new(chain.smiles, n, mn))
    }

    /// Generates a random (statistical) copolymer.
    ///
    /// `fractions` — weight fraction of each repeat unit (must sum to 1.0).
    /// The BigSMILES must contain exactly `fractions.len()` repeat units.
    ///
    /// Uses an optional seed (set via [`Self::seed`]) for reproducibility.
    pub fn random_copolymer(&self, fractions: &[f64]) -> Result<PolymerChain, PolySimError> {
        let sum: f64 = fractions.iter().sum();
        if (sum - 1.0).abs() > 1e-6 {
            return Err(PolySimError::InvalidFractions { sum });
        }

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

        let mut rng: Box<dyn RngCore> = match self.seed {
            Some(s) => Box::new(StdRng::seed_from_u64(s)),
            None => Box::new(rand::rng()),
        };

        let dist = WeightedIndex::new(fractions)
            .map_err(|e| PolySimError::BuildStrategy(format!("invalid weight fractions: {e}")))?;

        let sequence = match &self.strategy {
            BuildStrategy::ByRepeatCount(n) => {
                let n = *n;
                if n == 0 {
                    return Err(PolySimError::BuildStrategy(
                        "repeat count must be ≥ 1".to_string(),
                    ));
                }
                (0..n).map(|_| dist.sample(&mut *rng)).collect::<Vec<_>>()
            }
            BuildStrategy::ByTargetMn(target) => {
                build_incremental_sequence(&units, *target, average_mass, &mut *rng, &dist)?
            }
            BuildStrategy::ByExactMass(target) => {
                build_incremental_sequence(&units, *target, monoisotopic_mass, &mut *rng, &dist)?
            }
        };

        let smiles_seq: Vec<&str> = sequence.iter().map(|&i| units[i]).collect();
        let body = build_copolymer_smiles(&smiles_seq)?;
        let smiles = self.with_end_groups(&body);
        let n = sequence.len();
        let chain = PolymerChain::new(smiles, n, 0.0);
        let mn = average_mass(&chain);
        Ok(PolymerChain::new(chain.smiles, n, mn))
    }

    /// Generates an alternating copolymer (–A–B–A–B– or –A–B–C–A–B–C–).
    ///
    /// The BigSMILES must contain at least 2 repeat units.
    pub fn alternating_copolymer(&self) -> Result<PolymerChain, PolySimError> {
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

        let sequence: Vec<usize> = match &self.strategy {
            BuildStrategy::ByRepeatCount(n) => {
                let n = *n;
                if n == 0 {
                    return Err(PolySimError::BuildStrategy(
                        "repeat count must be ≥ 1".to_string(),
                    ));
                }
                (0..n).map(|i| i % k).collect()
            }
            BuildStrategy::ByTargetMn(target) => {
                build_incremental_alternating(&units, *target, average_mass)?
            }
            BuildStrategy::ByExactMass(target) => {
                build_incremental_alternating(&units, *target, monoisotopic_mass)?
            }
        };

        let smiles_seq: Vec<&str> = sequence.iter().map(|&i| units[i]).collect();
        let body = build_copolymer_smiles(&smiles_seq)?;
        let smiles = self.with_end_groups(&body);
        let n = sequence.len();
        let chain = PolymerChain::new(smiles, n, 0.0);
        let mn = average_mass(&chain);
        Ok(PolymerChain::new(chain.smiles, n, mn))
    }

    /// Generates a block copolymer (–AAAA–BBBB–).
    ///
    /// `block_lengths` — number of repeat units per block, in order.
    /// The BigSMILES must contain exactly `block_lengths.len()` repeat units.
    ///
    /// The `BuildStrategy` is ignored — `block_lengths` fully determines the chain.
    pub fn block_copolymer(&self, block_lengths: &[usize]) -> Result<PolymerChain, PolySimError> {
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

        if block_lengths.len() != stoch.repeat_units.len() {
            return Err(PolySimError::RepeatUnitCount {
                architecture: "block copolymer (block_lengths count mismatch)",
                got: block_lengths.len(),
                need_min: stoch.repeat_units.len(),
            });
        }

        let units: Vec<&str> = stoch
            .repeat_units
            .iter()
            .map(|f| f.smiles_raw.as_str())
            .collect();

        let smiles_seq: Vec<&str> = block_lengths
            .iter()
            .zip(units.iter())
            .flat_map(|(&len, &unit)| std::iter::repeat_n(unit, len))
            .collect();

        let n = smiles_seq.len();
        if n == 0 {
            return Err(PolySimError::BuildStrategy(
                "total block length must be ≥ 1".to_string(),
            ));
        }

        let body = build_copolymer_smiles(&smiles_seq)?;
        let smiles = self.with_end_groups(&body);
        let chain = PolymerChain::new(smiles, n, 0.0);
        let mn = average_mass(&chain);
        Ok(PolymerChain::new(chain.smiles, n, mn))
    }

    /// Generates a gradient copolymer where the composition of monomer A varies
    /// along the chain according to the given [`GradientProfile`].
    ///
    /// The BigSMILES must contain exactly 2 repeat units (A and B).
    pub fn gradient_copolymer(
        &self,
        profile: &GradientProfile,
    ) -> Result<PolymerChain, PolySimError> {
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

        // Resolve chain length using unit A
        let n = self.resolve_n(units[0])?;
        if n == 0 {
            return Err(PolySimError::BuildStrategy(
                "repeat count must be >= 1".to_string(),
            ));
        }

        let mut rng: Box<dyn RngCore> = match self.seed {
            Some(s) => Box::new(StdRng::seed_from_u64(s)),
            None => Box::new(rand::rng()),
        };

        let mut count_a: usize = 0;
        let mut sequence = Vec::with_capacity(n);

        for i in 0..n {
            let f_a = gradient_fraction(profile, i, n);
            let pick: f64 = rng.random();
            let idx = if pick < f_a { 0 } else { 1 };
            if idx == 0 {
                count_a += 1;
            }
            sequence.push(idx);
        }

        let smiles_seq: Vec<&str> = sequence.iter().map(|&i| units[i]).collect();
        let body = build_copolymer_smiles(&smiles_seq)?;
        let smiles = self.with_end_groups(&body);
        let chain = PolymerChain::new(smiles, n, 0.0);
        let mn = average_mass(&chain);

        let frac_a = count_a as f64 / n as f64;
        let composition = vec![
            MonomerUnit::new(units[0], frac_a),
            MonomerUnit::new(units[1], 1.0 - frac_a),
        ];

        Ok(PolymerChain::new(chain.smiles, n, mn)
            .with_composition(composition)
            .with_architecture(Architecture::Gradient))
    }

    /// Generates a cyclic homopolymer (ring closure connecting first and last atom).
    ///
    /// The BigSMILES must contain exactly 1 repeat unit.
    pub fn cyclic_homopolymer(&self) -> Result<PolymerChain, PolySimError> {
        let stoch = self
            .bigsmiles
            .first_stochastic()
            .ok_or(PolySimError::NoStochasticObject)?;

        if stoch.repeat_units.len() != 1 {
            return Err(PolySimError::RepeatUnitCount {
                architecture: "cyclic homopolymer",
                got: stoch.repeat_units.len(),
                need_min: 1,
            });
        }

        let fragment = &stoch.repeat_units[0];
        let n = self.resolve_n(&fragment.smiles_raw)?;

        if n == 0 {
            return Err(PolySimError::BuildStrategy(
                "repeat count must be >= 1".to_string(),
            ));
        }

        let linear = build_linear_smiles(&fragment.smiles_raw, n)?;
        let smiles = make_cyclic_smiles(&linear);
        let chain = PolymerChain::new(smiles, n, 0.0);
        let mn = average_mass(&chain);
        Ok(PolymerChain::new(chain.smiles, n, mn).with_architecture(Architecture::Cyclic))
    }

    /// Prepends prefix and appends suffix SMILES segments from the BigSMILES.
    fn with_end_groups(&self, body: &str) -> String {
        let prefix = collect_smiles_segments(self.bigsmiles.prefix_segments());
        let suffix = collect_smiles_segments(self.bigsmiles.suffix_segments());
        let mut result = String::with_capacity(prefix.len() + body.len() + suffix.len());
        result.push_str(&prefix);
        result.push_str(body);
        result.push_str(&suffix);
        result
    }

    fn resolve_n(&self, smiles_raw: &str) -> Result<usize, PolySimError> {
        match &self.strategy {
            BuildStrategy::ByRepeatCount(n) => Ok(*n),
            BuildStrategy::ByTargetMn(target) => {
                resolve_n_by_mass(smiles_raw, *target, average_mass)
            }
            BuildStrategy::ByExactMass(target) => {
                resolve_n_by_mass(smiles_raw, *target, monoisotopic_mass)
            }
        }
    }
}

// --- internal helpers -------------------------------------------------------

/// Déduit le nombre de répétitions à partir d'une masse cible.
///
/// Construit deux chaînes d'essai (n=1 et n=2) pour déterminer la masse par
/// unité et la masse des groupements terminaux, puis résout par extrapolation
/// linéaire : MW(n) = n × mw_per_unit + mw_end.
///
/// `mass_fn` peut être [`average_mass`] (pour [`BuildStrategy::ByTargetMn`]) ou
/// [`monoisotopic_mass`] (pour [`BuildStrategy::ByExactMass`]).
pub(crate) fn resolve_n_by_mass(
    smiles_raw: &str,
    target: f64,
    mass_fn: fn(&PolymerChain) -> f64,
) -> Result<usize, PolySimError> {
    let mw1 = mass_fn(&PolymerChain::new(
        build_linear_smiles(smiles_raw, 1)?,
        1,
        0.0,
    ));
    let mw2 = mass_fn(&PolymerChain::new(
        build_linear_smiles(smiles_raw, 2)?,
        2,
        0.0,
    ));
    let mw_per_unit = mw2 - mw1;
    let mw_end = mw1 - mw_per_unit;
    let n = ((target - mw_end) / mw_per_unit).round().max(1.0) as usize;
    Ok(n)
}

/// Calibrates per-unit masses for each distinct repeat unit via 2-point method.
///
/// Returns `(unit_masses, m_end)` where:
/// - `unit_masses[i]` is the mass contribution of one copy of unit i
/// - `m_end` is the end-group mass (constant across all units)
fn calibrate_unit_masses(
    units: &[&str],
    mass_fn: fn(&PolymerChain) -> f64,
) -> Result<(Vec<f64>, f64), PolySimError> {
    let mut unit_masses = Vec::with_capacity(units.len());
    let mut m_end_sum = 0.0;

    for &unit in units {
        let mw1 = mass_fn(&PolymerChain::new(build_linear_smiles(unit, 1)?, 1, 0.0));
        let mw2 = mass_fn(&PolymerChain::new(build_linear_smiles(unit, 2)?, 2, 0.0));
        let m0 = mw2 - mw1;
        unit_masses.push(m0);
        m_end_sum += mw1 - m0;
    }

    // Average end-group mass (should be ~identical for all units, but average for safety)
    let m_end = m_end_sum / units.len() as f64;
    Ok((unit_masses, m_end))
}

/// Builds a copolymer unit sequence incrementally for random copolymers.
///
/// Adds units one at a time (sampled from weighted distribution), tracking
/// accumulated mass. Stops when closest to the target mass.
fn build_incremental_sequence(
    units: &[&str],
    target: f64,
    mass_fn: fn(&PolymerChain) -> f64,
    rng: &mut dyn RngCore,
    dist: &WeightedIndex<f64>,
) -> Result<Vec<usize>, PolySimError> {
    let (unit_masses, m_end) = calibrate_unit_masses(units, mass_fn)?;

    let mut sequence = Vec::new();
    let mut running_mass = m_end;

    loop {
        let idx = dist.sample(rng);
        running_mass += unit_masses[idx];
        sequence.push(idx);

        if running_mass >= target {
            // Check if removing last unit gets closer
            if sequence.len() > 1 {
                let mass_without = running_mass - unit_masses[idx];
                if (mass_without - target).abs() < (running_mass - target).abs() {
                    sequence.pop();
                }
            }
            break;
        }
    }

    Ok(sequence)
}

/// Builds a copolymer unit sequence incrementally for alternating copolymers.
fn build_incremental_alternating(
    units: &[&str],
    target: f64,
    mass_fn: fn(&PolymerChain) -> f64,
) -> Result<Vec<usize>, PolySimError> {
    let (unit_masses, m_end) = calibrate_unit_masses(units, mass_fn)?;
    let k = units.len();

    let mut sequence = Vec::new();
    let mut running_mass = m_end;

    loop {
        let idx = sequence.len() % k;
        running_mass += unit_masses[idx];
        sequence.push(idx);

        if running_mass >= target {
            if sequence.len() > 1 {
                let mass_without = running_mass - unit_masses[idx];
                if (mass_without - target).abs() < (running_mass - target).abs() {
                    sequence.pop();
                }
            }
            break;
        }
    }

    Ok(sequence)
}

/// Builds the SMILES string for a linear chain of `n` repeat units.
///
/// Ring closure numbers are renumbered for each copy. Because each copy is
/// self-contained (every ring opened within a copy is also closed within that
/// copy), the offsets cycle over 1..=99, allowing chains of arbitrary length.
///
/// # Errors
///
/// Returns [`PolySimError::RingNumberOverflow`] if the repeat unit itself uses
/// more than 99 distinct ring-closure numbers (already invalid SMILES).
pub(crate) fn build_linear_smiles(smiles_raw: &str, n: usize) -> Result<String, PolySimError> {
    let max_ring = max_ring_number(smiles_raw);

    // Pathological case: the repeat unit alone already overflows SMILES ring numbers.
    if max_ring > 99 {
        return Err(PolySimError::RingNumberOverflow {
            max_ring,
            max_supported: 99,
        });
    }

    // Number of distinct copies before ring numbers must be recycled.
    // Since each copy closes its own rings before the next copy starts,
    // the same numbers can be safely reused.
    let cycle_length: usize = if max_ring == 0 {
        usize::MAX // no ring closures — no cycling needed
    } else {
        99 / max_ring as usize
    };

    let mut result = String::with_capacity(smiles_raw.len() * n);
    for i in 0..n {
        let slot = i % cycle_length;
        let offset = slot as u32 * max_ring;
        result.push_str(&renumber_ring_closures(smiles_raw, offset));
    }
    Ok(result)
}

/// Builds the SMILES string for a copolymer from a heterogeneous sequence of
/// repeat-unit SMILES fragments.
///
/// Ring closure numbers are renumbered globally so they never collide across
/// consecutive units, regardless of which unit type follows which.
pub(crate) fn build_copolymer_smiles(unit_sequence: &[&str]) -> Result<String, PolySimError> {
    // Compute max ring number across ALL distinct units.
    let global_max_ring = unit_sequence
        .iter()
        .map(|u| max_ring_number(u))
        .max()
        .unwrap_or(0);

    if global_max_ring > 99 {
        return Err(PolySimError::RingNumberOverflow {
            max_ring: global_max_ring,
            max_supported: 99,
        });
    }

    let cycle_length: usize = if global_max_ring == 0 {
        usize::MAX
    } else {
        99 / global_max_ring as usize
    };

    let total_len: usize = unit_sequence.iter().map(|u| u.len()).sum();
    let mut result = String::with_capacity(total_len + unit_sequence.len() * 4);

    for (i, &unit) in unit_sequence.iter().enumerate() {
        let slot = i % cycle_length;
        let offset = slot as u32 * global_max_ring;
        result.push_str(&renumber_ring_closures(unit, offset));
    }

    Ok(result)
}

/// Returns the highest ring-closure number used in a SMILES string.
///
/// Digits inside `[...]` (isotopes, hydrogen counts, charges, atom classes)
/// are ignored.
pub(crate) fn max_ring_number(smiles: &str) -> u32 {
    let mut max = 0u32;
    let mut in_bracket = false;
    let mut chars = smiles.chars().peekable();

    while let Some(c) = chars.next() {
        match c {
            '[' => in_bracket = true,
            ']' => in_bracket = false,
            _ if in_bracket => {}
            '%' => {
                // Two-digit notation: %dd
                let d1 = chars.next().unwrap_or('0');
                let d2 = chars.next().unwrap_or('0');
                if d1.is_ascii_digit() && d2.is_ascii_digit() {
                    let n = (d1 as u32 - '0' as u32) * 10 + (d2 as u32 - '0' as u32);
                    max = max.max(n);
                }
            }
            c if c.is_ascii_digit() => {
                max = max.max(c as u32 - '0' as u32);
            }
            _ => {}
        }
    }
    max
}

/// Returns a copy of `smiles` with every ring-closure number incremented by `offset`.
///
/// When `offset` is 0 the string is returned unchanged.
/// Digits inside `[...]` are never modified.
pub(crate) fn renumber_ring_closures(smiles: &str, offset: u32) -> String {
    if offset == 0 {
        return smiles.to_string();
    }
    let mut result = String::with_capacity(smiles.len() + 4);
    let mut in_bracket = false;
    let mut chars = smiles.chars().peekable();

    while let Some(c) = chars.next() {
        match c {
            '[' => {
                in_bracket = true;
                result.push(c);
            }
            ']' => {
                in_bracket = false;
                result.push(c);
            }
            _ if in_bracket => result.push(c),
            '%' => {
                let d1 = chars.next().unwrap_or('0');
                let d2 = chars.next().unwrap_or('0');
                if d1.is_ascii_digit() && d2.is_ascii_digit() {
                    let n = (d1 as u32 - '0' as u32) * 10 + (d2 as u32 - '0' as u32);
                    let new_n = n + offset;
                    result.push('%');
                    result.push_str(&format!("{new_n:02}"));
                } else {
                    result.push('%');
                    result.push(d1);
                    result.push(d2);
                }
            }
            c if c.is_ascii_digit() => {
                let n = c as u32 - '0' as u32;
                let new_n = n + offset;
                if new_n <= 9 {
                    result.push(char::from_digit(new_n, 10).unwrap());
                } else {
                    result.push('%');
                    result.push_str(&format!("{new_n:02}"));
                }
            }
            _ => result.push(c),
        }
    }
    result
}

/// Extracts plain SMILES text from a slice of BigSMILES segments,
/// ignoring any stochastic objects.
pub(crate) fn collect_smiles_segments(segs: &[BigSmilesSegment]) -> String {
    segs.iter()
        .filter_map(|s| match s {
            BigSmilesSegment::Smiles(mol) => Some(format!("{mol}")),
            _ => None,
        })
        .collect()
}

/// Computes the fraction of monomer A at position `i` in a chain of length `n`.
pub(crate) fn gradient_fraction(profile: &GradientProfile, i: usize, n: usize) -> f64 {
    match profile {
        GradientProfile::Linear { f_start, f_end } => {
            if n <= 1 {
                *f_start
            } else {
                f_start + (f_end - f_start) * i as f64 / (n - 1) as f64
            }
        }
        GradientProfile::Sigmoid { f_start, f_end } => {
            if n <= 1 {
                *f_start
            } else {
                let x = (i as f64 / n as f64 - 0.5) * 10.0;
                let sigma = 1.0 / (1.0 + (-x).exp());
                f_start + (f_end - f_start) * sigma
            }
        }
    }
}

/// Converts a linear SMILES into a cyclic one by inserting ring closure label "1"
/// after the first atom and appending "1" at the end.
///
/// Handles bracket atoms (`[...]`) and two-letter organic atoms (`Cl`, `Br`).
fn make_cyclic_smiles(linear: &str) -> String {
    let mut result = String::with_capacity(linear.len() + 2);
    let mut chars = linear.chars().peekable();

    // Find and copy the first atom, then insert "1"
    if let Some(c) = chars.next() {
        if c == '[' {
            // Bracket atom: copy up to and including ']'
            result.push(c);
            for ch in chars.by_ref() {
                result.push(ch);
                if ch == ']' {
                    break;
                }
            }
        } else {
            result.push(c);
            // Check for two-letter organic atoms (Cl, Br, Si, etc.)
            if c.is_ascii_uppercase() {
                if let Some(&next) = chars.peek() {
                    if next.is_ascii_lowercase() && next != 'c'
                        || matches!(
                            (c, next),
                            ('C', 'l')
                                | ('B', 'r')
                                | ('S', 'i')
                                | ('S', 'e')
                                | ('A', 'l')
                                | ('A', 's')
                                | ('A', 'r')
                                | ('A', 't')
                                | ('M', 'g')
                                | ('N', 'a')
                                | ('G', 'e')
                        )
                    {
                        result.push(chars.next().unwrap());
                    }
                }
            }
        }
        result.push('1');
    }

    // Copy the rest
    for ch in chars {
        result.push(ch);
    }

    // Append ring closure at end
    result.push('1');
    result
}
