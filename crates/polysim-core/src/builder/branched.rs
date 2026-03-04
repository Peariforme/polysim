use bigsmiles::BigSmiles;
use rand::prelude::*;
use rand::rngs::StdRng;

use crate::{
    error::PolySimError,
    polymer::{Architecture, MonomerUnit, PolymerChain},
    properties::molecular_weight::{average_mass, monoisotopic_mass},
};

use super::linear::{
    build_linear_smiles, collect_smiles_segments, max_ring_number, renumber_ring_closures,
    resolve_n_by_mass,
};
use super::strategy::BuildStrategy;

/// Builder for non-linear polymer architectures (comb, graft, star, dendrimer).
///
/// Unlike [`LinearBuilder`](super::linear::LinearBuilder), this builder takes
/// two BigSMILES strings: one for the **backbone** and one for the **branch**.
pub struct BranchedBuilder {
    /// BigSMILES of the backbone chain.
    backbone: BigSmiles,
    /// BigSMILES of the branch / side chain.
    branch: BigSmiles,
    /// Strategy that controls backbone chain length.
    strategy: BuildStrategy,
    /// Optional seed for reproducible random placement.
    seed: Option<u64>,
}

impl BranchedBuilder {
    /// Creates a new builder from backbone and branch BigSMILES strings plus a
    /// build strategy that governs the backbone length.
    pub fn new(backbone: BigSmiles, branch: BigSmiles, strategy: BuildStrategy) -> Self {
        Self {
            backbone,
            branch,
            strategy,
            seed: None,
        }
    }

    /// Set a random seed for reproducible generation.
    pub fn seed(mut self, seed: u64) -> Self {
        self.seed = Some(seed);
        self
    }

    /// Generates a comb (regularly branched) polymer.
    ///
    /// `branch_every` -- attach one branch every N backbone repeat units.
    pub fn comb_polymer(&self, branch_every: usize) -> Result<PolymerChain, PolySimError> {
        let backbone_raw = self.first_repeat_unit(&self.backbone, "comb_polymer backbone")?;
        let branch_raw = self.first_repeat_unit(&self.branch, "comb_polymer branch")?;
        let n = self.resolve_n(&backbone_raw)?;

        if n == 0 {
            return Err(PolySimError::BuildStrategy(
                "repeat count must be >= 1".to_string(),
            ));
        }

        let smiles = build_comb_smiles(&backbone_raw, &branch_raw, n, branch_every)?;
        let smiles = self.with_backbone_end_groups(&smiles);

        let branch_count = (1..=n).filter(|i| i % branch_every == 0).count();
        let total_units = n + branch_count;

        let chain = PolymerChain::new(smiles, total_units, 0.0);
        let mn = average_mass(&chain);

        let backbone_frac = n as f64 / total_units as f64;
        let branch_frac = branch_count as f64 / total_units as f64;
        let composition = vec![
            MonomerUnit::new(&backbone_raw, backbone_frac),
            MonomerUnit::new(&branch_raw, branch_frac),
        ];

        Ok(PolymerChain::new(chain.smiles, total_units, mn)
            .with_composition(composition)
            .with_architecture(Architecture::Comb {
                branch_spacing: branch_every,
            }))
    }

    /// Generates a graft copolymer (random branch-point placement).
    ///
    /// `graft_fraction` -- fraction of backbone repeat units that carry a branch
    /// (0.0 = no grafting, 1.0 = every backbone unit is grafted).
    ///
    /// `seed` -- optional random seed for reproducibility (overrides builder seed).
    pub fn graft_copolymer(
        &self,
        graft_fraction: f64,
        seed: Option<u64>,
    ) -> Result<PolymerChain, PolySimError> {
        let backbone_raw = self.first_repeat_unit(&self.backbone, "graft_copolymer backbone")?;
        let branch_raw = self.first_repeat_unit(&self.branch, "graft_copolymer branch")?;
        let n = self.resolve_n(&backbone_raw)?;

        if n == 0 {
            return Err(PolySimError::BuildStrategy(
                "repeat count must be >= 1".to_string(),
            ));
        }

        let effective_seed = seed.or(self.seed);
        let mut rng: Box<dyn RngCore> = match effective_seed {
            Some(s) => Box::new(StdRng::seed_from_u64(s)),
            None => Box::new(rand::rng()),
        };

        let max_ring_bb = max_ring_number(&backbone_raw);
        if max_ring_bb > 99 {
            return Err(PolySimError::RingNumberOverflow {
                max_ring: max_ring_bb,
                max_supported: 99,
            });
        }
        let cycle_length: usize = if max_ring_bb == 0 {
            usize::MAX
        } else {
            99 / max_ring_bb as usize
        };

        let mut result = String::new();
        let mut branch_count = 0usize;

        for i in 0..n {
            let slot = i % cycle_length;
            let offset = slot as u32 * max_ring_bb;
            let unit = renumber_ring_closures(&backbone_raw, offset);
            result.push_str(&unit);

            let roll: f64 = rng.random();
            if roll < graft_fraction {
                result.push('(');
                result.push_str(&branch_raw);
                result.push(')');
                branch_count += 1;
            }
        }

        let smiles = self.with_backbone_end_groups(&result);
        let total_units = n + branch_count;

        let chain = PolymerChain::new(smiles, total_units, 0.0);
        let mn = average_mass(&chain);

        let backbone_frac = n as f64 / total_units as f64;
        let branch_frac = branch_count as f64 / total_units as f64;
        let composition = vec![
            MonomerUnit::new(&backbone_raw, backbone_frac),
            MonomerUnit::new(&branch_raw, branch_frac),
        ];

        Ok(PolymerChain::new(chain.smiles, total_units, mn)
            .with_composition(composition)
            .with_architecture(Architecture::Graft { graft_fraction }))
    }

    /// Generates a star polymer with `arms` arms radiating from a central atom.
    ///
    /// Each arm is a homopolymer of the backbone BigSMILES repeat unit.
    /// The arm length is determined by the build strategy.
    ///
    /// Supports 3 to 12 arms.
    pub fn star_polymer(&self, arms: usize) -> Result<PolymerChain, PolySimError> {
        if !(3..=12).contains(&arms) {
            return Err(PolySimError::BuildStrategy(format!(
                "star polymer requires 3..=12 arms, got {arms}"
            )));
        }

        let unit_raw = self.first_repeat_unit(&self.backbone, "star_polymer")?;
        let arm_length = self.resolve_n(&unit_raw)?;

        if arm_length == 0 {
            return Err(PolySimError::BuildStrategy(
                "arm length must be >= 1".to_string(),
            ));
        }

        // Build each arm SMILES
        let arm_smiles = build_linear_smiles(&unit_raw, arm_length)?;

        // Star SMILES: C(ARM1)(ARM2)...(ARM_{n-1})ARM_n
        // The central "C" is the hub atom.
        let mut result = String::from("C");
        for i in 0..arms {
            if i < arms - 1 {
                result.push('(');
                result.push_str(&arm_smiles);
                result.push(')');
            } else {
                result.push_str(&arm_smiles);
            }
        }

        let total_units = arms * arm_length;
        let chain = PolymerChain::new(result, total_units, 0.0);
        let mn = average_mass(&chain);

        Ok(PolymerChain::new(chain.smiles, total_units, mn)
            .with_architecture(Architecture::Star { arms }))
    }

    /// Generates a dendrimer of the given `generation` with `branching_factor`
    /// sub-branches per node.
    ///
    /// - `generation`: 1..=6
    /// - `branching_factor`: number of sub-branches per terminal node (typically 2)
    pub fn dendrimer(
        &self,
        generation: usize,
        branching_factor: usize,
    ) -> Result<PolymerChain, PolySimError> {
        if generation == 0 || generation > 6 {
            return Err(PolySimError::BuildStrategy(format!(
                "dendrimer generation must be 1..=6, got {generation}"
            )));
        }

        let unit_raw = self.first_repeat_unit(&self.backbone, "dendrimer")?;

        let smiles = build_dendrimer_smiles(&unit_raw, generation, branching_factor);

        // Total units in a dendrimer:
        // generation g with branching_factor b has:
        //   1 (core) + b + b^2 + ... + b^g = (b^(g+1) - 1) / (b - 1) if b > 1
        //   or g+1 if b == 1
        let total_units = if branching_factor <= 1 {
            generation + 1
        } else {
            (branching_factor.pow(generation as u32 + 1) - 1) / (branching_factor - 1)
        };

        let chain = PolymerChain::new(smiles, total_units, 0.0);
        let mn = average_mass(&chain);

        Ok(PolymerChain::new(chain.smiles, total_units, mn)
            .with_architecture(Architecture::Dendrimer { generation }))
    }

    // --- private helpers -------------------------------------------------------

    /// Extracts the first repeat unit SMILES from a BigSMILES string.
    fn first_repeat_unit(&self, bs: &BigSmiles, context: &str) -> Result<String, PolySimError> {
        let stoch = bs
            .first_stochastic()
            .ok_or(PolySimError::NoStochasticObject)?;

        if stoch.repeat_units.is_empty() {
            return Err(PolySimError::RepeatUnitCount {
                architecture: "branched polymer",
                got: 0,
                need_min: 1,
            });
        }

        let _ = context;
        Ok(stoch.repeat_units[0].smiles_raw.clone())
    }

    /// Resolves repeat count from the build strategy.
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

    /// Prepends prefix and appends suffix SMILES segments from the backbone BigSMILES.
    fn with_backbone_end_groups(&self, body: &str) -> String {
        let prefix = collect_smiles_segments(self.backbone.prefix_segments());
        let suffix = collect_smiles_segments(self.backbone.suffix_segments());
        let mut result = String::with_capacity(prefix.len() + body.len() + suffix.len());
        result.push_str(&prefix);
        result.push_str(body);
        result.push_str(&suffix);
        result
    }
}

// --- free functions -----------------------------------------------------------

/// Builds a comb polymer SMILES by inserting `(branch)` after every
/// `branch_every`-th backbone unit.
fn build_comb_smiles(
    backbone_raw: &str,
    branch_raw: &str,
    n: usize,
    branch_every: usize,
) -> Result<String, PolySimError> {
    let max_ring_bb = max_ring_number(backbone_raw);
    if max_ring_bb > 99 {
        return Err(PolySimError::RingNumberOverflow {
            max_ring: max_ring_bb,
            max_supported: 99,
        });
    }

    let cycle_length: usize = if max_ring_bb == 0 {
        usize::MAX
    } else {
        99 / max_ring_bb as usize
    };

    let mut result = String::new();
    for i in 0..n {
        let slot = i % cycle_length;
        let offset = slot as u32 * max_ring_bb;
        let unit = renumber_ring_closures(backbone_raw, offset);
        result.push_str(&unit);

        if branch_every > 0 && (i + 1) % branch_every == 0 {
            result.push('(');
            result.push_str(branch_raw);
            result.push(')');
        }
    }
    Ok(result)
}

/// Recursively builds the SMILES for a dendrimer.
///
/// At generation 0, returns a single copy of `unit`.
/// At generation g, returns `unit` followed by `branching_factor` branches,
/// each being a sub-dendrimer of generation g-1.
fn build_dendrimer_smiles(unit: &str, generation: usize, branching_factor: usize) -> String {
    if generation == 0 {
        return unit.to_string();
    }

    let sub = build_dendrimer_smiles(unit, generation - 1, branching_factor);
    let mut result = unit.to_string();

    for i in 0..branching_factor {
        if i < branching_factor - 1 {
            result.push('(');
            result.push_str(&sub);
            result.push(')');
        } else {
            // Last branch without wrapping parens (SMILES convention)
            result.push_str(&sub);
        }
    }

    result
}
