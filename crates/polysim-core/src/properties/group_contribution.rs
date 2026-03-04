//! Group contribution method infrastructure for predicting polymer properties.
//!
//! This module implements the Van Krevelen group-contribution approach, which
//! decomposes a polymer repeat unit into functional groups and sums their
//! additive contributions to estimate bulk properties.
//!
//! # Reference
//!
//! Van Krevelen, D. W. & te Nijenhuis, K. (2009).
//! *Properties of Polymers*, 4th ed., Elsevier.

use opensmiles::{
    ast::{BondType, Molecule},
    parse as parse_smiles,
};

use crate::error::PolySimError;
use crate::polymer::PolymerChain;

// ---------------------------------------------------------------------------
// Core types
// ---------------------------------------------------------------------------

/// A functional group with its additive property contributions.
///
/// Each group carries Van Krevelen contributions that are summed to predict
/// macroscopic properties of the polymer.
#[derive(Debug, Clone, PartialEq)]
pub struct Group {
    /// Human-readable name (e.g. "-CH2-").
    pub name: &'static str,
    /// Molar contribution to Tg via Van Krevelen (K * g/mol).
    pub yg: f64,
    /// Molar contribution to Tm via Van Krevelen (K * g/mol). 0.0 for amorphous groups.
    pub ym: f64,
    /// Van der Waals volume (cm^3/mol).
    pub vw: f64,
    /// Cohesive energy (J/mol).
    pub ecoh: f64,
    /// Molar refraction Lorentz-Lorenz (cm^3/mol).
    pub ri: f64,
}

/// Result of matching a single group in the decomposition.
#[derive(Debug, Clone)]
pub struct GroupMatch {
    /// The matched group.
    pub group: &'static Group,
    /// Number of occurrences found.
    pub count: usize,
}

/// Trait for group-contribution based property prediction.
///
/// Implementors consume a set of `GroupMatch` results and return the predicted
/// property value.
pub trait GroupContributionMethod {
    /// Predicts a property value from group-contribution sums.
    fn predict(&self, groups: &[GroupMatch]) -> f64;
}

// ---------------------------------------------------------------------------
// Static group database (Van Krevelen, 1990 / 2009)
// ---------------------------------------------------------------------------

/// Methyl group -CH3.
static GROUP_CH3: Group = Group {
    name: "-CH3",
    yg: 4.0,
    ym: 2.0,
    vw: 13.67,
    ecoh: 4500.0,
    ri: 5.67,
};

/// Methylene group -CH2-.
static GROUP_CH2: Group = Group {
    name: "-CH2-",
    yg: 2.7,
    ym: 4.7,
    vw: 10.23,
    ecoh: 4100.0,
    ri: 4.65,
};

/// Methine group -CH<.
static GROUP_CH: Group = Group {
    name: "-CH<",
    yg: 1.0,
    ym: 3.0,
    vw: 6.78,
    ecoh: 3400.0,
    ri: 3.63,
};

/// Quaternary carbon >C<.
static GROUP_C: Group = Group {
    name: ">C<",
    yg: -1.0,
    ym: 2.2,
    vw: 3.33,
    ecoh: 2100.0,
    ri: 2.61,
};

/// Phenyl group -C6H5 (pendant aromatic ring).
static GROUP_PHENYL: Group = Group {
    name: "-C6H5",
    yg: 31.0,
    ym: 35.0,
    vw: 71.6,
    ecoh: 31900.0,
    ri: 25.93,
};

/// Para-phenylene group -C6H4- (in-chain aromatic ring).
static GROUP_PHENYLENE: Group = Group {
    name: "-C6H4-",
    yg: 28.5,
    ym: 33.0,
    vw: 67.0,
    ecoh: 31500.0,
    ri: 24.5,
};

/// Ether group -O-.
static GROUP_ETHER: Group = Group {
    name: "-O-",
    yg: 5.3,
    ym: 5.0,
    vw: 5.0,
    ecoh: 4200.0,
    ri: 1.64,
};

/// Ester group -COO-.
static GROUP_ESTER: Group = Group {
    name: "-COO-",
    yg: 25.0,
    ym: 18.0,
    vw: 22.0,
    ecoh: 18000.0,
    ri: 6.38,
};

/// Ketone group -CO-.
static GROUP_KETONE: Group = Group {
    name: "-CO-",
    yg: 12.0,
    ym: 15.0,
    vw: 17.0,
    ecoh: 17400.0,
    ri: 4.6,
};

/// Hydroxyl group -OH.
static GROUP_OH: Group = Group {
    name: "-OH",
    yg: 23.6,
    ym: 35.0,
    vw: 10.0,
    ecoh: 29800.0,
    ri: 2.75,
};

/// Carboxylic acid group -COOH.
static GROUP_COOH: Group = Group {
    name: "-COOH",
    yg: 45.0,
    ym: 45.0,
    vw: 28.5,
    ecoh: 27000.0,
    ri: 6.42,
};

/// Secondary amide group -CONH-.
static GROUP_AMIDE: Group = Group {
    name: "-CONH-",
    yg: 40.0,
    ym: 60.0,
    vw: 23.6,
    ecoh: 36000.0,
    ri: 7.35,
};

/// Primary amide group -CONH2.
static GROUP_AMIDE_PRIMARY: Group = Group {
    name: "-CONH2",
    yg: 50.0,
    ym: 70.0,
    vw: 23.6,
    ecoh: 50000.0,
    ri: 7.35,
};

/// Nitrile group -CN.
static GROUP_CN: Group = Group {
    name: "-CN",
    yg: 15.0,
    ym: 30.0,
    vw: 15.0,
    ecoh: 24000.0,
    ri: 5.55,
};

/// Chloro group -Cl.
static GROUP_CL: Group = Group {
    name: "-Cl",
    yg: 16.0,
    ym: 15.0,
    vw: 12.0,
    ecoh: 12800.0,
    ri: 5.84,
};

/// Fluoro group -F.
static GROUP_F: Group = Group {
    name: "-F",
    yg: 5.0,
    ym: 8.0,
    vw: 5.8,
    ecoh: 4200.0,
    ri: 0.81,
};

/// Siloxane group -Si-O-.
static GROUP_SILOXANE: Group = Group {
    name: "-SiO-",
    yg: 0.3,
    ym: 0.5,
    vw: 21.0,
    ecoh: 4200.0,
    ri: 6.5,
};

// ---------------------------------------------------------------------------
// Group database & SMILES decomposition
// ---------------------------------------------------------------------------

/// Database of Van Krevelen functional groups with SMILES decomposition logic.
///
/// The database decomposes a SMILES string into its constituent functional
/// groups by analysing atom types, bond connectivity, and ring membership.
pub struct GroupDatabase;

impl GroupDatabase {
    /// Decomposes a polymer chain into functional groups using the opensmiles
    /// molecular graph.
    ///
    /// The algorithm:
    /// 1. Identify aromatic 6-membered carbon rings (phenyl / phenylene).
    /// 2. Identify multi-atom functional groups (-COO-, -CONH-, -CONH2, -COOH, -CO-, -CN).
    /// 3. Identify heteroatom pendant groups (-OH, -Cl, -F, -O-, -SiO-).
    /// 4. Classify remaining aliphatic carbons by hydrogen count (CH3, CH2, CH, C).
    ///
    /// # Errors
    ///
    /// Returns `PolySimError::GroupDecomposition` if the SMILES cannot be parsed.
    pub fn decompose(chain: &PolymerChain) -> Result<Vec<GroupMatch>, PolySimError> {
        let mol = parse_smiles(&chain.smiles)
            .map_err(|e| PolySimError::GroupDecomposition(format!("invalid SMILES: {e}")))?;

        let num_atoms = mol.nodes().len();

        // Track which atoms have been consumed by a multi-atom group.
        let mut consumed = vec![false; num_atoms];

        // Build adjacency list for connectivity queries.
        let adj = build_adjacency(&mol);

        // Counters for each group.
        let mut ch3 = 0usize;
        let mut ch2 = 0usize;
        let mut ch = 0usize;
        let mut c_quat = 0usize;
        let mut phenyl = 0usize;
        let mut phenylene = 0usize;
        let mut ether = 0usize;
        let mut ester = 0usize;
        let mut ketone = 0usize;
        let mut oh = 0usize;
        let mut cooh = 0usize;
        let mut amide = 0usize;
        let mut amide_primary = 0usize;
        let mut cn = 0usize;
        let mut cl = 0usize;
        let mut f = 0usize;
        let mut siloxane = 0usize;

        // --- Pass 1: Aromatic rings ---
        let rings = mol.aromatic_rings();
        for ring in &rings {
            if ring.size() != 6 {
                continue;
            }
            // Check all ring atoms are carbon.
            let all_carbon = ring
                .nodes
                .iter()
                .all(|&idx| mol.nodes()[idx as usize].atom().element().atomic_number() == 6);
            if !all_carbon {
                continue;
            }

            // Count heavy-atom neighbours outside the ring to determine
            // whether this is a pendant phenyl (-C6H5) or in-chain phenylene (-C6H4-).
            let mut external_heavy_bonds = 0usize;
            for &idx in &ring.nodes {
                for &(neighbour, _bond_type) in &adj[idx as usize] {
                    if !ring.nodes.contains(&(neighbour as u16)) {
                        let n_atomic = mol.nodes()[neighbour].atom().element().atomic_number();
                        if n_atomic != 0 && n_atomic != 1 {
                            external_heavy_bonds += 1;
                        }
                    }
                }
            }

            if external_heavy_bonds >= 2 {
                phenylene += 1;
            } else {
                phenyl += 1;
            }

            // Mark ring atoms as consumed.
            for &idx in &ring.nodes {
                consumed[idx as usize] = true;
            }
        }

        // --- Pass 2: Multi-atom functional groups on non-consumed atoms ---

        // Helper: check if atom at idx is element with given atomic number.
        let is_element =
            |idx: usize, z: u8| -> bool { mol.nodes()[idx].atom().element().atomic_number() == z };

        // Detect -SiO- (siloxane): Si bonded to O.
        for i in 0..num_atoms {
            if consumed[i] || !is_element(i, 14) {
                continue; // not Si
            }
            // Find an adjacent O that is not consumed.
            let mut found_o = None;
            for &(nb, _) in &adj[i] {
                if !consumed[nb] && is_element(nb, 8) {
                    found_o = Some(nb);
                    break;
                }
            }
            if let Some(o_idx) = found_o {
                siloxane += 1;
                consumed[i] = true;
                consumed[o_idx] = true;
            }
        }

        // Detect -COO- (ester), -COOH, -CONH-, -CONH2, -CO- (ketone), -CN (nitrile).
        // We iterate over carbon atoms bonded to O or N via double/single bonds.
        for i in 0..num_atoms {
            if consumed[i] || !is_element(i, 6) {
                continue;
            }

            // Find double-bonded O neighbour (C=O).
            let mut double_o: Option<usize> = None;
            // Find single-bonded O neighbour.
            let mut single_o: Vec<usize> = Vec::new();
            // Find N neighbours.
            let mut n_neighbours: Vec<usize> = Vec::new();
            // Find triple-bonded N (C#N / nitrile).
            let mut triple_n: Option<usize> = None;

            for &(nb, bt) in &adj[i] {
                if consumed[nb] {
                    continue;
                }
                let z = mol.nodes()[nb].atom().element().atomic_number();
                match (z, bt) {
                    (8, BondType::Double) => {
                        double_o = Some(nb);
                    }
                    (8, _) => {
                        single_o.push(nb);
                    }
                    (7, BondType::Triple) => {
                        triple_n = Some(nb);
                    }
                    (7, _) => {
                        n_neighbours.push(nb);
                    }
                    _ => {}
                }
            }

            // -CN (nitrile): C#N
            if let Some(n_idx) = triple_n {
                cn += 1;
                consumed[i] = true;
                consumed[n_idx] = true;
                continue;
            }

            if let Some(o_dbl) = double_o {
                // We have C=O.

                // -COOH: C(=O)(OH) where OH has 1 hydrogen.
                let oh_idx = single_o
                    .iter()
                    .find(|&&idx| mol.nodes()[idx].hydrogens() >= 1);

                // -COO- (ester): C(=O)(O-R) where O is bonded to another heavy atom.
                let ester_o = single_o
                    .iter()
                    .find(|&&idx| mol.nodes()[idx].hydrogens() == 0);

                if !n_neighbours.is_empty() {
                    let n_idx = n_neighbours[0];
                    let n_h = mol.nodes()[n_idx].hydrogens();
                    if n_h >= 2 {
                        // -CONH2 (primary amide)
                        amide_primary += 1;
                    } else {
                        // -CONH- (secondary amide)
                        amide += 1;
                    }
                    consumed[i] = true;
                    consumed[o_dbl] = true;
                    consumed[n_idx] = true;
                } else if let Some(&oh_i) = oh_idx {
                    // -COOH
                    cooh += 1;
                    consumed[i] = true;
                    consumed[o_dbl] = true;
                    consumed[oh_i] = true;
                } else if let Some(&ester_i) = ester_o {
                    // -COO- (ester)
                    ester += 1;
                    consumed[i] = true;
                    consumed[o_dbl] = true;
                    consumed[ester_i] = true;
                } else {
                    // -CO- (ketone / aldehyde)
                    ketone += 1;
                    consumed[i] = true;
                    consumed[o_dbl] = true;
                }
            }
        }

        // --- Pass 3: Remaining heteroatom pendant groups ---
        for (flag, node) in consumed.iter_mut().zip(mol.nodes().iter()) {
            if *flag {
                continue;
            }
            let z = node.atom().element().atomic_number();
            match z {
                // Oxygen: -OH if has hydrogen, otherwise ether -O-.
                8 => {
                    if node.hydrogens() >= 1 {
                        oh += 1;
                    } else {
                        ether += 1;
                    }
                    *flag = true;
                }
                // Chlorine
                17 => {
                    cl += 1;
                    *flag = true;
                }
                // Fluorine
                9 => {
                    f += 1;
                    *flag = true;
                }
                // Nitrogen not consumed by amide/nitrile detection: skip for now
                // (amine groups would be an extension)
                _ => {}
            }
        }

        // --- Pass 4: Remaining aliphatic carbons ---
        for (flag, node) in consumed.iter_mut().zip(mol.nodes().iter()) {
            if *flag {
                continue;
            }
            let z = node.atom().element().atomic_number();
            if z != 6 {
                continue;
            }
            match node.hydrogens() {
                3 => ch3 += 1,
                2 => ch2 += 1,
                1 => ch += 1,
                0 => c_quat += 1,
                _ => ch3 += 1, // CH4 counted as CH3 (terminal methane-like)
            }
            *flag = true;
        }

        // --- Build result ---
        let mut matches = Vec::new();
        let mut push = |group: &'static Group, count: usize| {
            if count > 0 {
                matches.push(GroupMatch { group, count });
            }
        };
        push(&GROUP_CH3, ch3);
        push(&GROUP_CH2, ch2);
        push(&GROUP_CH, ch);
        push(&GROUP_C, c_quat);
        push(&GROUP_PHENYL, phenyl);
        push(&GROUP_PHENYLENE, phenylene);
        push(&GROUP_ETHER, ether);
        push(&GROUP_ESTER, ester);
        push(&GROUP_KETONE, ketone);
        push(&GROUP_OH, oh);
        push(&GROUP_COOH, cooh);
        push(&GROUP_AMIDE, amide);
        push(&GROUP_AMIDE_PRIMARY, amide_primary);
        push(&GROUP_CN, cn);
        push(&GROUP_CL, cl);
        push(&GROUP_F, f);
        push(&GROUP_SILOXANE, siloxane);

        Ok(matches)
    }
}

// ---------------------------------------------------------------------------
// Convenience summation helpers
// ---------------------------------------------------------------------------

/// Sums a given property contribution across all matched groups.
///
/// `extract` selects which field of `Group` to use (e.g. `|g| g.yg`).
pub fn sum_contribution(groups: &[GroupMatch], extract: fn(&Group) -> f64) -> f64 {
    groups
        .iter()
        .map(|gm| extract(gm.group) * gm.count as f64)
        .sum()
}

/// Total Van der Waals volume (cm^3/mol) from group contributions.
pub fn total_vw(groups: &[GroupMatch]) -> f64 {
    sum_contribution(groups, |g| g.vw)
}

/// Total cohesive energy (J/mol) from group contributions.
pub fn total_ecoh(groups: &[GroupMatch]) -> f64 {
    sum_contribution(groups, |g| g.ecoh)
}

/// Total molar refraction (cm^3/mol) from group contributions.
pub fn total_ri(groups: &[GroupMatch]) -> f64 {
    sum_contribution(groups, |g| g.ri)
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Adjacency list entry: (neighbour_index, bond_type).
type AdjList = Vec<Vec<(usize, BondType)>>;

/// Build an adjacency list from the opensmiles `Molecule`.
fn build_adjacency(mol: &Molecule) -> AdjList {
    let n = mol.nodes().len();
    let mut adj: AdjList = vec![Vec::new(); n];
    for bond in mol.bonds() {
        let s = bond.source() as usize;
        let t = bond.target() as usize;
        let k = bond.kind();
        adj[s].push((t, k));
        adj[t].push((s, k));
    }
    adj
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::{linear::LinearBuilder, BuildStrategy};

    /// Helper to build a homopolymer chain from a BigSMILES string.
    fn build_chain(bigsmiles: &str, n: usize) -> PolymerChain {
        let bs = crate::parse(bigsmiles).expect("valid BigSMILES");
        LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(n))
            .homopolymer()
            .expect("build should succeed")
    }

    #[test]
    fn polyethylene_groups() {
        // PE: {[]CC[]} => CCCC...CC, each repeat = 1 CH2 + 1 CH2
        // For n=5: CCCCCCCCCC = 10 carbons
        // Terminal CH3 + internal CH2s + terminal CH3
        let chain = build_chain("{[]CC[]}", 5);
        let groups = GroupDatabase::decompose(&chain).unwrap();

        let ch3_count: usize = groups
            .iter()
            .filter(|gm| gm.group.name == "-CH3")
            .map(|gm| gm.count)
            .sum();
        let ch2_count: usize = groups
            .iter()
            .filter(|gm| gm.group.name == "-CH2-")
            .map(|gm| gm.count)
            .sum();

        // 10 C atoms in a linear chain: 2 terminal CH3, 8 internal CH2
        assert_eq!(ch3_count, 2, "PE should have 2 terminal CH3 groups");
        assert_eq!(ch2_count, 8, "PE n=5 should have 8 CH2 groups");
    }

    #[test]
    fn pvc_groups() {
        // PVC: {[]C(Cl)C[]} => C(Cl)C repeated
        // Each repeat unit has 1 CHCl + 1 CH2
        let chain = build_chain("{[]C(Cl)C[]}", 3);
        let groups = GroupDatabase::decompose(&chain).unwrap();

        let cl_count: usize = groups
            .iter()
            .filter(|gm| gm.group.name == "-Cl")
            .map(|gm| gm.count)
            .sum();
        assert_eq!(cl_count, 3, "PVC n=3 should have 3 Cl groups");
    }

    #[test]
    fn sum_vw_polyethylene() {
        let chain = build_chain("{[]CC[]}", 5);
        let groups = GroupDatabase::decompose(&chain).unwrap();
        let vw = total_vw(&groups);
        // 2 * CH3(13.67) + 8 * CH2(10.23) = 27.34 + 81.84 = 109.18
        assert!((vw - 109.18).abs() < 0.1, "Vw = {vw}");
    }
}
