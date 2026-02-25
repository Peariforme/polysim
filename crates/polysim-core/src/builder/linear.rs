use bigsmiles::{BigSmiles, BigSmilesSegment, StochasticObject};

use crate::{error::PolySimError, polymer::PolymerChain};

use super::strategy::BuildStrategy;

/// Builder for linear polymer architectures.
///
/// Supports homopolymers, random/alternating/block copolymers — all derived
/// from a single BigSMILES string.
pub struct LinearBuilder {
    bigsmiles: BigSmiles,
    strategy: BuildStrategy,
}

impl LinearBuilder {
    /// Creates a new builder from a parsed BigSMILES and a build strategy.
    pub fn new(bigsmiles: BigSmiles, strategy: BuildStrategy) -> Self {
        Self {
            bigsmiles,
            strategy,
        }
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
        let stoch =
            find_first_stochastic(&self.bigsmiles).ok_or(PolySimError::NoStochasticObject)?;

        if stoch.repeat_units.len() != 1 {
            return Err(PolySimError::RepeatUnitCount {
                architecture: "homopolymer",
                got: stoch.repeat_units.len(),
                need: 1,
            });
        }

        let fragment = &stoch.repeat_units[0];
        let n = self.resolve_n()?;

        if n == 0 {
            return Err(PolySimError::BuildStrategy(
                "repeat count must be ≥ 1".to_string(),
            ));
        }

        let smiles = build_linear_smiles(&fragment.smiles_raw, n)?;
        Ok(PolymerChain::new(smiles, n, 0.0)) // Mn = 0.0 — MW calculation not yet implemented
    }

    /// Generates a random (statistical) copolymer.
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

    /// Generates an alternating copolymer (–A–B–A–B–).
    ///
    /// The BigSMILES must contain exactly 2 repeat units.
    pub fn alternating_copolymer(&self) -> Result<PolymerChain, PolySimError> {
        todo!("implement alternating copolymer generation")
    }

    /// Generates a block copolymer (–AAAA–BBBB–).
    ///
    /// `block_lengths` — number of repeat units per block, in order.
    /// The BigSMILES must contain exactly `block_lengths.len()` repeat units.
    pub fn block_copolymer(&self, _block_lengths: &[usize]) -> Result<PolymerChain, PolySimError> {
        todo!("implement block copolymer generation")
    }

    fn resolve_n(&self) -> Result<usize, PolySimError> {
        match &self.strategy {
            BuildStrategy::ByRepeatCount(n) => Ok(*n),
            BuildStrategy::ByTargetMn(_) | BuildStrategy::ByExactMass(_) => {
                Err(PolySimError::BuildStrategy(
                    "ByTargetMn / ByExactMass require molecular weight calculation \
                     (not yet implemented)"
                        .to_string(),
                ))
            }
        }
    }
}

// --- internal helpers -------------------------------------------------------

fn find_first_stochastic(bs: &BigSmiles) -> Option<&StochasticObject> {
    bs.segments.iter().find_map(|seg| match seg {
        BigSmilesSegment::Stochastic(obj) => Some(obj),
        _ => None,
    })
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
fn build_linear_smiles(smiles_raw: &str, n: usize) -> Result<String, PolySimError> {
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

/// Returns the highest ring-closure number used in a SMILES string.
///
/// Digits inside `[...]` (isotopes, hydrogen counts, charges, atom classes)
/// are ignored.
fn max_ring_number(smiles: &str) -> u32 {
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
fn renumber_ring_closures(smiles: &str, offset: u32) -> String {
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
