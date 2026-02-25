use bigsmiles::{BigSmiles, BigSmilesSegment, StochasticObject};

use crate::{error::PolySimError, polymer::PolymerChain};

use super::strategy::BuildStrategy;

/// Builder pour les architectures polymères linéaires.
pub struct LinearBuilder {
    bigsmiles: BigSmiles,
    strategy: BuildStrategy,
}

impl LinearBuilder {
    pub fn new(bigsmiles: BigSmiles, strategy: BuildStrategy) -> Self {
        Self { bigsmiles, strategy }
    }

    /// Génère un homopolymère linéaire (une seule unité de répétition).
    ///
    /// # Erreurs
    /// - [`PolySimError::NoStochasticObject`] si le BigSMILES ne contient pas d'objet stochastique.
    /// - [`PolySimError::RepeatUnitCount`] si l'objet stochastique contient ≠ 1 unité de répétition.
    /// - [`PolySimError::BuildStrategy`] si la stratégie produit n = 0.
    pub fn homopolymer(&self) -> Result<PolymerChain, PolySimError> {
        let stoch = find_first_stochastic(&self.bigsmiles)
            .ok_or(PolySimError::NoStochasticObject)?;

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
                "le nombre de répétitions doit être ≥ 1".to_string(),
            ));
        }

        let smiles = build_linear_smiles(&fragment.smiles_raw, n)?;
        Ok(PolymerChain::new(smiles, n, 0.0)) // Mn = 0.0 — calcul MW non encore implémenté
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
    pub fn block_copolymer(&self, _block_lengths: &[usize]) -> Result<PolymerChain, PolySimError> {
        todo!("implement block copolymer generation")
    }

    fn resolve_n(&self) -> Result<usize, PolySimError> {
        match &self.strategy {
            BuildStrategy::ByRepeatCount(n) => Ok(*n),
            BuildStrategy::ByTargetMn(_) | BuildStrategy::ByExactMass(_) => {
                Err(PolySimError::BuildStrategy(
                    "ByTargetMn / ByExactMass nécessitent le calcul de masse moléculaire (non encore implémenté)".to_string(),
                ))
            }
        }
    }
}

// ── Helpers internes ─────────────────────────────────────────────────────────

fn find_first_stochastic(bs: &BigSmiles) -> Option<&StochasticObject> {
    bs.segments.iter().find_map(|seg| match seg {
        BigSmilesSegment::Stochastic(obj) => Some(obj),
        _ => None,
    })
}

/// Construit le SMILES d'une chaîne linéaire de `n` unités de répétition.
///
/// Les numéros de ring closure sont renommés pour chaque copie. En SMILES,
/// un numéro de ring peut être réutilisé une fois qu'il est fermé — les copies
/// étant auto-contenues, on cycle les offsets sur 1..99, ce qui permet des
/// chaînes de longueur arbitraire.
///
/// # Erreurs
/// - [`PolySimError::RingNumberOverflow`] si le repeat unit lui-même contient > 99 ring closures
///   (ce qui est déjà un SMILES invalide).
fn build_linear_smiles(smiles_raw: &str, n: usize) -> Result<String, PolySimError> {
    let max_ring = max_ring_number(smiles_raw);

    // Cas pathologique : le repeat unit seul dépasse déjà 99 ring numbers
    if max_ring > 99 {
        return Err(PolySimError::RingNumberOverflow {
            max_ring,
            max_supported: 99,
        });
    }

    // Nombre de copies différentes avant de devoir recycler les ring numbers.
    // Comme chaque copie ferme ses propres rings avant que la suivante commence,
    // les numéros peuvent être réutilisés en toute sécurité.
    let cycle_length: usize = if max_ring == 0 {
        usize::MAX // pas de rings → pas de cycle nécessaire
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

/// Retourne le numéro de ring closure maximal utilisé dans un SMILES.
/// Les chiffres à l'intérieur de `[...]` (isotopes, charges, etc.) sont ignorés.
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
                // Notation %dd
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

/// Renumérote tous les ring closures d'un SMILES en ajoutant `offset`.
///
/// `offset = 0` retourne le SMILES inchangé.
/// Les chiffres à l'intérieur de `[...]` ne sont jamais modifiés.
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
