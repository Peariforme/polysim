use std::collections::BTreeMap;

use opensmiles::parse as parse_smiles;

use crate::polymer::PolymerChain;

/// Calcule la formule moléculaire brute d'une chaîne en notation Hill.
///
/// La notation Hill place **C** en premier, puis **H**, puis les autres éléments
/// par ordre alphabétique du symbole. Les hydrogènes implicites sont inclus.
///
/// # Exemple
///
/// ```rust
/// use polysim_core::{parse, builder::{linear::LinearBuilder, BuildStrategy},
///                    properties::formula::molecular_formula};
///
/// let bs = parse("{[]CC[]}").unwrap();
/// let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(10))
///     .homopolymer()
///     .unwrap();
/// // Polyéthylène n=10 → C₂₀H₄₂
/// assert_eq!(molecular_formula(&chain), "C20H42");
/// ```
pub fn molecular_formula(chain: &PolymerChain) -> String {
    let mol = parse_smiles(&chain.smiles).expect("chain SMILES must be valid SMILES");
    let mut counts: BTreeMap<&'static str, usize> = BTreeMap::new();

    for node in mol.nodes() {
        let atomic_num = node.atom().element().atomic_number();
        if atomic_num == 0 {
            continue; // wildcard (*)
        }
        if let Some(sym) = element_symbol(atomic_num) {
            *counts.entry(sym).or_insert(0) += 1;
        }
        let h = node.hydrogens() as usize;
        if h > 0 {
            *counts.entry("H").or_insert(0) += h;
        }
    }

    hill_notation(&counts)
}

/// Nombre total d'atomes dans la chaîne (atomes lourds + hydrogènes implicites/explicites).
///
/// # Exemple
///
/// ```rust
/// use polysim_core::{parse, builder::{linear::LinearBuilder, BuildStrategy},
///                    properties::formula::total_atom_count};
///
/// let bs = parse("{[]CC[]}").unwrap();
/// let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(1))
///     .homopolymer()
///     .unwrap();
/// // CC = éthane C₂H₆ → 2 + 6 = 8 atomes
/// assert_eq!(total_atom_count(&chain), 8);
/// ```
pub fn total_atom_count(chain: &PolymerChain) -> usize {
    let mol = parse_smiles(&chain.smiles).expect("chain SMILES must be valid SMILES");
    mol.nodes()
        .iter()
        .map(|node| 1 + node.hydrogens() as usize)
        .sum()
}

/// Formate les counts en notation Hill : C en premier, H en second,
/// puis les autres éléments par ordre alphabétique de symbole.
fn hill_notation(counts: &BTreeMap<&'static str, usize>) -> String {
    let mut result = String::new();
    let has_carbon = counts.contains_key("C");

    if has_carbon {
        // C et H en premier
        for sym in ["C", "H"] {
            if let Some(&n) = counts.get(sym) {
                result.push_str(sym);
                if n > 1 {
                    result.push_str(&n.to_string());
                }
            }
        }
        // Reste par ordre alphabétique (BTreeMap est déjà trié)
        for (&sym, &n) in counts {
            if sym == "C" || sym == "H" {
                continue;
            }
            result.push_str(sym);
            if n > 1 {
                result.push_str(&n.to_string());
            }
        }
    } else {
        // Pas de carbone → tout par ordre alphabétique
        for (&sym, &n) in counts {
            result.push_str(sym);
            if n > 1 {
                result.push_str(&n.to_string());
            }
        }
    }
    result
}

/// Retourne le symbole IUPAC de l'élément pour le numéro atomique donné.
///
/// Couvre les éléments courants en chimie des polymères.
/// Retourne `None` pour les éléments inconnus ou rares.
fn element_symbol(atomic_number: u8) -> Option<&'static str> {
    match atomic_number {
        1 => Some("H"),
        5 => Some("B"),
        6 => Some("C"),
        7 => Some("N"),
        8 => Some("O"),
        9 => Some("F"),
        14 => Some("Si"),
        15 => Some("P"),
        16 => Some("S"),
        17 => Some("Cl"),
        35 => Some("Br"),
        53 => Some("I"),
        _ => None,
    }
}
