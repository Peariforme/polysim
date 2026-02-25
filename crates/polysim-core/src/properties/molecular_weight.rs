use opensmiles::{parse as parse_smiles, AtomSymbol};

use crate::polymer::PolymerChain;

/// Masse standard de l'hydrogène (IUPAC 2021), en g/mol.
const H_AVERAGE_MASS: f64 = 1.008;

/// Masse du proton (¹H), en g/mol.
const H_MONO_MASS: f64 = 1.00782503207;

/// Calcule la masse moléculaire moyenne (poids atomiques IUPAC) de la chaîne, en g/mol.
///
/// Chaque atome lourd contribue par sa masse standard (moyenne isotopique), et les
/// hydrogènes implicites/explicites sont ajoutés avec la masse standard de l'hydrogène.
///
/// # Exemple
///
/// ```rust
/// use polysim_core::{parse, builder::{linear::LinearBuilder, BuildStrategy},
///                    properties::molecular_weight::average_mass};
///
/// let bs = parse("{[]CC[]}").unwrap();
/// let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(1))
///     .homopolymer()
///     .unwrap();
/// // CC = éthane C₂H₆ ≈ 30.07 g/mol
/// let mw = average_mass(&chain);
/// assert!((mw - 30.070).abs() < 0.01, "got {mw}");
/// ```
pub fn average_mass(chain: &PolymerChain) -> f64 {
    let mol = parse_smiles(&chain.smiles).expect("chain SMILES must be valid SMILES");
    mol.nodes().iter().fold(0.0, |acc, node| {
        // atom.mass() renvoie la masse standard (ou la masse isotopique si explicite [¹³C])
        acc + node.atom().mass() + node.hydrogens() as f64 * H_AVERAGE_MASS
    })
}

/// Calcule la masse monoisotopique de la chaîne (nucléide le plus abondant), en g/mol.
///
/// Pour les atomes sans isotope explicite, utilise le nucléide le plus abondant de chaque
/// élément (ex. ¹²C = 12.000, ¹⁶O = 15.9949…). Pour les atomes avec isotope explicite
/// (`[13C]`), respecte l'isotope spécifié.
///
/// # Exemple
///
/// ```rust
/// use polysim_core::{parse, builder::{linear::LinearBuilder, BuildStrategy},
///                    properties::molecular_weight::monoisotopic_mass};
///
/// let bs = parse("{[]CC[]}").unwrap();
/// let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(1))
///     .homopolymer()
///     .unwrap();
/// // CC = éthane C₂H₆, masse monoisotopique ≈ 30.047 g/mol
/// let m = monoisotopic_mass(&chain);
/// assert!((m - 30.047).abs() < 0.01, "got {m}");
/// ```
pub fn monoisotopic_mass(chain: &PolymerChain) -> f64 {
    let mol = parse_smiles(&chain.smiles).expect("chain SMILES must be valid SMILES");
    mol.nodes().iter().fold(0.0, |acc, node| {
        let atom = node.atom();
        let heavy_mass = if atom.isotope().is_some() {
            // Isotope explicitement spécifié → respecter (ex. [13C])
            atom.mass()
        } else {
            most_abundant_isotope_mass(atom.element())
        };
        acc + heavy_mass + node.hydrogens() as f64 * H_MONO_MASS
    })
}

/// Retourne la masse du nucléide le plus abondant pour chaque élément.
///
/// Pour les éléments organiques courants en chimie des polymères, les valeurs exactes
/// sont codées en dur. Pour les éléments rares, la masse standard IUPAC est utilisée
/// comme approximation.
fn most_abundant_isotope_mass(element: &AtomSymbol) -> f64 {
    match element.atomic_number() {
        0 => 0.0,                     // Wildcard (*)
        1 => H_MONO_MASS,             // ¹H (99.985 %)
        5 => 11.0093054,              // ¹¹B (80.1 %)
        6 => 12.0,                    // ¹²C (98.89 %)
        7 => 14.0030740048,           // ¹⁴N (99.63 %)
        8 => 15.9949146221,           // ¹⁶O (99.76 %)
        9 => 18.9984032,              // ¹⁹F (100 %)
        14 => 27.9769265325,          // ²⁸Si (92.23 %)
        15 => 30.97376163,            // ³¹P (100 %)
        16 => 31.97207100,            // ³²S (95.02 %)
        17 => 34.96885268,            // ³⁵Cl (75.77 %)
        35 => 78.9183371,             // ⁷⁹Br (50.69 %)
        53 => 126.904468,             // ¹²⁷I (100 %)
        _ => element.standard_mass(), // fallback : masse IUPAC pour éléments rares
    }
}
