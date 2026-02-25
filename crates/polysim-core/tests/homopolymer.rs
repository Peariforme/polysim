use bigsmiles::parse;
use polysim_core::{
    builder::{linear::LinearBuilder, BuildStrategy},
    error::PolySimError,
};

// ── ByRepeatCount — cas nominaux ────────────────────────────────────────────

#[test]
fn polyethylene_n1() {
    let bs = parse("{[]CC[]}").unwrap();
    let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(1))
        .homopolymer()
        .unwrap();
    assert_eq!(chain.smiles, "CC");
    assert_eq!(chain.repeat_count, 1);
}

#[test]
fn polyethylene_n3() {
    let bs = parse("{[]CC[]}").unwrap();
    let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(3))
        .homopolymer()
        .unwrap();
    assert_eq!(chain.smiles, "CCCCCC");
    assert_eq!(chain.repeat_count, 3);
}

#[test]
fn polypropylene_n2() {
    let bs = parse("{[]CC(C)[]}").unwrap();
    let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(2))
        .homopolymer()
        .unwrap();
    assert_eq!(chain.smiles, "CC(C)CC(C)");
    assert_eq!(chain.repeat_count, 2);
}

// ── Renommage des anneaux ────────────────────────────────────────────────────

#[test]
fn polystyrene_ring_renumbering_n2() {
    let bs = parse("{[]CC(c1ccccc1)[]}").unwrap();
    let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(2))
        .homopolymer()
        .unwrap();
    // Chaque copie doit utiliser un numéro de ring unique
    assert_eq!(chain.smiles, "CC(c1ccccc1)CC(c2ccccc2)");
    assert_eq!(chain.repeat_count, 2);
}

#[test]
fn polystyrene_ring_renumbering_n3() {
    let bs = parse("{[]CC(c1ccccc1)[]}").unwrap();
    let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(3))
        .homopolymer()
        .unwrap();
    assert_eq!(chain.smiles, "CC(c1ccccc1)CC(c2ccccc2)CC(c3ccccc3)");
    assert_eq!(chain.repeat_count, 3);
}

#[test]
fn bracket_atom_digits_not_renumbered() {
    // [13C] contient le digit 1 et 3, ils NE doivent PAS être traités comme des ring closures
    let bs = parse("{[][13C][13C][]}").unwrap();
    let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(2))
        .homopolymer()
        .unwrap();
    assert_eq!(chain.smiles, "[13C][13C][13C][13C]");
}

// ── Recyclage des ring numbers (SMILES autorise la réutilisation d'un ring fermé) ─

#[test]
fn polystyrene_ring_cycling_n100() {
    // max_ring = 1, cycle_length = 99
    // copie 0 → ring 1, copie 1 → ring 2, …, copie 98 → ring 99,
    // copie 99 → ring 1 à nouveau (ring 1 de la copie 0 est déjà fermé ✓)
    let bs = parse("{[]CC(c1ccccc1)[]}").unwrap();
    let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(100))
        .homopolymer()
        .unwrap();
    assert_eq!(chain.repeat_count, 100);
    // La 100e copie (index 99) doit recycler ring 1
    assert!(
        chain.smiles.ends_with("CC(c1ccccc1)"),
        "la dernière copie doit utiliser le ring 1 recyclé, smiles={}",
        &chain.smiles[chain.smiles.len().saturating_sub(40)..]
    );
}

// ── Cas d'erreur ─────────────────────────────────────────────────────────────

#[test]
fn repeat_count_zero_is_error() {
    let bs = parse("{[]CC[]}").unwrap();
    let result = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(0))
        .homopolymer();
    assert!(
        matches!(result, Err(PolySimError::BuildStrategy(_))),
        "got: {result:?}"
    );
}

#[test]
fn no_stochastic_object_is_error() {
    // BigSMILES pur SMILES — pas d'objet stochastique
    let bs = parse("CCO").unwrap();
    let result = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(3))
        .homopolymer();
    assert!(
        matches!(result, Err(PolySimError::NoStochasticObject)),
        "got: {result:?}"
    );
}

#[test]
fn multiple_repeat_units_is_error() {
    // Copolymère → pas un homopolymère
    let bs = parse("{[$]CC[$],[$]CC(C)[$]}").unwrap();
    let result = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(3))
        .homopolymer();
    assert!(
        matches!(result, Err(PolySimError::RepeatUnitCount { .. })),
        "got: {result:?}"
    );
}
