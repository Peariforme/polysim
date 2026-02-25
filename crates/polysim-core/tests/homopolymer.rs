use bigsmiles::parse;
use polysim_core::{
    builder::{linear::LinearBuilder, BuildStrategy},
    error::PolySimError,
};

// ── ByRepeatCount — nominal cases ────────────────────────────────────────────

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

// ── Ring renumbering ─────────────────────────────────────────────────────────

#[test]
fn polystyrene_ring_renumbering_n2() {
    let bs = parse("{[]CC(c1ccccc1)[]}").unwrap();
    let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(2))
        .homopolymer()
        .unwrap();
    // Each copy must use a unique ring number
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
    // [13C] contains digits 1 and 3, which must NOT be treated as ring closures
    let bs = parse("{[][13C][13C][]}").unwrap();
    let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(2))
        .homopolymer()
        .unwrap();
    assert_eq!(chain.smiles, "[13C][13C][13C][13C]");
}

// ── Ring number cycling (SMILES allows reuse of a closed ring number) ────────

#[test]
fn polystyrene_ring_cycling_n100() {
    // max_ring = 1, cycle_length = 99
    // copy 0  → ring 1, copy 1 → ring 2, …, copy 98 → ring 99,
    // copy 99 → ring 1 again (ring 1 from copy 0 is already closed ✓)
    let bs = parse("{[]CC(c1ccccc1)[]}").unwrap();
    let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(100))
        .homopolymer()
        .unwrap();
    assert_eq!(chain.repeat_count, 100);
    // Copy 99 (0-indexed) must recycle ring 1
    assert!(
        chain.smiles.ends_with("CC(c1ccccc1)"),
        "last copy must use recycled ring 1, tail={}",
        &chain.smiles[chain.smiles.len().saturating_sub(40)..]
    );
}

// ── Error cases ───────────────────────────────────────────────────────────────

#[test]
fn repeat_count_zero_is_error() {
    let bs = parse("{[]CC[]}").unwrap();
    let result = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(0)).homopolymer();
    assert!(
        matches!(result, Err(PolySimError::BuildStrategy(_))),
        "got: {result:?}"
    );
}

#[test]
fn no_stochastic_object_is_error() {
    // Plain SMILES — no stochastic object
    let bs = parse("CCO").unwrap();
    let result = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(3)).homopolymer();
    assert!(
        matches!(result, Err(PolySimError::NoStochasticObject)),
        "got: {result:?}"
    );
}

#[test]
fn multiple_repeat_units_is_error() {
    // Copolymer → not a homopolymer
    let bs = parse("{[$]CC[$],[$]CC(C)[$]}").unwrap();
    let result = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(3)).homopolymer();
    assert!(
        matches!(result, Err(PolySimError::RepeatUnitCount { .. })),
        "got: {result:?}"
    );
}
