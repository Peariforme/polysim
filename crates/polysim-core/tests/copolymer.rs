use polysim_core::{
    builder::{linear::LinearBuilder, BuildStrategy, EnsembleBuilder},
    distribution::SchulzZimm,
    parse, PolySimError,
};

// ═══ Alternating copolymer ══════════════════════════════════════════════════

#[test]
fn alternating_ethylene_propylene_n6() {
    let bs = parse("{[$]CC[$],[$]CC(C)[$]}").unwrap();
    let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(6))
        .alternating_copolymer()
        .unwrap();
    assert_eq!(chain.repeat_count, 6);
    // Pattern: CC CC(C) CC CC(C) CC CC(C)
    assert_eq!(chain.smiles, "CCCC(C)CCCC(C)CCCC(C)");
    assert!(chain.mn > 0.0);
}

#[test]
fn alternating_3_units_cycles() {
    let bs = parse("{[$]CC[$],[$]CC(C)[$],[$]CCO[$]}").unwrap();
    let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(6))
        .alternating_copolymer()
        .unwrap();
    assert_eq!(chain.repeat_count, 6);
    // Pattern: CC CC(C) CCO CC CC(C) CCO
    assert_eq!(chain.smiles, "CCCC(C)CCOCCCC(C)CCO");
}

#[test]
fn alternating_needs_at_least_2_units() {
    let bs = parse("{[$]CC[$]}").unwrap();
    let result = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(6)).alternating_copolymer();
    assert!(matches!(
        result,
        Err(PolySimError::RepeatUnitCount {
            architecture: "alternating copolymer",
            ..
        })
    ));
}

#[test]
fn alternating_by_target_mn() {
    let bs = parse("{[$]CC[$],[$]CC(C)[$]}").unwrap();
    let chain = LinearBuilder::new(bs, BuildStrategy::ByTargetMn(500.0))
        .alternating_copolymer()
        .unwrap();
    // Should be close to 500 g/mol
    let relative_error = (chain.mn - 500.0).abs() / 500.0;
    assert!(
        relative_error < 0.15,
        "Mn = {:.1}, expected ~500, relative error = {:.3}",
        chain.mn,
        relative_error
    );
}

// ═══ Block copolymer ════════════════════════════════════════════════════════

#[test]
fn block_3_3() {
    let bs = parse("{[$]CC[$],[$]CC(C)[$]}").unwrap();
    let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(6))
        .block_copolymer(&[3, 3])
        .unwrap();
    assert_eq!(chain.repeat_count, 6);
    // Pattern: CC CC CC CC(C) CC(C) CC(C)
    assert_eq!(chain.smiles, "CCCCCCCC(C)CC(C)CC(C)");
    assert!(chain.mn > 0.0);
}

#[test]
fn block_wrong_count() {
    let bs = parse("{[$]CC[$],[$]CC(C)[$]}").unwrap();
    let result =
        LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(6)).block_copolymer(&[3, 2, 1]); // 3 blocks but only 2 units
    assert!(matches!(result, Err(PolySimError::RepeatUnitCount { .. })));
}

#[test]
fn block_with_1_unit_is_error() {
    let bs = parse("{[$]CC[$]}").unwrap();
    let result = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(6)).block_copolymer(&[6]);
    assert!(matches!(
        result,
        Err(PolySimError::RepeatUnitCount {
            architecture: "block copolymer",
            ..
        })
    ));
}

#[test]
fn block_ring_renumbering() {
    // Styrene (6 ring closures) + ethylene (0 ring closures)
    let bs = parse("{[$]CC(c1ccccc1)[$],[$]CC[$]}").unwrap();
    let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(4))
        .block_copolymer(&[2, 2])
        .unwrap();
    assert_eq!(chain.repeat_count, 4);
    // Ring numbers for the two styrene copies must differ
    let smiles = &chain.smiles;
    // First styrene: c1ccccc1, second: c2ccccc2
    assert!(smiles.contains("c1ccccc1"), "first styrene ring: {smiles}");
    assert!(smiles.contains("c2ccccc2"), "second styrene ring: {smiles}");
}

// ═══ Random copolymer ═══════════════════════════════════════════════════════

#[test]
fn random_n10_seeded() {
    let bs = parse("{[$]CC[$],[$]CC(C)[$]}").unwrap();
    let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(10))
        .seed(42)
        .random_copolymer(&[0.5, 0.5])
        .unwrap();
    assert_eq!(chain.repeat_count, 10);
    assert!(chain.mn > 0.0);
}

#[test]
fn random_seed_reproducibility() {
    let bs1 = parse("{[$]CC[$],[$]CC(C)[$]}").unwrap();
    let bs2 = parse("{[$]CC[$],[$]CC(C)[$]}").unwrap();
    let c1 = LinearBuilder::new(bs1, BuildStrategy::ByRepeatCount(20))
        .seed(42)
        .random_copolymer(&[0.6, 0.4])
        .unwrap();
    let c2 = LinearBuilder::new(bs2, BuildStrategy::ByRepeatCount(20))
        .seed(42)
        .random_copolymer(&[0.6, 0.4])
        .unwrap();
    assert_eq!(c1.smiles, c2.smiles);
}

#[test]
fn random_fractions_sum_error() {
    let bs = parse("{[$]CC[$],[$]CC(C)[$]}").unwrap();
    let result =
        LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(10)).random_copolymer(&[0.3, 0.3]);
    assert!(matches!(result, Err(PolySimError::InvalidFractions { .. })));
}

#[test]
fn random_fractions_count_mismatch() {
    let bs = parse("{[$]CC[$],[$]CC(C)[$]}").unwrap();
    let result = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(10)).random_copolymer(&[1.0]); // 1 fraction but 2 units
    assert!(matches!(result, Err(PolySimError::RepeatUnitCount { .. })));
}

#[test]
fn random_by_target_mn() {
    let bs = parse("{[$]CC[$],[$]CC(C)[$]}").unwrap();
    let chain = LinearBuilder::new(bs, BuildStrategy::ByTargetMn(1000.0))
        .seed(42)
        .random_copolymer(&[0.5, 0.5])
        .unwrap();
    let relative_error = (chain.mn - 1000.0).abs() / 1000.0;
    assert!(
        relative_error < 0.15,
        "Mn = {:.1}, expected ~1000, relative error = {:.3}",
        chain.mn,
        relative_error
    );
}

// ═══ Validation: homopolymer rejects >1 unit ════════════════════════════════

#[test]
fn homo_with_2_units_is_error() {
    let bs = parse("{[$]CC[$],[$]CC(C)[$]}").unwrap();
    let result = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(10)).homopolymer();
    assert!(matches!(
        result,
        Err(PolySimError::RepeatUnitCount {
            architecture: "homopolymer",
            ..
        })
    ));
}

// ═══ Ensemble copolymers ════════════════════════════════════════════════════

#[test]
fn random_ensemble_basic() {
    let bs = parse("{[$]CC[$],[$]CC(C)[$]}").unwrap();
    let ensemble = EnsembleBuilder::new(bs, SchulzZimm, 2000.0, 2.0)
        .num_chains(50)
        .seed(42)
        .random_copolymer_ensemble(&[0.5, 0.5])
        .unwrap();
    assert_eq!(ensemble.len(), 50);
}

#[test]
fn alternating_ensemble_basic() {
    let bs = parse("{[$]CC[$],[$]CC(C)[$]}").unwrap();
    let ensemble = EnsembleBuilder::new(bs, SchulzZimm, 2000.0, 2.0)
        .num_chains(50)
        .seed(42)
        .alternating_copolymer_ensemble()
        .unwrap();
    assert_eq!(ensemble.len(), 50);
}

#[test]
fn block_ensemble_basic() {
    let bs = parse("{[$]CC[$],[$]CC(C)[$]}").unwrap();
    let ensemble = EnsembleBuilder::new(bs, SchulzZimm, 2000.0, 2.0)
        .num_chains(50)
        .seed(42)
        .block_copolymer_ensemble(&[0.5, 0.5])
        .unwrap();
    assert_eq!(ensemble.len(), 50);
}
