use polysim_core::{
    builder::{linear::LinearBuilder, BuildStrategy, GradientProfile},
    parse, Architecture, PolySimError,
};

// ═══ Gradient copolymer ═════════════════════════════════════════════════════

#[test]
fn gradient_linear_n20_seeded() {
    let bs = parse("{[$]CC[$],[$]CC(C)[$]}").unwrap();
    let profile = GradientProfile::Linear {
        f_start: 1.0,
        f_end: 0.0,
    };
    let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(20))
        .seed(42)
        .gradient_copolymer(&profile)
        .unwrap();
    assert_eq!(chain.repeat_count, 20);
    assert!(chain.mn > 0.0);
    assert_eq!(chain.architecture, Architecture::Gradient);
    // composition should have 2 entries summing to ~1.0
    assert_eq!(chain.composition.len(), 2);
    let total: f64 = chain.composition.iter().map(|m| m.fraction).sum();
    assert!((total - 1.0).abs() < 1e-9);
}

#[test]
fn gradient_sigmoid_architecture() {
    let bs = parse("{[$]CC[$],[$]CC(C)[$]}").unwrap();
    let profile = GradientProfile::Sigmoid {
        f_start: 0.9,
        f_end: 0.1,
    };
    let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(30))
        .seed(7)
        .gradient_copolymer(&profile)
        .unwrap();
    assert_eq!(chain.architecture, Architecture::Gradient);
    assert_eq!(chain.repeat_count, 30);
}

#[test]
fn gradient_needs_exactly_2_units() {
    let bs = parse("{[$]CC[$]}").unwrap();
    let profile = GradientProfile::Linear {
        f_start: 1.0,
        f_end: 0.0,
    };
    let result =
        LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(10)).gradient_copolymer(&profile);
    assert!(matches!(result, Err(PolySimError::RepeatUnitCount { .. })));
}

#[test]
fn gradient_seed_reproducibility() {
    let profile = GradientProfile::Linear {
        f_start: 0.8,
        f_end: 0.2,
    };
    let bs1 = parse("{[$]CC[$],[$]CC(C)[$]}").unwrap();
    let bs2 = parse("{[$]CC[$],[$]CC(C)[$]}").unwrap();
    let c1 = LinearBuilder::new(bs1, BuildStrategy::ByRepeatCount(15))
        .seed(99)
        .gradient_copolymer(&profile)
        .unwrap();
    let c2 = LinearBuilder::new(bs2, BuildStrategy::ByRepeatCount(15))
        .seed(99)
        .gradient_copolymer(&profile)
        .unwrap();
    assert_eq!(c1.smiles, c2.smiles);
}

// ═══ Cyclic polymer ══════════════════════════════════════════════════════════

#[test]
fn cyclic_pe_n4_has_ring_closure() {
    let bs = parse("{[]CC[]}").unwrap();
    let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(4))
        .cyclic_homopolymer()
        .unwrap();
    assert_eq!(chain.repeat_count, 4);
    assert_eq!(chain.architecture, Architecture::Cyclic);
    // SMILES must contain ring closure digit "1"
    assert!(
        chain.smiles.contains('1'),
        "cyclic smiles must have ring closure: {}",
        chain.smiles
    );
    assert!(chain.mn > 0.0);
}

#[test]
fn cyclic_architecture_field() {
    let bs = parse("{[$]CC(C)[$]}").unwrap();
    let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(3))
        .cyclic_homopolymer()
        .unwrap();
    assert_eq!(chain.architecture, Architecture::Cyclic);
}

#[test]
fn cyclic_needs_exactly_1_unit() {
    let bs = parse("{[$]CC[$],[$]CC(C)[$]}").unwrap();
    let result = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(4)).cyclic_homopolymer();
    assert!(matches!(result, Err(PolySimError::RepeatUnitCount { .. })));
}

// ═══ End groups ══════════════════════════════════════════════════════════════

#[test]
fn end_groups_prepended_and_appended() {
    // BigSMILES with prefix "CC" and suffix "CC"
    let bs = parse("CC{[$]CC[$]}CC").unwrap();
    let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(3))
        .homopolymer()
        .unwrap();
    assert!(
        chain.smiles.starts_with("CC"),
        "should start with end group: {}",
        chain.smiles
    );
    assert!(
        chain.smiles.ends_with("CC"),
        "should end with end group: {}",
        chain.smiles
    );
}

#[test]
fn end_groups_included_in_mn() {
    let bs_no_eg = parse("{[$]CC[$]}").unwrap();
    let bs_with_eg = parse("CC{[$]CC[$]}CC").unwrap();
    let chain_no = LinearBuilder::new(bs_no_eg, BuildStrategy::ByRepeatCount(5))
        .homopolymer()
        .unwrap();
    let chain_with = LinearBuilder::new(bs_with_eg, BuildStrategy::ByRepeatCount(5))
        .homopolymer()
        .unwrap();
    // With end groups, Mn must be higher
    assert!(
        chain_with.mn > chain_no.mn,
        "Mn with end groups ({:.2}) should exceed Mn without ({:.2})",
        chain_with.mn,
        chain_no.mn
    );
}
