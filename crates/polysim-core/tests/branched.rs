use polysim_core::{
    builder::{branched::BranchedBuilder, BuildStrategy},
    parse,
};

// --- Comb ---

#[test]
fn comb_pe_n4_branch_every2() {
    let backbone = parse("{[]CC[]}").unwrap(); // polyethylene
    let branch = parse("{[]CC(C)[]}").unwrap(); // polypropylene
    let chain = BranchedBuilder::new(backbone, branch, BuildStrategy::ByRepeatCount(4))
        .comb_polymer(2)
        .unwrap();

    // 4 backbone units, branch every 2 => branches at positions 2 and 4
    assert!(
        chain.smiles.contains('('),
        "comb smiles must contain branches"
    );

    // 2 branches expected
    let open_parens = chain.smiles.matches('(').count();
    assert!(
        open_parens >= 2,
        "expected at least 2 branches, smiles={}",
        chain.smiles
    );

    // Architecture check
    assert!(matches!(
        chain.architecture,
        polysim_core::Architecture::Comb { branch_spacing: 2 }
    ));
}

// --- Graft ---

#[test]
fn graft_fraction_respected() {
    let backbone = parse("{[]CC[]}").unwrap();
    let branch = parse("{[]CC(C)[]}").unwrap();

    // With seed and 100 backbone units, observed fraction should be close to 0.5
    let chain = BranchedBuilder::new(backbone, branch, BuildStrategy::ByRepeatCount(100))
        .seed(42)
        .graft_copolymer(0.5, None)
        .unwrap();

    // Count branches: each branch adds '(' to the SMILES
    let branch_count = chain.smiles.matches("(CC(C))").count();
    let fraction = branch_count as f64 / 100.0;

    // Should be roughly 0.5 (within [0.3, 0.7] for 100 units)
    assert!(
        (0.3..=0.7).contains(&fraction),
        "graft fraction should be near 0.5, got {fraction} (branches={branch_count})"
    );

    assert!(matches!(
        chain.architecture,
        polysim_core::Architecture::Graft { .. }
    ));
}

// --- Star ---

#[test]
fn star_3arms_pe() {
    let backbone = parse("{[]CC[]}").unwrap();
    let branch = parse("{[]CC[]}").unwrap(); // branch not used for star
    let chain = BranchedBuilder::new(backbone, branch, BuildStrategy::ByRepeatCount(2))
        .star_polymer(3)
        .unwrap();

    // Star SMILES: C(CCCC)(CCCC)CCCC
    // Should start with C(
    assert!(
        chain.smiles.starts_with("C("),
        "star smiles should start with C(, got: {}",
        chain.smiles
    );

    // Total repeat_count = 3 arms * 2 units = 6
    assert_eq!(chain.repeat_count, 6);

    assert!(matches!(
        chain.architecture,
        polysim_core::Architecture::Star { arms: 3 }
    ));
}

#[test]
fn star_arms_out_of_range_low() {
    let backbone = parse("{[]CC[]}").unwrap();
    let branch = parse("{[]CC[]}").unwrap();
    let result =
        BranchedBuilder::new(backbone, branch, BuildStrategy::ByRepeatCount(2)).star_polymer(2);
    assert!(result.is_err());
}

#[test]
fn star_arms_out_of_range_high() {
    let backbone = parse("{[]CC[]}").unwrap();
    let branch = parse("{[]CC[]}").unwrap();
    let result =
        BranchedBuilder::new(backbone, branch, BuildStrategy::ByRepeatCount(2)).star_polymer(13);
    assert!(result.is_err());
}

// --- Dendrimer ---

#[test]
fn dendrimer_g1_bf2() {
    let backbone = parse("{[]CC[]}").unwrap();
    let branch = parse("{[]CC[]}").unwrap();
    let chain = BranchedBuilder::new(backbone, branch, BuildStrategy::ByRepeatCount(1))
        .dendrimer(1, 2)
        .unwrap();

    // G1 with bf=2: core + 2 branches = 3 units
    // SMILES: CC(CC)CC
    assert_eq!(chain.repeat_count, 3, "g1 bf2 should have 3 units");

    // Count occurrences of the repeat unit "CC"
    let cc_count = chain.smiles.matches("CC").count();
    assert!(
        cc_count >= 3,
        "expected at least 3 CC units in dendrimer g1 bf2, smiles={}",
        chain.smiles
    );

    assert!(matches!(
        chain.architecture,
        polysim_core::Architecture::Dendrimer { generation: 1 }
    ));
}

#[test]
fn dendrimer_g2_bf2() {
    let backbone = parse("{[]CC[]}").unwrap();
    let branch = parse("{[]CC[]}").unwrap();
    let chain = BranchedBuilder::new(backbone, branch, BuildStrategy::ByRepeatCount(1))
        .dendrimer(2, 2)
        .unwrap();

    // G2 with bf=2: 1 + 2 + 4 = 7 units
    assert_eq!(chain.repeat_count, 7, "g2 bf2 should have 7 units");

    assert!(matches!(
        chain.architecture,
        polysim_core::Architecture::Dendrimer { generation: 2 }
    ));
}

#[test]
fn dendrimer_generation_too_high() {
    let backbone = parse("{[]CC[]}").unwrap();
    let branch = parse("{[]CC[]}").unwrap();
    let result =
        BranchedBuilder::new(backbone, branch, BuildStrategy::ByRepeatCount(1)).dendrimer(7, 2);
    assert!(result.is_err());
}

#[test]
fn dendrimer_generation_zero() {
    let backbone = parse("{[]CC[]}").unwrap();
    let branch = parse("{[]CC[]}").unwrap();
    let result =
        BranchedBuilder::new(backbone, branch, BuildStrategy::ByRepeatCount(1)).dendrimer(0, 2);
    assert!(result.is_err());
}
