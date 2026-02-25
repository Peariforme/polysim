use bigsmiles::parse;
use polysim_core::{
    builder::{linear::LinearBuilder, BuildStrategy},
    properties::molecular_weight::{average_mass, monoisotopic_mass},
};

// ─── Helpers ────────────────────────────────────────────────────────────────

fn build_pe(n: usize) -> polysim_core::PolymerChain {
    let bs = parse("{[]CC[]}").unwrap();
    LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(n))
        .homopolymer()
        .unwrap()
}

fn build_pp(n: usize) -> polysim_core::PolymerChain {
    let bs = parse("{[]CC(C)[]}").unwrap();
    LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(n))
        .homopolymer()
        .unwrap()
}

fn build_ps(n: usize) -> polysim_core::PolymerChain {
    let bs = parse("{[]CC(c1ccccc1)[]}").unwrap();
    LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(n))
        .homopolymer()
        .unwrap()
}

/// Vérifie que `got` est dans ±`tol` de `expected`.
fn assert_close(got: f64, expected: f64, tol: f64, label: &str) {
    assert!(
        (got - expected).abs() < tol,
        "{label}: got {got:.4}, expected {expected:.4} (±{tol})"
    );
}

// ─── average_mass — polyéthylène ────────────────────────────────────────────

#[test]
fn average_mass_pe_n1() {
    // CC = éthane, C₂H₆ = 2×12.011 + 6×1.008 = 30.070 g/mol
    assert_close(average_mass(&build_pe(1)), 30.070, 0.01, "PE n=1");
}

#[test]
fn average_mass_pe_n3() {
    // C₆H₁₄ = hexane = 86.175 g/mol
    assert_close(average_mass(&build_pe(3)), 86.175, 0.01, "PE n=3");
}

#[test]
fn average_mass_pe_n10() {
    // C₂₀H₄₂ = icosane = 282.554 g/mol
    assert_close(average_mass(&build_pe(10)), 282.554, 0.01, "PE n=10");
}

// ─── average_mass — polypropylène ───────────────────────────────────────────

#[test]
fn average_mass_pp_n1() {
    // CC(C) = propane, C₃H₈ = 3×12.011 + 8×1.008 = 44.097 g/mol
    assert_close(average_mass(&build_pp(1)), 44.097, 0.01, "PP n=1");
}

#[test]
fn average_mass_pp_n3() {
    // C₉H₂₀ = 9×12.011 + 20×1.008 = 128.255 g/mol
    assert_close(average_mass(&build_pp(3)), 128.255, 0.01, "PP n=3");
}

// ─── average_mass — polystyrène ─────────────────────────────────────────────

#[test]
fn average_mass_ps_n1() {
    // CC(c1ccccc1) = éthylbenzène (sans -CH₃ terminal : styrène hydrogéné)
    // C₈H₁₀ = 8×12.011 + 10×1.008 = 96.088 + 10.080 = 106.168 g/mol
    assert_close(average_mass(&build_ps(1)), 106.168, 0.01, "PS n=1");
}

// ─── average_mass est linéaire en n ─────────────────────────────────────────

#[test]
fn average_mass_is_linear_in_n() {
    // MW(n) doit être linéaire : MW(3) - MW(2) ≈ MW(2) - MW(1)
    let mw1 = average_mass(&build_pe(1));
    let mw2 = average_mass(&build_pe(2));
    let mw3 = average_mass(&build_pe(3));
    let delta12 = mw2 - mw1;
    let delta23 = mw3 - mw2;
    assert_close(delta12, delta23, 0.001, "linéarité PE");
}

// ─── monoisotopic_mass ───────────────────────────────────────────────────────

#[test]
fn monoisotopic_mass_pe_n1() {
    // C₂H₆ monoisotopique : 2×12.0 + 6×1.00782503207 = 30.047 g/mol
    assert_close(monoisotopic_mass(&build_pe(1)), 30.047, 0.01, "PE mono n=1");
}

#[test]
fn monoisotopic_mass_pe_n10() {
    // C₂₀H₄₂ mono : 20×12.0 + 42×1.00782503207 = 282.329 g/mol
    assert_close(
        monoisotopic_mass(&build_pe(10)),
        282.329,
        0.01,
        "PE mono n=10",
    );
}

#[test]
fn monoisotopic_mass_less_than_average() {
    // Masse monoisotopique < masse moyenne pour tous les éléments lourds organiques
    let chain = build_pe(10);
    assert!(
        monoisotopic_mass(&chain) < average_mass(&chain),
        "masse monoisotopique doit être < masse moyenne"
    );
}

// ─── chain.mn renseigné à la construction ───────────────────────────────────

#[test]
fn mn_populated_for_repeat_count() {
    let chain = build_pe(10);
    assert!(
        chain.mn > 0.0,
        "chain.mn doit être renseigné après construction, got {}",
        chain.mn
    );
    assert_close(
        chain.mn,
        average_mass(&chain),
        1e-9,
        "chain.mn == average_mass",
    );
}

#[test]
fn mn_populated_for_pp() {
    let chain = build_pp(5);
    assert!(chain.mn > 0.0);
    assert_close(chain.mn, average_mass(&chain), 1e-9, "PP chain.mn");
}

// ─── BuildStrategy::ByTargetMn ──────────────────────────────────────────────

#[test]
fn by_target_mn_pe_n10() {
    // Cible ≈ MW de PE avec 10 unités répétées
    let bs = parse("{[]CC[]}").unwrap();
    let chain = LinearBuilder::new(bs, BuildStrategy::ByTargetMn(282.554))
        .homopolymer()
        .unwrap();
    assert_eq!(chain.repeat_count, 10, "ByTargetMn doit donner n=10");
    assert_close(chain.mn, 282.554, 1.0, "MW de la chaîne construite");
}

#[test]
fn by_target_mn_pe_n1() {
    let bs = parse("{[]CC[]}").unwrap();
    let chain = LinearBuilder::new(bs, BuildStrategy::ByTargetMn(30.070))
        .homopolymer()
        .unwrap();
    assert_eq!(chain.repeat_count, 1);
}

#[test]
fn by_target_mn_pp_n5() {
    // PP n=5 : C₁₅H₃₂ = 15×12.011 + 32×1.008 = 212.421 g/mol
    let bs = parse("{[]CC(C)[]}").unwrap();
    let chain = LinearBuilder::new(bs, BuildStrategy::ByTargetMn(212.421))
        .homopolymer()
        .unwrap();
    assert_eq!(chain.repeat_count, 5);
}

#[test]
fn by_target_mn_rounds_to_nearest() {
    // PE: MW(n=1)=30.070 (éthane), MW(n=2)=58.122 (butane), mw_per_unit≈28.052
    // Midpoint entre n=1 et n=2 ≈ 44.1 → en dessous = n=1, au dessus = n=2
    let bs = parse("{[]CC[]}").unwrap();
    // Cible entre n=1 et midpoint → doit donner n=1
    let chain = LinearBuilder::new(bs.clone(), BuildStrategy::ByTargetMn(35.0))
        .homopolymer()
        .unwrap();
    assert_eq!(chain.repeat_count, 1);
    // Cible au-dessus du midpoint → doit donner n=2
    let chain2 = LinearBuilder::new(bs, BuildStrategy::ByTargetMn(50.0))
        .homopolymer()
        .unwrap();
    assert_eq!(chain2.repeat_count, 2);
}

// ─── BuildStrategy::ByExactMass ─────────────────────────────────────────────

#[test]
fn by_exact_mass_pe_n10() {
    // C₂₀H₄₂ monoisotopique ≈ 282.329 g/mol
    let bs = parse("{[]CC[]}").unwrap();
    let chain = LinearBuilder::new(bs, BuildStrategy::ByExactMass(282.329))
        .homopolymer()
        .unwrap();
    assert_eq!(chain.repeat_count, 10);
}

#[test]
fn by_exact_mass_pe_n1() {
    // C₂H₆ monoisotopique ≈ 30.047 g/mol
    let bs = parse("{[]CC[]}").unwrap();
    let chain = LinearBuilder::new(bs, BuildStrategy::ByExactMass(30.047))
        .homopolymer()
        .unwrap();
    assert_eq!(chain.repeat_count, 1);
}
