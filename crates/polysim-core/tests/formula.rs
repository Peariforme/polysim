use bigsmiles::parse;
use polysim_core::{
    builder::{linear::LinearBuilder, BuildStrategy},
    properties::formula::{molecular_formula, total_atom_count},
};

// ─── Helpers ────────────────────────────────────────────────────────────────

fn build(bigsmiles: &str, n: usize) -> polysim_core::PolymerChain {
    let bs = parse(bigsmiles).unwrap();
    LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(n))
        .homopolymer()
        .unwrap()
}

fn build_pe(n: usize) -> polysim_core::PolymerChain {
    build("{[]CC[]}", n)
}

fn build_pp(n: usize) -> polysim_core::PolymerChain {
    build("{[]CC(C)[]}", n)
}

fn build_ps(n: usize) -> polysim_core::PolymerChain {
    build("{[]CC(c1ccccc1)[]}", n)
}

// ─── molecular_formula — polyéthylène ───────────────────────────────────────

#[test]
fn formula_pe_n1() {
    // CC = éthane C₂H₆
    assert_eq!(molecular_formula(&build_pe(1)), "C2H6");
}

#[test]
fn formula_pe_n3() {
    // C₆H₁₄ = hexane
    assert_eq!(molecular_formula(&build_pe(3)), "C6H14");
}

#[test]
fn formula_pe_n10() {
    // C₂₀H₄₂ = icosane
    assert_eq!(molecular_formula(&build_pe(10)), "C20H42");
}

// ─── molecular_formula — polypropylène ──────────────────────────────────────

#[test]
fn formula_pp_n1() {
    // CC(C) = propane C₃H₈
    assert_eq!(molecular_formula(&build_pp(1)), "C3H8");
}

#[test]
fn formula_pp_n3() {
    // C₉H₂₀ : 9 C, 20 H
    assert_eq!(molecular_formula(&build_pp(3)), "C9H20");
}

// ─── molecular_formula — polystyrène ────────────────────────────────────────

#[test]
fn formula_ps_n1() {
    // CC(c1ccccc1) = éthylbenzène C₈H₁₀
    // C1(3H) + C2(2H) + c1(0H) + 5×c(1H) = 10 H, 8 C
    assert_eq!(molecular_formula(&build_ps(1)), "C8H10");
}

#[test]
fn formula_ps_n2() {
    // C₁₆H₁₈ : 2 unités styrène
    assert_eq!(molecular_formula(&build_ps(2)), "C16H18");
}

// ─── Notation Hill — C en premier ───────────────────────────────────────────

#[test]
fn formula_hill_carbon_is_first() {
    // La formule doit toujours commencer par C pour les polymères organiques
    for n in [1, 5, 10] {
        let f = molecular_formula(&build_pe(n));
        assert!(
            f.starts_with('C'),
            "Hill notation: 'C' must come first, got '{f}' for PE n={n}"
        );
        let f = molecular_formula(&build_pp(n));
        assert!(
            f.starts_with('C'),
            "Hill notation: 'C' must come first, got '{f}' for PP n={n}"
        );
    }
}

#[test]
fn formula_hill_hydrogen_is_second() {
    // 'H' doit être juste après le premier groupe 'C{digits}'
    let f = molecular_formula(&build_pe(10)); // "C20H42"
    let after_c = f.trim_start_matches(|c: char| c == 'C' || c.is_ascii_digit());
    assert!(
        after_c.starts_with('H'),
        "Hill notation: 'H' must follow 'C', got '{f}'"
    );
}

// ─── Cohérence formule ↔ masse ───────────────────────────────────────────────

#[test]
fn formula_pe_matches_mass_rule() {
    // PE : formule CₙH₂ₙ₊₂ (alcane linéaire)
    // Pour n unités répétées (–CH₂CH₂–)ₙ avec groupes terminaux H:
    // C : 2n, H : 4n+2
    for n in [1_usize, 2, 5, 10] {
        let f = molecular_formula(&build_pe(n));
        let expected = format!("C{}H{}", 2 * n, 4 * n + 2);
        assert_eq!(f, expected, "PE n={n}");
    }
}

#[test]
fn formula_pp_matches_mass_rule() {
    // PP : formule C₃ₙH₆ₙ₊₂ (polypropylène, alcane ramifié)
    for n in [1_usize, 2, 3, 5] {
        let f = molecular_formula(&build_pp(n));
        let expected = format!("C{}H{}", 3 * n, 6 * n + 2);
        assert_eq!(f, expected, "PP n={n}");
    }
}

// ─── total_atom_count — polyéthylène ────────────────────────────────────────

#[test]
fn atom_count_pe_n1() {
    // C₂H₆ → 2 + 6 = 8 atomes
    assert_eq!(total_atom_count(&build_pe(1)), 8);
}

#[test]
fn atom_count_pe_n3() {
    // C₆H₁₄ → 6 + 14 = 20 atomes
    assert_eq!(total_atom_count(&build_pe(3)), 20);
}

#[test]
fn atom_count_pe_n10() {
    // C₂₀H₄₂ → 20 + 42 = 62 atomes
    assert_eq!(total_atom_count(&build_pe(10)), 62);
}

// ─── total_atom_count — polypropylène ───────────────────────────────────────

#[test]
fn atom_count_pp_n1() {
    // C₃H₈ → 3 + 8 = 11 atomes
    assert_eq!(total_atom_count(&build_pp(1)), 11);
}

#[test]
fn atom_count_pp_n3() {
    // C₉H₂₀ → 9 + 20 = 29 atomes
    assert_eq!(total_atom_count(&build_pp(3)), 29);
}

// ─── total_atom_count — polystyrène ─────────────────────────────────────────

#[test]
fn atom_count_ps_n1() {
    // C₈H₁₀ → 8 + 10 = 18 atomes
    assert_eq!(total_atom_count(&build_ps(1)), 18);
}

// ─── Linéarité de total_atom_count ──────────────────────────────────────────

#[test]
fn atom_count_is_linear_in_n_for_pe() {
    // PE : CₙH₂ₙ₊₂ → total = 2n + (4n+2) = 6n+2
    // Donc Δatoms entre n et n+1 = 6 (constant)
    let counts: Vec<usize> = (1..=5).map(|n| total_atom_count(&build_pe(n))).collect();
    let deltas: Vec<usize> = counts.windows(2).map(|w| w[1] - w[0]).collect();
    let first_delta = deltas[0];
    for (i, &d) in deltas.iter().enumerate() {
        assert_eq!(
            d, first_delta,
            "atom count doit croître linéairement: Δ[{i}]={d} ≠ Δ[0]={first_delta}"
        );
    }
}

#[test]
fn atom_count_is_linear_in_n_for_pp() {
    // PP : C₃ₙH₆ₙ₊₂ → total = 3n + (6n+2) = 9n+2
    // Δatoms = 9 entre chaque n consécutif
    let counts: Vec<usize> = (1..=4).map(|n| total_atom_count(&build_pp(n))).collect();
    let deltas: Vec<usize> = counts.windows(2).map(|w| w[1] - w[0]).collect();
    let first_delta = deltas[0];
    for (i, &d) in deltas.iter().enumerate() {
        assert_eq!(
            d, first_delta,
            "PP atom count doit croître linéairement: Δ[{i}]={d} ≠ Δ[0]={first_delta}"
        );
    }
}

// ─── Cohérence formula ↔ total_atom_count ────────────────────────────────────

#[test]
fn atom_count_equals_sum_of_formula_counts_pe() {
    // Pour PE, on connaît la formule → on peut vérifier la cohérence
    for n in [1_usize, 5, 10] {
        let chain = build_pe(n);
        let formula = molecular_formula(&chain);
        let total = total_atom_count(&chain);

        // Formule PE : "C{2n}H{4n+2}"
        let c_count = 2 * n;
        let h_count = 4 * n + 2;
        let expected_total = c_count + h_count;

        assert_eq!(
            total, expected_total,
            "PE n={n}: formula='{formula}', total={total}, expected={expected_total}"
        );
    }
}
