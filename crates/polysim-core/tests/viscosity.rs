//! Tests d'intégration pour la viscosité intrinsèque (Mark-Houwink).
//!
//! Référence : Brandrup & Immergut, *Polymer Handbook*, 4th ed., Section VII.

use bigsmiles::parse;
use polysim_core::{
    builder::{linear::LinearBuilder, BuildStrategy},
    properties::viscosity::{intrinsic_viscosity, MarkHouwinkParams},
};

// ── Helper ────────────────────────────────────────────────────────────────────

fn build_homo(bigsmiles: &str, n: usize) -> polysim_core::PolymerChain {
    let bs = parse(bigsmiles).expect("valid BigSMILES");
    LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(n))
        .homopolymer()
        .expect("build should succeed")
}

// ── Plage physique ────────────────────────────────────────────────────────────

/// [η] doit être positif pour tout polymère et toute longueur de chaîne.
#[test]
fn viscosity_result_is_positive() {
    let chain = build_homo("{[]CC[]}", 100);
    let eta = intrinsic_viscosity(&chain, &MarkHouwinkParams::pe_decalin_135c()).unwrap();
    assert!(eta > 0.0, "[η] PE doit etre positif, obtenu: {eta:.4}");
}

/// PS dans toluène, n=100 : [η] ∈ [0.1, 5.0] dL/g.
///
/// Mn(PS, n=100) ≈ 10 415 g/mol.
/// [η] = 1.16e-2 × 10415^0.73 / 100 ≈ 0.19 dL/g.
#[test]
fn viscosity_ps_toluene_in_physical_range_n100() {
    let chain = build_homo("{[]CC(c1ccccc1)[]}", 100);
    let eta = intrinsic_viscosity(&chain, &MarkHouwinkParams::ps_toluene_25c()).unwrap();
    assert!(
        eta > 0.05 && eta < 5.0,
        "PS [η](n=100) = {eta:.4} dL/g, attendu dans [0.05, 5.0]"
    );
}

/// PMMA dans acétone, n=100 : [η] ∈ [0.01, 5.0] dL/g.
///
/// Mn(PMMA, n=100) ≈ 10 012 g/mol.
/// [η] = 7.5e-3 × 10012^0.70 / 100 ≈ 0.047 dL/g.
#[test]
fn viscosity_pmma_acetone_in_physical_range_n100() {
    let chain = build_homo("{[]CC(C)(C(=O)OC)[]}", 100);
    let eta = intrinsic_viscosity(&chain, &MarkHouwinkParams::pmma_acetone_25c()).unwrap();
    assert!(
        eta > 0.01 && eta < 5.0,
        "PMMA [η](n=100) = {eta:.4} dL/g, attendu dans [0.01, 5.0]"
    );
}

// ── Monotonie ─────────────────────────────────────────────────────────────────

/// [η] croît avec Mn : n=100 < n=1000.
#[test]
fn viscosity_increases_with_mn() {
    let chain_100 = build_homo("{[]CC(c1ccccc1)[]}", 100);
    let chain_1000 = build_homo("{[]CC(c1ccccc1)[]}", 1000);
    let eta_100 = intrinsic_viscosity(&chain_100, &MarkHouwinkParams::ps_toluene_25c()).unwrap();
    let eta_1000 = intrinsic_viscosity(&chain_1000, &MarkHouwinkParams::ps_toluene_25c()).unwrap();
    assert!(
        eta_1000 > eta_100,
        "[η] PS : n=100 → {eta_100:.4}, n=1000 → {eta_1000:.4} — doit croitre"
    );
}

/// [η] PE croît avec Mn : n=50 < n=500.
#[test]
fn viscosity_pe_increases_with_mn() {
    let eta_50 = intrinsic_viscosity(
        &build_homo("{[]CC[]}", 50),
        &MarkHouwinkParams::pe_decalin_135c(),
    )
    .unwrap();
    let eta_500 = intrinsic_viscosity(
        &build_homo("{[]CC[]}", 500),
        &MarkHouwinkParams::pe_decalin_135c(),
    )
    .unwrap();
    assert!(
        eta_500 > eta_50,
        "[η] PE : n=50 → {eta_50:.4}, n=500 → {eta_500:.4} — doit croitre"
    );
}

// ── Preset Mark-Houwink ───────────────────────────────────────────────────────

/// Le preset PS/toluène donne un résultat cohérent avec la littérature.
///
/// Brandrup & Immergut (1999) Section VII : K = 1.16e-2, a = 0.73.
/// Pour Mn ≈ 10^5 g/mol → [η] ≈ 1.1 dL/g.
#[test]
fn viscosity_ps_toluene_preset_works() {
    let params = MarkHouwinkParams::ps_toluene_25c();
    assert_eq!(params.solvent, "toluene");
    assert!(
        (params.k - 1.16e-2).abs() < 1e-5,
        "K PS/toluene = {}",
        params.k
    );
    assert!(
        (params.a - 0.73).abs() < 1e-4,
        "a PS/toluene = {}",
        params.a
    );
    assert_eq!(params.temperature, 298.15);

    let chain = build_homo("{[]CC(c1ccccc1)[]}", 1000);
    let eta = intrinsic_viscosity(&chain, &params).unwrap();
    assert!(
        eta > 0.3 && eta < 2.0,
        "PS [η](n=1000) = {eta:.4} dL/g, attendu dans [0.3, 2.0]"
    );
}

/// Le preset PE/décaline donne un résultat cohérent avec la littérature.
#[test]
fn viscosity_pe_decalin_preset_works() {
    let params = MarkHouwinkParams::pe_decalin_135c();
    assert_eq!(params.solvent, "decalin");
    assert!(
        (params.k - 6.2e-2).abs() < 1e-5,
        "K PE/decalin = {}",
        params.k
    );
    assert!(
        (params.a - 0.70).abs() < 1e-4,
        "a PE/decalin = {}",
        params.a
    );

    let chain = build_homo("{[]CC[]}", 1000);
    let eta = intrinsic_viscosity(&chain, &params).unwrap();
    assert!(
        eta > 0.3 && eta < 3.0,
        "PE [η](n=1000) = {eta:.4} dL/g, attendu dans [0.3, 3.0]"
    );
}

/// Le preset PMMA/acétone donne un résultat cohérent avec la littérature.
#[test]
fn viscosity_pmma_acetone_preset_works() {
    let params = MarkHouwinkParams::pmma_acetone_25c();
    assert_eq!(params.solvent, "acetone");

    let chain = build_homo("{[]CC(C)(C(=O)OC)[]}", 1000);
    let eta = intrinsic_viscosity(&chain, &params).unwrap();
    assert!(
        eta > 0.1 && eta < 2.0,
        "PMMA [η](n=1000) = {eta:.4} dL/g, attendu dans [0.1, 2.0]"
    );
}

// ── Cas limites ───────────────────────────────────────────────────────────────

/// n=1 : ne doit pas paniquer.
#[test]
fn viscosity_n1_no_panic() {
    let chain = build_homo("{[]CC[]}", 1);
    let result = intrinsic_viscosity(&chain, &MarkHouwinkParams::pe_decalin_135c());
    assert!(result.is_ok(), "n=1 ne doit pas paniquer");
    assert!(result.unwrap() > 0.0);
}

/// K négatif retourne une erreur.
#[test]
fn viscosity_negative_k_returns_error() {
    let chain = build_homo("{[]CC[]}", 50);
    let bad_params = MarkHouwinkParams {
        k: -1.0,
        a: 0.70,
        solvent: String::new(),
        temperature: 298.15,
    };
    assert!(
        intrinsic_viscosity(&chain, &bad_params).is_err(),
        "K negatif doit retourner une erreur"
    );
}

/// Exposant `a` hors plage [0, 2] retourne une erreur.
#[test]
fn viscosity_exponent_out_of_range_returns_error() {
    let chain = build_homo("{[]CC[]}", 50);
    let bad_params = MarkHouwinkParams {
        k: 1.0e-2,
        a: 3.0,
        solvent: String::new(),
        temperature: 298.15,
    };
    assert!(
        intrinsic_viscosity(&chain, &bad_params).is_err(),
        "a=3.0 hors [0,2] doit retourner une erreur"
    );
}

// ── Cohérence entre polymères ─────────────────────────────────────────────────

/// Pour des paramètres identiques (K, a), le polymère le plus lourd a [η] plus élevé.
/// PE (28 g/mol/u) vs PS (104 g/mol/u) — même n → PS a Mn bien plus élevé → [η] plus grand.
#[test]
fn viscosity_higher_mn_polymer_has_higher_eta_same_params() {
    let params = MarkHouwinkParams {
        k: 1.0e-2,
        a: 0.70,
        solvent: "generic".into(),
        temperature: 298.15,
    };
    let pe_chain = build_homo("{[]CC[]}", 100);
    let ps_chain = build_homo("{[]CC(c1ccccc1)[]}", 100);
    let eta_pe = intrinsic_viscosity(&pe_chain, &params).unwrap();
    let eta_ps = intrinsic_viscosity(&ps_chain, &params).unwrap();
    assert!(
        eta_ps > eta_pe,
        "PS (Mn plus eleve) doit avoir [η] > PE pour memes parametres MH: eta_pe={eta_pe:.4}, eta_ps={eta_ps:.4}"
    );
}
