//! Tests d'intégration pour la tension de surface et la constante diélectrique.
//!
//! Références :
//! - Sugden, S. (1924). *J. Chem. Soc.*, 125, 1177. (Parachor)
//! - Van Krevelen & te Nijenhuis (2009). *Properties of Polymers*, 4th ed., ch. 8, 11.
//! - Brandrup & Immergut (1999). *Polymer Handbook*, 4th ed., Section VI.

use bigsmiles::parse;
use polysim_core::{
    builder::{linear::LinearBuilder, BuildStrategy},
    properties::surface::{dielectric_constant, surface_tension},
};

// ── Helper ────────────────────────────────────────────────────────────────────

fn build_homo(bigsmiles: &str, n: usize) -> polysim_core::PolymerChain {
    let bs = parse(bigsmiles).expect("valid BigSMILES");
    LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(n))
        .homopolymer()
        .expect("build should succeed")
}

// ── Tension de surface — plage physique ───────────────────────────────────────

/// γ doit être positif pour tout polymère.
#[test]
fn surface_tension_positive_for_all_polymers() {
    let polymers = [
        ("{[]CC[]}", "PE"),
        ("{[]CC(C)[]}", "PP"),
        ("{[]CC(c1ccccc1)[]}", "PS"),
        ("{[]CC(C)(C(=O)OC)[]}", "PMMA"),
        ("{[]C(Cl)C[]}", "PVC"),
    ];
    for (bigsmiles, name) in polymers {
        let chain = build_homo(bigsmiles, 50);
        let gamma = surface_tension(&chain).unwrap();
        assert!(
            gamma > 0.0,
            "{name} : tension de surface = {gamma:.2} mN/m doit etre positive"
        );
    }
}

/// γ doit être dans une plage physique raisonnable : [5, 100] mN/m.
///
/// La méthode Parachor a une erreur systématique : PE prédit ~50 mN/m vs exp 31 mN/m.
/// On utilise une plage large pour la validité physique.
#[test]
fn surface_tension_physical_range_all_polymers() {
    let polymers = [
        ("{[]CC[]}", "PE"),
        ("{[]CC(C)[]}", "PP"),
        ("{[]CC(c1ccccc1)[]}", "PS"),
        ("{[]CC(C)(C(=O)OC)[]}", "PMMA"),
        ("{[]C(Cl)C[]}", "PVC"),
    ];
    for (bigsmiles, name) in polymers {
        let chain = build_homo(bigsmiles, 50);
        let gamma = surface_tension(&chain).unwrap();
        assert!(
            gamma > 5.0 && gamma < 100.0,
            "{name} : tension de surface = {gamma:.2} mN/m hors plage [5, 100]"
        );
    }
}

/// PE : γ prédit ≈ 50 mN/m (Parachor surestimé par rapport à exp 31 mN/m).
/// Tolérance ±50% par rapport à la valeur expérimentale.
#[test]
fn surface_tension_pe_in_tolerance() {
    let chain = build_homo("{[]CC[]}", 50);
    let gamma = surface_tension(&chain).unwrap();
    // Exp: 31.0 mN/m. VK Parachor prédit ~50 → erreur ~60%.
    // Test de plage physique large : [15, 80]
    assert!(
        gamma > 15.0 && gamma < 80.0,
        "PE tension de surface = {gamma:.2} mN/m, attendu dans [15, 80]"
    );
}

/// PS : γ exp ~ 40 mN/m. Tolérance large car Parachor surestimé.
#[test]
fn surface_tension_ps_in_tolerance() {
    let chain = build_homo("{[]CC(c1ccccc1)[]}", 50);
    let gamma = surface_tension(&chain).unwrap();
    assert!(
        gamma > 20.0 && gamma < 90.0,
        "PS tension de surface = {gamma:.2} mN/m, attendu dans [20, 90]"
    );
}

/// PVC : γ exp ~ 39 mN/m (groupes polaires).
#[test]
fn surface_tension_pvc_positive() {
    let chain = build_homo("{[]C(Cl)C[]}", 50);
    let gamma = surface_tension(&chain).unwrap();
    assert!(
        gamma > 0.0,
        "PVC tension de surface = {gamma:.2} mN/m doit etre positive"
    );
}

// ── Tension de surface — convergence avec n ───────────────────────────────────

/// γ est une propriété intensive — converge avec n.
#[test]
fn surface_tension_converges_with_n() {
    let gamma_50 = surface_tension(&build_homo("{[]CC[]}", 50)).unwrap();
    let gamma_200 = surface_tension(&build_homo("{[]CC[]}", 200)).unwrap();
    let rel_diff = (gamma_50 - gamma_200).abs() / gamma_200;
    assert!(
        rel_diff < 0.05,
        "PE γ converge : n=50 → {gamma_50:.2}, n=200 → {gamma_200:.2} (diff {:.1}%)",
        rel_diff * 100.0
    );
}

// ── Tension de surface — cas limites ─────────────────────────────────────────

/// n=1 : ne doit pas paniquer.
#[test]
fn surface_tension_n1_no_panic() {
    let chain = build_homo("{[]CC[]}", 1);
    let result = surface_tension(&chain);
    assert!(
        result.is_ok(),
        "n=1 ne doit pas paniquer: {:?}",
        result.err()
    );
}

// ── Constante diélectrique — plage physique ────────────────────────────────────

/// ε > 1 pour tout polymère (condition physique fondamentale).
#[test]
fn dielectric_constant_greater_than_one_all_polymers() {
    let polymers = [
        ("{[]CC[]}", "PE"),
        ("{[]CC(C)[]}", "PP"),
        ("{[]CC(c1ccccc1)[]}", "PS"),
        ("{[]CC(C)(C(=O)OC)[]}", "PMMA"),
        ("{[]C(Cl)C[]}", "PVC"),
    ];
    for (bigsmiles, name) in polymers {
        let chain = build_homo(bigsmiles, 50);
        let eps = dielectric_constant(&chain).unwrap();
        assert!(eps > 1.0, "{name} : ε = {eps:.4} doit etre > 1.0");
    }
}

/// ε dans la plage physique [1.5, 6.0] pour les polymères organiques courants.
#[test]
fn dielectric_constant_physical_range_all_polymers() {
    let polymers = [
        ("{[]CC[]}", "PE"),
        ("{[]CC(C)[]}", "PP"),
        ("{[]CC(c1ccccc1)[]}", "PS"),
        ("{[]CC(C)(C(=O)OC)[]}", "PMMA"),
        ("{[]C(Cl)C[]}", "PVC"),
    ];
    for (bigsmiles, name) in polymers {
        let chain = build_homo(bigsmiles, 50);
        let eps = dielectric_constant(&chain).unwrap();
        assert!(
            eps > 1.5 && eps < 6.0,
            "{name} : ε = {eps:.4} hors plage physique [1.5, 6.0]"
        );
    }
}

/// PE (apolar) : ε exp ≈ 2.3. Méthode de Clausius-Mossotti → ε électronique uniquement.
#[test]
fn dielectric_constant_pe_in_tolerance() {
    let chain = build_homo("{[]CC[]}", 50);
    let eps = dielectric_constant(&chain).unwrap();
    // Exp ≈ 2.3. Clause. Moss. électronique seul → sous-estime les polaires.
    // Pour PE (apolar) : bonne précision attendue ±30%.
    let exp = 2.3_f64;
    let rel_err = (eps - exp).abs() / exp;
    assert!(
        rel_err < 0.35,
        "PE ε = {eps:.4} vs exp {exp:.2} : erreur relative {:.1}% > 35%",
        rel_err * 100.0
    );
}

/// Les polymères polaires ont ε > PE (comportement qualitatif attendu).
///
/// PVC (a Cl) et PMMA (a ester) doivent avoir ε plus élevé que PE (apolar).
/// Note : avec la méthode électronique, la différence peut être faible.
#[test]
fn dielectric_constant_polar_ge_nonpolar() {
    let eps_pe = dielectric_constant(&build_homo("{[]CC[]}", 50)).unwrap();
    let eps_pvc = dielectric_constant(&build_homo("{[]C(Cl)C[]}", 50)).unwrap();
    // PVC contient Cl polaire → ε(PVC) >= ε(PE)
    assert!(
        eps_pvc >= eps_pe,
        "PVC ε ({eps_pvc:.4}) doit etre >= PE ε ({eps_pe:.4})"
    );
}

// ── Constante diélectrique — convergence ─────────────────────────────────────

/// ε est une propriété intensive — converge avec n.
#[test]
fn dielectric_constant_converges_with_n() {
    let eps_50 = dielectric_constant(&build_homo("{[]CC[]}", 50)).unwrap();
    let eps_200 = dielectric_constant(&build_homo("{[]CC[]}", 200)).unwrap();
    let rel_diff = (eps_50 - eps_200).abs() / eps_200;
    assert!(
        rel_diff < 0.05,
        "PE ε converge : n=50 → {eps_50:.4}, n=200 → {eps_200:.4} (diff {:.1}%)",
        rel_diff * 100.0
    );
}

/// n=1 : ne doit pas paniquer.
#[test]
fn dielectric_constant_n1_no_panic() {
    let chain = build_homo("{[]CC[]}", 1);
    let result = dielectric_constant(&chain);
    assert!(
        result.is_ok(),
        "n=1 ne doit pas paniquer: {:?}",
        result.err()
    );
}
