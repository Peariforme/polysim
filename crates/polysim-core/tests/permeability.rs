//! Tests d'intégration pour la perméabilité aux gaz (méthode Permachor) — US-2.5.1.
//!
//! La méthode Permachor donne log₁₀(P) = Σ τᵢ · countᵢ sur toute la chaîne.
//! L'unité de P est le Barrer (1 Barrer = 10⁻¹⁰ cm³(STP)·cm/cm²·s·cmHg).
//!
//! Référence : Van Krevelen & te Nijenhuis (2009), ch. 22.

use bigsmiles::parse;
use polysim_core::{
    builder::{linear::LinearBuilder, BuildStrategy},
    properties::permeability::{gas_permeability, Gas},
};

// ── Helper ────────────────────────────────────────────────────────────────────

fn build_homo(bigsmiles: &str, n: usize) -> polysim_core::PolymerChain {
    let bs = parse(bigsmiles).expect("valid BigSMILES");
    LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(n))
        .homopolymer()
        .expect("build should succeed")
}

// ── Valeurs positives ─────────────────────────────────────────────────────────

/// Tous les gaz retournent une valeur positive pour PE.
#[test]
fn all_gases_return_positive_for_pe() {
    let chain = build_homo("{[]CC[]}", 1);
    for gas in [Gas::O2, Gas::CO2, Gas::N2, Gas::He, Gas::H2] {
        let p = gas_permeability(&chain, gas).unwrap();
        assert!(
            p > 0.0,
            "{gas:?} : perméabilité = {p:.4e} doit etre positive"
        );
    }
}

/// Tous les gaz retournent Ok (pas d'erreur) pour tous les polymères courants.
#[test]
fn all_gases_return_ok_for_common_polymers() {
    let polymers = [
        ("{[]CC[]}", "PE"),
        ("{[]CC(C)[]}", "PP"),
        ("{[]CC(c1ccccc1)[]}", "PS"),
        ("{[]CC(C)(C(=O)OC)[]}", "PMMA"),
        ("{[]C(Cl)C[]}", "PVC"),
    ];
    for (bigsmiles, name) in polymers {
        let chain = build_homo(bigsmiles, 1);
        for gas in [Gas::O2, Gas::CO2, Gas::N2, Gas::He, Gas::H2] {
            assert!(
                gas_permeability(&chain, gas).is_ok(),
                "{name} {gas:?} doit retourner Ok"
            );
        }
    }
}

// ── Relations qualitatives ────────────────────────────────────────────────────

/// PE est plus perméable à O₂ que PS (PE est plus flexible, moins encombré).
///
/// n=1 : PE O₂ ≈ 4.17 Barrer vs PS O₂ ≈ 3.89 Barrer.
#[test]
fn pe_more_permeable_to_o2_than_ps() {
    let pe = build_homo("{[]CC[]}", 1);
    let ps = build_homo("{[]CC(c1ccccc1)[]}", 1);
    let p_pe = gas_permeability(&pe, Gas::O2).unwrap();
    let p_ps = gas_permeability(&ps, Gas::O2).unwrap();
    assert!(
        p_pe > p_ps,
        "PE O₂ ({p_pe:.4}) doit etre > PS O₂ ({p_ps:.4})"
    );
}

/// CO₂ > O₂ pour PE (règle générale : CO₂ est plus soluble que O₂).
///
/// n=1 : PE CO₂ ≈ 17.4 Barrer vs PE O₂ ≈ 4.17 Barrer.
#[test]
fn co2_more_permeable_than_o2_in_pe() {
    let pe = build_homo("{[]CC[]}", 1);
    let p_o2 = gas_permeability(&pe, Gas::O2).unwrap();
    let p_co2 = gas_permeability(&pe, Gas::CO2).unwrap();
    assert!(
        p_co2 > p_o2,
        "PE CO₂ ({p_co2:.4}) doit etre > PE O₂ ({p_o2:.4})"
    );
}

/// PVC est bien moins perméable à O₂ que PE (groupes Cl bloquants).
///
/// n=1 : PVC O₂ ≈ 0.55 vs PE O₂ ≈ 4.17 Barrer.
#[test]
fn pvc_less_permeable_to_o2_than_pe() {
    let pe = build_homo("{[]CC[]}", 1);
    let pvc = build_homo("{[]C(Cl)C[]}", 1);
    let p_pe = gas_permeability(&pe, Gas::O2).unwrap();
    let p_pvc = gas_permeability(&pvc, Gas::O2).unwrap();
    assert!(
        p_pvc < p_pe,
        "PVC O₂ ({p_pvc:.4}) doit etre < PE O₂ ({p_pe:.4})"
    );
}

/// N₂ < O₂ pour PE (N₂ est toujours moins perméable que O₂).
#[test]
fn n2_less_permeable_than_o2_in_pe() {
    let pe = build_homo("{[]CC[]}", 1);
    let p_o2 = gas_permeability(&pe, Gas::O2).unwrap();
    let p_n2 = gas_permeability(&pe, Gas::N2).unwrap();
    assert!(
        p_n2 < p_o2,
        "PE N₂ ({p_n2:.4}) doit etre < PE O₂ ({p_o2:.4})"
    );
}

/// He > O₂ > N₂ pour PE (ordre de perméabilité classique pour les gaz légers).
#[test]
fn pe_permeability_order_he_gt_o2_gt_n2() {
    let pe = build_homo("{[]CC[]}", 1);
    let p_he = gas_permeability(&pe, Gas::He).unwrap();
    let p_o2 = gas_permeability(&pe, Gas::O2).unwrap();
    let p_n2 = gas_permeability(&pe, Gas::N2).unwrap();
    assert!(
        p_he > p_o2 && p_o2 > p_n2,
        "Ordre attendu He > O₂ > N₂ pour PE : He={p_he:.4}, O₂={p_o2:.4}, N₂={p_n2:.4}"
    );
}

// ── Unité Barrer ──────────────────────────────────────────────────────────────

/// Pour n=1, la perméabilité O₂ de PE est dans une plage de Barrer réaliste.
///
/// n=1 : PE O₂ ≈ 4.17 Barrer. Exp PE amorphe ≈ 2–7 Barrer.
#[test]
fn pe_o2_permeability_realistic_for_n1() {
    let pe = build_homo("{[]CC[]}", 1);
    let p = gas_permeability(&pe, Gas::O2).unwrap();
    assert!(
        p > 0.5 && p < 50.0,
        "PE O₂ (n=1) = {p:.4} Barrer, attendu dans [0.5, 50]"
    );
}

/// La perméabilité croît exponentiellement avec n (somme log sur la chaîne entière).
#[test]
fn permeability_increases_exponentially_with_n() {
    let p1 = gas_permeability(&build_homo("{[]CC[]}", 1), Gas::O2).unwrap();
    let p5 = gas_permeability(&build_homo("{[]CC[]}", 5), Gas::O2).unwrap();
    // log10(P5) = 5 × log10(P1) → P5 = P1^5
    assert!(
        p5 > p1,
        "Perméabilité doit croitre avec n : n=1 → {p1:.4}, n=5 → {p5:.4}"
    );
}

// ── Cas limites ───────────────────────────────────────────────────────────────

/// n=1 : ne doit pas paniquer.
#[test]
fn permeability_n1_no_panic() {
    let chain = build_homo("{[]CC[]}", 1);
    let result = gas_permeability(&chain, Gas::O2);
    assert!(
        result.is_ok(),
        "n=1 ne doit pas paniquer: {:?}",
        result.err()
    );
}

/// H₂O : PE et PVC retournent tous deux une valeur positive.
/// Note : τ(H₂O) est négatif pour -CH2-, -CH3, et -Cl → les deux polymères ont P < 1 Barrer.
#[test]
fn h2o_permeability_positive_for_pe_and_pvc() {
    let pe = build_homo("{[]CC[]}", 1);
    let pvc = build_homo("{[]C(Cl)C[]}", 1);
    let p_pe = gas_permeability(&pe, Gas::H2O).unwrap();
    let p_pvc = gas_permeability(&pvc, Gas::H2O).unwrap();
    assert!(p_pe > 0.0, "PE H₂O ({p_pe:.4}) doit etre positif");
    assert!(p_pvc > 0.0, "PVC H₂O ({p_pvc:.4}) doit etre positif");
}
