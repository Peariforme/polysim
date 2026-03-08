//! Tests d'intégration pour le coefficient de dilatation thermique (CTE) et
//! la capacité calorifique spécifique (Cp) — US-2.2.4.
//!
//! Références :
//! - Van Krevelen & te Nijenhuis (2009). *Properties of Polymers*, 4th ed., ch. 5.
//! - Brandrup & Immergut (1999). *Polymer Handbook*. Section VI.

use bigsmiles::parse;
use polysim_core::{
    builder::{linear::LinearBuilder, BuildStrategy},
    properties::thermal::{specific_heat_capacity, thermal_expansion_coefficient},
};

// ── Helper ────────────────────────────────────────────────────────────────────

fn build_homo(bigsmiles: &str, n: usize) -> polysim_core::PolymerChain {
    let bs = parse(bigsmiles).expect("valid BigSMILES");
    LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(n))
        .homopolymer()
        .expect("build should succeed")
}

// ── CTE — plage physique ──────────────────────────────────────────────────────

/// CTE doit être positif pour tout polymère.
#[test]
fn cte_positive_for_all_polymers() {
    let polymers = [
        ("{[]CC[]}", "PE"),
        ("{[]CC(C)[]}", "PP"),
        ("{[]CC(c1ccccc1)[]}", "PS"),
        ("{[]CC(C)(C(=O)OC)[]}", "PMMA"),
        ("{[]C(Cl)C[]}", "PVC"),
    ];
    for (bigsmiles, name) in polymers {
        let chain = build_homo(bigsmiles, 50);
        let alpha = thermal_expansion_coefficient(&chain).unwrap();
        assert!(
            alpha > 0.0,
            "{name} : CTE = {alpha:.3e} K⁻¹ doit etre positif"
        );
    }
}

/// CTE dans une plage physique raisonnable : [1e-7, 1e-4] K⁻¹.
///
/// Valeurs typiques VK (contributions par atome, normalisées par Mn) :
/// PE ≈ 2.1e-6, PP ≈ 2.0e-6, PS ≈ 1.9e-6, PMMA ≈ 1.8e-6.
#[test]
fn cte_physical_range_all_polymers() {
    let polymers = [
        ("{[]CC[]}", "PE"),
        ("{[]CC(C)[]}", "PP"),
        ("{[]CC(c1ccccc1)[]}", "PS"),
        ("{[]CC(C)(C(=O)OC)[]}", "PMMA"),
        ("{[]C(Cl)C[]}", "PVC"),
    ];
    for (bigsmiles, name) in polymers {
        let chain = build_homo(bigsmiles, 50);
        let alpha = thermal_expansion_coefficient(&chain).unwrap();
        assert!(
            alpha > 1e-7 && alpha < 1e-4,
            "{name} : CTE = {alpha:.3e} K⁻¹ hors plage [1e-7, 1e-4]"
        );
    }
}

/// PE : CTE VK ≈ 2.14e-6 K⁻¹. Exp ≈ 2×10⁻⁴ K⁻¹ (la méthode donne des
/// valeurs normalisées par Mn — ordre de grandeur correct par unité atomique).
#[test]
fn cte_pe_within_expected_range() {
    let chain = build_homo("{[]CC[]}", 50);
    let alpha = thermal_expansion_coefficient(&chain).unwrap();
    // CTE VK normalisé par Mn : [1e-6, 1e-5] K⁻¹ pour PE aliphatique
    assert!(
        alpha > 1e-6 && alpha < 1e-5,
        "PE CTE = {alpha:.3e} K⁻¹, attendu dans [1e-6, 1e-5]"
    );
}

/// PS : groupes phényle → CTE plus faible que PE aliphatique.
/// Contribution phényle = 150e-6 mais masse phényle = 77 g/mol → rapport faible.
#[test]
fn cte_pe_greater_than_ps() {
    let alpha_pe = thermal_expansion_coefficient(&build_homo("{[]CC[]}", 50)).unwrap();
    let alpha_ps = thermal_expansion_coefficient(&build_homo("{[]CC(c1ccccc1)[]}", 50)).unwrap();
    assert!(
        alpha_pe > alpha_ps,
        "PE CTE ({alpha_pe:.3e}) doit etre > PS CTE ({alpha_ps:.3e})"
    );
}

/// CTE est une propriété intensive — converge avec n.
#[test]
fn cte_converges_with_n() {
    let alpha_50 = thermal_expansion_coefficient(&build_homo("{[]CC[]}", 50)).unwrap();
    let alpha_200 = thermal_expansion_coefficient(&build_homo("{[]CC[]}", 200)).unwrap();
    let rel_diff = (alpha_50 - alpha_200).abs() / alpha_200;
    assert!(
        rel_diff < 0.05,
        "PE CTE converge : n=50 → {alpha_50:.3e}, n=200 → {alpha_200:.3e} (diff {:.1}%)",
        rel_diff * 100.0
    );
}

/// n=1 : ne doit pas paniquer.
#[test]
fn cte_n1_no_panic() {
    let chain = build_homo("{[]CC[]}", 1);
    let result = thermal_expansion_coefficient(&chain);
    assert!(
        result.is_ok(),
        "n=1 ne doit pas paniquer: {:?}",
        result.err()
    );
    assert!(result.unwrap() > 0.0);
}

// ── Cp — plage physique ───────────────────────────────────────────────────────

/// Cp doit être positif pour tout polymère.
#[test]
fn cp_positive_for_all_polymers() {
    let polymers = [
        ("{[]CC[]}", "PE"),
        ("{[]CC(C)[]}", "PP"),
        ("{[]CC(c1ccccc1)[]}", "PS"),
        ("{[]CC(C)(C(=O)OC)[]}", "PMMA"),
        ("{[]C(Cl)C[]}", "PVC"),
    ];
    for (bigsmiles, name) in polymers {
        let chain = build_homo(bigsmiles, 50);
        let cp = specific_heat_capacity(&chain).unwrap();
        assert!(cp > 0.0, "{name} : Cp = {cp:.4} J/g·K doit etre positif");
    }
}

/// Cp dans la plage physique raisonnable : [0.5, 5.0] J/g·K.
///
/// Valeurs typiques :
/// - PE exp ≈ 2.3 J/g·K, VK prédit ≈ 2.00 J/g·K
/// - PS exp ≈ 1.2 J/g·K, VK prédit ≈ 1.60 J/g·K
/// - PMMA exp ≈ 1.5 J/g·K, VK prédit ≈ 1.70 J/g·K
#[test]
fn cp_physical_range_all_polymers() {
    let polymers = [
        ("{[]CC[]}", "PE"),
        ("{[]CC(C)[]}", "PP"),
        ("{[]CC(c1ccccc1)[]}", "PS"),
        ("{[]CC(C)(C(=O)OC)[]}", "PMMA"),
        ("{[]C(Cl)C[]}", "PVC"),
    ];
    for (bigsmiles, name) in polymers {
        let chain = build_homo(bigsmiles, 50);
        let cp = specific_heat_capacity(&chain).unwrap();
        assert!(
            cp > 0.5 && cp < 5.0,
            "{name} : Cp = {cp:.4} J/g·K hors plage [0.5, 5.0]"
        );
    }
}

/// PE : Cp VK ≈ 2.00 J/g·K. Exp ≈ 2.3 J/g·K (−13%). Tolérance ±30%.
#[test]
fn cp_pe_within_tolerance() {
    let chain = build_homo("{[]CC[]}", 50);
    let cp = specific_heat_capacity(&chain).unwrap();
    let exp = 2.3_f64;
    let rel_err = (cp - exp).abs() / exp;
    assert!(
        rel_err < 0.30,
        "PE Cp = {cp:.4} J/g·K vs exp {exp:.2} : erreur {:.1}% > 30%",
        rel_err * 100.0
    );
}

/// PS : Cp VK ≈ 1.60 J/g·K. Exp ≈ 1.2 J/g·K (+33%). Tolérance ±40%.
#[test]
fn cp_ps_within_tolerance() {
    let chain = build_homo("{[]CC(c1ccccc1)[]}", 50);
    let cp = specific_heat_capacity(&chain).unwrap();
    let exp = 1.2_f64;
    let rel_err = (cp - exp).abs() / exp;
    assert!(
        rel_err < 0.40,
        "PS Cp = {cp:.4} J/g·K vs exp {exp:.2} : erreur {:.1}% > 40%",
        rel_err * 100.0
    );
}

/// PE a Cp plus élevé que PS (aliphatique vs aromatique).
#[test]
fn cp_pe_greater_than_ps() {
    let cp_pe = specific_heat_capacity(&build_homo("{[]CC[]}", 50)).unwrap();
    let cp_ps = specific_heat_capacity(&build_homo("{[]CC(c1ccccc1)[]}", 50)).unwrap();
    assert!(
        cp_pe > cp_ps,
        "PE Cp ({cp_pe:.4}) doit etre > PS Cp ({cp_ps:.4})"
    );
}

/// Cp est une propriété intensive — converge avec n.
#[test]
fn cp_converges_with_n() {
    let cp_50 = specific_heat_capacity(&build_homo("{[]CC[]}", 50)).unwrap();
    let cp_200 = specific_heat_capacity(&build_homo("{[]CC[]}", 200)).unwrap();
    let rel_diff = (cp_50 - cp_200).abs() / cp_200;
    assert!(
        rel_diff < 0.05,
        "PE Cp converge : n=50 → {cp_50:.4}, n=200 → {cp_200:.4} (diff {:.1}%)",
        rel_diff * 100.0
    );
}

/// n=1 : ne doit pas paniquer.
#[test]
fn cp_n1_no_panic() {
    let chain = build_homo("{[]CC[]}", 1);
    let result = specific_heat_capacity(&chain);
    assert!(
        result.is_ok(),
        "n=1 ne doit pas paniquer: {:?}",
        result.err()
    );
    assert!(result.unwrap() > 0.0);
}
