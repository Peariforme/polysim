use bigsmiles::parse;
use polysim_core::{
    builder::{linear::LinearBuilder, BuildStrategy},
    properties::solubility::flory_huggins_chi,
};

// ── Helpers ───────────────────────────────────────────────────────────────────

fn build_homo(bigsmiles: &str, n: usize) -> polysim_core::PolymerChain {
    let bs = parse(bigsmiles).expect("valid BigSMILES");
    LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(n))
        .homopolymer()
        .expect("build should succeed")
}

// BigSMILES des polymères courants
const PE_SMILES: &str = "{[]CC[]}";
const PP_SMILES: &str = "{[]CC(C)[]}";
const PS_SMILES: &str = "{[]CC(c1ccccc1)[]}";
const PVC_SMILES: &str = "{[]C(Cl)C[]}";

// ── Tests fondamentaux ────────────────────────────────────────────────────────

/// chi(PS, PS, 298K) = 0 — un polymère avec lui-même ne donne aucune
/// incompatibilité par définition de la formule Flory-Huggins.
#[test]
fn chi_self_is_zero() {
    let ps = build_homo(PS_SMILES, 50);
    let chi = flory_huggins_chi(&ps, &ps, 298.0).unwrap();
    assert!(chi < 1e-10, "chi(PS, PS) doit être ~0, obtenu {chi:.6}");
}

/// chi(PE, PVC, 298K) > 0 — deux polymères différents ont une interaction > 0.
#[test]
fn chi_positive() {
    let pe = build_homo(PE_SMILES, 50);
    let pvc = build_homo(PVC_SMILES, 50);
    let chi = flory_huggins_chi(&pe, &pvc, 298.0).unwrap();
    assert!(chi > 0.0, "chi(PE, PVC) doit être positif, obtenu {chi:.4}");
}

/// chi(A, B) == chi(B, A) — la formule Hildebrand est symétrique par construction
/// car (delta_A - delta_B)^2 = (delta_B - delta_A)^2.
#[test]
fn chi_symmetric() {
    let pe = build_homo(PE_SMILES, 50);
    let pvc = build_homo(PVC_SMILES, 50);
    let chi_ab = flory_huggins_chi(&pe, &pvc, 298.0).unwrap();
    let chi_ba = flory_huggins_chi(&pvc, &pe, 298.0).unwrap();
    assert!(
        (chi_ab - chi_ba).abs() < 1e-10,
        "chi doit être symétrique : chi(PE,PVC)={chi_ab:.6}, chi(PVC,PE)={chi_ba:.6}"
    );
}

/// chi ~ 1/T : augmenter T doit diminuer chi (PE/PVC).
#[test]
fn chi_decreases_with_temp() {
    let pe = build_homo(PE_SMILES, 50);
    let pvc = build_homo(PVC_SMILES, 50);
    let chi_low_t = flory_huggins_chi(&pe, &pvc, 298.0).unwrap();
    let chi_high_t = flory_huggins_chi(&pe, &pvc, 398.0).unwrap();
    assert!(
        chi_high_t < chi_low_t,
        "chi doit diminuer avec T : chi(298K)={chi_low_t:.4}, chi(398K)={chi_high_t:.4}"
    );
}

/// chi(PE, PVC) > chi(PE, PP) — PVC est plus différent de PE que PP (qui est
/// simplement PE avec un méthyle pendant, delta PP est proche de delta PE).
#[test]
fn chi_pe_pvc_ordering() {
    let pe = build_homo(PE_SMILES, 50);
    let pp = build_homo(PP_SMILES, 50);
    let pvc = build_homo(PVC_SMILES, 50);
    let chi_pe_pp = flory_huggins_chi(&pe, &pp, 298.0).unwrap();
    let chi_pe_pvc = flory_huggins_chi(&pe, &pvc, 298.0).unwrap();
    assert!(
        chi_pe_pvc > chi_pe_pp,
        "chi(PE,PVC)={chi_pe_pvc:.4} doit être > chi(PE,PP)={chi_pe_pp:.4}"
    );
}

// ── Invariants de la formule ──────────────────────────────────────────────────

/// chi est proportionnel à 1/T (loi de Flory-Huggins) :
/// chi(T1) / chi(T2) doit être égal à T2 / T1.
#[test]
fn chi_inverse_temperature_law() {
    let pe = build_homo(PE_SMILES, 50);
    let pvc = build_homo(PVC_SMILES, 50);
    let t1 = 298.0_f64;
    let t2 = 596.0_f64; // 2 * T1
    let chi_t1 = flory_huggins_chi(&pe, &pvc, t1).unwrap();
    let chi_t2 = flory_huggins_chi(&pe, &pvc, t2).unwrap();
    // chi_t1 / chi_t2 doit être ≈ t2 / t1 = 2.0
    let ratio = chi_t1 / chi_t2;
    assert!(
        (ratio - 2.0).abs() < 0.001,
        "chi inversement proportionnel à T : ratio={ratio:.4}, attendu 2.0"
    );
}

/// chi ne doit pas dépendre de la longueur de chaîne (propriété intensive).
#[test]
fn chi_independent_of_chain_length() {
    let pe_short = build_homo(PE_SMILES, 10);
    let pe_long = build_homo(PE_SMILES, 100);
    let pvc_short = build_homo(PVC_SMILES, 10);
    let pvc_long = build_homo(PVC_SMILES, 100);
    let chi_short = flory_huggins_chi(&pe_short, &pvc_short, 298.0).unwrap();
    let chi_long = flory_huggins_chi(&pe_long, &pvc_long, 298.0).unwrap();
    let rel_diff = (chi_short - chi_long).abs() / chi_long;
    assert!(
        rel_diff < 0.10,
        "chi doit être quasi-indépendant de n : chi(n=10)={chi_short:.4}, chi(n=100)={chi_long:.4} (diff={:.1}%)",
        rel_diff * 100.0
    );
}

// ── Cas limites ───────────────────────────────────────────────────────────────

/// n=1 ne doit pas paniquer.
#[test]
fn chi_n1_no_panic() {
    let pe = build_homo(PE_SMILES, 1);
    let pvc = build_homo(PVC_SMILES, 1);
    let result = flory_huggins_chi(&pe, &pvc, 298.0);
    assert!(result.is_ok(), "chi(n=1) ne doit pas paniquer");
    assert!(result.unwrap() >= 0.0);
}

// ── Limitation documentée — méthode Hildebrand apolaire ──────────────────────

/// PE/PS : chi VK prédit ~0.033 alors que le seuil d'immiscibilité est chi > 0.5.
/// La méthode Hildebrand ne distingue pas suffisamment deux polymères apolaires
/// (PE : delta~16.3, PS : delta~17.2 (MPa)^0.5 après correction Vmol).
/// En pratique PE/PS sont immiscibles mais VK ne le prédit pas.
#[test]
#[ignore = "Hildebrand ne distingue pas les paires apolaires/apolaires — chi(PE/PS)≈0.033 vs immiscibilité réelle (chi > 0.5). Limitation connue de la méthode VK."]
fn chi_pe_ps_immiscible() {
    let pe = build_homo(PE_SMILES, 50);
    let ps = build_homo(PS_SMILES, 50);
    let chi = flory_huggins_chi(&pe, &ps, 298.0).unwrap();
    assert!(
        chi > 0.5,
        "chi(PE, PS) prédit {chi:.4} mais l'expérience indique une immiscibilité (chi > 0.5)"
    );
}
