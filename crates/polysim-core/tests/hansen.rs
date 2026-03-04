use bigsmiles::parse;
use polysim_core::{
    builder::{linear::LinearBuilder, BuildStrategy},
    properties::solubility::{
        hansen_solubility_parameters, hildebrand_solubility_parameter, red_distance, HansenParams,
    },
};

// ── Helpers ───────────────────────────────────────────────────────────────────

fn build_homo(bigsmiles: &str, n: usize) -> polysim_core::PolymerChain {
    let bs = parse(bigsmiles).expect("valid BigSMILES");
    LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(n))
        .homopolymer()
        .expect("build should succeed")
}

/// Construit un HansenParams solvant manuellement.
fn solvent(delta_d: f64, delta_p: f64, delta_h: f64) -> HansenParams {
    let delta_t = (delta_d.powi(2) + delta_p.powi(2) + delta_h.powi(2)).sqrt();
    HansenParams {
        delta_d,
        delta_p,
        delta_h,
        delta_t,
    }
}

// ── Tests de cohérence mathématique ──────────────────────────────────────────

/// delta_t² = delta_d² + delta_p² + delta_h² (identité de Hansen)
#[test]
fn hansen_components_sum_to_total() {
    let polymers = [
        ("{[]CC[]}", "PE"),
        ("{[]CC(c1ccccc1)[]}", "PS"),
        ("{[]CC(C)(C(=O)OC)[]}", "PMMA"),
        ("{[]C(Cl)C[]}", "PVC"),
        ("{[]C(=O)CCCCC(=O)NCCCCCCN[]}", "Nylon-6,6"),
    ];
    for (bigsmiles, name) in polymers {
        let chain = build_homo(bigsmiles, 50);
        let h = hansen_solubility_parameters(&chain).unwrap();
        let expected = (h.delta_d.powi(2) + h.delta_p.powi(2) + h.delta_h.powi(2)).sqrt();
        assert!(
            (h.delta_t - expected).abs() < 0.01,
            "{name}: delta_t={:.4} doit etre sqrt(d²+p²+h²)={expected:.4}",
            h.delta_t
        );
    }
}

/// delta_t (Hansen) doit etre numeriquement egal a delta (Hildebrand)
/// car Ecoh = Ed + Ep + Eh par construction
#[test]
fn hansen_total_equals_hildebrand() {
    let polymers = [
        ("{[]CC[]}", "PE"),
        ("{[]CC(c1ccccc1)[]}", "PS"),
        ("{[]CC(C)(C(=O)OC)[]}", "PMMA"),
        ("{[]C(Cl)C[]}", "PVC"),
        ("{[]C(=O)CCCCC(=O)NCCCCCCN[]}", "Nylon-6,6"),
    ];
    for (bigsmiles, name) in polymers {
        let chain = build_homo(bigsmiles, 50);
        let h = hansen_solubility_parameters(&chain).unwrap();
        let hild = hildebrand_solubility_parameter(&chain).unwrap();
        assert!(
            (h.delta_t - hild).abs() < 0.01,
            "{name}: Hansen delta_t={:.4} doit etre egal a Hildebrand={:.4}",
            h.delta_t,
            hild
        );
    }
}

/// RED d'un polymere avec lui-meme doit etre 0
#[test]
fn hansen_red_self_is_zero() {
    let chain = build_homo("{[]CC[]}", 50);
    let h = hansen_solubility_parameters(&chain).unwrap();
    let red = red_distance(&h, &h, 5.0);
    assert!(
        red < 0.001,
        "RED(polymer, polymer) doit etre ~0, got {red:.4}"
    );
}

// ── Tests de composantes caractéristiques ────────────────────────────────────

/// PE : purement dispersif — delta_p ≈ 0 et delta_h ≈ 0
/// (CH2 et CH3 n'ont que ed, pas de ep/eh)
#[test]
fn hansen_pe_purely_dispersive() {
    let chain = build_homo("{[]CC[]}", 50);
    let h = hansen_solubility_parameters(&chain).unwrap();
    assert!(
        h.delta_d > 15.0,
        "PE delta_d doit etre > 15, got {:.2}",
        h.delta_d
    );
    assert!(
        h.delta_p < 0.1,
        "PE delta_p doit etre ~0 (apolare), got {:.2}",
        h.delta_p
    );
    assert!(
        h.delta_h < 0.1,
        "PE delta_h doit etre ~0 (pas de H-bond), got {:.2}",
        h.delta_h
    );
}

/// PS : purement dispersif — phényle n'a que Ed (ep=eh=0)
#[test]
fn hansen_ps_purely_dispersive() {
    let chain = build_homo("{[]CC(c1ccccc1)[]}", 50);
    let h = hansen_solubility_parameters(&chain).unwrap();
    assert!(
        h.delta_d > 18.0,
        "PS delta_d doit etre > 18, got {:.2}",
        h.delta_d
    );
    assert!(
        h.delta_p < 1.0,
        "PS delta_p doit etre ~0 (arene apolaire), got {:.2}",
        h.delta_p
    );
    assert!(
        h.delta_h < 1.0,
        "PS delta_h doit etre ~0, got {:.2}",
        h.delta_h
    );
}

/// PMMA : composante polaire et H-bonding non nulles (groupe ester)
#[test]
fn hansen_pmma_has_polar_and_hbond() {
    let chain = build_homo("{[]CC(C)(C(=O)OC)[]}", 50);
    let h = hansen_solubility_parameters(&chain).unwrap();
    assert!(
        h.delta_p > 3.0,
        "PMMA delta_p doit etre > 3 (ester polaire), got {:.2}",
        h.delta_p
    );
    assert!(
        h.delta_h > 3.0,
        "PMMA delta_h doit etre > 3 (ester H-bond accepteur), got {:.2}",
        h.delta_h
    );
}

/// Nylon-6,6 : forte composante H-bonding (amide)
#[test]
fn hansen_nylon66_strong_hbond() {
    let chain = build_homo("{[]C(=O)CCCCC(=O)NCCCCCCN[]}", 50);
    let h = hansen_solubility_parameters(&chain).unwrap();
    assert!(
        h.delta_h > 10.0,
        "Nylon-6,6 delta_h doit etre > 10 (amide H-bond fort), got {:.2}",
        h.delta_h
    );
    // Nylon delta_h > PS delta_h (pas de H-bond)
    let ps = hansen_solubility_parameters(&build_homo("{[]CC(c1ccccc1)[]}", 50)).unwrap();
    assert!(
        h.delta_h > ps.delta_h,
        "Nylon delta_h ({:.2}) doit etre > PS delta_h ({:.2})",
        h.delta_h,
        ps.delta_h
    );
}

// ── Tests RED (miscibilité) ───────────────────────────────────────────────────

/// PS / toluène : miscibles — RED < 1 avec R0=8.8 (rayon sphère Hansen PS, Hansen 2007)
/// Toluène Hansen : delta_d=18.0, delta_p=1.4, delta_h=2.0 (Hansen 2007)
#[test]
fn hansen_red_ps_toluene_miscible() {
    let ps_chain = build_homo("{[]CC(c1ccccc1)[]}", 50);
    let ps = hansen_solubility_parameters(&ps_chain).unwrap();
    let toluene = solvent(18.0, 1.4, 2.0);
    // R0 = 8.8 (MPa)^0.5 : rayon de la sphere de solubilite de PS (Hansen 2007 Table A.1)
    let red = red_distance(&ps, &toluene, 8.8);
    assert!(
        red < 1.0,
        "PS/toluene doit etre miscible (RED < 1), got RED={red:.3}"
    );
}

/// PS / eau : immiscibles — RED > 1 avec R0=8.8
/// Eau Hansen : delta_d=15.6, delta_p=16.0, delta_h=42.3 (Hansen 2007)
#[test]
fn hansen_red_ps_water_immiscible() {
    let ps_chain = build_homo("{[]CC(c1ccccc1)[]}", 50);
    let ps = hansen_solubility_parameters(&ps_chain).unwrap();
    let water = solvent(15.6, 16.0, 42.3);
    let red = red_distance(&ps, &water, 8.8);
    assert!(
        red > 1.0,
        "PS/eau doit etre immiscible (RED > 1), got RED={red:.3}"
    );
}

/// PE / eau : immiscibles — PE est dispersif pur, eau est très polaire/H-bond
/// Eau Hansen : delta_d=15.6, delta_p=16.0, delta_h=42.3 (Hansen 2007)
#[test]
fn hansen_red_pe_water_immiscible() {
    let pe_chain = build_homo("{[]CC[]}", 50);
    let pe = hansen_solubility_parameters(&pe_chain).unwrap();
    let water = solvent(15.6, 16.0, 42.3);
    // R0 PE ~ 8.6 (MPa)^0.5
    let red = red_distance(&pe, &water, 8.6);
    assert!(
        red > 1.0,
        "PE/eau doit etre immiscible (RED > 1), got RED={red:.3}"
    );
}

/// Nylon-6,6 vs PE : Nylon est plus proche de l'eau que PE (affinite H-bond)
/// Ce test verifie la tendance qualitative sans exiger RED < 1 (VK surestime delta(PE))
#[test]
fn hansen_nylon66_closer_to_water_than_pe() {
    let nylon_chain = build_homo("{[]C(=O)CCCCC(=O)NCCCCCCN[]}", 50);
    let pe_chain = build_homo("{[]CC[]}", 50);
    let nylon = hansen_solubility_parameters(&nylon_chain).unwrap();
    let pe = hansen_solubility_parameters(&pe_chain).unwrap();
    let water = solvent(15.6, 16.0, 42.3);
    let r0 = 10.0;
    let red_nylon = red_distance(&nylon, &water, r0);
    let red_pe = red_distance(&pe, &water, r0);
    assert!(
        red_nylon < red_pe,
        "Nylon doit etre plus proche de l'eau que PE : RED(Nylon)={red_nylon:.3}, RED(PE)={red_pe:.3}"
    );
}

// ── Tests de plage physique ───────────────────────────────────────────────────

/// Toutes les composantes doivent etre dans [0, 50] (MPa)^0.5
#[test]
fn hansen_physical_range_all_polymers() {
    let polymers = [
        ("{[]CC[]}", "PE"),
        ("{[]CC(C)[]}", "PP"),
        ("{[]CC(c1ccccc1)[]}", "PS"),
        ("{[]CC(C)(C(=O)OC)[]}", "PMMA"),
        ("{[]C(Cl)C[]}", "PVC"),
    ];
    for (bigsmiles, name) in polymers {
        let chain = build_homo(bigsmiles, 50);
        let h = hansen_solubility_parameters(&chain).unwrap();
        assert!(
            h.delta_d > 0.0 && h.delta_d < 50.0,
            "{name}: delta_d={:.2} hors plage",
            h.delta_d
        );
        assert!(
            h.delta_p >= 0.0 && h.delta_p < 50.0,
            "{name}: delta_p={:.2} hors plage",
            h.delta_p
        );
        assert!(
            h.delta_h >= 0.0 && h.delta_h < 50.0,
            "{name}: delta_h={:.2} hors plage",
            h.delta_h
        );
        assert!(
            h.delta_t > 0.0 && h.delta_t < 50.0,
            "{name}: delta_t={:.2} hors plage",
            h.delta_t
        );
    }
}

// ── Tests de cas limites ──────────────────────────────────────────────────────

/// n=1 : ne doit pas paniquer
#[test]
fn hansen_n1_no_panic() {
    let chain = build_homo("{[]CC[]}", 1);
    let result = hansen_solubility_parameters(&chain);
    assert!(result.is_ok(), "n=1 ne doit pas paniquer");
}

/// Les composantes convergent avec n (propriete intensive)
#[test]
fn hansen_converges_with_n() {
    let h50 = hansen_solubility_parameters(&build_homo("{[]CC(C)(C(=O)OC)[]}", 50)).unwrap();
    let h200 = hansen_solubility_parameters(&build_homo("{[]CC(C)(C(=O)OC)[]}", 200)).unwrap();
    assert!(
        (h50.delta_t - h200.delta_t).abs() / h200.delta_t < 0.02,
        "PMMA delta_t converge : n=50 -> {:.3}, n=200 -> {:.3}",
        h50.delta_t,
        h200.delta_t
    );
}
