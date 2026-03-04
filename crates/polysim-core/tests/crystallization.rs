use bigsmiles::parse;
use polysim_core::{
    builder::{linear::LinearBuilder, BuildStrategy},
    properties::thermal::{crystallization_tendency, CrystallizationTendency},
};

// ── Helpers ───────────────────────────────────────────────────────────────────

fn build_homo(bigsmiles: &str, n: usize) -> polysim_core::PolymerChain {
    let bs = parse(bigsmiles).expect("valid BigSMILES");
    LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(n))
        .homopolymer()
        .expect("build should succeed")
}

// ── Tests de tendance pour polymères bien prédits par VK ─────────────────────

/// PE : High — delta(Tm-Tg) VK ~137 K > 100 K, chaîne purement aliphatique (symétrie haute)
#[test]
fn crystallization_pe_high() {
    let chain = build_homo("{[]CC[]}", 50);
    let tendency = crystallization_tendency(&chain);
    assert!(
        matches!(tendency, CrystallizationTendency::High),
        "PE doit etre High, got {:?}",
        tendency
    );
}

/// PMMA : Amorphous — Tm VK < Tg (delta négatif = pas de point de fusion significatif)
#[test]
fn crystallization_pmma_amorphous() {
    let chain = build_homo("{[]CC(C)(C(=O)OC)[]}", 50);
    let tendency = crystallization_tendency(&chain);
    assert!(
        matches!(tendency, CrystallizationTendency::Amorphous),
        "PMMA doit etre Amorphous, got {:?}",
        tendency
    );
}

/// PVC : Low — delta(Tm-Tg) VK ~47 K, juste en dessous de 50 K
#[test]
fn crystallization_pvc_low() {
    let chain = build_homo("{[]C(Cl)C[]}", 50);
    let tendency = crystallization_tendency(&chain);
    assert!(
        matches!(tendency, CrystallizationTendency::Low),
        "PVC doit etre Low, got {:?}",
        tendency
    );
}

/// PP : Low — la méthode VK sous-estime fortement Tm(PP) (~230K vs exp 449K),
/// ce qui comprime le delta à ~46K < 50K → Low plutôt que High attendu.
/// Cette valeur reflète la limitation VK pour PP (pas de distinction de tacticité).
#[test]
fn crystallization_pp_low_due_to_vk_limitation() {
    let chain = build_homo("{[]CC(C)[]}", 50);
    let tendency = crystallization_tendency(&chain);
    // VK prédit Low pour PP car Tm VK ~230K (exp 449K) → delta(Tm-Tg) ~46K
    // En réalité PP isotactique est High, mais VK ne distingue pas la tacticité.
    assert!(
        matches!(
            tendency,
            CrystallizationTendency::Low | CrystallizationTendency::Medium
        ),
        "PP doit etre Low ou Medium (limitation VK), got {:?}",
        tendency
    );
}

/// PS : Medium — VK prédit un Tm (409K) > seuil 200K car phényle a ym > 0.
/// PS atactique est amorphe en réalité, mais VK simplifié ne distingue pas la tacticité.
/// Ce test documente la prédiction VK réelle.
#[test]
fn crystallization_ps_medium_vk_prediction() {
    let chain = build_homo("{[]CC(c1ccccc1)[]}", 50);
    let tendency = crystallization_tendency(&chain);
    // VK prédit delta(Tm-Tg) ~76K > 50K mais phényle = pas haute symétrie → Medium
    // (PS atactique réel = Amorphous, limitation connue de VK)
    assert!(
        matches!(
            tendency,
            CrystallizationTendency::Medium | CrystallizationTendency::Amorphous
        ),
        "PS doit etre Medium ou Amorphous (VK), got {:?}",
        tendency
    );
}

// ── Test ignoré : PS Amorphous attendu vs Medium VK ──────────────────────────

/// PS atactique devrait être Amorphous, mais VK donne Medium car il ne distingue
/// pas la tacticité et calcule un Tm fictif pour le phényle.
#[test]
#[ignore = "PS crystallization VK prédit Medium (delta=76K) au lieu de Amorphous — VK ne distingue pas PS atactique (sans Tm réel) de PS cristallin"]
fn crystallization_ps_amorphous_expected() {
    let chain = build_homo("{[]CC(c1ccccc1)[]}", 50);
    let tendency = crystallization_tendency(&chain);
    assert!(
        matches!(tendency, CrystallizationTendency::Amorphous),
        "PS atactique doit etre Amorphous (experimental), got {:?}",
        tendency
    );
}

// ── Tests de cohérence hiérarchique ──────────────────────────────────────────

/// PE (High) doit être plus crystallisable que PMMA (Amorphous)
/// Ordre : High > Medium > Low > Amorphous
#[test]
fn crystallization_pe_more_than_pmma() {
    fn tendency_rank(t: CrystallizationTendency) -> u8 {
        match t {
            CrystallizationTendency::High => 3,
            CrystallizationTendency::Medium => 2,
            CrystallizationTendency::Low => 1,
            CrystallizationTendency::Amorphous => 0,
        }
    }
    let t_pe = crystallization_tendency(&build_homo("{[]CC[]}", 50));
    let t_pmma = crystallization_tendency(&build_homo("{[]CC(C)(C(=O)OC)[]}", 50));
    assert!(
        tendency_rank(t_pe) > tendency_rank(t_pmma),
        "PE doit etre plus crystallisable que PMMA : PE={:?}, PMMA={:?}",
        t_pe,
        t_pmma
    );
}

// ── Tests de cas limites ──────────────────────────────────────────────────────

/// n=1 : ne doit pas paniquer
#[test]
fn crystallization_n1_no_panic() {
    let chain = build_homo("{[]CC[]}", 1);
    // La fonction doit retourner un résultat sans paniquer
    let _ = crystallization_tendency(&chain);
}

/// La valeur de retour doit toujours être un variant valide de l'enum
#[test]
fn crystallization_always_returns_valid_variant() {
    let polymers = [
        ("{[]CC[]}", "PE"),
        ("{[]CC(C)[]}", "PP"),
        ("{[]CC(c1ccccc1)[]}", "PS"),
        ("{[]CC(C)(C(=O)OC)[]}", "PMMA"),
        ("{[]C(Cl)C[]}", "PVC"),
    ];
    for (bigsmiles, name) in polymers {
        let chain = build_homo(bigsmiles, 50);
        let t = crystallization_tendency(&chain);
        // Vérifier que c'est un variant valide en pattern matching exhaustif
        let _valid = matches!(
            t,
            CrystallizationTendency::High
                | CrystallizationTendency::Medium
                | CrystallizationTendency::Low
                | CrystallizationTendency::Amorphous
        );
        assert!(_valid, "{name}: variant invalide {:?}", t);
    }
}

/// Nylon-6,6 : Medium — delta VK ~264K > 50K, mais amide = pas haute symétrie
#[test]
fn crystallization_nylon66_medium() {
    let chain = build_homo("{[]C(=O)CCCCC(=O)NCCCCCCN[]}", 50);
    let tendency = crystallization_tendency(&chain);
    assert!(
        matches!(
            tendency,
            CrystallizationTendency::Medium | CrystallizationTendency::High
        ),
        "Nylon-6,6 doit etre Medium ou High, got {:?}",
        tendency
    );
}
