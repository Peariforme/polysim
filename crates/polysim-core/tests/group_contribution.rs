use bigsmiles::parse;
use polysim_core::{
    builder::{linear::LinearBuilder, BuildStrategy},
    properties::group_contribution::{total_ecoh, total_ri, total_vw, GroupDatabase},
};

// ── Helpers ───────────────────────────────────────────────────────────────────

fn build_homo(bigsmiles: &str, n: usize) -> polysim_core::PolymerChain {
    let bs = parse(bigsmiles).expect("valid BigSMILES");
    LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(n))
        .homopolymer()
        .expect("build should succeed")
}

fn group_count(
    groups: &[polysim_core::properties::group_contribution::GroupMatch],
    name: &str,
) -> usize {
    groups
        .iter()
        .filter(|gm| gm.group.name == name)
        .map(|gm| gm.count)
        .sum()
}

// ── Tests de decomposition : PE ───────────────────────────────────────────────

/// PE n=1 : CC = ethane = 2xCH3
#[test]
fn decompose_pe_n1() {
    let chain = build_homo("{[]CC[]}", 1);
    let groups = GroupDatabase::decompose(&chain).unwrap();
    assert_eq!(group_count(&groups, "-CH3"), 2, "PE n=1 : 2 CH3");
    assert_eq!(group_count(&groups, "-CH2-"), 0, "PE n=1 : 0 CH2");
}

/// PE n=5 : CCCCCCCCCC = 2xCH3 + 8xCH2
#[test]
fn decompose_pe_n5() {
    let chain = build_homo("{[]CC[]}", 5);
    let groups = GroupDatabase::decompose(&chain).unwrap();
    assert_eq!(group_count(&groups, "-CH3"), 2, "PE n=5 : 2 CH3 terminaux");
    assert_eq!(group_count(&groups, "-CH2-"), 8, "PE n=5 : 8 CH2");
    let total: usize = groups.iter().map(|gm| gm.count).sum();
    assert_eq!(total, 10, "PE n=5 : 10 atomes C au total");
}

/// PE n=10 : 2xCH3 + 18xCH2
#[test]
fn decompose_pe_n10() {
    let chain = build_homo("{[]CC[]}", 10);
    let groups = GroupDatabase::decompose(&chain).unwrap();
    assert_eq!(group_count(&groups, "-CH3"), 2);
    assert_eq!(group_count(&groups, "-CH2-"), 18);
}

// ── Tests de decomposition : PP ───────────────────────────────────────────────

/// PP n=1 : CC(C) = propane.
/// C central a 2 voisins C (C1 et branche CH3) -> 2H -> CH2.
/// Le groupe CH< n'apparait qu'a partir de n=2.
#[test]
fn decompose_pp_n1() {
    let chain = build_homo("{[]CC(C)[]}", 1);
    let groups = GroupDatabase::decompose(&chain).unwrap();
    // CC(C) : C1(CH3) - C2(CH2, 2 voisins C) - C3(CH3)
    assert_eq!(group_count(&groups, "-CH3"), 2, "PP n=1 : 2 CH3");
    assert_eq!(
        group_count(&groups, "-CH2-"),
        1,
        "PP n=1 : 1 CH2 (carbone central 2 voisins C)"
    );
    assert_eq!(
        group_count(&groups, "-CH<"),
        0,
        "PP n=1 : pas de CH (apparait a n>=2)"
    );
}

/// PP n=2 : CC(C)CC(C) = 3xCH3 + 2xCH2 + 1xCH
#[test]
fn decompose_pp_n2() {
    let chain = build_homo("{[]CC(C)[]}", 2);
    let groups = GroupDatabase::decompose(&chain).unwrap();
    assert_eq!(group_count(&groups, "-CH3"), 3, "PP n=2 : 3 CH3");
    assert_eq!(group_count(&groups, "-CH2-"), 2, "PP n=2 : 2 CH2");
    assert_eq!(group_count(&groups, "-CH<"), 1, "PP n=2 : 1 CH");
}

/// PP n=3 : aucun groupe polaire
#[test]
fn decompose_pp_n3_no_polar() {
    let chain = build_homo("{[]CC(C)[]}", 3);
    let groups = GroupDatabase::decompose(&chain).unwrap();
    assert_eq!(group_count(&groups, "-O-"), 0, "PP : pas d'oxygene");
    assert_eq!(group_count(&groups, "-Cl"), 0, "PP : pas de Cl");
    assert_eq!(group_count(&groups, "-COO-"), 0, "PP : pas d'ester");
    assert_eq!(group_count(&groups, "-C6H5"), 0, "PP : pas de phenyle");
    let ch3 = group_count(&groups, "-CH3");
    let ch2 = group_count(&groups, "-CH2-");
    let ch = group_count(&groups, "-CH<");
    assert!(ch3 + ch2 + ch > 0, "PP n=3 : doit avoir des groupes CH");
}

// ── Tests de decomposition : PS ───────────────────────────────────────────────

/// PS n=1 : 1 phenyle pendant
#[test]
fn decompose_ps_n1() {
    let chain = build_homo("{[]CC(c1ccccc1)[]}", 1);
    let groups = GroupDatabase::decompose(&chain).unwrap();
    assert_eq!(
        group_count(&groups, "-C6H5"),
        1,
        "PS n=1 : 1 phenyle pendant"
    );
    assert_eq!(
        group_count(&groups, "-C6H4-"),
        0,
        "PS n=1 : pas de phenylene"
    );
}

/// PS n=3 : 3 phenyles
#[test]
fn decompose_ps_n3() {
    let chain = build_homo("{[]CC(c1ccccc1)[]}", 3);
    let groups = GroupDatabase::decompose(&chain).unwrap();
    assert_eq!(group_count(&groups, "-C6H5"), 3, "PS n=3 : 3 phenyles");
    assert_eq!(
        group_count(&groups, "-C6H4-"),
        0,
        "PS n=3 : pas de phenylene"
    );
}

/// PS n=10 : 10 phenyles
#[test]
fn decompose_ps_n10() {
    let chain = build_homo("{[]CC(c1ccccc1)[]}", 10);
    let groups = GroupDatabase::decompose(&chain).unwrap();
    assert_eq!(group_count(&groups, "-C6H5"), 10, "PS n=10 : 10 phenyles");
}

// ── Tests de decomposition : PMMA ─────────────────────────────────────────────

/// PMMA n=1 : CC(C)(C(=O)OC) = 3xCH3 + 1xCH + 1xCOO.
/// Pour n=1 le C2 a 3 voisins C -> 1H -> CH (pas quaternaire).
#[test]
fn decompose_pmma_n1() {
    let chain = build_homo("{[]CC(C)(C(=O)OC)[]}", 1);
    let groups = GroupDatabase::decompose(&chain).unwrap();
    let ester = group_count(&groups, "-COO-");
    let ch3 = group_count(&groups, "-CH3");
    assert!(ester >= 1, "PMMA n=1 : au moins 1 ester, got {ester}");
    assert!(ch3 >= 2, "PMMA n=1 : au moins 2 CH3, got {ch3}");
    assert_eq!(group_count(&groups, "-Cl"), 0);
    assert_eq!(group_count(&groups, "-C6H5"), 0);
}

/// PMMA n=2 : le carbone quaternaire apparait
#[test]
fn decompose_pmma_n2_has_quaternary() {
    let chain = build_homo("{[]CC(C)(C(=O)OC)[]}", 2);
    let groups = GroupDatabase::decompose(&chain).unwrap();
    let c_quat = group_count(&groups, ">C<");
    let ester = group_count(&groups, "-COO-");
    assert!(
        c_quat >= 1,
        "PMMA n=2 : au moins 1 quaternaire, got {c_quat}"
    );
    assert_eq!(ester, 2, "PMMA n=2 : 2 esters");
}

/// PMMA n=3 : 3 esters
#[test]
fn decompose_pmma_n3() {
    let chain = build_homo("{[]CC(C)(C(=O)OC)[]}", 3);
    let groups = GroupDatabase::decompose(&chain).unwrap();
    assert_eq!(group_count(&groups, "-COO-"), 3, "PMMA n=3 : 3 esters");
}

// ── Tests de decomposition : PET ──────────────────────────────────────────────

/// PET n=1 : contient au moins 1 ester
#[test]
fn decompose_pet_n1() {
    let chain = build_homo("{[]C(=O)c1ccc(cc1)C(=O)OCCO[]}", 1);
    let groups = GroupDatabase::decompose(&chain).unwrap();
    let ester = group_count(&groups, "-COO-");
    assert!(ester >= 1, "PET n=1 : au moins 1 ester, got {ester}");
    let total: usize = groups.iter().map(|gm| gm.count).sum();
    assert!(
        total > 0,
        "PET n=1 : la decomposition ne doit pas etre vide"
    );
}

/// PET n=1 : le benzene disubstitue contribue (phenylene ou CH aromatiques)
#[test]
fn decompose_pet_contains_aromatic_contribution() {
    let chain = build_homo("{[]C(=O)c1ccc(cc1)C(=O)OCCO[]}", 1);
    let groups = GroupDatabase::decompose(&chain).unwrap();
    let phenylene = group_count(&groups, "-C6H4-");
    // Si le phénylène n'est pas détecté, les carbones aromatiques tombent en CH/C_quat
    let ch_aromatic = group_count(&groups, "-CH<");
    assert!(
        phenylene >= 1 || ch_aromatic >= 4,
        "PET : phenylene={phenylene}, CH aromatic={ch_aromatic} -- l'anneau benzene doit etre present"
    );
}

// ── Tests de decomposition : PVC ──────────────────────────────────────────────

/// PVC n=1 : 1 Cl
#[test]
fn decompose_pvc_n1() {
    let chain = build_homo("{[]C(Cl)C[]}", 1);
    let groups = GroupDatabase::decompose(&chain).unwrap();
    assert_eq!(group_count(&groups, "-Cl"), 1, "PVC n=1 : 1 Cl");
    assert_eq!(group_count(&groups, "-O-"), 0);
    assert_eq!(group_count(&groups, "-COO-"), 0);
}

/// PVC n=3 : 3 chlores
#[test]
fn decompose_pvc_n3() {
    let chain = build_homo("{[]C(Cl)C[]}", 3);
    let groups = GroupDatabase::decompose(&chain).unwrap();
    assert_eq!(group_count(&groups, "-Cl"), 3, "PVC n=3 : 3 groupes Cl");
}

/// PVC n=10 : 10 chlores
#[test]
fn decompose_pvc_n10() {
    let chain = build_homo("{[]C(Cl)C[]}", 10);
    let groups = GroupDatabase::decompose(&chain).unwrap();
    assert_eq!(group_count(&groups, "-Cl"), 10, "PVC n=10 : 10 groupes Cl");
}

// ── Tests de decomposition : PEO ──────────────────────────────────────────────

/// PEO n=3 : CCOCCOCCO = CH3(1) + CH2(5) + ether(2) + OH(1)
/// Le dernier O est terminal (1H -> OH, pas ether)
#[test]
fn decompose_peo_n3() {
    let chain = build_homo("{[]CCO[]}", 3);
    let groups = GroupDatabase::decompose(&chain).unwrap();
    let ether = group_count(&groups, "-O-");
    let oh = group_count(&groups, "-OH");
    // n=3 : 2 ethers internes + 1 OH terminal
    assert_eq!(ether, 2, "PEO n=3 : 2 groupes ether, got {ether}");
    assert_eq!(oh, 1, "PEO n=3 : 1 OH terminal, got {oh}");
    assert_eq!(group_count(&groups, "-Cl"), 0);
    assert_eq!(group_count(&groups, "-C6H5"), 0);
}

/// PEO n=5 : 4 ethers + 1 OH terminal
#[test]
fn decompose_peo_n5() {
    let chain = build_homo("{[]CCO[]}", 5);
    let groups = GroupDatabase::decompose(&chain).unwrap();
    assert_eq!(
        group_count(&groups, "-O-"),
        4,
        "PEO n=5 : 4 ethers internes"
    );
    assert_eq!(group_count(&groups, "-OH"), 1, "PEO n=5 : 1 OH terminal");
}

// ── Tests de couverture : tous les polymeres de reference decomposables ────────

#[test]
fn all_reference_polymers_decomposable() {
    let polymers = [
        ("{[]CC[]}", "PE"),
        ("{[]CC(C)[]}", "PP"),
        ("{[]CC(c1ccccc1)[]}", "PS"),
        ("{[]CC(C)(C(=O)OC)[]}", "PMMA"),
        ("{[]C(Cl)C[]}", "PVC"),
        ("{[]CCO[]}", "PEO"),
    ];

    for (bigsmiles, name) in polymers {
        let chain = build_homo(bigsmiles, 5);
        let result = GroupDatabase::decompose(&chain);
        assert!(
            result.is_ok(),
            "{name} : decomposition doit reussir, got {:?}",
            result.err()
        );
        let groups = result.unwrap();
        assert!(
            !groups.is_empty(),
            "{name} : doit contenir au moins un groupe"
        );
    }
}

/// SMILES valide -> decomposition OK (verification de compilation)
#[test]
fn valid_smiles_returns_ok() {
    use polysim_core::error::PolySimError;
    let chain = build_homo("{[]CC[]}", 1);
    let result = GroupDatabase::decompose(&chain);
    assert!(result.is_ok(), "SMILES valide doit reussir");
    let _ = PolySimError::GroupDecomposition("test".to_string());
}

// ── Tests de coherence physique ───────────────────────────────────────────────

/// Hierarchie d'energie de cohesion : PE < PS (Van Krevelen)
#[test]
fn polar_polymers_higher_ecoh() {
    let pe_chain = build_homo("{[]CC[]}", 10);
    let ecoh_pe = total_ecoh(&GroupDatabase::decompose(&pe_chain).unwrap());

    let ps_chain = build_homo("{[]CC(c1ccccc1)[]}", 10);
    let ecoh_ps = total_ecoh(&GroupDatabase::decompose(&ps_chain).unwrap());

    let pvc_chain = build_homo("{[]C(Cl)C[]}", 10);
    let ecoh_pvc = total_ecoh(&GroupDatabase::decompose(&pvc_chain).unwrap());

    assert!(
        ecoh_pe < ecoh_ps,
        "Ecoh(PS) doit etre > Ecoh(PE) : PE={ecoh_pe:.0}, PS={ecoh_ps:.0}"
    );
    assert!(
        ecoh_pe < ecoh_pvc,
        "Ecoh(PVC) doit etre > Ecoh(PE) : PE={ecoh_pe:.0}, PVC={ecoh_pvc:.0}"
    );
}

/// PS a un volume de Van der Waals plus grand que PE (n=1)
#[test]
fn larger_groups_higher_vw() {
    let ps_chain = build_homo("{[]CC(c1ccccc1)[]}", 1);
    let vw_ps = total_vw(&GroupDatabase::decompose(&ps_chain).unwrap());

    let pe_chain = build_homo("{[]CC[]}", 1);
    let vw_pe = total_vw(&GroupDatabase::decompose(&pe_chain).unwrap());

    assert!(
        vw_ps > vw_pe,
        "Vw(PS n=1) > Vw(PE n=1) : PS={vw_ps:.2}, PE={vw_pe:.2}"
    );
}

/// La refraction molaire est positive pour tous les polymeres
#[test]
fn molar_refraction_positive() {
    let polymers = [
        ("{[]CC[]}", "PE"),
        ("{[]CC(c1ccccc1)[]}", "PS"),
        ("{[]CC(C)(C(=O)OC)[]}", "PMMA"),
        ("{[]C(Cl)C[]}", "PVC"),
    ];

    for (bigsmiles, name) in polymers {
        let chain = build_homo(bigsmiles, 5);
        let groups = GroupDatabase::decompose(&chain).unwrap();
        let ri = total_ri(&groups);
        assert!(ri > 0.0, "{name} : refraction molaire > 0, got {ri}");
    }
}

// ── Tests de reference numerique VK ──────────────────────────────────────────

/// PE n=5 : Vw = 2xCH3(13.67) + 8xCH2(10.23) = 109.18 cm3/mol
#[test]
fn pe_vw_matches_vk_reference() {
    let chain = build_homo("{[]CC[]}", 5);
    let groups = GroupDatabase::decompose(&chain).unwrap();
    let vw = total_vw(&groups);
    // 2*13.67 + 8*10.23 = 27.34 + 81.84 = 109.18
    assert!(
        (vw - 109.18).abs() < 0.1,
        "PE n=5 : Vw attendu 109.18 cm3/mol, got {vw:.4}"
    );
}

/// PE n=10 : Ecoh = 2xCH3(4500) + 18xCH2(4100) = 82800 J/mol
#[test]
fn pe_ecoh_matches_vk_reference() {
    let chain = build_homo("{[]CC[]}", 10);
    let groups = GroupDatabase::decompose(&chain).unwrap();
    let ecoh = total_ecoh(&groups);
    let expected = 2.0 * 4500.0 + 18.0 * 4100.0; // 82800
    assert!(
        (ecoh - expected).abs() < 1.0,
        "PE n=10 : Ecoh attendu {expected:.0} J/mol, got {ecoh:.0}"
    );
}

/// PS n=1 : Vw doit etre > 80 cm3/mol (phenyle seul = 71.6)
#[test]
fn ps_vw_positive_and_large() {
    let chain = build_homo("{[]CC(c1ccccc1)[]}", 1);
    let groups = GroupDatabase::decompose(&chain).unwrap();
    let vw = total_vw(&groups);
    assert!(vw > 80.0, "PS n=1 : Vw > 80 cm3/mol, got {vw:.2}");
}

/// PVC n=3 : Ecoh contient la contribution des 3 Cl (3x12800 = 38400)
#[test]
fn pvc_ecoh_contains_cl_contribution() {
    let chain = build_homo("{[]C(Cl)C[]}", 3);
    let groups = GroupDatabase::decompose(&chain).unwrap();
    let ecoh = total_ecoh(&groups);
    // 3 Cl x 12800 = 38400, plus carbones -> ecoh > 38400
    assert!(
        ecoh > 38400.0,
        "PVC n=3 : Ecoh > 38400 J/mol, got {ecoh:.0}"
    );
}

/// Ether corrige (yg=5.3) : PEO n=1 a une contribution Yg > 5.0
#[test]
fn ether_group_uses_corrected_yg() {
    use polysim_core::properties::group_contribution::sum_contribution;
    let chain = build_homo("{[]CCO[]}", 1);
    let groups = GroupDatabase::decompose(&chain).unwrap();
    let yg_sum = sum_contribution(&groups, |g| g.yg);
    // yg(ether) = 5.3 + yg(CH2 ou CH3) -> yg_sum > 5.0
    assert!(
        yg_sum > 5.0,
        "PEO n=1 : yg_sum > 5.0 (ether corrige a 5.3), got {yg_sum:.2}"
    );
}

// ── Tests de croissance lineaire ──────────────────────────────────────────────

/// Vw de PE croit lineairement avec n
#[test]
fn vw_grows_linearly_with_n() {
    let vw5 = total_vw(&GroupDatabase::decompose(&build_homo("{[]CC[]}", 5)).unwrap());
    let vw10 = total_vw(&GroupDatabase::decompose(&build_homo("{[]CC[]}", 10)).unwrap());
    let vw15 = total_vw(&GroupDatabase::decompose(&build_homo("{[]CC[]}", 15)).unwrap());

    let delta1 = vw10 - vw5;
    let delta2 = vw15 - vw10;
    assert!(
        (delta1 - delta2).abs() < 0.01,
        "Vw croit lineairement : delta1={delta1:.4}, delta2={delta2:.4}"
    );
}

/// Ecoh de PE croit lineairement avec n
#[test]
fn ecoh_grows_linearly_with_n() {
    let ecoh5 = total_ecoh(&GroupDatabase::decompose(&build_homo("{[]CC[]}", 5)).unwrap());
    let ecoh10 = total_ecoh(&GroupDatabase::decompose(&build_homo("{[]CC[]}", 10)).unwrap());
    let ecoh20 = total_ecoh(&GroupDatabase::decompose(&build_homo("{[]CC[]}", 20)).unwrap());

    let rate_5_10 = (ecoh10 - ecoh5) / 5.0;
    let rate_10_20 = (ecoh20 - ecoh10) / 10.0;
    assert!(
        (rate_5_10 - rate_10_20).abs() < 1.0,
        "Ecoh croit a taux constant : {rate_5_10:.2} vs {rate_10_20:.2}"
    );
}

// ── Tests cas limites ─────────────────────────────────────────────────────────

/// PE n=1 : decomposition ne doit pas etre vide
#[test]
fn decompose_minimal_pe_n1_not_empty() {
    let chain = build_homo("{[]CC[]}", 1);
    let groups = GroupDatabase::decompose(&chain).unwrap();
    assert!(
        !groups.is_empty(),
        "PE n=1 : la decomposition ne doit pas etre vide"
    );
}

/// PE n=1000 : doit reussir sans panique
#[test]
fn decompose_pe_n1000_succeeds() {
    let chain = build_homo("{[]CC[]}", 1000);
    let result = GroupDatabase::decompose(&chain);
    assert!(result.is_ok(), "PE n=1000 : decomposition doit reussir");
    let groups = result.unwrap();
    assert_eq!(
        group_count(&groups, "-CH3"),
        2,
        "PE n=1000 : toujours 2 terminaux"
    );
    assert_eq!(group_count(&groups, "-CH2-"), 1998, "PE n=1000 : 1998 CH2");
}

/// PS n=100 : decomposition correcte meme avec le cycling des numeros de ring
#[test]
fn decompose_ps_n100_ring_cycling() {
    let chain = build_homo("{[]CC(c1ccccc1)[]}", 100);
    let result = GroupDatabase::decompose(&chain);
    assert!(result.is_ok(), "PS n=100 : decomposition doit reussir");
    let groups = result.unwrap();
    assert_eq!(
        group_count(&groups, "-C6H5"),
        100,
        "PS n=100 : 100 phenyles"
    );
}
