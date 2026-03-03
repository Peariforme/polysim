use polysim_core::{
    builder::EnsembleBuilder,
    distribution::{Flory, SchulzZimm},
    parse,
    polymer::{PolymerChain, PolymerEnsemble},
    properties::ensemble::EnsembleStats,
};

#[test]
fn ensemble_mn_mw_pdi_manual() {
    // Two chains: 100 g/mol and 300 g/mol
    // Mn = (100 + 300) / 2 = 200
    // Mw = (100² + 300²) / (100 + 300) = 100000 / 400 = 250
    // PDI = 250 / 200 = 1.25
    let chains = vec![
        PolymerChain::new("CC".to_string(), 1, 100.0),
        PolymerChain::new("CCCCCC".to_string(), 3, 300.0),
    ];
    let ensemble = PolymerEnsemble::new(chains).unwrap();
    assert!((ensemble.mn() - 200.0).abs() < 0.01);
    assert!((ensemble.mw() - 250.0).abs() < 0.01);
    assert!((ensemble.pdi() - 1.25).abs() < 0.01);
}

#[test]
fn ensemble_len() {
    let chains = vec![
        PolymerChain::new("CC".to_string(), 1, 100.0),
        PolymerChain::new("CCCC".to_string(), 2, 200.0),
        PolymerChain::new("CCCCCC".to_string(), 3, 300.0),
    ];
    let ensemble = PolymerEnsemble::new(chains).unwrap();
    assert_eq!(ensemble.len(), 3);
}

#[test]
fn ensemble_empty_returns_error() {
    let result = PolymerEnsemble::new(vec![]);
    assert!(result.is_err());
}

#[test]
fn ensemble_stats_display() {
    let chains = vec![
        PolymerChain::new("CC".to_string(), 1, 100.0),
        PolymerChain::new("CCCCCC".to_string(), 3, 300.0),
    ];
    let ensemble = PolymerEnsemble::new(chains).unwrap();
    let stats = EnsembleStats::from_ensemble(&ensemble);
    let display = format!("{stats}");
    assert!(display.contains("Mn:"));
    assert!(display.contains("Mw:"));
    assert!(display.contains("PDI:"));
    assert!(display.contains("Median:"));
    assert!(display.contains("Std dev:"));
}

#[test]
fn ensemble_stats_median_and_std_dev() {
    let chains = vec![
        PolymerChain::new("CC".to_string(), 1, 100.0),
        PolymerChain::new("CCCC".to_string(), 2, 200.0),
        PolymerChain::new("CCCCCC".to_string(), 3, 300.0),
    ];
    let ensemble = PolymerEnsemble::new(chains).unwrap();
    let stats = EnsembleStats::from_ensemble(&ensemble);
    assert!((stats.mn_median - 200.0).abs() < 0.01);
    // std_dev of [100, 200, 300] = sqrt(((−100)² + 0² + 100²) / 3) = sqrt(20000/3) ≈ 81.65
    assert!((stats.mn_std_dev - 81.65).abs() < 0.1);
}

#[test]
fn ensemble_builder_polyethylene_flory() {
    let bs = parse("{[]CC[]}").unwrap();
    let ensemble = EnsembleBuilder::new(bs, Flory, 2805.0, 2.0)
        .num_chains(200)
        .seed(42)
        .homopolymer_ensemble()
        .unwrap();

    assert_eq!(ensemble.len(), 200);
    assert!(
        (ensemble.mn() - 2805.0).abs() / 2805.0 < 0.2,
        "Mn = {:.1}, expected ~2805",
        ensemble.mn()
    );
}

#[test]
fn ensemble_builder_polyethylene_schulz_zimm() {
    let bs = parse("{[]CC[]}").unwrap();
    let ensemble = EnsembleBuilder::new(bs, SchulzZimm, 2805.0, 1.5)
        .num_chains(200)
        .seed(42)
        .homopolymer_ensemble()
        .unwrap();

    assert_eq!(ensemble.len(), 200);
    let pdi = ensemble.pdi();
    assert!((pdi - 1.5).abs() < 0.5, "PDI = {pdi:.3}, expected ~1.5");
}

#[test]
fn ensemble_builder_seed_reproducibility() {
    let bs1 = parse("{[]CC[]}").unwrap();
    let bs2 = parse("{[]CC[]}").unwrap();

    let e1 = EnsembleBuilder::new(bs1, Flory, 2805.0, 2.0)
        .num_chains(50)
        .seed(99)
        .homopolymer_ensemble()
        .unwrap();
    let e2 = EnsembleBuilder::new(bs2, Flory, 2805.0, 2.0)
        .num_chains(50)
        .seed(99)
        .homopolymer_ensemble()
        .unwrap();

    let mns1: Vec<f64> = e1.chains().iter().map(|c| c.mn).collect();
    let mns2: Vec<f64> = e2.chains().iter().map(|c| c.mn).collect();
    assert_eq!(mns1, mns2, "Same seed should produce identical ensembles");
}

#[test]
fn ensemble_builder_no_stochastic_object() {
    let bs = parse("CC").unwrap();
    let result = EnsembleBuilder::new(bs, Flory, 2805.0, 2.0).homopolymer_ensemble();
    assert!(result.is_err());
}
