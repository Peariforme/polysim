use polysim_core::distribution::{ChainLengthDistribution, Flory, LogNormal, SchulzZimm};
use rand::{rngs::StdRng, SeedableRng};

const NUM_SAMPLES: usize = 10_000;
const M0_PE: f64 = 28.05; // approximate repeat-unit mass for polyethylene
const SEED: u64 = 42;

fn seeded_rng() -> StdRng {
    StdRng::seed_from_u64(SEED)
}

#[test]
fn flory_mean_xn_near_target() {
    let mut rng = seeded_rng();
    let lengths = Flory.sample(28050.0, 2.0, M0_PE, NUM_SAMPLES, &mut rng);
    let mean: f64 = lengths.iter().map(|&n| n as f64).sum::<f64>() / lengths.len() as f64;
    let target_xn = 28050.0 / M0_PE; // ~1000
    assert!(
        (mean - target_xn).abs() / target_xn < 0.1,
        "Flory mean Xn = {mean:.1}, expected ~{target_xn:.1}"
    );
}

#[test]
fn flory_pdi_near_two() {
    let mut rng = seeded_rng();
    let lengths = Flory.sample(28050.0, 2.0, M0_PE, NUM_SAMPLES, &mut rng);
    let (mn, mw) = compute_mn_mw(&lengths, M0_PE);
    let pdi = mw / mn;
    assert!(
        (pdi - 2.0).abs() < 0.2,
        "Flory PDI = {pdi:.3}, expected ~2.0"
    );
}

#[test]
fn log_normal_matches_target_mn() {
    let target_mn = 28050.0;
    let mut rng = seeded_rng();
    let lengths = LogNormal.sample(target_mn, 1.5, M0_PE, NUM_SAMPLES, &mut rng);
    let (mn, _) = compute_mn_mw(&lengths, M0_PE);
    assert!(
        (mn - target_mn).abs() / target_mn < 0.15,
        "LogNormal Mn = {mn:.1}, expected ~{target_mn:.1}"
    );
}

#[test]
fn log_normal_matches_target_pdi() {
    let target_pdi = 1.5;
    let mut rng = seeded_rng();
    let lengths = LogNormal.sample(28050.0, target_pdi, M0_PE, NUM_SAMPLES, &mut rng);
    let (mn, mw) = compute_mn_mw(&lengths, M0_PE);
    let pdi = mw / mn;
    assert!(
        (pdi - target_pdi).abs() / target_pdi < 0.2,
        "LogNormal PDI = {pdi:.3}, expected ~{target_pdi}"
    );
}

#[test]
fn schulz_zimm_matches_target_mn() {
    let target_mn = 28050.0;
    let mut rng = seeded_rng();
    let lengths = SchulzZimm.sample(target_mn, 2.0, M0_PE, NUM_SAMPLES, &mut rng);
    let (mn, _) = compute_mn_mw(&lengths, M0_PE);
    assert!(
        (mn - target_mn).abs() / target_mn < 0.15,
        "SchulzZimm Mn = {mn:.1}, expected ~{target_mn:.1}"
    );
}

#[test]
fn schulz_zimm_matches_target_pdi() {
    let target_pdi = 2.0;
    let mut rng = seeded_rng();
    let lengths = SchulzZimm.sample(28050.0, target_pdi, M0_PE, NUM_SAMPLES, &mut rng);
    let (mn, mw) = compute_mn_mw(&lengths, M0_PE);
    let pdi = mw / mn;
    assert!(
        (pdi - target_pdi).abs() / target_pdi < 0.2,
        "SchulzZimm PDI = {pdi:.3}, expected ~{target_pdi}"
    );
}

#[test]
fn all_lengths_at_least_one() {
    for dist in [
        &Flory as &dyn ChainLengthDistribution,
        &LogNormal,
        &SchulzZimm,
    ] {
        let mut rng = seeded_rng();
        let lengths = dist.sample(2805.0, 2.0, M0_PE, 1_000, &mut rng);
        assert!(
            lengths.iter().all(|&n| n >= 1),
            "{}: found chain with n < 1",
            dist.name()
        );
    }
}

#[test]
fn seeded_sampling_is_reproducible() {
    let mut rng1 = StdRng::seed_from_u64(123);
    let mut rng2 = StdRng::seed_from_u64(123);
    let a = Flory.sample(2805.0, 2.0, M0_PE, 50, &mut rng1);
    let b = Flory.sample(2805.0, 2.0, M0_PE, 50, &mut rng2);
    assert_eq!(a, b, "Same seed should produce identical results");
}

/// Helper: compute Mn and Mw from repeat-unit counts (approximate, using M0).
fn compute_mn_mw(lengths: &[usize], m0: f64) -> (f64, f64) {
    let n = lengths.len() as f64;
    let sum_mi: f64 = lengths.iter().map(|&l| l as f64 * m0).sum();
    let sum_mi2: f64 = lengths.iter().map(|&l| (l as f64 * m0).powi(2)).sum();
    let mn = sum_mi / n;
    let mw = sum_mi2 / sum_mi;
    (mn, mw)
}
