use colored::Colorize;
use polysim_core::{
    builder::EnsembleBuilder,
    distribution::{ChainLengthDistribution, Flory, LogNormal, SchulzZimm},
    parse,
    polymer::PolymerEnsemble,
    properties::ensemble::EnsembleStats,
    PolySimError,
};

use crate::display;
use crate::DistributionKind;

/// Entry point for the `generate` subcommand.
pub fn run(
    bigsmiles_str: &str,
    mn: f64,
    pdi: f64,
    distribution: &DistributionKind,
    num_chains: usize,
    seed: Option<u64>,
) -> Result<(), i32> {
    let bs = parse(bigsmiles_str).map_err(report_err)?;

    if matches!(distribution, DistributionKind::Flory) && (pdi - 2.0).abs() > 0.01 {
        eprintln!(
            "  {} Flory distribution has an inherent PDI of ~2.0; \
             --pdi {} is ignored.",
            "warning:".yellow().bold(),
            pdi
        );
        eprintln!(
            "  {}",
            "Use --distribution schulz-zimm or log-normal to target a specific PDI.".dimmed()
        );
        eprintln!();
    }

    let ensemble =
        build_ensemble(distribution, bs, mn, pdi, num_chains, seed).map_err(report_err)?;

    let stats = EnsembleStats::from_ensemble(&ensemble);
    display::print_ensemble_report(bigsmiles_str, distribution, mn, pdi, &stats);
    Ok(())
}

fn build_ensemble(
    distribution: &DistributionKind,
    bs: polysim_core::BigSmiles,
    mn: f64,
    pdi: f64,
    num_chains: usize,
    seed: Option<u64>,
) -> Result<PolymerEnsemble, PolySimError> {
    match distribution {
        DistributionKind::Flory => build(Flory, bs, mn, pdi, num_chains, seed),
        DistributionKind::LogNormal => build(LogNormal, bs, mn, pdi, num_chains, seed),
        DistributionKind::SchulzZimm => build(SchulzZimm, bs, mn, pdi, num_chains, seed),
    }
}

fn build<D: ChainLengthDistribution>(
    dist: D,
    bs: polysim_core::BigSmiles,
    mn: f64,
    pdi: f64,
    num_chains: usize,
    seed: Option<u64>,
) -> Result<PolymerEnsemble, PolySimError> {
    let mut builder = EnsembleBuilder::new(bs, dist, mn, pdi).num_chains(num_chains);
    if let Some(s) = seed {
        builder = builder.seed(s);
    }
    builder.homopolymer_ensemble()
}

fn report_err(e: impl std::fmt::Display) -> i32 {
    eprintln!("{} {e}", "error:".red().bold());
    1
}
