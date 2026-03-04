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
use crate::{Architecture, ArchitectureArgs, DistributionKind};

/// Entry point for the `generate` subcommand.
pub fn run(
    bigsmiles_str: &str,
    mn: f64,
    pdi: f64,
    distribution: &DistributionKind,
    num_chains: usize,
    seed: Option<u64>,
    arch_args: &ArchitectureArgs,
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

    let ensemble = build_ensemble(distribution, bs, mn, pdi, num_chains, seed, arch_args)
        .map_err(report_err)?;

    let stats = EnsembleStats::from_ensemble(&ensemble);
    display::print_ensemble_report(
        bigsmiles_str,
        distribution,
        &arch_args.arch,
        mn,
        pdi,
        &stats,
    );
    Ok(())
}

fn build_ensemble(
    distribution: &DistributionKind,
    bs: polysim_core::BigSmiles,
    mn: f64,
    pdi: f64,
    num_chains: usize,
    seed: Option<u64>,
    arch_args: &ArchitectureArgs,
) -> Result<PolymerEnsemble, PolySimError> {
    match distribution {
        DistributionKind::Flory => build(Flory, bs, mn, pdi, num_chains, seed, arch_args),
        DistributionKind::LogNormal => build(LogNormal, bs, mn, pdi, num_chains, seed, arch_args),
        DistributionKind::SchulzZimm => build(SchulzZimm, bs, mn, pdi, num_chains, seed, arch_args),
    }
}

fn build<D: ChainLengthDistribution>(
    dist: D,
    bs: polysim_core::BigSmiles,
    mn: f64,
    pdi: f64,
    num_chains: usize,
    seed: Option<u64>,
    arch_args: &ArchitectureArgs,
) -> Result<PolymerEnsemble, PolySimError> {
    let mut builder = EnsembleBuilder::new(bs, dist, mn, pdi).num_chains(num_chains);
    if let Some(s) = seed {
        builder = builder.seed(s);
    }

    match arch_args.arch {
        Architecture::Homo => builder.homopolymer_ensemble(),
        Architecture::Random => {
            let fractions = arch_args.fractions.as_deref().unwrap_or(&[]);
            builder.random_copolymer_ensemble(fractions)
        }
        Architecture::Alternating => builder.alternating_copolymer_ensemble(),
        Architecture::Block => {
            let ratios = arch_args.block_ratios.as_deref().unwrap_or(&[]);
            builder.block_copolymer_ensemble(ratios)
        }
        Architecture::Gradient => {
            let profile = arch_args.gradient_profile();
            builder.gradient_copolymer_ensemble(&profile)
        }
    }
}

fn report_err(e: impl std::fmt::Display) -> i32 {
    eprintln!("{} {e}", "error:".red().bold());
    1
}
