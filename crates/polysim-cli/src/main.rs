mod commands;
mod display;
mod report;
mod utils;

use clap::{Args, Parser, Subcommand, ValueEnum};
use polysim_core::BuildStrategy;

/// Polymer structure generator and property simulator.
#[derive(Parser)]
#[command(
    name = "polysim",
    version,
    about = "Polymer structure generator and property simulator",
    long_about = "Generates concrete polymer chains from BigSMILES notation\n\
                  and computes physical/chemical properties on them."
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Analyze a polymer chain from a BigSMILES string.
    ///
    /// Generates a single ideal chain and computes its properties:
    /// Mn, Mw, dispersity, molecular formula, monoisotopic mass, and atom count.
    Analyze {
        /// BigSMILES string, e.g. "{[]CC[]}" for polyethylene.
        ///
        /// The stochastic object {…} must contain exactly one repeat unit
        /// (homopolymer). Copolymers will be supported in a future release.
        bigsmiles: String,

        #[command(flatten)]
        strategy: StrategyArgs,
    },

    /// Generate a polydisperse ensemble of polymer chains.
    ///
    /// Samples chain lengths from a statistical distribution and reports
    /// ensemble-averaged properties (Mn, Mw, PDI).
    Generate {
        /// BigSMILES string, e.g. "{[]CC[]}" for polyethylene.
        bigsmiles: String,

        /// Target number-average molecular weight (g/mol).
        #[arg(long)]
        mn: f64,

        /// Target polydispersity index (Mw/Mn).
        #[arg(long, default_value = "2.0")]
        pdi: f64,

        /// Chain length distribution model.
        #[arg(long, value_enum, default_value = "schulz-zimm")]
        distribution: DistributionKind,

        /// Number of chains to generate.
        #[arg(long, env = "POLYSIM_NUM_CHAINS", default_value = "100")]
        num_chains: usize,

        /// Random seed for reproducible results.
        #[arg(long)]
        seed: Option<u64>,
    },
}

/// Build strategy — exactly one of the three flags must be provided.
#[derive(Args)]
#[group(required = true, multiple = false)]
pub(crate) struct StrategyArgs {
    /// Build chain with exactly N repeat units.
    #[arg(long, value_name = "N", help_heading = "Build strategy")]
    pub(crate) by_repeat: Option<usize>,

    /// Build chain targeting the given number-average molecular weight (g/mol).
    #[arg(long, value_name = "MN", help_heading = "Build strategy")]
    pub(crate) by_mn: Option<f64>,

    /// Build chain targeting the given exact monoisotopic mass (g/mol).
    #[arg(long, value_name = "MASS", help_heading = "Build strategy")]
    pub(crate) by_mass: Option<f64>,
}

impl StrategyArgs {
    pub(crate) fn build_strategy(&self) -> BuildStrategy {
        self.by_repeat
            .map(BuildStrategy::ByRepeatCount)
            .or_else(|| self.by_mn.map(BuildStrategy::ByTargetMn))
            .or_else(|| self.by_mass.map(BuildStrategy::ByExactMass))
            .expect("clap enforces required group")
    }

    pub(crate) fn label(&self) -> String {
        self.by_repeat
            .map(|n| format!("By repeat count  ·  n = {n}"))
            .or_else(|| {
                self.by_mn
                    .map(|mn| format!("By target Mn  ·  Mn = {mn:.3} g/mol"))
            })
            .or_else(|| {
                self.by_mass
                    .map(|mass| format!("By exact monoisotopic mass  ·  m = {mass:.3} g/mol"))
            })
            .expect("clap enforces required group")
    }
}

#[derive(Clone, ValueEnum)]
pub(crate) enum DistributionKind {
    Flory,
    LogNormal,
    SchulzZimm,
}

impl DistributionKind {
    pub(crate) fn label(&self) -> &'static str {
        match self {
            Self::Flory => "Flory (most probable)",
            Self::LogNormal => "Log-normal",
            Self::SchulzZimm => "Schulz-Zimm",
        }
    }
}

fn main() {
    let cli = Cli::parse();
    match cli.command {
        Commands::Analyze {
            bigsmiles,
            strategy,
        } => {
            if let Err(code) = commands::analyze::run(&bigsmiles, &strategy) {
                std::process::exit(code);
            }
        }
        Commands::Generate {
            bigsmiles,
            mn,
            pdi,
            distribution,
            num_chains,
            seed,
        } => {
            if let Err(code) =
                commands::generate::run(&bigsmiles, mn, pdi, &distribution, num_chains, seed)
            {
                std::process::exit(code);
            }
        }
    }
}
