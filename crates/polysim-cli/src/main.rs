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
    /// Predict physical properties of a polymer chain.
    ///
    /// Generates a single ideal chain and computes physical/chemical properties
    /// using group-contribution methods (Van Krevelen).
    Properties {
        /// BigSMILES string, e.g. "{[]CC[]}" for polyethylene.
        bigsmiles: String,

        #[command(flatten)]
        strategy: StrategyArgs,

        #[command(flatten)]
        arch: ArchitectureArgs,

        /// Output format.
        #[arg(long, value_enum, default_value = "table")]
        format: OutputFormat,
    },

    /// Analyze a polymer chain from a BigSMILES string.
    ///
    /// Generates a single ideal chain and computes its properties:
    /// Mn, Mw, dispersity, molecular formula, monoisotopic mass, and atom count.
    Analyze {
        /// BigSMILES string, e.g. "{[]CC[]}" for polyethylene.
        bigsmiles: String,

        #[command(flatten)]
        strategy: StrategyArgs,

        #[command(flatten)]
        arch: ArchitectureArgs,
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

        #[command(flatten)]
        arch: ArchitectureArgs,
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

/// Polymer architecture and copolymer parameters.
#[derive(Args)]
pub(crate) struct ArchitectureArgs {
    /// Polymer architecture.
    #[arg(
        long,
        value_enum,
        default_value = "homo",
        help_heading = "Architecture"
    )]
    pub(crate) arch: Architecture,

    /// Weight fractions for random copolymer (comma-separated, e.g. "0.6,0.4").
    #[arg(long, value_delimiter = ',', help_heading = "Architecture")]
    pub(crate) fractions: Option<Vec<f64>>,

    /// Block lengths for block copolymer analysis (comma-separated, e.g. "50,30").
    #[arg(long, value_delimiter = ',', help_heading = "Architecture")]
    pub(crate) block_lengths: Option<Vec<usize>>,

    /// Block ratios for block copolymer ensemble (comma-separated, e.g. "0.5,0.5").
    #[arg(long, value_delimiter = ',', help_heading = "Architecture")]
    pub(crate) block_ratios: Option<Vec<f64>>,

    /// Random seed for reproducible random/gradient copolymers.
    #[arg(long, help_heading = "Architecture")]
    pub(crate) copolymer_seed: Option<u64>,

    /// Gradient profile shape (for --arch gradient).
    #[arg(
        long,
        value_enum,
        default_value = "linear",
        help_heading = "Architecture"
    )]
    pub(crate) gradient_profile: GradientProfileKind,

    /// Starting fraction of monomer A (for --arch gradient).
    #[arg(long, default_value = "1.0", help_heading = "Architecture")]
    pub(crate) gradient_f_start: f64,

    /// Ending fraction of monomer A (for --arch gradient).
    #[arg(long, default_value = "0.0", help_heading = "Architecture")]
    pub(crate) gradient_f_end: f64,
}

impl ArchitectureArgs {
    pub(crate) fn gradient_profile(&self) -> polysim_core::GradientProfile {
        match self.gradient_profile {
            GradientProfileKind::Linear => polysim_core::GradientProfile::Linear {
                f_start: self.gradient_f_start,
                f_end: self.gradient_f_end,
            },
            GradientProfileKind::Sigmoid => polysim_core::GradientProfile::Sigmoid {
                f_start: self.gradient_f_start,
                f_end: self.gradient_f_end,
            },
        }
    }
}

#[derive(Clone, ValueEnum)]
pub(crate) enum GradientProfileKind {
    Linear,
    Sigmoid,
}

#[derive(Clone, ValueEnum)]
pub(crate) enum Architecture {
    Homo,
    Random,
    Alternating,
    Block,
    Gradient,
}

impl Architecture {
    pub(crate) fn label(&self) -> &'static str {
        match self {
            Self::Homo => "Homopolymer",
            Self::Random => "Random copolymer",
            Self::Alternating => "Alternating copolymer",
            Self::Block => "Block copolymer",
            Self::Gradient => "Gradient copolymer",
        }
    }
}

#[derive(Clone, ValueEnum)]
pub(crate) enum OutputFormat {
    Table,
    Json,
    Csv,
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
        Commands::Properties {
            bigsmiles,
            strategy,
            arch,
            format,
        } => {
            if let Err(code) = commands::properties::run(&bigsmiles, &strategy, &arch, &format) {
                std::process::exit(code);
            }
        }
        Commands::Analyze {
            bigsmiles,
            strategy,
            arch,
        } => {
            if let Err(code) = commands::analyze::run(&bigsmiles, &strategy, &arch) {
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
            arch,
        } => {
            if let Err(code) =
                commands::generate::run(&bigsmiles, mn, pdi, &distribution, num_chains, seed, &arch)
            {
                std::process::exit(code);
            }
        }
    }
}
