mod commands;
mod display;
mod utils;

use clap::{Args, Parser, Subcommand};

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
}

/// Build strategy — exactly one of the three flags must be provided.
#[derive(Args)]
#[group(required = true, multiple = false)]
pub struct StrategyArgs {
    /// Build chain with exactly N repeat units.
    #[arg(long, value_name = "N", help_heading = "Build strategy")]
    pub by_repeat: Option<usize>,

    /// Build chain targeting the given number-average molecular weight (g/mol).
    #[arg(long, value_name = "MN", help_heading = "Build strategy")]
    pub by_mn: Option<f64>,

    /// Build chain targeting the given exact monoisotopic mass (g/mol).
    #[arg(long, value_name = "MASS", help_heading = "Build strategy")]
    pub by_mass: Option<f64>,
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
    }
}
