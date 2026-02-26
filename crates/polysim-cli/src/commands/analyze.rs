use colored::Colorize;
use polysim_core::{
    builder::{linear::LinearBuilder, BuildStrategy},
    parse,
    properties::{
        formula::{molecular_formula, total_atom_count},
        molecular_weight::monoisotopic_mass,
    },
};

use crate::display;
use crate::utils::bigsmiles_ext;
use crate::StrategyArgs;

/// All data needed to render one analysis report.
pub struct AnalysisResult {
    pub bigsmiles_str: String,
    pub strategy_label: String,
    pub begin_block: Option<String>,
    pub end_block: Option<String>,
    pub smiles: String,
    pub repeat_count: usize,
    pub mn: f64,
    pub mono_mass: f64,
    /// Raw (ASCII) molecular formula, subscript conversion is done at render time.
    pub formula_raw: String,
    pub n_atoms: usize,
    /// Mn − target, present only when `--by-mn` was used.
    pub delta_mn: Option<f64>,
    /// monoisotopic mass − target, present only when `--by-mass` was used.
    pub delta_mass: Option<f64>,
}

/// Entry point for the `analyze` subcommand.
pub fn run(bigsmiles_str: &str, args: &StrategyArgs) -> Result<(), i32> {
    let (strategy, strategy_label) = resolve_strategy(args);

    let bigsmiles = parse(bigsmiles_str).map_err(|e| {
        eprintln!("{} {e}", "error:".red().bold());
        1_i32
    })?;

    let chain = LinearBuilder::new(bigsmiles.clone(), strategy)
        .homopolymer()
        .map_err(|e| {
            eprintln!("{} {e}", "error:".red().bold());
            1_i32
        })?;

    let mono_mass = monoisotopic_mass(&chain);

    let result = AnalysisResult {
        bigsmiles_str: bigsmiles_str.to_owned(),
        strategy_label,
        begin_block: bigsmiles_ext::before_stochastic(&bigsmiles),
        end_block: bigsmiles_ext::after_stochastic(&bigsmiles),
        smiles: chain.smiles.clone(),
        repeat_count: chain.repeat_count,
        mn: chain.mn,
        mono_mass,
        formula_raw: molecular_formula(&chain),
        n_atoms: total_atom_count(&chain),
        delta_mn: args.by_mn.map(|t| chain.mn - t),
        delta_mass: args.by_mass.map(|t| mono_mass - t),
    };

    display::print_report(&result);
    Ok(())
}

// ─── Helpers ─────────────────────────────────────────────────────────────────

fn resolve_strategy(args: &StrategyArgs) -> (BuildStrategy, String) {
    if let Some(n) = args.by_repeat {
        (
            BuildStrategy::ByRepeatCount(n),
            format!("By repeat count  ·  n = {n}"),
        )
    } else if let Some(mn) = args.by_mn {
        (
            BuildStrategy::ByTargetMn(mn),
            format!("By target Mn  ·  Mn = {mn:.3} g/mol"),
        )
    } else if let Some(mass) = args.by_mass {
        (
            BuildStrategy::ByExactMass(mass),
            format!("By exact monoisotopic mass  ·  m = {mass:.3} g/mol"),
        )
    } else {
        unreachable!("clap enforces required group")
    }
}
