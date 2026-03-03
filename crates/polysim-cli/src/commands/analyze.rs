use bigsmiles::BigSmilesSegment;
use colored::Colorize;
use polysim_core::{
    builder::linear::LinearBuilder,
    parse,
    properties::{
        formula::{molecular_formula, total_atom_count},
        molecular_weight::monoisotopic_mass,
    },
};

use crate::display;
use crate::report::AnalysisResult;
use crate::StrategyArgs;

/// Entry point for the `analyze` subcommand.
pub fn run(bigsmiles_str: &str, args: &StrategyArgs) -> Result<(), i32> {
    let bigsmiles = parse(bigsmiles_str).map_err(report_err)?;

    let chain = LinearBuilder::new(bigsmiles.clone(), args.build_strategy())
        .homopolymer()
        .map_err(report_err)?;

    let mono_mass = monoisotopic_mass(&chain);

    let result = AnalysisResult {
        bigsmiles_str: bigsmiles_str.to_owned(),
        strategy_label: args.label(),
        begin_block: segments_to_smiles(bigsmiles.prefix_segments()),
        end_block: segments_to_smiles(bigsmiles.suffix_segments()),
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

fn report_err(e: impl std::fmt::Display) -> i32 {
    eprintln!("{} {e}", "error:".red().bold());
    1
}

fn segments_to_smiles(segs: &[BigSmilesSegment]) -> Option<String> {
    let s: String = segs
        .iter()
        .filter_map(|seg| match seg {
            BigSmilesSegment::Smiles(mol) => Some(format!("{mol}")),
            BigSmilesSegment::Stochastic(_) => None,
        })
        .collect();
    (!s.is_empty()).then_some(s)
}
