use clap::{Args, Parser, Subcommand};
use colored::Colorize;
use comfy_table::{Attribute, Cell, Color as TableColor, ContentArrangement, Table};
use polysim_core::{
    builder::{linear::LinearBuilder, BuildStrategy},
    parse,
    properties::{
        formula::{molecular_formula, total_atom_count},
        molecular_weight::monoisotopic_mass,
    },
};

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
struct StrategyArgs {
    /// Build chain with exactly N repeat units.
    #[arg(long, value_name = "N", help_heading = "Build strategy")]
    by_repeat: Option<usize>,

    /// Build chain targeting the given number-average molecular weight (g/mol).
    #[arg(long, value_name = "MN", help_heading = "Build strategy")]
    by_mn: Option<f64>,

    /// Build chain targeting the given exact monoisotopic mass (g/mol).
    #[arg(long, value_name = "MASS", help_heading = "Build strategy")]
    by_mass: Option<f64>,
}

fn main() {
    let cli = Cli::parse();
    match cli.command {
        Commands::Analyze { bigsmiles, strategy } => {
            if let Err(code) = run_analyze(&bigsmiles, &strategy) {
                std::process::exit(code);
            }
        }
    }
}

fn run_analyze(bigsmiles_str: &str, args: &StrategyArgs) -> Result<(), i32> {
    // --- Determine build strategy ----------------------------------------
    let (strategy, strategy_label) = if let Some(n) = args.by_repeat {
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
    };

    // --- Parse BigSMILES --------------------------------------------------
    let bigsmiles = parse(bigsmiles_str).map_err(|e| {
        eprintln!("{} {e}", "error:".red().bold());
        1
    })?;

    // --- Build chain ------------------------------------------------------
    let chain = LinearBuilder::new(bigsmiles, strategy)
        .homopolymer()
        .map_err(|e| {
            eprintln!("{} {e}", "error:".red().bold());
            1
        })?;

    // --- Compute properties -----------------------------------------------
    let mono_mass = monoisotopic_mass(&chain);
    let formula_raw = molecular_formula(&chain);
    let n_atoms = total_atom_count(&chain);
    let formula_display = subscript_digits(&formula_raw);

    // --- Banner -----------------------------------------------------------
    println!();
    let title = "  polysim — Polymer Chain Analysis  ";
    let bar = "─".repeat(title.len());
    println!("  ╭{bar}╮");
    println!("  │{}│", title.bold().cyan());
    println!("  ╰{bar}╯");
    println!();

    // --- Input summary ----------------------------------------------------
    println!(
        "  {:<11}{}",
        "BigSMILES".bold(),
        bigsmiles_str.yellow()
    );
    println!("  {:<11}{}", "Strategy".bold(), strategy_label);
    println!(
        "  {:<11}{}",
        "SMILES".bold(),
        truncate(&chain.smiles, 60).dimmed()
    );
    println!();

    // --- Results table ----------------------------------------------------
    let mut table = Table::new();
    table.load_preset(comfy_table::presets::UTF8_FULL);
    table.set_content_arrangement(ContentArrangement::Dynamic);

    table.set_header(vec![
        Cell::new("Property").add_attribute(Attribute::Bold),
        Cell::new("Value").add_attribute(Attribute::Bold),
    ]);

    table.add_row(vec![
        Cell::new("Repeat units (n)"),
        Cell::new(chain.repeat_count.to_string()).fg(TableColor::Cyan),
    ]);
    table.add_row(vec![
        Cell::new("Mn  (number-average Mw)"),
        Cell::new(format!("{:.3} g/mol", chain.mn)).fg(TableColor::Green),
    ]);
    table.add_row(vec![
        Cell::new("Mw ¹"),
        Cell::new(format!("{:.3} g/mol", chain.mn)).fg(TableColor::Green),
    ]);
    table.add_row(vec![
        Cell::new("Dispersity  Đ ¹"),
        Cell::new("1.000").fg(TableColor::Green),
    ]);
    table.add_row(vec![
        Cell::new("Monoisotopic mass"),
        Cell::new(format!("{mono_mass:.3} g/mol")).fg(TableColor::Yellow),
    ]);
    table.add_row(vec![
        Cell::new("Molecular formula"),
        Cell::new(&formula_display).fg(TableColor::Magenta),
    ]);
    table.add_row(vec![
        Cell::new("Total atoms"),
        Cell::new(n_atoms.to_string()).fg(TableColor::Cyan),
    ]);

    for line in table.to_string().lines() {
        println!("  {line}");
    }

    // --- Footnote ---------------------------------------------------------
    println!();
    println!(
        "  {} Single ideal chain — Mw = Mn, Đ = 1.000",
        "¹".dimmed()
    );
    println!("    {}",
        "Material simulation (real distributions) will be available in a future release."
            .dimmed()
            .italic()
    );
    println!();

    Ok(())
}

/// Remplace les chiffres ASCII par leurs équivalents Unicode en exposant inférieur.
fn subscript_digits(s: &str) -> String {
    s.chars()
        .map(|c| match c {
            '0' => '₀',
            '1' => '₁',
            '2' => '₂',
            '3' => '₃',
            '4' => '₄',
            '5' => '₅',
            '6' => '₆',
            '7' => '₇',
            '8' => '₈',
            '9' => '₉',
            _ => c,
        })
        .collect()
}

/// Tronque une chaîne longue avec « … » au milieu.
fn truncate(s: &str, max_len: usize) -> String {
    if s.len() <= max_len {
        s.to_string()
    } else {
        let half = (max_len.saturating_sub(1)) / 2;
        format!("{}…{}", &s[..half], &s[s.len() - half..])
    }
}

// ─── Tests unitaires ─────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // subscript_digits --------------------------------------------------------

    #[test]
    fn subscript_digits_converts_all_ten() {
        assert_eq!(subscript_digits("0123456789"), "₀₁₂₃₄₅₆₇₈₉");
    }

    #[test]
    fn subscript_digits_typical_formula() {
        assert_eq!(subscript_digits("C20H42"), "C₂₀H₄₂");
    }

    #[test]
    fn subscript_digits_formula_with_heteroatoms() {
        assert_eq!(subscript_digits("C8H8O2"), "C₈H₈O₂");
    }

    #[test]
    fn subscript_digits_no_digits_unchanged() {
        assert_eq!(subscript_digits("CHONSFClBrI"), "CHONSFClBrI");
    }

    #[test]
    fn subscript_digits_empty_string() {
        assert_eq!(subscript_digits(""), "");
    }

    // truncate ----------------------------------------------------------------

    #[test]
    fn truncate_short_string_returned_unchanged() {
        assert_eq!(truncate("CCCC", 10), "CCCC");
    }

    #[test]
    fn truncate_exact_max_len_returned_unchanged() {
        let s = "A".repeat(20);
        assert_eq!(truncate(&s, 20), s);
    }

    #[test]
    fn truncate_long_string_contains_ellipsis() {
        let long = "A".repeat(100);
        let result = truncate(&long, 20);
        assert!(
            result.contains('…'),
            "truncated string must contain '…', got: '{result}'"
        );
    }

    #[test]
    fn truncate_long_string_is_shorter_than_original() {
        let long = "A".repeat(100);
        let result = truncate(&long, 20);
        // The result should be significantly shorter than 100 chars
        assert!(result.len() < 100, "result len={}", result.len());
    }

    #[test]
    fn truncate_preserves_start_and_end() {
        // "AAAA...BBBB" → the first and last chars should be from the original
        let s = format!("{}{}", "A".repeat(50), "B".repeat(50));
        let result = truncate(&s, 20);
        assert!(result.starts_with('A'), "start should be preserved");
        assert!(result.ends_with('B'), "end should be preserved");
    }
}
