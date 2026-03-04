use colored::Colorize;
use comfy_table::{Attribute, Cell, Color as TableColor, ContentArrangement, Table};

use crate::report::AnalysisResult;
use crate::utils::format::{delta_style, subscript_digits, truncate};

/// Prints the full analysis report to stdout.
pub fn print_report(r: &AnalysisResult) {
    print_banner();
    print_summary(r);
    print_table(r);
    print_footnote();
}

// ─── Sections ────────────────────────────────────────────────────────────────

fn print_banner() {
    println!();
    let title = "  polysim — Polymer Chain Analysis  ";
    let bar = "─".repeat(title.chars().count());
    println!("  ╭{bar}╮");
    println!("  │{}│", title.bold().cyan());
    println!("  ╰{bar}╯");
    println!();
}

fn print_summary(r: &AnalysisResult) {
    println!("  {:<11}{}", "BigSMILES".bold(), r.bigsmiles_str.yellow());
    println!("  {:<11}{}", "Arch".bold(), r.architecture_label.cyan());
    println!("  {:<11}{}", "Strategy".bold(), r.strategy_label);
    if let Some(ref bb) = r.begin_block {
        println!("  {:<11}{}", "Begin".bold(), bb.yellow());
    }
    if let Some(ref eb) = r.end_block {
        println!("  {:<11}{}", "End".bold(), eb.yellow());
    }
    println!(
        "  {:<11}{}",
        "SMILES".bold(),
        truncate(&r.smiles, 60).dimmed()
    );
    println!();
}

fn print_table(r: &AnalysisResult) {
    let table = build_table(r);
    for line in table.to_string().lines() {
        println!("  {line}");
    }
}

fn print_footnote() {
    println!();
    println!("  {} Single ideal chain — Mw = Mn, Đ = 1.000", "¹".dimmed());
    println!(
        "    {}",
        "Material simulation (real distributions) will be available in a future release."
            .dimmed()
            .italic()
    );
    println!();
}

// ─── Table construction ──────────────────────────────────────────────────────

fn build_table(r: &AnalysisResult) -> Table {
    let mut table = Table::new();
    table.load_preset(comfy_table::presets::UTF8_FULL);
    table.set_content_arrangement(ContentArrangement::Dynamic);
    table.set_header(vec![
        Cell::new("Property").add_attribute(Attribute::Bold),
        Cell::new("Value").add_attribute(Attribute::Bold),
    ]);

    table.add_row(vec![
        Cell::new("Repeat units (n)"),
        Cell::new(r.repeat_count.to_string()).fg(TableColor::Cyan),
    ]);

    add_mn_rows(&mut table, r);

    table.add_row(vec![
        Cell::new("Mw¹"),
        Cell::new(format!("{:.3} g/mol", r.mn)).fg(TableColor::Green),
    ]);
    table.add_row(vec![
        Cell::new("Dispersity Đ¹"),
        Cell::new("1.000").fg(TableColor::Green),
    ]);

    add_mono_rows(&mut table, r);

    table.add_row(vec![
        Cell::new("Molecular formula"),
        Cell::new(subscript_digits(&r.formula_raw)).fg(TableColor::Magenta),
    ]);
    table.add_row(vec![
        Cell::new("Total atoms"),
        Cell::new(r.n_atoms.to_string()).fg(TableColor::Cyan),
    ]);

    table
}

fn add_mn_rows(table: &mut Table, r: &AnalysisResult) {
    table.add_row(vec![
        Cell::new("Mn (number-average)"),
        Cell::new(format!("{:.3} g/mol", r.mn)).fg(TableColor::Green),
    ]);
    if let Some(d) = r.delta_mn {
        let (sign, color) = delta_style(d, r.mn);
        table.add_row(vec![
            Cell::new("Δ Mn (achieved − target)").fg(TableColor::DarkGrey),
            Cell::new(format!("{sign}{d:.3} g/mol")).fg(color),
        ]);
    }
}

fn add_mono_rows(table: &mut Table, r: &AnalysisResult) {
    table.add_row(vec![
        Cell::new("Monoisotopic mass"),
        Cell::new(format!("{:.3} g/mol", r.mono_mass)).fg(TableColor::Yellow),
    ]);
    if let Some(d) = r.delta_mass {
        let (sign, color) = delta_style(d, r.mono_mass);
        table.add_row(vec![
            Cell::new("Δ mono (achieved − target)").fg(TableColor::DarkGrey),
            Cell::new(format!("{sign}{d:.3} g/mol")).fg(color),
        ]);
    }
}

// ═══ Ensemble report ═════════════════════════════════════════════════════════

use crate::{Architecture, DistributionKind};
use polysim_core::properties::ensemble::EnsembleStats;

/// Prints the ensemble generation report to stdout.
pub fn print_ensemble_report(
    bigsmiles_str: &str,
    distribution: &DistributionKind,
    architecture: &Architecture,
    target_mn: f64,
    target_pdi: f64,
    stats: &EnsembleStats,
) {
    print_ensemble_banner();
    print_ensemble_summary(
        bigsmiles_str,
        distribution,
        architecture,
        target_mn,
        target_pdi,
        stats,
    );
    print_ensemble_table(stats, target_mn, target_pdi);
    println!();
}

fn print_ensemble_banner() {
    println!();
    let title = "  polysim — Polydisperse Ensemble  ";
    let bar = "─".repeat(title.chars().count());
    println!("  ╭{bar}╮");
    println!("  │{}│", title.bold().cyan());
    println!("  ╰{bar}╯");
    println!();
}

fn print_ensemble_summary(
    bigsmiles_str: &str,
    distribution: &DistributionKind,
    architecture: &Architecture,
    target_mn: f64,
    target_pdi: f64,
    stats: &EnsembleStats,
) {
    println!("  {:<15}{}", "BigSMILES".bold(), bigsmiles_str.yellow());
    println!("  {:<15}{}", "Arch".bold(), architecture.label().cyan());
    println!(
        "  {:<15}{}",
        "Distribution".bold(),
        distribution.label().cyan()
    );
    println!("  {:<15}{:.1} g/mol", "Target Mn".bold(), target_mn);
    println!("  {:<15}{:.3}", "Target PDI".bold(), target_pdi);
    println!("  {:<15}{}", "Chains".bold(), stats.num_chains);
    println!();
}

fn print_ensemble_table(stats: &EnsembleStats, target_mn: f64, target_pdi: f64) {
    let mut table = Table::new();
    table.load_preset(comfy_table::presets::UTF8_FULL);
    table.set_content_arrangement(ContentArrangement::Dynamic);
    table.set_header(vec![
        Cell::new("Property").add_attribute(Attribute::Bold),
        Cell::new("Value").add_attribute(Attribute::Bold),
    ]);

    // Mn
    table.add_row(vec![
        Cell::new("Mn (number-average)"),
        Cell::new(format!("{:.1} g/mol", stats.mn)).fg(TableColor::Green),
    ]);
    let delta_mn = stats.mn - target_mn;
    let (sign, color) = delta_style(delta_mn, target_mn);
    table.add_row(vec![
        Cell::new("Δ Mn (achieved − target)").fg(TableColor::DarkGrey),
        Cell::new(format!("{sign}{delta_mn:.1} g/mol")).fg(color),
    ]);

    // Mw
    table.add_row(vec![
        Cell::new("Mw (weight-average)"),
        Cell::new(format!("{:.1} g/mol", stats.mw)).fg(TableColor::Green),
    ]);

    // PDI
    let pdi_color = if (stats.pdi - target_pdi).abs() / target_pdi < 0.1 {
        TableColor::Green
    } else {
        TableColor::Yellow
    };
    table.add_row(vec![
        Cell::new("Dispersity Đ (Mw/Mn)"),
        Cell::new(format!("{:.3}", stats.pdi)).fg(pdi_color),
    ]);

    // Median
    table.add_row(vec![
        Cell::new("Mn median"),
        Cell::new(format!("{:.1} g/mol", stats.mn_median)).fg(TableColor::Green),
    ]);

    // Std dev
    table.add_row(vec![
        Cell::new("Mn std dev (σ)"),
        Cell::new(format!("{:.1} g/mol", stats.mn_std_dev)).fg(TableColor::Yellow),
    ]);

    // Range
    table.add_row(vec![
        Cell::new("Mn range"),
        Cell::new(format!("{:.1} — {:.1} g/mol", stats.mn_min, stats.mn_max)).fg(TableColor::Cyan),
    ]);

    for line in table.to_string().lines() {
        println!("  {line}");
    }
}
