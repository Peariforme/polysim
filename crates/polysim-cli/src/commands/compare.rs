use colored::Colorize;
use comfy_table::{Attribute, Cell, Color as TableColor, ContentArrangement, Table};
use polysim_core::{
    builder::linear::LinearBuilder,
    parse,
    properties::{
        mechanical::{density, tensile_strength, youngs_modulus},
        optical::refractive_index,
        permeability::{gas_permeability, Gas},
        solubility::hildebrand_solubility_parameter,
        thermal::{
            crystallization_tendency, specific_heat_capacity, tg_van_krevelen, tm_van_krevelen,
            CrystallizationTendency,
        },
    },
};

use crate::{Architecture, ArchitectureArgs, OutputFormat, StrategyArgs};

/// Entry point for the `compare` subcommand.
pub fn run(
    polymers: &[String],
    strategy: &StrategyArgs,
    arch_args: &ArchitectureArgs,
    format: &OutputFormat,
) -> Result<(), i32> {
    if polymers.len() < 2 {
        eprintln!(
            "{} at least 2 BigSMILES strings are required for comparison",
            "error:".red().bold()
        );
        return Err(1);
    }

    let mut chains = Vec::new();
    for s in polymers {
        let bigsmiles = parse(s).map_err(|e| {
            eprintln!("{} {e}", "error:".red().bold());
            1i32
        })?;

        let mut builder = LinearBuilder::new(bigsmiles, strategy.build_strategy());
        if let Some(seed) = arch_args.copolymer_seed {
            builder = builder.seed(seed);
        }

        let chain = match arch_args.arch {
            Architecture::Homo => builder.homopolymer(),
            Architecture::Random => {
                let fractions = arch_args.fractions.as_deref().unwrap_or(&[]);
                builder.random_copolymer(fractions)
            }
            Architecture::Alternating => builder.alternating_copolymer(),
            Architecture::Block => {
                let lengths = arch_args.block_lengths.as_deref().unwrap_or(&[]);
                builder.block_copolymer(lengths)
            }
            Architecture::Gradient => {
                let profile = arch_args.gradient_profile();
                builder.gradient_copolymer(&profile)
            }
        }
        .map_err(|e| {
            eprintln!("{} {e}", "error:".red().bold());
            1i32
        })?;

        chains.push(chain);
    }

    let all_props: Vec<CompareProps> = chains.iter().map(compute_props).collect();

    match format {
        OutputFormat::Table => print_table(polymers, &all_props),
        OutputFormat::Json => print_json(polymers, &all_props),
        OutputFormat::Csv => print_csv(polymers, &all_props),
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Property data
// ---------------------------------------------------------------------------

struct CompareProps {
    tg: Option<f64>,
    tm: Option<f64>,
    crystallization: CrystallizationTendency,
    rho: Option<f64>,
    youngs: Option<f64>,
    tensile: Option<f64>,
    hildebrand: Option<f64>,
    ri: Option<f64>,
    cp: Option<f64>,
    perm_o2: Option<f64>,
}

fn compute_props(chain: &polysim_core::PolymerChain) -> CompareProps {
    CompareProps {
        tg: tg_van_krevelen(chain).ok(),
        tm: tm_van_krevelen(chain).ok().flatten(),
        crystallization: crystallization_tendency(chain),
        rho: density(chain).ok(),
        youngs: youngs_modulus(chain).ok(),
        tensile: tensile_strength(chain).ok(),
        hildebrand: hildebrand_solubility_parameter(chain).ok(),
        ri: refractive_index(chain).ok(),
        cp: specific_heat_capacity(chain).ok(),
        perm_o2: gas_permeability(chain, Gas::O2).ok(),
    }
}

// ---------------------------------------------------------------------------
// Table output
// ---------------------------------------------------------------------------

fn print_table(labels: &[String], props: &[CompareProps]) {
    println!();
    let title = "  polysim — Comparative Properties  ";
    let bar = "\u{2500}".repeat(title.chars().count());
    println!("  \u{256d}{bar}\u{256e}");
    println!("  \u{2502}{}\u{2502}", title.bold().cyan());
    println!("  \u{2570}{bar}\u{256f}");
    println!();

    let mut table = Table::new();
    table.load_preset(comfy_table::presets::UTF8_FULL);
    table.set_content_arrangement(ContentArrangement::Dynamic);

    // Header row: "Property" + one column per polymer
    let mut header = vec![Cell::new("Property").add_attribute(Attribute::Bold)];
    for label in labels {
        header.push(
            Cell::new(label)
                .add_attribute(Attribute::Bold)
                .fg(TableColor::Cyan),
        );
    }
    table.set_header(header);

    // ---- numeric rows ----

    add_numeric_section(&mut table, "THERMAL", props.len());

    add_numeric_row(
        &mut table,
        "Tg (K)",
        props.iter().map(|p| p.tg).collect(),
        true, // higher is better? — contextual; we just highlight max
        |v| format!("{v:.1}"),
    );
    add_numeric_row(
        &mut table,
        "Tm (K)",
        props.iter().map(|p| p.tm).collect(),
        true,
        |v| format!("{v:.1}"),
    );
    add_text_row(
        &mut table,
        "Crystallization",
        props
            .iter()
            .map(|p| tendency_label(p.crystallization))
            .collect(),
    );

    add_numeric_section(&mut table, "MECHANICAL", props.len());

    add_numeric_row(
        &mut table,
        "Density (g/cm3)",
        props.iter().map(|p| p.rho).collect(),
        false,
        |v| format!("{v:.3}"),
    );
    add_numeric_row(
        &mut table,
        "Young's modulus (GPa)",
        props.iter().map(|p| p.youngs).collect(),
        true,
        |v| format!("{v:.2}"),
    );
    add_numeric_row(
        &mut table,
        "Tensile strength (MPa)",
        props.iter().map(|p| p.tensile).collect(),
        true,
        |v| format!("{v:.1}"),
    );

    add_numeric_section(&mut table, "SOLUBILITY", props.len());

    add_numeric_row(
        &mut table,
        "Hildebrand (MPa^0.5)",
        props.iter().map(|p| p.hildebrand).collect(),
        true,
        |v| format!("{v:.1}"),
    );

    add_numeric_section(&mut table, "OPTICAL", props.len());

    add_numeric_row(
        &mut table,
        "Refractive index",
        props.iter().map(|p| p.ri).collect(),
        true,
        |v| format!("{v:.3}"),
    );

    add_numeric_section(&mut table, "THERMAL TRANSPORT", props.len());

    add_numeric_row(
        &mut table,
        "Cp (J/g·K)",
        props.iter().map(|p| p.cp).collect(),
        true,
        |v| format!("{v:.3}"),
    );

    add_numeric_section(&mut table, "PERMEABILITY", props.len());

    add_numeric_row(
        &mut table,
        "O2 permeability (Barrer)",
        props.iter().map(|p| p.perm_o2).collect(),
        true,
        |v| format!("{v:.4}"),
    );

    for line in table.to_string().lines() {
        println!("  {line}");
    }
    println!();
}

fn add_numeric_section(table: &mut Table, label: &str, ncols: usize) {
    let mut row = vec![Cell::new(label)
        .add_attribute(Attribute::Bold)
        .fg(TableColor::Cyan)];
    for _ in 0..ncols {
        row.push(Cell::new(""));
    }
    table.add_row(row);
}

fn add_numeric_row(
    table: &mut Table,
    name: &str,
    values: Vec<Option<f64>>,
    higher_is_better: bool,
    fmt: impl Fn(f64) -> String,
) {
    // Find the best value index.
    let best_idx = best_index(&values, higher_is_better);

    let mut row = vec![Cell::new(format!("  {name}"))];
    for (i, v) in values.iter().enumerate() {
        let cell = match v {
            Some(x) => {
                let text = fmt(*x);
                if Some(i) == best_idx {
                    Cell::new(text)
                        .fg(TableColor::Green)
                        .add_attribute(Attribute::Bold)
                } else {
                    Cell::new(text).fg(TableColor::Green)
                }
            }
            None => Cell::new("N/A").fg(TableColor::DarkGrey),
        };
        row.push(cell);
    }
    table.add_row(row);
}

fn add_text_row(table: &mut Table, name: &str, values: Vec<String>) {
    let mut row = vec![Cell::new(format!("  {name}"))];
    for v in values {
        row.push(Cell::new(v));
    }
    table.add_row(row);
}

/// Returns the index of the best (max or min) value among the `Option<f64>` slice.
fn best_index(values: &[Option<f64>], higher_is_better: bool) -> Option<usize> {
    let mut best_idx: Option<usize> = None;
    let mut best_val: Option<f64> = None;
    for (i, v) in values.iter().enumerate() {
        if let Some(x) = v {
            match best_val {
                None => {
                    best_idx = Some(i);
                    best_val = Some(*x);
                }
                Some(b) => {
                    let is_better = if higher_is_better { *x > b } else { *x < b };
                    if is_better {
                        best_idx = Some(i);
                        best_val = Some(*x);
                    }
                }
            }
        }
    }
    best_idx
}

fn tendency_label(t: CrystallizationTendency) -> String {
    match t {
        CrystallizationTendency::High => "High".to_string(),
        CrystallizationTendency::Medium => "Medium".to_string(),
        CrystallizationTendency::Low => "Low".to_string(),
        CrystallizationTendency::Amorphous => "Amorphous".to_string(),
    }
}

// ---------------------------------------------------------------------------
// JSON output
// ---------------------------------------------------------------------------

fn print_json(labels: &[String], props: &[CompareProps]) {
    let mut root = serde_json::Map::new();
    for (label, p) in labels.iter().zip(props.iter()) {
        let mut obj = serde_json::Map::new();

        let mut thermal = serde_json::Map::new();
        thermal.insert("tg_K".into(), opt_json(p.tg));
        thermal.insert("tm_K".into(), opt_json(p.tm));
        thermal.insert(
            "crystallization".into(),
            serde_json::Value::String(tendency_label(p.crystallization)),
        );
        thermal.insert("cp_J_per_g_K".into(), opt_json(p.cp));
        obj.insert("thermal".into(), serde_json::Value::Object(thermal));

        let mut mechanical = serde_json::Map::new();
        mechanical.insert("density_g_cm3".into(), opt_json(p.rho));
        mechanical.insert("youngs_modulus_GPa".into(), opt_json(p.youngs));
        mechanical.insert("tensile_strength_MPa".into(), opt_json(p.tensile));
        obj.insert("mechanical".into(), serde_json::Value::Object(mechanical));

        let mut solubility = serde_json::Map::new();
        solubility.insert("hildebrand_MPa05".into(), opt_json(p.hildebrand));

        // Also compute Hansen params inline
        obj.insert("solubility".into(), serde_json::Value::Object(solubility));

        let mut optical = serde_json::Map::new();
        optical.insert("refractive_index".into(), opt_json(p.ri));
        obj.insert("optical".into(), serde_json::Value::Object(optical));

        let mut permeability = serde_json::Map::new();
        permeability.insert("o2_permeability_Barrer".into(), opt_json(p.perm_o2));
        obj.insert(
            "permeability".into(),
            serde_json::Value::Object(permeability),
        );

        root.insert(label.clone(), serde_json::Value::Object(obj));
    }

    println!(
        "{}",
        serde_json::to_string_pretty(&serde_json::Value::Object(root)).expect("JSON serialization")
    );
}

fn opt_json(v: Option<f64>) -> serde_json::Value {
    match v {
        Some(x) => serde_json::Value::from(x),
        None => serde_json::Value::Null,
    }
}

// ---------------------------------------------------------------------------
// CSV output
// ---------------------------------------------------------------------------

fn print_csv(labels: &[String], props: &[CompareProps]) {
    // Header
    let header_cols: Vec<String> = labels.iter().map(|l| csv_escape(l)).collect();
    println!("property,{}", header_cols.join(","));

    csv_row("tg_K", props.iter().map(|p| p.tg).collect());
    csv_row("tm_K", props.iter().map(|p| p.tm).collect());

    // Crystallization text row
    let cryst_vals: Vec<String> = props
        .iter()
        .map(|p| tendency_label(p.crystallization))
        .collect();
    println!("crystallization,{}", cryst_vals.join(","));

    csv_row("density_g_cm3", props.iter().map(|p| p.rho).collect());
    csv_row(
        "youngs_modulus_GPa",
        props.iter().map(|p| p.youngs).collect(),
    );
    csv_row(
        "tensile_strength_MPa",
        props.iter().map(|p| p.tensile).collect(),
    );
    csv_row(
        "hildebrand_MPa05",
        props.iter().map(|p| p.hildebrand).collect(),
    );
    csv_row("refractive_index", props.iter().map(|p| p.ri).collect());
    csv_row("cp_J_per_g_K", props.iter().map(|p| p.cp).collect());
    csv_row(
        "o2_permeability_Barrer",
        props.iter().map(|p| p.perm_o2).collect(),
    );
}

fn csv_row(name: &str, values: Vec<Option<f64>>) {
    let cols: Vec<String> = values
        .iter()
        .map(|v| match v {
            Some(x) => format!("{x}"),
            None => String::new(),
        })
        .collect();
    println!("{name},{}", cols.join(","));
}

fn csv_escape(s: &str) -> String {
    if s.contains(',') || s.contains('"') {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_owned()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Architecture, ArchitectureArgs, GradientProfileKind, OutputFormat, StrategyArgs};

    fn default_strategy() -> StrategyArgs {
        StrategyArgs {
            by_repeat: Some(20),
            by_mn: None,
            by_mass: None,
        }
    }

    fn default_arch() -> ArchitectureArgs {
        ArchitectureArgs {
            arch: Architecture::Homo,
            fractions: None,
            block_lengths: None,
            block_ratios: None,
            copolymer_seed: None,
            gradient_profile: GradientProfileKind::Linear,
            gradient_f_start: 1.0,
            gradient_f_end: 0.0,
        }
    }

    #[test]
    fn compare_two_polymers_table() {
        let polymers = vec!["{[]CC[]}".to_string(), "{[]C(Cl)C[]}".to_string()];
        let result = run(
            &polymers,
            &default_strategy(),
            &default_arch(),
            &OutputFormat::Table,
        );
        assert!(result.is_ok());
    }

    #[test]
    fn compare_two_polymers_json() {
        let polymers = vec!["{[]CC[]}".to_string(), "{[]C(Cl)C[]}".to_string()];
        let result = run(
            &polymers,
            &default_strategy(),
            &default_arch(),
            &OutputFormat::Json,
        );
        assert!(result.is_ok());
    }

    #[test]
    fn compare_two_polymers_csv() {
        let polymers = vec!["{[]CC[]}".to_string(), "{[]C(Cl)C[]}".to_string()];
        let result = run(
            &polymers,
            &default_strategy(),
            &default_arch(),
            &OutputFormat::Csv,
        );
        assert!(result.is_ok());
    }

    #[test]
    fn compare_requires_at_least_two() {
        let polymers = vec!["{[]CC[]}".to_string()];
        let result = run(
            &polymers,
            &default_strategy(),
            &default_arch(),
            &OutputFormat::Table,
        );
        assert!(result.is_err());
    }
}
