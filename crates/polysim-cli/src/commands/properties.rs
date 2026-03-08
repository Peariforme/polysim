use colored::Colorize;
use comfy_table::{Attribute, Cell, Color as TableColor, ContentArrangement, Table};
use polysim_core::{
    builder::linear::LinearBuilder,
    parse,
    properties::{
        mechanical::{density, tensile_strength, youngs_modulus},
        optical::refractive_index,
        solubility::{hansen_solubility_parameters, hildebrand_solubility_parameter, HansenParams},
        thermal::{
            crystallization_tendency, tg_van_krevelen, tm_van_krevelen, CrystallizationTendency,
        },
    },
    PolymerChain,
};

use crate::{Architecture, ArchitectureArgs, OutputFormat, StrategyArgs};

/// Entry point for the `properties` subcommand.
pub fn run(
    bigsmiles_str: &str,
    args: &StrategyArgs,
    arch_args: &ArchitectureArgs,
    format: &OutputFormat,
) -> Result<(), i32> {
    let bigsmiles = parse(bigsmiles_str).map_err(report_err)?;

    let mut builder = LinearBuilder::new(bigsmiles, args.build_strategy());
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
    .map_err(report_err)?;

    let props = compute_properties(&chain);

    match format {
        OutputFormat::Table => print_table(bigsmiles_str, args, arch_args, &props),
        OutputFormat::Json => print_json(&props),
        OutputFormat::Csv => print_csv(&props),
    }

    Ok(())
}

// ---- Property data --------------------------------------------------------

struct PropertySet {
    tg: Option<f64>,
    tm: Option<f64>,
    crystallization: CrystallizationTendency,
    rho: Option<f64>,
    youngs: Option<f64>,
    tensile: Option<f64>,
    hildebrand: Option<f64>,
    hansen: Option<HansenParams>,
    ri: Option<f64>,
}

fn compute_properties(chain: &PolymerChain) -> PropertySet {
    PropertySet {
        tg: tg_van_krevelen(chain).ok(),
        tm: tm_van_krevelen(chain).ok().flatten(),
        crystallization: crystallization_tendency(chain),
        rho: density(chain).ok(),
        youngs: youngs_modulus(chain).ok(),
        tensile: tensile_strength(chain).ok(),
        hildebrand: hildebrand_solubility_parameter(chain).ok(),
        hansen: hansen_solubility_parameters(chain).ok(),
        ri: refractive_index(chain).ok(),
    }
}

// ---- Table output ---------------------------------------------------------

fn print_table(
    bigsmiles_str: &str,
    args: &StrategyArgs,
    arch_args: &ArchitectureArgs,
    p: &PropertySet,
) {
    println!();
    let title = "  polysim — Physical Properties  ";
    let bar = "\u{2500}".repeat(title.chars().count());
    println!("  \u{256d}{bar}\u{256e}");
    println!("  \u{2502}{}\u{2502}", title.bold().cyan());
    println!("  \u{2570}{bar}\u{256f}");
    println!();

    println!("  {:<11}{}", "BigSMILES".bold(), bigsmiles_str.yellow());
    println!("  {:<11}{}", "Arch".bold(), arch_args.arch.label().cyan());
    println!("  {:<11}{}", "Strategy".bold(), args.label());
    println!();

    let mut table = Table::new();
    table.load_preset(comfy_table::presets::UTF8_FULL);
    table.set_content_arrangement(ContentArrangement::Dynamic);
    table.set_header(vec![
        Cell::new("Property").add_attribute(Attribute::Bold),
        Cell::new("Value").add_attribute(Attribute::Bold),
        Cell::new("Unit").add_attribute(Attribute::Bold),
        Cell::new("Confidence").add_attribute(Attribute::Bold),
    ]);

    // Thermal
    add_section(&mut table, "THERMAL");
    add_row(
        &mut table,
        "Tg (Van Krevelen)",
        p.tg.map(|v| format!("{v:.1}")),
        "K",
        stars(3),
    );
    add_row(
        &mut table,
        "Tm (Van Krevelen)",
        p.tm.map(|v| format!("{v:.1}")),
        "K",
        stars(2),
    );
    add_row(
        &mut table,
        "Crystallization",
        Some(tendency_label(p.crystallization)),
        "",
        stars(1),
    );

    // Mechanical
    add_section(&mut table, "MECHANICAL");
    add_row(
        &mut table,
        "Density",
        p.rho.map(|v| format!("{v:.3}")),
        "g/cm3",
        stars(3),
    );
    add_row(
        &mut table,
        "Young's modulus",
        p.youngs.map(|v| format!("{v:.2}")),
        "GPa",
        stars(1),
    );
    add_row(
        &mut table,
        "Tensile strength",
        p.tensile.map(|v| format!("{v:.1}")),
        "MPa",
        stars(1),
    );

    // Solubility
    add_section(&mut table, "SOLUBILITY");
    add_row(
        &mut table,
        "Hildebrand delta",
        p.hildebrand.map(|v| format!("{v:.1}")),
        "MPa^0.5",
        stars(2),
    );
    if let Some(ref h) = p.hansen {
        add_row(
            &mut table,
            "Hansen delta_d",
            Some(format!("{:.1}", h.delta_d)),
            "MPa^0.5",
            stars(2),
        );
        add_row(
            &mut table,
            "Hansen delta_p",
            Some(format!("{:.1}", h.delta_p)),
            "MPa^0.5",
            stars(2),
        );
        add_row(
            &mut table,
            "Hansen delta_h",
            Some(format!("{:.1}", h.delta_h)),
            "MPa^0.5",
            stars(2),
        );
    } else {
        add_row(&mut table, "Hansen delta_d", None, "MPa^0.5", stars(2));
        add_row(&mut table, "Hansen delta_p", None, "MPa^0.5", stars(2));
        add_row(&mut table, "Hansen delta_h", None, "MPa^0.5", stars(2));
    }

    // Optical
    add_section(&mut table, "OPTICAL");
    add_row(
        &mut table,
        "Refractive index",
        p.ri.map(|v| format!("{v:.3}")),
        "",
        stars(2),
    );

    for line in table.to_string().lines() {
        println!("  {line}");
    }
    println!();
}

fn add_section(table: &mut Table, label: &str) {
    table.add_row(vec![
        Cell::new(label)
            .add_attribute(Attribute::Bold)
            .fg(TableColor::Cyan),
        Cell::new(""),
        Cell::new(""),
        Cell::new(""),
    ]);
}

fn add_row(table: &mut Table, name: &str, value: Option<String>, unit: &str, confidence: &str) {
    let val_cell = match value {
        Some(v) => Cell::new(v).fg(TableColor::Green),
        None => Cell::new("N/A").fg(TableColor::DarkGrey),
    };
    table.add_row(vec![
        Cell::new(format!("  {name}")),
        val_cell,
        Cell::new(unit).fg(TableColor::DarkGrey),
        Cell::new(confidence).fg(TableColor::Yellow),
    ]);
}

fn stars(n: u8) -> &'static str {
    match n {
        3 => "\u{2605}\u{2605}\u{2605}",
        2 => "\u{2605}\u{2605}\u{2606}",
        1 => "\u{2605}\u{2606}\u{2606}",
        _ => "\u{2606}\u{2606}\u{2606}",
    }
}

fn tendency_label(t: CrystallizationTendency) -> String {
    match t {
        CrystallizationTendency::High => "High".to_string(),
        CrystallizationTendency::Medium => "Medium".to_string(),
        CrystallizationTendency::Low => "Low".to_string(),
        CrystallizationTendency::Amorphous => "Amorphous".to_string(),
    }
}

// ---- JSON output ----------------------------------------------------------

fn print_json(p: &PropertySet) {
    let obj = to_json_value(p);
    println!(
        "{}",
        serde_json::to_string_pretty(&obj).expect("JSON serialization")
    );
}

fn to_json_value(p: &PropertySet) -> serde_json::Value {
    let mut thermal = serde_json::Map::new();
    thermal.insert("tg_K".into(), opt_json(p.tg));
    thermal.insert("tm_K".into(), opt_json(p.tm));
    thermal.insert(
        "crystallization".into(),
        serde_json::Value::String(tendency_label(p.crystallization)),
    );

    let mut mechanical = serde_json::Map::new();
    mechanical.insert("density_g_cm3".into(), opt_json(p.rho));
    mechanical.insert("youngs_modulus_GPa".into(), opt_json(p.youngs));
    mechanical.insert("tensile_strength_MPa".into(), opt_json(p.tensile));

    let mut solubility = serde_json::Map::new();
    solubility.insert("hildebrand_MPa05".into(), opt_json(p.hildebrand));
    solubility.insert(
        "hansen_delta_d_MPa05".into(),
        opt_json(p.hansen.as_ref().map(|h| h.delta_d)),
    );
    solubility.insert(
        "hansen_delta_p_MPa05".into(),
        opt_json(p.hansen.as_ref().map(|h| h.delta_p)),
    );
    solubility.insert(
        "hansen_delta_h_MPa05".into(),
        opt_json(p.hansen.as_ref().map(|h| h.delta_h)),
    );

    let mut optical = serde_json::Map::new();
    optical.insert("refractive_index".into(), opt_json(p.ri));

    let mut root = serde_json::Map::new();
    root.insert("thermal".into(), serde_json::Value::Object(thermal));
    root.insert("mechanical".into(), serde_json::Value::Object(mechanical));
    root.insert("solubility".into(), serde_json::Value::Object(solubility));
    root.insert("optical".into(), serde_json::Value::Object(optical));

    serde_json::Value::Object(root)
}

fn opt_json(v: Option<f64>) -> serde_json::Value {
    match v {
        Some(x) => serde_json::Value::from(x),
        None => serde_json::Value::Null,
    }
}

// ---- CSV output -----------------------------------------------------------

fn print_csv(p: &PropertySet) {
    println!("property,value,unit,confidence");
    csv_row("tg", p.tg, "K", 3);
    csv_row("tm", p.tm, "K", 2);
    println!("crystallization,{},, 1", tendency_label(p.crystallization));
    csv_row("density", p.rho, "g/cm3", 3);
    csv_row("youngs_modulus", p.youngs, "GPa", 1);
    csv_row("tensile_strength", p.tensile, "MPa", 1);
    csv_row("hildebrand", p.hildebrand, "MPa^0.5", 2);
    csv_row(
        "hansen_delta_d",
        p.hansen.as_ref().map(|h| h.delta_d),
        "MPa^0.5",
        2,
    );
    csv_row(
        "hansen_delta_p",
        p.hansen.as_ref().map(|h| h.delta_p),
        "MPa^0.5",
        2,
    );
    csv_row(
        "hansen_delta_h",
        p.hansen.as_ref().map(|h| h.delta_h),
        "MPa^0.5",
        2,
    );
    csv_row("refractive_index", p.ri, "", 2);
}

fn csv_row(name: &str, value: Option<f64>, unit: &str, confidence: u8) {
    match value {
        Some(v) => println!("{name},{v},{unit},{confidence}"),
        None => println!("{name},,{unit},{confidence}"),
    }
}

fn report_err(e: impl std::fmt::Display) -> i32 {
    eprintln!("{} {e}", "error:".red().bold());
    1
}
