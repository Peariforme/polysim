//! Tests d'intégration pour la sous-commande `polysim compare` — US-2.7.2.
//!
//! Chaque test lance le binaire en subprocess via `assert_cmd`.
//! `NO_COLOR=1` est positionné pour supprimer les codes ANSI.

use assert_cmd::Command;
use predicates::str::contains;
use serde_json::Value;

// ─── Helper ──────────────────────────────────────────────────────────────────

fn polysim() -> Command {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_polysim"));
    cmd.env("NO_COLOR", "1");
    cmd
}

// ─── Sortie table (par défaut) ────────────────────────────────────────────────

/// La commande compare réussit avec deux polymères.
#[test]
fn compare_two_polymers_exits_ok() {
    polysim()
        .args(["compare", "{[]CC[]}", "{[]CC(C)[]}", "--by-repeat", "50"])
        .assert()
        .success();
}

/// La sortie table contient la section THERMAL.
#[test]
fn compare_table_contains_thermal_section() {
    polysim()
        .args(["compare", "{[]CC[]}", "{[]CC(C)[]}", "--by-repeat", "50"])
        .assert()
        .success()
        .stdout(contains("THERMAL"));
}

/// La sortie table contient la section MECHANICAL.
#[test]
fn compare_table_contains_mechanical_section() {
    polysim()
        .args(["compare", "{[]CC[]}", "{[]CC(C)[]}", "--by-repeat", "50"])
        .assert()
        .success()
        .stdout(contains("MECHANICAL"));
}

/// La sortie table contient la section PERMEABILITY.
#[test]
fn compare_table_contains_permeability_section() {
    polysim()
        .args(["compare", "{[]CC[]}", "{[]CC(C)[]}", "--by-repeat", "50"])
        .assert()
        .success()
        .stdout(contains("PERMEABILITY"));
}

/// La sortie table contient Tg (K).
#[test]
fn compare_table_contains_tg_row() {
    polysim()
        .args(["compare", "{[]CC[]}", "{[]CC(C)[]}", "--by-repeat", "50"])
        .assert()
        .success()
        .stdout(contains("Tg (K)"));
}

/// Trois polymères simultanément — doit réussir.
#[test]
fn compare_three_polymers_exits_ok() {
    polysim()
        .args([
            "compare",
            "{[]CC[]}",
            "{[]CC(C)[]}",
            "{[]CC(c1ccccc1)[]}",
            "--by-repeat",
            "20",
        ])
        .assert()
        .success();
}

/// Un seul polymère — doit échouer (au moins 2 requis).
#[test]
fn compare_single_polymer_fails() {
    polysim()
        .args(["compare", "{[]CC[]}", "--by-repeat", "50"])
        .assert()
        .failure();
}

// ─── Sortie JSON ──────────────────────────────────────────────────────────────

/// La sortie JSON est un objet valide.
#[test]
fn compare_json_is_valid() {
    let output = polysim()
        .args([
            "compare",
            "{[]CC[]}",
            "{[]CC(C)[]}",
            "--by-repeat",
            "50",
            "--format",
            "json",
        ])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();

    let text = String::from_utf8(output).expect("UTF-8");
    let v: Value = serde_json::from_str(&text).expect("JSON valide");
    assert!(v.is_object(), "La racine doit être un objet JSON");
}

/// La sortie JSON contient les deux polymères comme clés racine.
#[test]
fn compare_json_contains_both_polymers() {
    let output = polysim()
        .args([
            "compare",
            "{[]CC[]}",
            "{[]CC(C)[]}",
            "--by-repeat",
            "50",
            "--format",
            "json",
        ])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();

    let text = String::from_utf8(output).expect("UTF-8");
    let v: Value = serde_json::from_str(&text).expect("JSON valide");
    assert!(v.get("{[]CC[]}").is_some(), "PE doit être dans le JSON");
    assert!(v.get("{[]CC(C)[]}").is_some(), "PP doit être dans le JSON");
}

/// Le JSON contient les clés thermiques attendues.
#[test]
fn compare_json_thermal_keys_present() {
    let output = polysim()
        .args([
            "compare",
            "{[]CC[]}",
            "{[]CC(C)[]}",
            "--by-repeat",
            "50",
            "--format",
            "json",
        ])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();

    let text = String::from_utf8(output).expect("UTF-8");
    let v: Value = serde_json::from_str(&text).expect("JSON valide");
    let pe = &v["{[]CC[]}"];
    assert!(
        pe["thermal"]["tg_K"].is_number(),
        "tg_K doit être un nombre"
    );
    assert!(
        pe["thermal"]["crystallization"].is_string(),
        "crystallization doit être une chaîne"
    );
    assert!(
        pe["thermal"]["cp_J_per_g_K"].is_number(),
        "cp_J_per_g_K doit être un nombre"
    );
}

/// Le JSON contient les clés mécaniques attendues.
#[test]
fn compare_json_mechanical_keys_present() {
    let output = polysim()
        .args([
            "compare",
            "{[]CC[]}",
            "{[]CC(C)[]}",
            "--by-repeat",
            "50",
            "--format",
            "json",
        ])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();

    let text = String::from_utf8(output).expect("UTF-8");
    let v: Value = serde_json::from_str(&text).expect("JSON valide");
    let pe = &v["{[]CC[]}"];
    assert!(
        pe["mechanical"]["density_g_cm3"].is_number(),
        "density_g_cm3 doit être un nombre"
    );
    assert!(
        pe["mechanical"]["youngs_modulus_GPa"].is_number()
            || pe["mechanical"]["youngs_modulus_GPa"].is_null(),
        "youngs_modulus_GPa doit être un nombre ou null"
    );
}

/// Le JSON contient la clé de perméabilité O₂.
#[test]
fn compare_json_permeability_o2_present() {
    let output = polysim()
        .args([
            "compare",
            "{[]CC[]}",
            "{[]CC(C)[]}",
            "--by-repeat",
            "50",
            "--format",
            "json",
        ])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();

    let text = String::from_utf8(output).expect("UTF-8");
    let v: Value = serde_json::from_str(&text).expect("JSON valide");
    let pe = &v["{[]CC[]}"];
    assert!(
        pe["permeability"]["o2_permeability_Barrer"].is_number(),
        "o2_permeability_Barrer doit être un nombre"
    );
}

/// PE et PP ont des Tg différentes (propriétés distinctes).
#[test]
fn compare_json_pe_and_pp_have_different_tg() {
    let output = polysim()
        .args([
            "compare",
            "{[]CC[]}",
            "{[]CC(C)[]}",
            "--by-repeat",
            "50",
            "--format",
            "json",
        ])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();

    let text = String::from_utf8(output).expect("UTF-8");
    let v: Value = serde_json::from_str(&text).expect("JSON valide");
    let tg_pe = v["{[]CC[]}"]["thermal"]["tg_K"].as_f64().unwrap();
    let tg_pp = v["{[]CC(C)[]}"]["thermal"]["tg_K"].as_f64().unwrap();
    assert!(
        (tg_pe - tg_pp).abs() > 1.0,
        "PE Tg ({tg_pe:.1} K) et PP Tg ({tg_pp:.1} K) doivent différer"
    );
}

/// PE a une densité différente de PP.
#[test]
fn compare_json_pe_and_pp_have_different_density() {
    let output = polysim()
        .args([
            "compare",
            "{[]CC[]}",
            "{[]CC(C)[]}",
            "--by-repeat",
            "50",
            "--format",
            "json",
        ])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();

    let text = String::from_utf8(output).expect("UTF-8");
    let v: Value = serde_json::from_str(&text).expect("JSON valide");
    let rho_pe = v["{[]CC[]}"]["mechanical"]["density_g_cm3"]
        .as_f64()
        .unwrap();
    let rho_pp = v["{[]CC(C)[]}"]["mechanical"]["density_g_cm3"]
        .as_f64()
        .unwrap();
    assert!(
        (rho_pe - rho_pp).abs() > 0.001,
        "PE densité ({rho_pe:.3}) et PP densité ({rho_pp:.3}) doivent différer"
    );
}

// ─── Sortie CSV ───────────────────────────────────────────────────────────────

/// La sortie CSV contient l'en-tête attendue.
#[test]
fn compare_csv_contains_header() {
    polysim()
        .args([
            "compare",
            "{[]CC[]}",
            "{[]CC(C)[]}",
            "--by-repeat",
            "50",
            "--format",
            "csv",
        ])
        .assert()
        .success()
        .stdout(contains("property"));
}

/// La sortie CSV contient la ligne tg_K.
#[test]
fn compare_csv_contains_tg_row() {
    polysim()
        .args([
            "compare",
            "{[]CC[]}",
            "{[]CC(C)[]}",
            "--by-repeat",
            "50",
            "--format",
            "csv",
        ])
        .assert()
        .success()
        .stdout(contains("tg_K"));
}

/// La sortie CSV contient la ligne de densité.
#[test]
fn compare_csv_contains_density_row() {
    polysim()
        .args([
            "compare",
            "{[]CC[]}",
            "{[]CC(C)[]}",
            "--by-repeat",
            "50",
            "--format",
            "csv",
        ])
        .assert()
        .success()
        .stdout(contains("density_g_cm3"));
}

/// La sortie CSV contient la ligne de perméabilité O₂.
#[test]
fn compare_csv_contains_permeability_row() {
    polysim()
        .args([
            "compare",
            "{[]CC[]}",
            "{[]CC(C)[]}",
            "--by-repeat",
            "50",
            "--format",
            "csv",
        ])
        .assert()
        .success()
        .stdout(contains("o2_permeability_Barrer"));
}

/// La sortie CSV contient la ligne de cristallisation.
#[test]
fn compare_csv_contains_crystallization_row() {
    polysim()
        .args([
            "compare",
            "{[]CC[]}",
            "{[]CC(C)[]}",
            "--by-repeat",
            "50",
            "--format",
            "csv",
        ])
        .assert()
        .success()
        .stdout(contains("crystallization"));
}
