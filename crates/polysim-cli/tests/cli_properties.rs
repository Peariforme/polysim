//! Tests d'intégration pour la sous-commande `polysim properties`.
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

#[test]
fn properties_ps_table_output_exits_ok() {
    polysim()
        .args(["properties", "{[]CC(c1ccccc1)[]}", "--by-repeat", "100"])
        .assert()
        .success();
}

#[test]
fn properties_ps_table_contains_thermal_section() {
    polysim()
        .args(["properties", "{[]CC(c1ccccc1)[]}", "--by-repeat", "100"])
        .assert()
        .success()
        .stdout(contains("THERMAL"));
}

#[test]
fn properties_ps_table_contains_tg_row() {
    polysim()
        .args(["properties", "{[]CC(c1ccccc1)[]}", "--by-repeat", "100"])
        .assert()
        .success()
        .stdout(contains("Tg"));
}

#[test]
fn properties_ps_table_contains_hildebrand_row() {
    polysim()
        .args(["properties", "{[]CC(c1ccccc1)[]}", "--by-repeat", "100"])
        .assert()
        .success()
        .stdout(contains("Hildebrand"));
}

#[test]
fn properties_ps_table_contains_mechanical_section() {
    polysim()
        .args(["properties", "{[]CC(c1ccccc1)[]}", "--by-repeat", "100"])
        .assert()
        .success()
        .stdout(contains("MECHANICAL"));
}

#[test]
fn properties_ps_table_contains_solubility_section() {
    polysim()
        .args(["properties", "{[]CC(c1ccccc1)[]}", "--by-repeat", "100"])
        .assert()
        .success()
        .stdout(contains("SOLUBILITY"));
}

#[test]
fn properties_ps_table_contains_optical_section() {
    polysim()
        .args(["properties", "{[]CC(c1ccccc1)[]}", "--by-repeat", "100"])
        .assert()
        .success()
        .stdout(contains("OPTICAL"));
}

#[test]
fn properties_table_echoes_bigsmiles() {
    polysim()
        .args(["properties", "{[]CC[]}", "--by-repeat", "50"])
        .assert()
        .success()
        .stdout(contains("{[]CC[]}"));
}

// ─── Sortie JSON ──────────────────────────────────────────────────────────────

#[test]
fn properties_pe_json_output_is_valid_json() {
    let output = polysim()
        .args([
            "properties",
            "{[]CC[]}",
            "--by-repeat",
            "50",
            "--format",
            "json",
        ])
        .assert()
        .success()
        .get_output()
        .clone();
    let _: Value = serde_json::from_slice(&output.stdout).expect("sortie doit etre du JSON valide");
}

#[test]
fn properties_pe_json_has_thermal_section() {
    let output = polysim()
        .args([
            "properties",
            "{[]CC[]}",
            "--by-repeat",
            "50",
            "--format",
            "json",
        ])
        .assert()
        .success()
        .get_output()
        .clone();
    let json: Value = serde_json::from_slice(&output.stdout).unwrap();
    assert!(json.get("thermal").is_some(), "JSON doit avoir 'thermal'");
}

#[test]
fn properties_pe_json_has_mechanical_section() {
    let output = polysim()
        .args([
            "properties",
            "{[]CC[]}",
            "--by-repeat",
            "50",
            "--format",
            "json",
        ])
        .assert()
        .success()
        .get_output()
        .clone();
    let json: Value = serde_json::from_slice(&output.stdout).unwrap();
    assert!(
        json.get("mechanical").is_some(),
        "JSON doit avoir 'mechanical'"
    );
}

#[test]
fn properties_pe_json_has_solubility_section() {
    let output = polysim()
        .args([
            "properties",
            "{[]CC[]}",
            "--by-repeat",
            "50",
            "--format",
            "json",
        ])
        .assert()
        .success()
        .get_output()
        .clone();
    let json: Value = serde_json::from_slice(&output.stdout).unwrap();
    assert!(
        json.get("solubility").is_some(),
        "JSON doit avoir 'solubility'"
    );
}

#[test]
fn properties_pe_json_has_optical_section() {
    let output = polysim()
        .args([
            "properties",
            "{[]CC[]}",
            "--by-repeat",
            "50",
            "--format",
            "json",
        ])
        .assert()
        .success()
        .get_output()
        .clone();
    let json: Value = serde_json::from_slice(&output.stdout).unwrap();
    assert!(json.get("optical").is_some(), "JSON doit avoir 'optical'");
}

#[test]
fn properties_pe_json_tg_is_positive_number() {
    let output = polysim()
        .args([
            "properties",
            "{[]CC[]}",
            "--by-repeat",
            "100",
            "--format",
            "json",
        ])
        .assert()
        .success()
        .get_output()
        .clone();
    let json: Value = serde_json::from_slice(&output.stdout).unwrap();
    let tg = json["thermal"]["tg_K"]
        .as_f64()
        .expect("tg_K doit etre un nombre");
    assert!(tg > 0.0, "Tg doit etre positif, obtenu: {tg}");
    assert!(tg < 1000.0, "Tg doit etre < 1000 K, obtenu: {tg}");
}

#[test]
fn properties_pe_json_tg_lt_tm_when_some() {
    let output = polysim()
        .args([
            "properties",
            "{[]CC[]}",
            "--by-repeat",
            "100",
            "--format",
            "json",
        ])
        .assert()
        .success()
        .get_output()
        .clone();
    let json: Value = serde_json::from_slice(&output.stdout).unwrap();
    if let (Some(tg), Some(tm)) = (
        json["thermal"]["tg_K"].as_f64(),
        json["thermal"]["tm_K"].as_f64(),
    ) {
        assert!(tg < tm, "Tg ({tg:.1} K) doit etre < Tm ({tm:.1} K) pour PE");
    }
}

#[test]
fn properties_pe_json_hildebrand_positive() {
    let output = polysim()
        .args([
            "properties",
            "{[]CC[]}",
            "--by-repeat",
            "100",
            "--format",
            "json",
        ])
        .assert()
        .success()
        .get_output()
        .clone();
    let json: Value = serde_json::from_slice(&output.stdout).unwrap();
    if let Some(hild) = json["solubility"]["hildebrand_MPa05"].as_f64() {
        assert!(hild > 0.0, "Hildebrand doit etre positif, obtenu: {hild}");
        assert!(
            hild < 50.0,
            "Hildebrand doit etre < 50 MPa^0.5, obtenu: {hild}"
        );
    }
}

#[test]
fn properties_pe_json_density_physical_range() {
    let output = polysim()
        .args([
            "properties",
            "{[]CC[]}",
            "--by-repeat",
            "100",
            "--format",
            "json",
        ])
        .assert()
        .success()
        .get_output()
        .clone();
    let json: Value = serde_json::from_slice(&output.stdout).unwrap();
    if let Some(rho) = json["mechanical"]["density_g_cm3"].as_f64() {
        assert!(
            rho > 0.5 && rho < 2.5,
            "Densite PE doit etre dans [0.5, 2.5], obtenu: {rho}"
        );
    }
}

// ─── Sortie CSV ───────────────────────────────────────────────────────────────

#[test]
fn properties_pe_csv_output_has_header() {
    polysim()
        .args([
            "properties",
            "{[]CC[]}",
            "--by-repeat",
            "50",
            "--format",
            "csv",
        ])
        .assert()
        .success()
        .stdout(contains("property,value,unit,confidence"));
}

#[test]
fn properties_pe_csv_has_tg_row() {
    polysim()
        .args([
            "properties",
            "{[]CC[]}",
            "--by-repeat",
            "50",
            "--format",
            "csv",
        ])
        .assert()
        .success()
        .stdout(contains("tg,"));
}

#[test]
fn properties_pe_csv_has_hildebrand_row() {
    polysim()
        .args([
            "properties",
            "{[]CC[]}",
            "--by-repeat",
            "50",
            "--format",
            "csv",
        ])
        .assert()
        .success()
        .stdout(contains("hildebrand,"));
}

// ─── Gestion d'erreurs ────────────────────────────────────────────────────────

#[test]
fn properties_invalid_smiles_returns_error() {
    polysim()
        .args(["properties", "NOT_VALID_SMILES", "--by-repeat", "10"])
        .assert()
        .failure();
}

#[test]
fn properties_invalid_smiles_reports_error_to_stderr() {
    polysim()
        .args(["properties", "NOT_VALID_SMILES", "--by-repeat", "10"])
        .assert()
        .failure()
        .stderr(contains("error:"));
}

#[test]
fn properties_no_strategy_flag_exits_failure() {
    polysim()
        .args(["properties", "{[]CC[]}"])
        .assert()
        .failure();
}

#[test]
fn properties_repeat_count_zero_exits_failure() {
    polysim()
        .args(["properties", "{[]CC[]}", "--by-repeat", "0"])
        .assert()
        .failure()
        .stderr(contains("error:"));
}

// ─── Sous-commande enregistrée dans l'aide globale ────────────────────────────

#[test]
fn help_mentions_properties_subcommand() {
    polysim()
        .arg("--help")
        .assert()
        .success()
        .stdout(contains("properties"));
}

// ─── Test PMMA (polymère avec ester) ─────────────────────────────────────────

#[test]
fn properties_pmma_json_has_tg() {
    // PMMA : {[]CC(C)(C(=O)OC)[]}
    let output = polysim()
        .args([
            "properties",
            "{[]CC(C)(C(=O)OC)[]}",
            "--by-repeat",
            "50",
            "--format",
            "json",
        ])
        .assert()
        .success()
        .get_output()
        .clone();
    let json: Value = serde_json::from_slice(&output.stdout).unwrap();
    let tg = json["thermal"]["tg_K"]
        .as_f64()
        .expect("tg_K doit etre present pour PMMA");
    // Tg VK ≈ 347 K (exp 378 K), tolérance ±40 K
    assert!(
        tg > 300.0 && tg < 420.0,
        "Tg PMMA doit etre dans [300, 420] K, obtenu: {tg:.1}"
    );
}

#[test]
fn properties_pmma_json_crystallization_is_amorphous() {
    // PMMA : Tg VK > Tm VK → cristallisation = Amorphous
    let output = polysim()
        .args([
            "properties",
            "{[]CC(C)(C(=O)OC)[]}",
            "--by-repeat",
            "50",
            "--format",
            "json",
        ])
        .assert()
        .success()
        .get_output()
        .clone();
    let json: Value = serde_json::from_slice(&output.stdout).unwrap();
    let cryst = json["thermal"]["crystallization"]
        .as_str()
        .expect("crystallization doit etre une chaine");
    assert_eq!(
        cryst, "Amorphous",
        "PMMA doit etre classe Amorphous, obtenu: {cryst}"
    );
}

// ─── Stratégie --by-mn ────────────────────────────────────────────────────────

#[test]
fn properties_pe_by_mn_exits_ok() {
    polysim()
        .args(["properties", "{[]CC[]}", "--by-mn", "5000"])
        .assert()
        .success();
}

// ─── Architecture -- test de la sortie cohérente ────────────────────────────

#[test]
fn properties_ps_json_crystallization_is_string() {
    let output = polysim()
        .args([
            "properties",
            "{[]CC(c1ccccc1)[]}",
            "--by-repeat",
            "100",
            "--format",
            "json",
        ])
        .assert()
        .success()
        .get_output()
        .clone();
    let json: Value = serde_json::from_slice(&output.stdout).unwrap();
    assert!(
        json["thermal"]["crystallization"].is_string(),
        "crystallization doit etre une chaine de caracteres"
    );
    let cryst = json["thermal"]["crystallization"].as_str().unwrap();
    assert!(
        ["High", "Medium", "Low", "Amorphous"].contains(&cryst),
        "crystallization doit etre l'un des 4 etats, obtenu: {cryst}"
    );
}
