//! Tests d'intégration pour le binaire `polysim`.
//!
//! Chaque test lance le binaire en subprocess via `assert_cmd`.
//! `NO_COLOR=1` est positionné sur toutes les commandes pour supprimer les
//! codes ANSI et simplifier les assertions sur le contenu texte.

use assert_cmd::Command;
use predicates::str::contains;

// ─── Helpers ────────────────────────────────────────────────────────────────

/// Retourne une commande `polysim` avec les couleurs désactivées.
///
/// Utilise `CARGO_BIN_EXE_polysim` (variable définie par cargo à la
/// compilation) pour localiser le binaire sans passer par `cargo_bin()`
/// (déprécié depuis assert_cmd 2.1).
fn polysim() -> Command {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_polysim"));
    cmd.env("NO_COLOR", "1");
    cmd
}

// ─── Help & version ──────────────────────────────────────────────────────────

#[test]
fn help_flag_exits_ok() {
    polysim().arg("--help").assert().success();
}

#[test]
fn help_mentions_analyze_subcommand() {
    polysim()
        .arg("--help")
        .assert()
        .success()
        .stdout(contains("analyze"));
}

#[test]
fn analyze_help_shows_strategy_flags() {
    polysim()
        .args(["analyze", "--help"])
        .assert()
        .success()
        .stdout(contains("--by-repeat"))
        .stdout(contains("--by-mn"))
        .stdout(contains("--by-mass"));
}

#[test]
fn version_flag_exits_ok() {
    polysim().arg("--version").assert().success();
}

#[test]
fn version_flag_shows_version_number() {
    polysim()
        .arg("--version")
        .assert()
        .success()
        // Le numéro de version commence toujours par "0."
        .stdout(contains("0."));
}

// ─── analyze — polyéthylène (--by-repeat) ────────────────────────────────────

#[test]
fn analyze_pe_by_repeat_exits_ok() {
    polysim()
        .args(["analyze", "{[]CC[]}", "--by-repeat", "10"])
        .assert()
        .success();
}

#[test]
fn analyze_pe_by_repeat_shows_correct_repeat_count() {
    polysim()
        .args(["analyze", "{[]CC[]}", "--by-repeat", "10"])
        .assert()
        .success()
        .stdout(contains("10"));
}

#[test]
fn analyze_pe_by_repeat_shows_mn() {
    // PE n=10 : C₂₀H₄₂, Mn ≈ 282.55 g/mol
    polysim()
        .args(["analyze", "{[]CC[]}", "--by-repeat", "10"])
        .assert()
        .success()
        .stdout(contains("282."));
}

#[test]
fn analyze_pe_by_repeat_shows_molecular_formula_with_subscripts() {
    // La formule doit être affichée avec des indices Unicode
    polysim()
        .args(["analyze", "{[]CC[]}", "--by-repeat", "10"])
        .assert()
        .success()
        .stdout(contains("C₂₀H₄₂"));
}

#[test]
fn analyze_pe_by_repeat_shows_total_atoms() {
    // C₂₀H₄₂ → 62 atomes
    polysim()
        .args(["analyze", "{[]CC[]}", "--by-repeat", "10"])
        .assert()
        .success()
        .stdout(contains("62"));
}

#[test]
fn analyze_pe_by_repeat_shows_monoisotopic_mass() {
    // C₂₀H₄₂ monoisotopique ≈ 282.329 g/mol
    polysim()
        .args(["analyze", "{[]CC[]}", "--by-repeat", "10"])
        .assert()
        .success()
        .stdout(contains("282.3"));
}

#[test]
fn analyze_pe_by_repeat_shows_dispersity_one() {
    polysim()
        .args(["analyze", "{[]CC[]}", "--by-repeat", "10"])
        .assert()
        .success()
        .stdout(contains("1.000"));
}

#[test]
fn analyze_pe_by_repeat_shows_single_chain_footnote() {
    polysim()
        .args(["analyze", "{[]CC[]}", "--by-repeat", "10"])
        .assert()
        .success()
        .stdout(contains("Single ideal chain"));
}

// ─── analyze — polypropylène (--by-repeat) ───────────────────────────────────

#[test]
fn analyze_pp_by_repeat_formula() {
    // PP n=3 → C₉H₂₀
    polysim()
        .args(["analyze", "{[]CC(C)[]}", "--by-repeat", "3"])
        .assert()
        .success()
        .stdout(contains("C₉H₂₀"));
}

#[test]
fn analyze_pp_by_repeat_atom_count() {
    // C₉H₂₀ → 29 atomes
    polysim()
        .args(["analyze", "{[]CC(C)[]}", "--by-repeat", "3"])
        .assert()
        .success()
        .stdout(contains("29"));
}

// ─── analyze — polystyrène (--by-repeat) ─────────────────────────────────────

#[test]
fn analyze_ps_by_repeat_formula() {
    // PS n=1 → C₈H₁₀
    polysim()
        .args(["analyze", "{[]CC(c1ccccc1)[]}", "--by-repeat", "1"])
        .assert()
        .success()
        .stdout(contains("C₈H₁₀"));
}

// ─── analyze — stratégie --by-mn ─────────────────────────────────────────────

#[test]
fn analyze_pe_by_mn_exits_ok() {
    // Cible Mn ≈ 282.554 → doit donner n=10
    polysim()
        .args(["analyze", "{[]CC[]}", "--by-mn", "282.554"])
        .assert()
        .success();
}

#[test]
fn analyze_pe_by_mn_resolves_correct_n() {
    // ByTargetMn(282.554) → n=10, formule C₂₀H₄₂
    polysim()
        .args(["analyze", "{[]CC[]}", "--by-mn", "282.554"])
        .assert()
        .success()
        .stdout(contains("C₂₀H₄₂"));
}

#[test]
fn analyze_ps_by_mn_large_target() {
    // PS Mn ≈ 5000 → n≈48, formule attendue commence par C₃
    polysim()
        .args(["analyze", "{[]CC(c1ccccc1)[]}", "--by-mn", "5000"])
        .assert()
        .success()
        // La formule de départ doit contenir un nombre à 3 chiffres de C
        .stdout(contains("C₃"));
}

// ─── analyze — stratégie --by-mass ───────────────────────────────────────────

#[test]
fn analyze_pe_by_mass_exits_ok() {
    // Masse monoisotopique de C₂₀H₄₂ ≈ 282.329 g/mol → doit donner n=10
    polysim()
        .args(["analyze", "{[]CC[]}", "--by-mass", "282.329"])
        .assert()
        .success();
}

#[test]
fn analyze_pe_by_mass_resolves_correct_formula() {
    polysim()
        .args(["analyze", "{[]CC[]}", "--by-mass", "282.329"])
        .assert()
        .success()
        .stdout(contains("C₂₀H₄₂"));
}

// ─── Gestion des erreurs ──────────────────────────────────────────────────────

#[test]
fn analyze_invalid_bigsmiles_exits_failure() {
    polysim()
        .args(["analyze", "not_a_bigsmiles", "--by-repeat", "5"])
        .assert()
        .failure();
}

#[test]
fn analyze_invalid_bigsmiles_reports_error_to_stderr() {
    polysim()
        .args(["analyze", "not_a_bigsmiles", "--by-repeat", "5"])
        .assert()
        .failure()
        .stderr(contains("error:"));
}

#[test]
fn analyze_no_strategy_flag_exits_failure() {
    // Aucune stratégie fournie → erreur clap (groupe requis)
    polysim().args(["analyze", "{[]CC[]}"]).assert().failure();
}

#[test]
fn analyze_two_strategy_flags_exits_failure() {
    // Les flags sont mutuellement exclusifs (group multiple = false)
    polysim()
        .args([
            "analyze",
            "{[]CC[]}",
            "--by-repeat",
            "10",
            "--by-mn",
            "282.0",
        ])
        .assert()
        .failure();
}

#[test]
fn analyze_repeat_count_zero_exits_failure() {
    // n=0 est rejeté par le builder
    polysim()
        .args(["analyze", "{[]CC[]}", "--by-repeat", "0"])
        .assert()
        .failure()
        .stderr(contains("error:"));
}

#[test]
fn analyze_no_stochastic_object_exits_failure() {
    // SMILES simple sans objet stochastique {…}
    polysim()
        .args(["analyze", "CCO", "--by-repeat", "5"])
        .assert()
        .failure()
        .stderr(contains("error:"));
}

#[test]
fn analyze_copolymer_bigsmiles_exits_failure() {
    // Copolymère → homopolymer() retourne une erreur (deux unités répétées)
    polysim()
        .args(["analyze", "{[$]CC[$],[$]CC(C)[$]}", "--by-repeat", "5"])
        .assert()
        .failure()
        .stderr(contains("error:"));
}

// ─── Contenu structurel de la sortie ─────────────────────────────────────────

#[test]
fn analyze_output_contains_all_property_labels() {
    polysim()
        .args(["analyze", "{[]CC[]}", "--by-repeat", "5"])
        .assert()
        .success()
        .stdout(contains("Repeat units"))
        .stdout(contains("Mn"))
        .stdout(contains("Mw"))
        .stdout(contains("Dispersity"))
        .stdout(contains("Monoisotopic mass"))
        .stdout(contains("Molecular formula"))
        .stdout(contains("Total atoms"));
}

#[test]
fn analyze_output_contains_strategy_label_by_repeat() {
    polysim()
        .args(["analyze", "{[]CC[]}", "--by-repeat", "5"])
        .assert()
        .success()
        .stdout(contains("By repeat count"));
}

#[test]
fn analyze_output_contains_strategy_label_by_mn() {
    polysim()
        .args(["analyze", "{[]CC[]}", "--by-mn", "100.0"])
        .assert()
        .success()
        .stdout(contains("By target Mn"));
}

#[test]
fn analyze_output_contains_strategy_label_by_mass() {
    polysim()
        .args(["analyze", "{[]CC[]}", "--by-mass", "100.0"])
        .assert()
        .success()
        .stdout(contains("By exact monoisotopic mass"));
}

#[test]
fn analyze_output_echoes_input_bigsmiles() {
    polysim()
        .args(["analyze", "{[]CC[]}", "--by-repeat", "1"])
        .assert()
        .success()
        .stdout(contains("{[]CC[]}"));
}

// ─── Libellé de stratégie dans la sortie ─────────────────────────────────────

#[test]
fn analyze_pe_by_repeat_strategy_label_shows_n() {
    // Le label doit contenir "n = 10" pour que la stratégie soit non-ambiguë.
    polysim()
        .args(["analyze", "{[]CC[]}", "--by-repeat", "10"])
        .assert()
        .success()
        .stdout(contains("n = 10"));
}

// ─── Lignes Δ (delta) ─────────────────────────────────────────────────────────

#[test]
fn analyze_pe_by_mn_shows_delta_mn_row() {
    polysim()
        .args(["analyze", "{[]CC[]}", "--by-mn", "200.0"])
        .assert()
        .success()
        .stdout(contains("Δ Mn"));
}

#[test]
fn analyze_pe_by_mass_shows_delta_mono_row() {
    polysim()
        .args(["analyze", "{[]CC[]}", "--by-mass", "200.0"])
        .assert()
        .success()
        .stdout(contains("Δ mono"));
}

// ─── Groupes terminaux (begin / end) ─────────────────────────────────────────

#[test]
fn analyze_output_shows_begin_and_end_when_terminal_groups_present() {
    // CC{[]CC[]}CC — "CC" en tête et en queue → lignes Begin et End visibles.
    polysim()
        .args(["analyze", "CC{[]CC[]}CC", "--by-repeat", "3"])
        .assert()
        .success()
        .stdout(contains("Begin"))
        .stdout(contains("End"));
}

// ═══════════════════════════════════════════════════════════════════════════════
// ─── generate ──────────────────────────────────────────────────────────────────
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn generate_pe_exits_ok() {
    polysim()
        .args(["generate", "{[]CC[]}", "--mn", "2805"])
        .assert()
        .success();
}

#[test]
fn generate_output_contains_ensemble_banner() {
    polysim()
        .args(["generate", "{[]CC[]}", "--mn", "2805", "--seed", "42"])
        .assert()
        .success()
        .stdout(contains("Polydisperse Ensemble"));
}

#[test]
fn generate_output_shows_mn_mw_pdi() {
    polysim()
        .args(["generate", "{[]CC[]}", "--mn", "2805", "--seed", "42"])
        .assert()
        .success()
        .stdout(contains("Mn (number-average)"))
        .stdout(contains("Mw (weight-average)"))
        .stdout(contains("Dispersity"));
}

#[test]
fn generate_output_shows_median_and_std_dev() {
    polysim()
        .args(["generate", "{[]CC[]}", "--mn", "2805", "--seed", "42"])
        .assert()
        .success()
        .stdout(contains("Mn median"))
        .stdout(contains("Mn std dev"));
}

#[test]
fn generate_output_shows_mn_range() {
    polysim()
        .args(["generate", "{[]CC[]}", "--mn", "2805", "--seed", "42"])
        .assert()
        .success()
        .stdout(contains("Mn range"));
}

#[test]
fn generate_seed_produces_reproducible_output() {
    let out1 = polysim()
        .args(["generate", "{[]CC[]}", "--mn", "2805", "--seed", "99"])
        .output()
        .unwrap();
    let out2 = polysim()
        .args(["generate", "{[]CC[]}", "--mn", "2805", "--seed", "99"])
        .output()
        .unwrap();
    assert_eq!(
        out1.stdout, out2.stdout,
        "Same seed should produce identical output"
    );
}

#[test]
fn generate_schulz_zimm_with_custom_pdi() {
    polysim()
        .args([
            "generate",
            "{[]CC[]}",
            "--mn",
            "2805",
            "--distribution",
            "schulz-zimm",
            "--pdi",
            "1.5",
            "--seed",
            "42",
        ])
        .assert()
        .success();
}

#[test]
fn generate_log_normal_distribution() {
    polysim()
        .args([
            "generate",
            "{[]CC[]}",
            "--mn",
            "2805",
            "--distribution",
            "log-normal",
            "--seed",
            "42",
        ])
        .assert()
        .success();
}

#[test]
fn generate_flory_distribution() {
    polysim()
        .args([
            "generate",
            "{[]CC[]}",
            "--mn",
            "2805",
            "--distribution",
            "flory",
            "--seed",
            "42",
        ])
        .assert()
        .success();
}

#[test]
fn generate_flory_pdi_warning_on_stderr() {
    polysim()
        .args([
            "generate",
            "{[]CC[]}",
            "--mn",
            "2805",
            "--distribution",
            "flory",
            "--pdi",
            "1.5",
            "--seed",
            "42",
        ])
        .assert()
        .success()
        .stderr(contains("warning"));
}

#[test]
fn generate_custom_num_chains() {
    polysim()
        .args([
            "generate",
            "{[]CC[]}",
            "--mn",
            "2805",
            "--num-chains",
            "50",
            "--seed",
            "42",
        ])
        .assert()
        .success()
        .stdout(contains("50"));
}

#[test]
fn generate_invalid_bigsmiles_exits_failure() {
    polysim()
        .args(["generate", "{{invalid", "--mn", "2805"])
        .assert()
        .failure();
}

#[test]
fn generate_no_stochastic_object_exits_failure() {
    polysim()
        .args(["generate", "CC", "--mn", "2805"])
        .assert()
        .failure();
}

#[test]
fn generate_missing_mn_exits_failure() {
    polysim().args(["generate", "{[]CC[]}"]).assert().failure();
}

#[test]
fn generate_echoes_input_bigsmiles() {
    polysim()
        .args(["generate", "{[]CC[]}", "--mn", "2805", "--seed", "42"])
        .assert()
        .success()
        .stdout(contains("{[]CC[]}"));
}

#[test]
fn generate_shows_distribution_label() {
    polysim()
        .args([
            "generate",
            "{[]CC[]}",
            "--mn",
            "2805",
            "--distribution",
            "schulz-zimm",
            "--seed",
            "42",
        ])
        .assert()
        .success()
        .stdout(contains("Schulz-Zimm"));
}

#[test]
fn generate_shows_target_mn() {
    polysim()
        .args(["generate", "{[]CC[]}", "--mn", "2805", "--seed", "42"])
        .assert()
        .success()
        .stdout(contains("2805"));
}
