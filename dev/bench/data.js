window.BENCHMARK_DATA = {
  "lastUpdate": 1772534189369,
  "repoUrl": "https://github.com/Peariforme/polysim",
  "entries": {
    "Polymer Builder Benchmarks": [
      {
        "commit": {
          "author": {
            "email": "47952322+Peariforme@users.noreply.github.com",
            "name": "Peariforme",
            "username": "Peariforme"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "5febffc0717777236fe50d78492ea37e629f00e2",
          "message": "fix: make benchmark workflow produce bencher-format output (#4)\n\n* fix: make benchmark workflow produce bencher-format output\n\nCriterion 0.5 removed the `--output-format bencher` CLI flag, so\n`benchmark-action/github-action-benchmark` (tool: cargo) could never\nfind the expected `test … bench: X ns/iter` lines in the output file.\n\n- Downgrade criterion from 0.5 to 0.4, which still supports\n  `--output-format bencher`\n- Switch the workflow step to `shell: bash` with `set -o pipefail` so\n  a cargo bench failure is no longer silently masked by `tee`'s exit\n  code\n\nhttps://claude.ai/code/session_01BdDTgdjMPWEbRxZmRpN7VM\n\n* fix: stop lib.rs bench target from receiving --output-format bencher\n\nRoot cause: `cargo bench -p polysim-core -- --output-format bencher`\npasses the flags to ALL bench targets, including the default `lib.rs`\ntarget which uses the libtest harness. libtest rejects `--output-format`\nas an unknown option, causing cargo bench to fail immediately. Because\nthe original step used a bare pipe (`| tee`), the non-zero exit code\nfrom cargo was masked by tee's success, leaving the output file with\nnothing but cargo's build progress; the benchmark action found no\nresults.\n\nFix:\n- Add `[lib] bench = false` in polysim-core/Cargo.toml so cargo never\n  invokes the libtest harness for lib.rs during `cargo bench`; only the\n  explicit `[[bench]]` criterion targets run\n- Restore criterion to 0.5 (was not the issue)\n- Keep the `set -o pipefail` + `shell: bash` from the previous commit\n  so any future cargo bench failure surfaces immediately\n\nhttps://claude.ai/code/session_01BdDTgdjMPWEbRxZmRpN7VM\n\n---------\n\nCo-authored-by: Claude <noreply@anthropic.com>",
          "timestamp": "2026-02-25T23:10:18+01:00",
          "tree_id": "c653c36d83126689448fe020d00c61d3d5606711",
          "url": "https://github.com/Peariforme/polysim/commit/5febffc0717777236fe50d78492ea37e629f00e2"
        },
        "date": 1772057675944,
        "tool": "cargo",
        "benches": [
          {
            "name": "homopolymer/polyethylene/10",
            "value": 1609,
            "range": "± 5",
            "unit": "ns/iter"
          },
          {
            "name": "homopolymer/polyethylene/100",
            "value": 12057,
            "range": "± 102",
            "unit": "ns/iter"
          },
          {
            "name": "homopolymer/polyethylene/1000",
            "value": 109827,
            "range": "± 618",
            "unit": "ns/iter"
          },
          {
            "name": "homopolymer/polyethylene/10000",
            "value": 1084324,
            "range": "± 10643",
            "unit": "ns/iter"
          },
          {
            "name": "homopolymer/polystyrene/10",
            "value": 8339,
            "range": "± 47",
            "unit": "ns/iter"
          },
          {
            "name": "homopolymer/polystyrene/100",
            "value": 86677,
            "range": "± 379",
            "unit": "ns/iter"
          },
          {
            "name": "homopolymer/polystyrene/1000",
            "value": 850850,
            "range": "± 5597",
            "unit": "ns/iter"
          },
          {
            "name": "molecular_weight/average_mass/polyethylene/10",
            "value": 1266,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "molecular_weight/average_mass/polyethylene/100",
            "value": 9756,
            "range": "± 59",
            "unit": "ns/iter"
          },
          {
            "name": "molecular_weight/average_mass/polyethylene/1000",
            "value": 87250,
            "range": "± 545",
            "unit": "ns/iter"
          },
          {
            "name": "molecular_weight/average_mass/polystyrene/10",
            "value": 7488,
            "range": "± 69",
            "unit": "ns/iter"
          },
          {
            "name": "molecular_weight/average_mass/polystyrene/100",
            "value": 72318,
            "range": "± 291",
            "unit": "ns/iter"
          },
          {
            "name": "molecular_weight/monoisotopic_mass/polyethylene/10",
            "value": 1405,
            "range": "± 11",
            "unit": "ns/iter"
          },
          {
            "name": "molecular_weight/monoisotopic_mass/polyethylene/100",
            "value": 10109,
            "range": "± 148",
            "unit": "ns/iter"
          },
          {
            "name": "molecular_weight/monoisotopic_mass/polyethylene/1000",
            "value": 88217,
            "range": "± 643",
            "unit": "ns/iter"
          },
          {
            "name": "molecular_weight/by_target_mn/polyethylene/282",
            "value": 2870,
            "range": "± 22",
            "unit": "ns/iter"
          },
          {
            "name": "molecular_weight/by_target_mn/polyethylene/2825",
            "value": 13313,
            "range": "± 159",
            "unit": "ns/iter"
          },
          {
            "name": "molecular_weight/by_target_mn/polyethylene/28255",
            "value": 111042,
            "range": "± 467",
            "unit": "ns/iter"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "47952322+Peariforme@users.noreply.github.com",
            "name": "Peariforme",
            "username": "Peariforme"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "53344d7848230ac187ddeee3cb61eb64b4fb7257",
          "message": "feat(cli): add polysim analyze command — V1 polymer chain analysis (#5)\n\n* feat(cli): add polysim analyze command — V1 polymer chain analysis\n\n- Add `polysim analyze <BIGSMILES> --by-repeat N | --by-mn MN | --by-mass MASS`\n  subcommand that generates a single ideal homopolymer chain and prints its\n  properties in a modern, coloured terminal table.\n\nProperties reported:\n  - Repeat units (n)\n  - Mn / Mw (= Mn for a single ideal chain, Đ = 1.000, clearly footnoted)\n  - Monoisotopic mass\n  - Molecular formula in Hill notation with Unicode subscripts (e.g. C₂₀H₄₂)\n  - Total atom count\n\npolysim-core additions:\n  - `properties::formula` module: `molecular_formula()` (Hill notation) and\n    `total_atom_count()`, both with doc-tests.\n\nCLI dependencies: clap 4 (derive), colored 2, comfy-table 7.\n\nCI/CD:\n  - `.github/workflows/dist.yml`: cross-platform release workflow triggered on\n    version tags; builds polysim for linux-x86_64, macos-arm64, macos-x86_64,\n    windows-x86_64, packages archives with SHA-256 checksums, and publishes a\n    GitHub Release via softprops/action-gh-release.\n\nhttps://claude.ai/code/session_016A62v72jssAzr24cNLGyNn\n\n* test(cli+core): add full test suite for CLI and formula module\n\npolysim-core / tests/formula.rs  (20 tests):\n  - molecular_formula: PE/PP/PS formulas, Hill notation order (C first, H second)\n  - algebraic rules: PE → C{2n}H{4n+2}, PP → C{3n}H{6n+2} verified for all n\n  - total_atom_count: PE/PP/PS values, linearity in n, consistency with formula\n\npolysim-cli / tests/cli.rs  (33 integration tests via assert_cmd):\n  - help / version flags\n  - analyze with all three strategies (--by-repeat, --by-mn, --by-mass)\n  - PE / PP / PS: Mn, formula with Unicode subscripts, total atoms, dispersity\n  - error paths: invalid BigSMILES, missing strategy, conflicting flags,\n    n=0, no stochastic object, copolymer → homopolymer mismatch\n  - output structure: all property labels, strategy echo, BigSMILES echo\n  Use CARGO_BIN_EXE_polysim (compile-time path) instead of deprecated cargo_bin().\n\npolysim-cli / src/main.rs  (10 unit tests in #[cfg(test)]):\n  - subscript_digits: all ten digits, typical formula, heteroatoms, no-op\n  - truncate: short/exact/long strings, ellipsis presence, start/end preserved\n\nhttps://claude.ai/code/session_016A62v72jssAzr24cNLGyNn\n\n* cli: add Δ rows for --by-mn/--by-mass and display begin/end blocks\n\n- For --by-mn and --by-mass, a \"Δ achieved − target\" row is inserted in\n  the table immediately below the targeted property (Mn or monoisotopic\n  mass).  The delta is coloured green when |Δ/target| < 0.5 %, yellow\n  otherwise, making quantization error immediately visible.\n\n- BigSMILES strings that contain SMILES fragments before/after the\n  stochastic object (initiator / terminator blocks, e.g. CC{[$]CC[$]}CC)\n  now show those fragments as \"Begin\" and \"End\" lines in the input\n  summary section.\n\n- Add bigsmiles as a direct dependency of polysim-cli (needed to match\n  on BigSmilesSegment in the helper functions).\n\nhttps://claude.ai/code/session_016A62v72jssAzr24cNLGyNn\n\n* chore: update Cargo.lock for bigsmiles dep in polysim-cli\n\nhttps://claude.ai/code/session_016A62v72jssAzr24cNLGyNn\n\n* refactor(cli): split monolithic main.rs into focused modules\n\nmain.rs was a single 385-line file mixing CLI parsing, business logic,\nterminal rendering, string utilities, and unit tests.  Splitting it into\ndedicated modules makes each responsibility independently readable,\ntestable, and extensible.\n\nNew layout:\n  src/main.rs                — CLI structs (Cli / Commands / StrategyArgs)\n                               + fn main(), nothing else\n  src/commands/analyze.rs    — run() orchestration: parse → build → report\n                               + AnalysisResult view-model struct\n  src/display.rs             — all terminal output (banner, summary, table,\n                               footnote); depends only on AnalysisResult\n  src/utils/format.rs        — pure functions: subscript_digits, truncate,\n                               delta_style; fully unit-tested\n  src/utils/bigsmiles_ext.rs — before/after_stochastic helpers; unit-tested\n                               with four cases (none / begin / end / both)\n\nAdditional changes:\n  - Install git hooks (make install-hooks) so fmt+clippy+tests run\n    automatically before every commit and push\n  - cargo fmt applied across the workspace\n  - delta_style now uses (!s.is_empty()).then_some(s) as per clippy idiom\n  - Four new delta_style unit tests added\n\nhttps://claude.ai/code/session_016A62v72jssAzr24cNLGyNn\n\n* chore: apply cargo fmt to cli tests and mark hooks as executable\n\n- cargo fmt reformatted two test functions in tests/cli.rs\n  (collapsed short .args([...]) chains to single lines)\n- make install-hooks set the execute bit on the three .githooks scripts\n\nhttps://claude.ai/code/session_016A62v72jssAzr24cNLGyNn\n\n* refactor: use bigsmiles 0.1.2 prefix/suffix_segments API\n\nReplace the local bigsmiles_ext module with the new methods\nBigSmiles::prefix_segments() and ::suffix_segments() added upstream\nin bigsmiles 0.1.2. Removes the duplicate logic and the workaround\nmodule entirely.\n\nhttps://claude.ai/code/session_016A62v72jssAzr24cNLGyNn\n\n* refactor(cli): code review fixes\n\n- display: fix banner len -> chars().count() (em dash alignment)\n- format: fix truncate byte-slicing -> char-based to avoid panic on multi-byte\n- format: delta_style uses > 0.0 so exact zero shows 0.000 without spurious +\n- report: extract AnalysisResult into src/report.rs, decouple display from analyze\n- main: add build_strategy()/label() methods on StrategyArgs; pub -> pub(crate)\n- analyze: extract fn report_err(), remove resolve_strategy, use StrategyArgs methods\n- display: fix number-average Mw -> number-average; remove double space before D\n- tests: add delta Mn, delta mono, begin/end, strategy label tests (37 total)\n\n* style: fmt\n\n---------\n\nCo-authored-by: Claude <noreply@anthropic.com>",
          "timestamp": "2026-03-03T11:32:07+01:00",
          "tree_id": "e52ad6526a5a782e881f6c0259217f3005184232",
          "url": "https://github.com/Peariforme/polysim/commit/53344d7848230ac187ddeee3cb61eb64b4fb7257"
        },
        "date": 1772534188479,
        "tool": "cargo",
        "benches": [
          {
            "name": "homopolymer/polyethylene/10",
            "value": 1635,
            "range": "± 28",
            "unit": "ns/iter"
          },
          {
            "name": "homopolymer/polyethylene/100",
            "value": 11740,
            "range": "± 248",
            "unit": "ns/iter"
          },
          {
            "name": "homopolymer/polyethylene/1000",
            "value": 109076,
            "range": "± 691",
            "unit": "ns/iter"
          },
          {
            "name": "homopolymer/polyethylene/10000",
            "value": 1064314,
            "range": "± 21247",
            "unit": "ns/iter"
          },
          {
            "name": "homopolymer/polystyrene/10",
            "value": 8423,
            "range": "± 50",
            "unit": "ns/iter"
          },
          {
            "name": "homopolymer/polystyrene/100",
            "value": 87760,
            "range": "± 481",
            "unit": "ns/iter"
          },
          {
            "name": "homopolymer/polystyrene/1000",
            "value": 862949,
            "range": "± 2161",
            "unit": "ns/iter"
          },
          {
            "name": "molecular_weight/average_mass/polyethylene/10",
            "value": 1297,
            "range": "± 8",
            "unit": "ns/iter"
          },
          {
            "name": "molecular_weight/average_mass/polyethylene/100",
            "value": 10149,
            "range": "± 32",
            "unit": "ns/iter"
          },
          {
            "name": "molecular_weight/average_mass/polyethylene/1000",
            "value": 91318,
            "range": "± 646",
            "unit": "ns/iter"
          },
          {
            "name": "molecular_weight/average_mass/polystyrene/10",
            "value": 7622,
            "range": "± 20",
            "unit": "ns/iter"
          },
          {
            "name": "molecular_weight/average_mass/polystyrene/100",
            "value": 74625,
            "range": "± 245",
            "unit": "ns/iter"
          },
          {
            "name": "molecular_weight/monoisotopic_mass/polyethylene/10",
            "value": 1434,
            "range": "± 10",
            "unit": "ns/iter"
          },
          {
            "name": "molecular_weight/monoisotopic_mass/polyethylene/100",
            "value": 10519,
            "range": "± 114",
            "unit": "ns/iter"
          },
          {
            "name": "molecular_weight/monoisotopic_mass/polyethylene/1000",
            "value": 89807,
            "range": "± 271",
            "unit": "ns/iter"
          },
          {
            "name": "molecular_weight/by_target_mn/polyethylene/282",
            "value": 2869,
            "range": "± 7",
            "unit": "ns/iter"
          },
          {
            "name": "molecular_weight/by_target_mn/polyethylene/2825",
            "value": 13598,
            "range": "± 39",
            "unit": "ns/iter"
          },
          {
            "name": "molecular_weight/by_target_mn/polyethylene/28255",
            "value": 114665,
            "range": "± 373",
            "unit": "ns/iter"
          }
        ]
      }
    ]
  }
}