window.BENCHMARK_DATA = {
  "lastUpdate": 1772057676623,
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
      }
    ]
  }
}