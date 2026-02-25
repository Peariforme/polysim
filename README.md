# polysim

[![CI](https://github.com/Peariforme/polysim/actions/workflows/ci.yml/badge.svg)](https://github.com/Peariforme/polysim/actions/workflows/ci.yml)
[![Crates.io](https://img.shields.io/crates/v/polysim-core.svg)](https://crates.io/crates/polysim-core)
[![docs.rs](https://docs.rs/polysim-core/badge.svg)](https://docs.rs/polysim-core)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

Polymer structure generator and physical property simulator written in Rust.

Given a [BigSMILES](https://olsenlabmit.github.io/BigSMILES/docs/line_notation.html) string,
polysim generates concrete polymer chains (as SMILES) and computes physical/chemical properties
such as glass transition temperature, crystallisation tendency, molecular weight, and more.

---

## Features

| Status | Feature |
|--------|---------|
| ✅ | Linear homopolymer generation |
| 🔜 | Random / alternating / block copolymers |
| 🔜 | Branched polymers, graft copolymers, macromonomers |
| ✅ | Chain length by repeat count |
| 🔜 | Chain length by target Mn or exact mass |
| ✅ | Tg estimation — Fox equation |
| 🔜 | Tg estimation — Van Krevelen group contributions |
| 🔜 | Crystallisation tendency |
| 🔜 | Monoisotopic mass & average molecular weight |
| 🔜 | Hildebrand solubility parameter |
| 🔜 | Melting temperature Tm |

---

## Quick start

```toml
# Cargo.toml
[dependencies]
polysim-core = "0.1"
```

```rust
use polysim_core::{parse, builder::{linear::LinearBuilder, BuildStrategy}};

// Generate a polystyrene chain with 50 repeat units
let bs = parse("{[]CC(c1ccccc1)[]}").unwrap();
let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(50))
    .homopolymer()
    .unwrap();

println!("{}", chain.smiles);        // full SMILES string
println!("{}", chain.repeat_count);  // 50
```

```rust
use polysim_core::properties::thermal::tg_fox;

// Tg of a 70/30 PS/PMMA blend  (PS: 373 K, PMMA: 378 K)
let tg = tg_fox(&[(0.70, 373.0), (0.30, 378.0)]);
println!("Tg ≈ {tg:.1} K");
```

---

## Polymer architectures

| Architecture | Builder | Key parameters |
|---|---|---|
| Homopolymer | `LinearBuilder::homopolymer` | repeat count or target Mn |
| Random copolymer | `LinearBuilder::random_copolymer` | weight fraction per monomer |
| Alternating copolymer | `LinearBuilder::alternating_copolymer` | exactly 2 repeat units |
| Block copolymer | `LinearBuilder::block_copolymer` | length of each block |
| Comb / branched | `BranchedBuilder::comb_polymer` | branch frequency |
| Graft copolymer | `BranchedBuilder::graft_copolymer` | graft fraction |
| Macromonomer | `BranchedBuilder::macromonomer` | side chain + end group |

### Build strategies

```rust
pub enum BuildStrategy {
    ByRepeatCount(usize),  // exact number of repeat units
    ByTargetMn(f64),       // target number-average Mn in g/mol  (🔜)
    ByExactMass(f64),      // target monoisotopic mass in g/mol  (🔜)
}
```

---

## Workspace layout

```
polysim/
├── crates/
│   ├── polysim-core/     # library (published to crates.io)
│   │   └── src/
│   │       ├── builder/      # chain generators
│   │       ├── polymer/      # PolymerChain type
│   │       └── properties/   # Tg, MW, …
│   └── polysim-cli/      # command-line tool (not yet published)
└── tests/
    └── integration.rs
```

---

## Development

### Prerequisites

- Rust stable (≥ 1.70)
- `make` (optional, for convenience targets)

### Setup

```bash
git clone https://github.com/Peariforme/polysim
cd polysim

# Install git hooks (pre-commit: fmt+clippy, pre-push: tests)
git config core.hooksPath .githooks
# or: make install-hooks
```

### Common commands

```bash
cargo test --all-features   # run tests
cargo fmt --all             # format code
cargo clippy --all-targets --all-features -- -D warnings  # lint
cargo bench                 # run benchmarks
cargo doc --open            # browse documentation
```

---

## References

- Lin, T.-S. *et al.* (2019). BigSMILES: A Structurally-Based Line Notation for Describing
  Macromolecules. *ACS Central Science* **5**, 1523–1531.
  [doi:10.1021/acscentsci.9b00476](https://doi.org/10.1021/acscentsci.9b00476)
- Fox, T. G. (1956). Influence of diluent and of copolymer composition on the glass
  temperature of a polymer system. *Bull. Am. Phys. Soc.* **1**, 123.
- Van Krevelen, D. W. & te Nijenhuis, K. (2009). *Properties of Polymers*, 4th ed.
  Elsevier.

---

## License

MIT — see [LICENSE](LICENSE).
