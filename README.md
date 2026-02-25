# polysim

[![CI](https://github.com/Peariforme/polysim/actions/workflows/ci.yml/badge.svg)](https://github.com/Peariforme/polysim/actions/workflows/ci.yml)
[![Crates.io](https://img.shields.io/crates/v/polysim-core.svg)](https://crates.io/crates/polysim-core)
[![docs.rs](https://docs.rs/polysim-core/badge.svg)](https://docs.rs/polysim-core)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

Polymer structure generator and physical property simulator written in Rust.

Given a [BigSMILES](https://olsenlabmit.github.io/BigSMILES/docs/line_notation.html) string,
polysim generates concrete polymer chains (as SMILES) and computes physical/chemical properties
such as glass transition temperature, molecular weight, and more.

---

## Features

| Status | Feature |
|--------|---------|
| ✅ | Linear homopolymer generation |
| 🔜 | Random / alternating / block copolymers |
| 🔜 | Branched polymers, graft copolymers, macromonomers |
| ✅ | Chain length by repeat count |
| ✅ | Chain length by target Mn (`ByTargetMn`) |
| ✅ | Chain length by monoisotopic mass (`ByExactMass`) |
| ✅ | Average molecular weight (IUPAC standard atomic weights) |
| ✅ | Monoisotopic mass (most abundant isotope per element) |
| ✅ | Tg estimation — Fox equation |
| 🔜 | Tg estimation — Van Krevelen group contributions |
| 🔜 | Crystallisation tendency |
| 🔜 | Hildebrand solubility parameter |
| 🔜 | Melting temperature Tm |

---

## Quick start

```toml
# Cargo.toml
[dependencies]
polysim-core = "0.1"
```

### Build a chain and read its molecular weight

```rust
use polysim_core::{parse, builder::{linear::LinearBuilder, BuildStrategy}};

// Polyéthylène — 100 unités répétées
let bs = parse("{[]CC[]}").unwrap();
let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(100))
    .homopolymer()
    .unwrap();

println!("Repeat units : {}", chain.repeat_count); // 100
println!("Mn           : {:.1} g/mol", chain.mn);  // 1410.7 g/mol
println!("SMILES       : {}…", &chain.smiles[..20]);
```

### Target a specific Mn

```rust
use polysim_core::{parse, builder::{linear::LinearBuilder, BuildStrategy}};

// Polypropylène — viser Mn ≈ 10 000 g/mol
let bs = parse("{[]CC(C)[]}").unwrap();
let chain = LinearBuilder::new(bs, BuildStrategy::ByTargetMn(10_000.0))
    .homopolymer()
    .unwrap();

println!("Repeat units : {}", chain.repeat_count); // ≈ 237
println!("Mn réel      : {:.1} g/mol", chain.mn);  // ≈ 9 996 g/mol
```

### Compute masses independently

```rust
use polysim_core::{parse, builder::{linear::LinearBuilder, BuildStrategy},
                   properties::molecular_weight::{average_mass, monoisotopic_mass}};

let bs = parse("{[]CC(c1ccccc1)[]}").unwrap(); // polystyrène
let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(10))
    .homopolymer()
    .unwrap();

println!("Masse moyenne       : {:.2} g/mol", average_mass(&chain));
println!("Masse monoisotopique: {:.2} g/mol", monoisotopic_mass(&chain));
```

### Glass transition temperature

```rust
use polysim_core::properties::thermal::tg_fox;

// Tg d'un mélange 70/30 PS/PMMA  (PS : 373 K, PMMA : 378 K)
let tg = tg_fox(&[(0.70, 373.0), (0.30, 378.0)]);
println!("Tg ≈ {tg:.1} K");
```

---

## Build strategies

```rust
pub enum BuildStrategy {
    /// Nombre exact d'unités répétées.
    ByRepeatCount(usize),

    /// Mn cible (masse moléculaire moyenne, g/mol).
    /// Le nombre de répétitions est déduit par extrapolation linéaire.
    ByTargetMn(f64),

    /// Masse monoisotopique cible (g/mol).
    /// Même logique que ByTargetMn mais avec les masses monoisotopiques.
    ByExactMass(f64),
}
```

Après construction, `chain.mn` contient toujours la masse moléculaire moyenne calculée,
quelle que soit la stratégie utilisée.

---

## Polymer architectures

| Architecture | Builder | Statut |
|---|---|---|
| Homopolymer | `LinearBuilder::homopolymer` | ✅ |
| Random copolymer | `LinearBuilder::random_copolymer` | 🔜 |
| Alternating copolymer | `LinearBuilder::alternating_copolymer` | 🔜 |
| Block copolymer | `LinearBuilder::block_copolymer` | 🔜 |
| Comb / branched | `BranchedBuilder::comb_polymer` | 🔜 |
| Graft copolymer | `BranchedBuilder::graft_copolymer` | 🔜 |
| Macromonomer | `BranchedBuilder::macromonomer` | 🔜 |

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
