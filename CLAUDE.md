# CLAUDE.md — polysim

## Ce que nous construisons

Un simulateur de polymères en Rust qui :
1. **Génère** des structures de chaînes polymères à partir d'une notation BigSMILES
2. **Calcule** des propriétés physico-chimiques de ces polymères

### Dépendance principale
- `bigsmiles-rs` : https://github.com/Peariforme/bigsmiles-rs
  - Crates utilisées : `bigsmiles` (parser BigSMILES) et `opensmiles` (parser SMILES)

---

## Structure du projet

```
polysim/
├── crates/
│   ├── polysim-core/         # Bibliothèque principale
│   │   └── src/
│   │       ├── builder/      # Générateurs de chaînes
│   │       │   ├── strategy.rs   # BuildStrategy enum
│   │       │   ├── linear.rs     # Homo, random, alternating, block
│   │       │   └── branched.rs   # Branché, graft, macromonomer
│   │       ├── polymer/
│   │       │   └── chain.rs      # PolymerChain (SMILES + metadata)
│   │       └── properties/
│   │           ├── molecular_weight.rs  # Mn, Mw, PDI
│   │           └── thermal.rs           # Tg (Fox, Van Krevelen), cristallisation
│   └── polysim-cli/          # Outil en ligne de commande
└── tests/
    └── integration.rs
```

---

## Architectures polymères supportées

| Architecture       | Builder            | Paramètres clés                        |
|--------------------|--------------------|----------------------------------------|
| Homopolymère       | `LinearBuilder`    | n répétitions ou Mn cible              |
| Copolymère random  | `LinearBuilder`    | fractions pondérales par monomère      |
| Copolymère alternant | `LinearBuilder`  | 2 unités de répétition                 |
| Copolymère bloc    | `LinearBuilder`    | longueur de chaque bloc                |
| Polymère branché   | `BranchedBuilder`  | 1 branche tous les N monomères         |
| Graft copolymère   | `BranchedBuilder`  | fraction de greffage                   |
| Macromonomère      | `BranchedBuilder`  | chaîne latérale + groupe polymérisable |

---

## Stratégie de construction (`BuildStrategy`)

```rust
pub enum BuildStrategy {
    ByRepeatCount(usize),  // Nombre exact d'unités de répétition
    ByTargetMn(f64),       // Masse moléculaire moyenne en nombre cible (g/mol)
    ByExactMass(f64),      // Masse exacte monoisotopique cible (g/mol)
}
```

---

## Propriétés à calculer

### Déjà scaffoldées
- `molecular_weight::monoisotopic_mass` — masse monoisotopique
- `molecular_weight::average_mass` — masse moléculaire moyenne
- `thermal::tg_fox` — Tg par équation de Fox (implémentée)
- `thermal::tg_van_krevelen` — Tg par contribution de groupes (TODO)
- `thermal::crystallization_tendency` — tendance à cristalliser (TODO)

### À ajouter
- Paramètre de solubilité (Hildebrand)
- Hydrophilicité / log P
- Module d'Young (contribution de groupes)
- Température de fusion Tm
- Distribution de masse (Mw, PDI) pour des ensembles de chaînes

---

## Références scientifiques

- Fox, T. G. (1956). *Bull. Am. Phys. Soc.* **1**, 123. — équation de Fox pour Tg
- Van Krevelen, D. W. (1990). *Properties of Polymers*, 3rd ed., Elsevier. — contributions de groupes
- Lin, T.-S. et al. (2019). *ACS Central Science* — BigSMILES

---

## Conventions de code

- Toutes les masses en **g/mol**
- Toutes les températures en **Kelvin**
- Les erreurs utilisent `thiserror` via `PolySimError`
- Les builders prennent un `BigSmiles` parsé (pas une `String`)
- Pas de `unwrap()` dans le code de bibliothèque
