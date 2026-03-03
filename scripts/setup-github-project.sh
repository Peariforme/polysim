#!/usr/bin/env bash
# =============================================================================
# setup-github-project.sh
#
# Creates all GitHub labels, milestones, issues, and project board for the
# polysim roadmap. Run this once after authenticating with `gh auth login`.
#
# Usage:
#   ./scripts/setup-github-project.sh
#
# Prerequisites:
#   - gh CLI installed and authenticated (gh auth login)
#   - jq installed
# =============================================================================
set -euo pipefail

REPO="Peariforme/polysim"

# ── Colors ────────────────────────────────────────────────────────────────────
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

info()  { echo -e "${CYAN}[INFO]${NC}  $*"; }
ok()    { echo -e "${GREEN}[OK]${NC}    $*"; }
warn()  { echo -e "${YELLOW}[WARN]${NC}  $*"; }
fail()  { echo -e "${RED}[FAIL]${NC}  $*"; exit 1; }

# ── Preflight ─────────────────────────────────────────────────────────────────
command -v gh  >/dev/null 2>&1 || fail "gh CLI not found. Install: https://cli.github.com/"
command -v jq  >/dev/null 2>&1 || fail "jq not found. Install: apt install jq / brew install jq"
gh auth status >/dev/null 2>&1 || fail "Not authenticated. Run: gh auth login"
info "Authenticated. Setting up polysim roadmap on $REPO..."

# ── Helper: create label if not exists ────────────────────────────────────────
create_label() {
  local name="$1" color="$2" description="$3"
  if gh label view "$name" --repo "$REPO" >/dev/null 2>&1; then
    warn "Label '$name' already exists, skipping."
  else
    gh label create "$name" --repo "$REPO" --color "$color" --description "$description"
    ok "Label '$name' created."
  fi
}

# ── Helper: create milestone if not exists ────────────────────────────────────
create_milestone() {
  local title="$1" description="$2"
  local existing
  existing=$(gh api "repos/$REPO/milestones" --paginate -q ".[] | select(.title == \"$title\") | .number" 2>/dev/null || true)
  if [[ -n "$existing" ]]; then
    warn "Milestone '$title' already exists (#$existing), skipping."
    echo "$existing"
  else
    local number
    number=$(gh api "repos/$REPO/milestones" -f title="$title" -f description="$description" -f state="open" -q '.number')
    ok "Milestone '$title' created (#$number)."
    echo "$number"
  fi
}

# ── Helper: create issue ──────────────────────────────────────────────────────
create_issue() {
  local title="$1" body="$2" milestone="$3"
  shift 3
  local labels=("$@")

  local label_args=()
  for l in "${labels[@]}"; do
    label_args+=(--label "$l")
  done

  local number
  number=$(gh issue create --repo "$REPO" \
    --title "$title" \
    --body "$body" \
    --milestone "$milestone" \
    "${label_args[@]}" \
    2>&1 | grep -oP '\d+$' || true)

  if [[ -n "$number" ]]; then
    ok "Issue #$number: $title"
  else
    warn "Failed to create issue: $title"
  fi
}

# =============================================================================
# 1. LABELS
# =============================================================================
info "Creating labels..."

# Epic labels
create_label "epic:copolymers"      "1d76db" "Epic 1.1: Linear copolymers"
create_label "epic:branched"        "1d76db" "Epic 1.2: Branched architectures"
create_label "epic:chain-repr"      "1d76db" "Epic 1.3: Chain representation"
create_label "epic:group-contrib"   "0e8a16" "Epic 2.1: Group contribution methods"
create_label "epic:thermal"         "0e8a16" "Epic 2.2: Thermal properties"
create_label "epic:solubility"      "0e8a16" "Epic 2.3: Solubility parameters"
create_label "epic:mechanical"      "0e8a16" "Epic 2.4: Mechanical properties"
create_label "epic:transport"       "0e8a16" "Epic 2.5: Transport properties"
create_label "epic:optical"         "0e8a16" "Epic 2.6: Optical & electrical"
create_label "epic:topology"        "d93f0b" "Epic 3.1: Topology & atoms"
create_label "epic:3d-gen"          "d93f0b" "Epic 3.2: 3D coordinate generation"
create_label "epic:io"              "d93f0b" "Epic 3.3: File format I/O"
create_label "epic:desktop"         "d93f0b" "Epic 3.4: Desktop application"
create_label "epic:forcefield"      "e11d48" "Epic 4.1: Force fields"
create_label "epic:md"              "e11d48" "Epic 4.2: Molecular dynamics"
create_label "epic:mc"              "e11d48" "Epic 4.3: Monte Carlo"
create_label "epic:dpd"             "e11d48" "Epic 4.4: DPD"
create_label "epic:analysis"        "c2e0c6" "Epic 5.x: Analysis & post-processing"
create_label "epic:performance"     "5319e7" "Epic 6.1: Performance"
create_label "epic:python"          "5319e7" "Epic 6.2: Python bindings"
create_label "epic:docs"            "5319e7" "Epic 6.3: Documentation"
create_label "epic:cli"             "bfdadc" "CLI commands & UX"

# Priority labels
create_label "priority:P0"  "b60205" "Critical — must-have for the phase"
create_label "priority:P1"  "fbca04" "Important — high value"
create_label "priority:P2"  "0e8a16" "Nice-to-have — can be deferred"

# Phase labels
create_label "phase:1"  "c5def5" "Phase 1: Polymer architectures (v0.3-v0.5)"
create_label "phase:2"  "bfd4f2" "Phase 2: Property predictions (v0.6-v1.0)"
create_label "phase:3"  "d4c5f9" "Phase 3: 3D generation + desktop (v1.1-v1.5)"
create_label "phase:4"  "f9d0c4" "Phase 4: Simulation engines (v2.0+)"
create_label "phase:5"  "c2e0c6" "Phase 5: Analysis (v2.1+)"
create_label "phase:6"  "e4e669" "Phase 6: Ecosystem (v2.5+)"

echo ""

# =============================================================================
# 2. MILESTONES
# =============================================================================
info "Creating milestones..."

M1=$(create_milestone "Phase 1: Polymer Architectures (v0.3-v0.5)" \
  "Complete all polymer architectures: copolymers (random, alternating, block, gradient), branched (comb, star, graft, dendrimer), cyclic. Enrich PolymerChain representation.")

M2=$(create_milestone "Phase 2: Property Predictions (v0.6-v1.0)" \
  "Group contribution methods (Van Krevelen, Bicerano). Predict: Tg, Tm, density, solubility (Hildebrand, Hansen), mechanical, optical, transport properties.")

M3=$(create_milestone "Phase 3: 3D Generation + Desktop App (v1.1-v1.5)" \
  "Topology data structures, 3D coordinate generation (random walk, RIS), simulation box with PBC, file I/O (LAMMPS, GROMACS, PDB), desktop application for visualization.")

M4=$(create_milestone "Phase 4: Simulation Engines (v2.0+)" \
  "Force fields (TraPPE-UA, OPLS-AA), MD (velocity Verlet, NVT/NPT), MC (Metropolis, CBMC, Wang-Landau), DPD.")

M5=$(create_milestone "Phase 5: Analysis (v2.1+)" \
  "Post-processing: RDF, MSD, Rg, S(q), autocorrelation, Green-Kubo viscosity, Tg by simulation.")

M6=$(create_milestone "Phase 6: Ecosystem (v2.5+)" \
  "Performance (Rayon, SIMD, neighbor lists), Python bindings (PyO3), tutorials, documentation, CLI run/init.")

echo ""

# =============================================================================
# 3. ISSUES — Phase 1
# =============================================================================
info "Creating Phase 1 issues..."

create_issue "US-1.3.1: Enrich PolymerChain with monomer composition" \
"## User Story
> As a developer, I want to enrich \`PolymerChain\` to store monomer composition (types and fractions), to enable copolymer property calculations.

## Acceptance Criteria
- [ ] New field: \`composition: Vec<MonomerUnit>\` with \`MonomerUnit { smiles: String, fraction: f64 }\`
- [ ] Backwards-compatible: homopolymers have a single-element composition
- [ ] Accessible via public API

## Files
- \`crates/polysim-core/src/polymer/chain.rs\`
- \`crates/polysim-core/src/polymer/ensemble.rs\`" \
  "Phase 1: Polymer Architectures (v0.3-v0.5)" \
  "epic:chain-repr" "priority:P0" "phase:1"

create_issue "US-1.3.2: Add Architecture enum to PolymerChain" \
"## User Story
> As a developer, I want to add an \`architecture\` field to \`PolymerChain\` (Linear, Star, Comb, Dendrimer, Cyclic) to condition property calculations.

## Acceptance Criteria
- [ ] Enum \`Architecture { Linear, Star { arms }, Comb { branch_spacing }, Dendrimer { generation }, Cyclic }\`
- [ ] Each builder populates this field automatically

## Files
- \`crates/polysim-core/src/polymer/chain.rs\`" \
  "Phase 1: Polymer Architectures (v0.3-v0.5)" \
  "epic:chain-repr" "priority:P0" "phase:1"

create_issue "US-1.1.1: Random (statistical) copolymer generation" \
"## User Story
> As a polymer chemist, I want to generate a random copolymer from a multi-monomer BigSMILES with molar fractions, to model real copolymers from radical polymerization.

## Acceptance Criteria
- [ ] \`LinearBuilder::random_copolymer(&[f64])\` generates SMILES with random monomer placement per fractions
- [ ] Fractions must sum to 1.0 (error otherwise)
- [ ] Seed support for reproducibility
- [ ] Correct Mn calculation
- [ ] Tests: monomer distribution converges to target fractions (N=1000)

## Files
- \`crates/polysim-core/src/builder/linear.rs\` (line ~85, currently \`todo!()\`)" \
  "Phase 1: Polymer Architectures (v0.3-v0.5)" \
  "epic:copolymers" "priority:P0" "phase:1"

create_issue "US-1.1.2: Alternating copolymer generation" \
"## User Story
> As a polymer chemist, I want to generate an alternating copolymer (A-B-A-B) to model polymers like SMA (styrene-maleic anhydride).

## Acceptance Criteria
- [ ] \`LinearBuilder::alternating_copolymer()\` generates a strictly alternating sequence
- [ ] BigSMILES must contain exactly 2 repeat units
- [ ] SMILES output verifies ABABAB pattern
- [ ] Tests on PE-alt-PP, PS-alt-PMA

## Files
- \`crates/polysim-core/src/builder/linear.rs\` (line ~96, currently \`todo!()\`)" \
  "Phase 1: Polymer Architectures (v0.3-v0.5)" \
  "epic:copolymers" "priority:P0" "phase:1"

create_issue "US-1.1.3: Block copolymer generation" \
"## User Story
> As a polymer chemist, I want to generate a block copolymer (AAAA-BBBB) with specified block lengths, to model systems like PS-b-PEO.

## Acceptance Criteria
- [ ] \`LinearBuilder::block_copolymer(&[usize])\` accepts N blocks with lengths
- [ ] Support 2+ blocks (diblock, triblock, multiblock)
- [ ] Correct ring closure handling between blocks
- [ ] Tests: PS-b-PMMA diblock, PS-b-PB-b-PS triblock (SBS)

## Files
- \`crates/polysim-core/src/builder/linear.rs\` (line ~106, currently \`todo!()\`)" \
  "Phase 1: Polymer Architectures (v0.3-v0.5)" \
  "epic:copolymers" "priority:P1" "phase:1"

create_issue "US-1.1.4: Gradient copolymer generation" \
"## User Story
> As a polymer chemist, I want to generate a gradient copolymer where composition varies progressively along the chain.

## Acceptance Criteria
- [ ] New method \`LinearBuilder::gradient_copolymer(profile: GradientProfile)\`
- [ ] Profiles: linear, sigmoid
- [ ] Monomer A fraction transitions from \`f_start\` to \`f_end\` along the chain
- [ ] Statistical verification of gradient monotonicity" \
  "Phase 1: Polymer Architectures (v0.3-v0.5)" \
  "epic:copolymers" "priority:P2" "phase:1"

create_issue "US-1.2.1: Comb polymer generation" \
"## User Story
> As a polymer chemist, I want to generate a comb polymer with regularly spaced branches, to model architectures like PEG-comb.

## Acceptance Criteria
- [ ] \`BranchedBuilder::comb_polymer(branch_every)\` works
- [ ] Backbone and branches can be different monomers
- [ ] Correct SMILES branching
- [ ] Mn includes backbone + branches

## Files
- \`crates/polysim-core/src/builder/branched.rs\` (currently \`todo!()\`)" \
  "Phase 1: Polymer Architectures (v0.3-v0.5)" \
  "epic:branched" "priority:P1" "phase:1"

create_issue "US-1.2.2: Graft copolymer generation" \
"## User Story
> As a polymer chemist, I want to generate a graft copolymer with randomly placed grafts, to model systems like HIPS.

## Acceptance Criteria
- [ ] \`BranchedBuilder::graft_copolymer(graft_fraction)\` places grafts randomly
- [ ] Seed for reproducibility
- [ ] Graft fraction is statistically respected

## Files
- \`crates/polysim-core/src/builder/branched.rs\`" \
  "Phase 1: Polymer Architectures (v0.3-v0.5)" \
  "epic:branched" "priority:P2" "phase:1"

create_issue "US-1.2.3: Star polymer generation" \
"## User Story
> As a polymer chemist, I want to generate star polymers with a defined number of arms.

## Acceptance Criteria
- [ ] \`BranchedBuilder::star_polymer(arms: usize)\`
- [ ] Support 3 to 12 arms
- [ ] Same monomer per arm (homostar) or different (miktostar)" \
  "Phase 1: Polymer Architectures (v0.3-v0.5)" \
  "epic:branched" "priority:P1" "phase:1"

create_issue "US-1.2.4: Dendrimer generation" \
"## User Story
> As a polymer chemist, I want to generate dendrimers of specified generation, to model PAMAM or poly(propylene imine).

## Acceptance Criteria
- [ ] \`BranchedBuilder::dendrimer(generation: usize, branching_factor: usize)\`
- [ ] Generations 1 to 6
- [ ] Correct SMILES with deep ring closures
- [ ] Monomer count per generation" \
  "Phase 1: Polymer Architectures (v0.3-v0.5)" \
  "epic:branched" "priority:P2" "phase:1"

create_issue "US-1.2.5: Cyclic (ring) polymer generation" \
"## User Story
> As a polymer chemist, I want to generate cyclic polymers without chain ends.

## Acceptance Criteria
- [ ] \`LinearBuilder::cyclic_homopolymer()\`
- [ ] SMILES forms a ring (ring closure between first and last atom)
- [ ] Mn does not include end group mass" \
  "Phase 1: Polymer Architectures (v0.3-v0.5)" \
  "epic:branched" "priority:P2" "phase:1"

create_issue "US-1.3.3: Explicit end groups support" \
"## User Story
> As a polymer chemist, I want to specify end groups explicitly, since they significantly influence properties at low Mn.

## Acceptance Criteria
- [ ] Support BigSMILES \`begin\`/\`end\` sections in builders
- [ ] End groups considered in MW, formula, and Tg calculations" \
  "Phase 1: Polymer Architectures (v0.3-v0.5)" \
  "epic:chain-repr" "priority:P1" "phase:1"

create_issue "US-1.4.1: CLI analyze — copolymer support" \
"## User Story
> As a CLI user, I want \`polysim analyze\` extended to support \`--copolymer random|alternating|block|gradient\`.

## Acceptance Criteria
- [ ] \`polysim analyze \"{[]CC[]; CC(C)[]}\" --copolymer random --fractions 0.7,0.3 --by-repeat 100\`
- [ ] \`polysim analyze \"{[]CC[]; CC(c1ccccc1)[]}\" --copolymer block --block-lengths 50,50\`
- [ ] Output: composition, architecture, properties" \
  "Phase 1: Polymer Architectures (v0.3-v0.5)" \
  "epic:cli" "priority:P1" "phase:1"

create_issue "US-1.4.2: CLI generate — copolymer ensemble support" \
"## User Story
> As a CLI user, I want \`polysim generate\` to support copolymer ensembles for compositional polydispersity studies.

## Acceptance Criteria
- [ ] Ensembles with composition and chain length variation" \
  "Phase 1: Polymer Architectures (v0.3-v0.5)" \
  "epic:cli" "priority:P2" "phase:1"

# =============================================================================
# 3. ISSUES — Phase 2
# =============================================================================
info "Creating Phase 2 issues..."

create_issue "US-2.1.1: Generic group contribution system" \
"## User Story
> As a developer, I want a generic group contribution system that decomposes SMILES into functional groups and sums contributions.

## Acceptance Criteria
- [ ] Trait \`GroupContributionMethod\` with \`fn predict(&self, groups: &[Group]) -> f64\`
- [ ] Van Krevelen group database (~40 groups for common polymers)
- [ ] SMILES -> functional groups decomposition (SMARTS or pattern matching)
- [ ] Tests: decomposition of PE, PP, PS, PMMA, PET, Nylon-6,6

## References
- Van Krevelen & te Nijenhuis (2009), *Properties of Polymers*, 4th ed., Elsevier" \
  "Phase 2: Property Predictions (v0.6-v1.0)" \
  "epic:group-contrib" "priority:P0" "phase:2"

create_issue "US-2.1.2: Extensible group database via TOML" \
"## User Story
> As a developer, I want to extend the group contribution database via an external TOML file.

## Acceptance Criteria
- [ ] TOML format for new groups and contributions
- [ ] Merge with default database
- [ ] Validation (physically reasonable intervals)" \
  "Phase 2: Property Predictions (v0.6-v1.0)" \
  "epic:group-contrib" "priority:P2" "phase:2"

create_issue "US-2.2.1: Tg prediction — Van Krevelen method" \
"## User Story
> As a polymer chemist, I want to predict Tg via Van Krevelen (Yg = sum contributions / repeat unit mass).

## Acceptance Criteria
- [ ] \`tg_van_krevelen(chain)\` returns Tg in K
- [ ] Precision: error < 15 K for PE, PP, PS, PMMA, PVC, PET
- [ ] Reference: Van Krevelen, Chapter 6

## Files
- \`crates/polysim-core/src/properties/thermal.rs\` (line 34, currently \`todo!()\`)" \
  "Phase 2: Property Predictions (v0.6-v1.0)" \
  "epic:thermal" "priority:P0" "phase:2"

create_issue "US-2.2.2: Melting temperature Tm prediction" \
"## User Story
> As a polymer chemist, I want to predict Tm by group contributions.

## Acceptance Criteria
- [ ] \`tm_van_krevelen(chain) -> Option<f64>\`
- [ ] Returns \`None\` for amorphous polymers
- [ ] Cross-validation: Tg/Tm ratio ~ 2/3 (Boyer-Beaman rule)
- [ ] Tests: PE (Tm ~ 411 K), PP iso (Tm ~ 449 K), PET (Tm ~ 538 K)" \
  "Phase 2: Property Predictions (v0.6-v1.0)" \
  "epic:thermal" "priority:P0" "phase:2"

create_issue "US-2.2.3: Crystallization tendency estimation" \
"## User Story
> As a polymer chemist, I want to estimate the crystallization tendency of a polymer.

## Acceptance Criteria
- [ ] Implement \`crystallization_tendency()\` (currently \`todo!()\` in thermal.rs)
- [ ] Based on chain symmetry and Tm - Tg difference
- [ ] Returns existing \`CrystallizationTendency\` enum" \
  "Phase 2: Property Predictions (v0.6-v1.0)" \
  "epic:thermal" "priority:P1" "phase:2"

create_issue "US-2.2.4: Thermal expansion coefficient (CTE) and heat capacity (Cp)" \
"## User Story
> As a polymer chemist, I want to predict CTE and Cp.

## Acceptance Criteria
- [ ] \`thermal_expansion_coefficient(chain) -> f64\` (1/K)
- [ ] \`specific_heat_capacity(chain) -> f64\` (J/g.K)
- [ ] Van Krevelen group contributions" \
  "Phase 2: Property Predictions (v0.6-v1.0)" \
  "epic:thermal" "priority:P2" "phase:2"

create_issue "US-2.3.1: Hildebrand solubility parameter" \
"## User Story
> As a polymer chemist, I want to calculate the Hildebrand solubility parameter (delta) to predict polymer-solvent miscibility.

## Acceptance Criteria
- [ ] \`hildebrand_solubility_parameter(chain) -> f64\` in (MPa)^0.5
- [ ] Based on Ecoh and Vw by group contributions
- [ ] Tests: PS (delta ~ 18.5), PE (delta ~ 16.2), PMMA (delta ~ 18.6)" \
  "Phase 2: Property Predictions (v0.6-v1.0)" \
  "epic:solubility" "priority:P0" "phase:2"

create_issue "US-2.3.2: Hansen solubility parameters (3D)" \
"## User Story
> As a polymer chemist, I want Hansen parameters (delta_d, delta_p, delta_h) for finer miscibility analysis.

## Acceptance Criteria
- [ ] \`hansen_solubility_parameters(chain) -> HansenParams { d, p, h }\`
- [ ] Decomposition: dispersive, polar, H-bond contributions
- [ ] RED (Relative Energy Difference) between polymer and solvent" \
  "Phase 2: Property Predictions (v0.6-v1.0)" \
  "epic:solubility" "priority:P1" "phase:2"

create_issue "US-2.3.3: Flory-Huggins chi interaction parameter" \
"## User Story
> As a polymer chemist, I want to evaluate polymer-polymer miscibility via Flory-Huggins chi.

## Acceptance Criteria
- [ ] \`flory_huggins_chi(polymer_a, polymer_b, temperature) -> f64\`
- [ ] Based on chi ~ (delta_A - delta_B)^2
- [ ] chi < 0.5: miscible, chi > 0.5: immiscible" \
  "Phase 2: Property Predictions (v0.6-v1.0)" \
  "epic:solubility" "priority:P2" "phase:2"

create_issue "US-2.4.1: Young's modulus and tensile strength prediction" \
"## User Story
> As a materials engineer, I want to predict Young's modulus and tensile strength.

## Acceptance Criteria
- [ ] \`youngs_modulus(chain) -> f64\` (GPa) — via Rao function (Van Krevelen)
- [ ] \`tensile_strength(chain) -> f64\` (MPa)
- [ ] Differentiation amorphous vs semi-crystalline
- [ ] Tests: PS (E ~ 3.0 GPa), PE-HD (E ~ 1.0 GPa), PMMA (E ~ 3.3 GPa)" \
  "Phase 2: Property Predictions (v0.6-v1.0)" \
  "epic:mechanical" "priority:P1" "phase:2"

create_issue "US-2.4.2: Density prediction" \
"## User Story
> As a materials engineer, I want to estimate polymer density at room temperature.

## Acceptance Criteria
- [ ] \`density(chain) -> f64\` (g/cm3) — Van der Waals volume + packing factor
- [ ] Tests: PE (0.95), PS (1.05), PMMA (1.18), PVC (1.40)" \
  "Phase 2: Property Predictions (v0.6-v1.0)" \
  "epic:mechanical" "priority:P0" "phase:2"

create_issue "US-2.5.1: Gas permeability prediction" \
"## User Story
> As a process engineer, I want to predict gas permeability (O2, CO2, N2, H2O) for packaging films.

## Acceptance Criteria
- [ ] \`gas_permeability(chain, gas: Gas) -> f64\` (Barrer)
- [ ] Permachor method (Van Krevelen) or Salame correlations
- [ ] Tests: PE (high O2 perm), PET (good barrier), EVOH (excellent barrier)" \
  "Phase 2: Property Predictions (v0.6-v1.0)" \
  "epic:transport" "priority:P2" "phase:2"

create_issue "US-2.5.2: Intrinsic viscosity prediction (Mark-Houwink)" \
"## User Story
> As a polymer chemist, I want to predict intrinsic viscosity [eta] via Mark-Houwink equation.

## Acceptance Criteria
- [ ] \`intrinsic_viscosity(chain, params: MarkHouwinkParams) -> f64\` (dL/g)
- [ ] Database of K, a constants for common polymer-solvent pairs" \
  "Phase 2: Property Predictions (v0.6-v1.0)" \
  "epic:transport" "priority:P2" "phase:2"

create_issue "US-2.6.1: Refractive index prediction" \
"## User Story
> As an optical engineer, I want to predict the refractive index of a polymer.

## Acceptance Criteria
- [ ] \`refractive_index(chain) -> f64\` — via Lorentz-Lorenz molar refraction
- [ ] Tests: PS (n ~ 1.59), PMMA (n ~ 1.49), PC (n ~ 1.585)" \
  "Phase 2: Property Predictions (v0.6-v1.0)" \
  "epic:optical" "priority:P1" "phase:2"

create_issue "US-2.6.2: Dielectric constant and surface tension prediction" \
"## User Story
> As an engineer, I want to predict the dielectric constant and surface tension.

## Acceptance Criteria
- [ ] \`dielectric_constant(chain) -> f64\`
- [ ] \`surface_tension(chain) -> f64\` (mN/m) — via Parachor
- [ ] Tests: PE (gamma ~ 31 mN/m), PS (gamma ~ 40 mN/m)" \
  "Phase 2: Property Predictions (v0.6-v1.0)" \
  "epic:optical" "priority:P2" "phase:2"

create_issue "US-2.7.1: CLI \`polysim properties\` command" \
"## User Story
> As a CLI user, I want a \`polysim properties\` command that generates a full property report.

## Acceptance Criteria
- [ ] \`polysim properties \"{[]CC(c1ccccc1)[]}\" --by-repeat 100\`
- [ ] Sections: Thermal, Mechanical, Solubility, Transport, Optical
- [ ] Output formats: console table, JSON (\`--format json\`), CSV (\`--format csv\`)
- [ ] Confidence indicator per property" \
  "Phase 2: Property Predictions (v0.6-v1.0)" \
  "epic:cli" "priority:P1" "phase:2"

create_issue "US-2.7.2: CLI \`polysim compare\` command" \
"## User Story
> As a CLI user, I want to compare properties of multiple polymers side by side.

## Acceptance Criteria
- [ ] \`polysim compare \"{[]CC[]}\" \"{[]CC(C)[]}\" \"{[]CC(c1ccccc1)[]}\" --by-repeat 100\`
- [ ] Comparative table with all properties" \
  "Phase 2: Property Predictions (v0.6-v1.0)" \
  "epic:cli" "priority:P2" "phase:2"

# =============================================================================
# 3. ISSUES — Phase 3
# =============================================================================
info "Creating Phase 3 issues..."

create_issue "US-3.1.1: Topology data structure" \
"## User Story
> As a developer, I want a Topology structure representing atoms, bonds, angles, dihedrals, and impropers.

## Acceptance Criteria
- [ ] \`Atom { id, element, mass, charge, position: [f64; 3] }\`
- [ ] \`Bond { i, j, bond_type }\`, \`Angle { i, j, k }\`, \`Dihedral { i, j, k, l }\`
- [ ] Incremental construction methods (add_atom, add_bond, etc.)
- [ ] Automatic angle and dihedral detection from bonds

## New Crate
\`crates/polysim-geom/\`" \
  "Phase 3: 3D Generation + Desktop App (v1.1-v1.5)" \
  "epic:topology" "priority:P0" "phase:3"

create_issue "US-3.1.2: PolymerChain (SMILES) to Topology conversion" \
"## User Story
> As a developer, I want to convert a PolymerChain to Topology with full connectivity.

## Acceptance Criteria
- [ ] SMILES parsing -> molecular graph (atoms + bonds)
- [ ] Automatic implicit hydrogen addition
- [ ] Bond, angle, dihedral detection
- [ ] Atom types assigned by local chemistry" \
  "Phase 3: 3D Generation + Desktop App (v1.1-v1.5)" \
  "epic:topology" "priority:P0" "phase:3"

create_issue "US-3.1.3: Simulation box with periodic boundary conditions" \
"## User Story
> As a developer, I want a simulation box with PBC support.

## Acceptance Criteria
- [ ] \`SimulationBox { origin, dimensions, periodicity: [bool; 3] }\`
- [ ] Support: orthogonal (orthorhombic) and triclinic
- [ ] \`wrap(position)\` and \`minimum_image(r_ij)\` methods" \
  "Phase 3: 3D Generation + Desktop App (v1.1-v1.5)" \
  "epic:topology" "priority:P0" "phase:3"

create_issue "US-3.2.1: Initial conformation by random walk" \
"## User Story
> As a computational chemist, I want to generate initial conformations by random walk.

## Acceptance Criteria
- [ ] Fixed bond lengths (C-C = 1.54 A, etc.)
- [ ] Tetrahedral bond angles (109.5 deg)
- [ ] Random dihedral rotation (or Boltzmann distribution)
- [ ] No atomic overlap (minimum distance check)
- [ ] Seed for reproducibility" \
  "Phase 3: 3D Generation + Desktop App (v1.1-v1.5)" \
  "epic:3d-gen" "priority:P0" "phase:3"

create_issue "US-3.2.2: RIS (Rotational Isomeric State) conformation model" \
"## User Story
> As a computational chemist, I want RIS conformations for more realistic structures.

## Acceptance Criteria
- [ ] Statistical weight matrices U for dihedral pairs
- [ ] RIS parameters for PE, PP, PS, PMMA (from literature)
- [ ] Monte Carlo generation with RIS weights
- [ ] C_inf (characteristic ratio) calculation as validation" \
  "Phase 3: 3D Generation + Desktop App (v1.1-v1.5)" \
  "epic:3d-gen" "priority:P1" "phase:3"

create_issue "US-3.2.3: Simulation box packing at target density" \
"## User Story
> As a simulator, I want to pack a simulation box with multiple chains at a target density.

## Acceptance Criteria
- [ ] \`PackingBuilder::new(box_size, density, chains)\`
- [ ] Random insertion with overlap rejection
- [ ] Alternative: grow-and-push method
- [ ] Resulting density within 5% of target" \
  "Phase 3: 3D Generation + Desktop App (v1.1-v1.5)" \
  "epic:3d-gen" "priority:P1" "phase:3"

create_issue "US-3.2.4: Energy minimization (steepest descent / CG)" \
"## User Story
> As a simulator, I want fast energy minimization to relax initial geometries.

## Acceptance Criteria
- [ ] Steepest descent with line search
- [ ] Convergence on max gradient or energy
- [ ] PBC support
- [ ] Configurable stopping criteria" \
  "Phase 3: 3D Generation + Desktop App (v1.1-v1.5)" \
  "epic:3d-gen" "priority:P2" "phase:3"

create_issue "US-3.3.1: LAMMPS data file export" \
"## User Story
> As a simulator, I want to export to LAMMPS data file format.

## Acceptance Criteria
- [ ] Complete header (atoms, bonds, angles, dihedrals)
- [ ] Sections: Atoms, Bonds, Angles, Dihedrals, Masses
- [ ] Styles: full, molecular
- [ ] Round-trip tests: write -> read -> compare

## New Crate
\`crates/polysim-io/\`" \
  "Phase 3: 3D Generation + Desktop App (v1.1-v1.5)" \
  "epic:io" "priority:P1" "phase:3"

create_issue "US-3.3.2: GROMACS export (.gro + .top)" \
"## User Story
> As a simulator, I want to export to GROMACS format.

## Acceptance Criteria
- [ ] .gro file (coordinates + velocities)
- [ ] .top file (topology, force field includes)
- [ ] PBC support" \
  "Phase 3: 3D Generation + Desktop App (v1.1-v1.5)" \
  "epic:io" "priority:P2" "phase:3"

create_issue "US-3.3.3: PDB and XYZ export/import" \
"## User Story
> As a chemist, I want to export/import PDB and XYZ formats.

## Acceptance Criteria
- [ ] PDB: ATOM/HETATM records, CONECT, CRYST1
- [ ] XYZ: standard format
- [ ] Read and write support" \
  "Phase 3: 3D Generation + Desktop App (v1.1-v1.5)" \
  "epic:io" "priority:P2" "phase:3"

create_issue "US-3.3.4: TOML simulation configuration format" \
"## User Story
> As a developer, I want a TOML configuration format to fully define a simulation.

## Acceptance Criteria
- [ ] System specification (polymer, density, box)
- [ ] Force field parameters
- [ ] Simulation parameters (integrator, dt, nsteps, T, P)
- [ ] Requested outputs (trajectory, energy, properties)
- [ ] Valid schema and clear error messages" \
  "Phase 3: 3D Generation + Desktop App (v1.1-v1.5)" \
  "epic:io" "priority:P1" "phase:3"

create_issue "US-3.4.1: Desktop app — 3D polymer chain visualization" \
"## User Story
> As a polymer chemist, I want to visualize generated polymer chains in 3D to verify architecture and conformation.

## Acceptance Criteria
- [ ] 3D rendering: ball-and-stick or space-filling
- [ ] Color by atom type, monomer, or charge
- [ ] Interactive rotation, zoom, pan
- [ ] Simulation box display (wireframe)

## Notes
Framework TBD (Tauri, egui, iced, etc.). New crate: \`crates/polysim-gui/\`" \
  "Phase 3: 3D Generation + Desktop App (v1.1-v1.5)" \
  "epic:desktop" "priority:P0" "phase:3"

create_issue "US-3.4.2: Desktop app — property panel" \
"## User Story
> As a chemist, I want a side panel displaying predicted properties next to the 3D view.

## Acceptance Criteria
- [ ] All Phase 2 properties displayed in a structured panel
- [ ] Real-time update when polymer changes
- [ ] Export report to JSON/CSV/PDF" \
  "Phase 3: 3D Generation + Desktop App (v1.1-v1.5)" \
  "epic:desktop" "priority:P1" "phase:3"

create_issue "US-3.4.3: Desktop app — interactive polymer builder" \
"## User Story
> As a chemist, I want to build polymers interactively by selecting monomers, architecture, length, etc.

## Acceptance Criteria
- [ ] Monomer selector (common monomers library)
- [ ] Architecture choice (linear, star, comb, block...)
- [ ] Parameters (length, fractions, arm count...)
- [ ] Instant generation with 3D preview" \
  "Phase 3: 3D Generation + Desktop App (v1.1-v1.5)" \
  "epic:desktop" "priority:P1" "phase:3"

create_issue "US-3.4.4: Desktop app — simulation box visualization" \
"## User Story
> As a simulator, I want to visualize a packed simulation box with PBC.

## Acceptance Criteria
- [ ] Box rendering and periodic images
- [ ] Individual chain selection
- [ ] Distance/angle/dihedral display on click
- [ ] Ghost atoms (periodic images) toggle" \
  "Phase 3: 3D Generation + Desktop App (v1.1-v1.5)" \
  "epic:desktop" "priority:P2" "phase:3"

create_issue "US-3.4.5: Desktop app — simulation configuration editor" \
"## User Story
> As a simulator, I want to configure a simulation via GUI instead of editing TOML manually.

## Acceptance Criteria
- [ ] Form for all simulation parameters
- [ ] Real-time parameter validation
- [ ] Generated TOML preview
- [ ] \"Launch simulation\" button" \
  "Phase 3: 3D Generation + Desktop App (v1.1-v1.5)" \
  "epic:desktop" "priority:P2" "phase:3"

create_issue "US-3.4.6: Desktop app — trajectory viewer" \
"## User Story
> As a simulator, I want to view simulation trajectories frame by frame with playback controls.

## Acceptance Criteria
- [ ] Load DCD, XTC, XYZ trajectories
- [ ] Controls: play, pause, step forward/backward, speed
- [ ] Timeline slider
- [ ] Property overlay (energy, temperature) on timeline" \
  "Phase 3: 3D Generation + Desktop App (v1.1-v1.5)" \
  "epic:desktop" "priority:P2" "phase:3"

# =============================================================================
# 3. ISSUES — Phase 4
# =============================================================================
info "Creating Phase 4 issues..."

create_issue "US-4.1.1: Generic Potential trait" \
"## User Story
> As a developer, I want a \`Potential\` trait for defining interaction potentials.

## Acceptance Criteria
- [ ] \`trait Potential { fn energy(&self, r: f64) -> f64; fn force(&self, r: f64) -> f64; }\`
- [ ] Implementations: LennardJones, Harmonic, CosineDihedral, FENE
- [ ] SIMD-friendly vectorized calculation

## New Crate
\`crates/polysim-ff/\`" \
  "Phase 4: Simulation Engines (v2.0+)" \
  "epic:forcefield" "priority:P0" "phase:4"

create_issue "US-4.1.2: TraPPE-UA force field" \
"## User Story
> As a simulator, I want TraPPE-UA for hydrocarbon polymer simulations.

## Acceptance Criteria
- [ ] Parameters for CH4, CH3, CH2, CH groups
- [ ] Bond, angle, dihedral parameters
- [ ] LJ with cutoff and tail corrections
- [ ] Tests: short PE energy vs LAMMPS reference" \
  "Phase 4: Simulation Engines (v2.0+)" \
  "epic:forcefield" "priority:P0" "phase:4"

create_issue "US-4.1.3: OPLS-AA force field" \
"## User Story
> As a simulator, I want OPLS-AA for detailed all-atom simulations.

## Acceptance Criteria
- [ ] Parameters for common polymer atoms (C, H, O, N)
- [ ] Combination rules (geometric epsilon, arithmetic sigma)
- [ ] Partial charge support" \
  "Phase 4: Simulation Engines (v2.0+)" \
  "epic:forcefield" "priority:P1" "phase:4"

create_issue "US-4.1.4: Custom force field via parameter file" \
"## User Story
> As a simulator, I want to define custom force fields via TOML.

## Acceptance Criteria
- [ ] TOML format for atom types, potentials, parameters
- [ ] Parameter validation
- [ ] Configurable mixing rules" \
  "Phase 4: Simulation Engines (v2.0+)" \
  "epic:forcefield" "priority:P2" "phase:4"

create_issue "US-4.2.1: Velocity Verlet integrator (NVE)" \
"## User Story
> As a simulator, I want a velocity Verlet integrator for NVE MD.

## Acceptance Criteria
- [ ] Total energy conservation (drift < 10^-4 kT over 10^6 steps)
- [ ] Configurable timestep (1-2 fs typical)
- [ ] PBC support
- [ ] Cell lists for O(N) scaling

## New Crate
\`crates/polysim-engine/\`" \
  "Phase 4: Simulation Engines (v2.0+)" \
  "epic:md" "priority:P0" "phase:4"

create_issue "US-4.2.2: Thermostats for NVT simulations" \
"## User Story
> As a simulator, I want thermostats for NVT simulations.

## Acceptance Criteria
- [ ] Nose-Hoover chain (deterministic, correct ensemble)
- [ ] Langevin (stochastic, good for relaxation)
- [ ] Velocity rescaling (Bussi-Donadio-Parrinello, correct canonical)
- [ ] Verification: Maxwell-Boltzmann velocity distribution" \
  "Phase 4: Simulation Engines (v2.0+)" \
  "epic:md" "priority:P0" "phase:4"

create_issue "US-4.2.3: Barostats for NPT simulations" \
"## User Story
> As a simulator, I want barostats for NPT simulations.

## Acceptance Criteria
- [ ] Berendsen (fast equilibration)
- [ ] Parrinello-Rahman (correct fluctuations)
- [ ] Isotropic and anisotropic support
- [ ] Verification: equilibrium density vs experimental (PE)" \
  "Phase 4: Simulation Engines (v2.0+)" \
  "epic:md" "priority:P1" "phase:4"

create_issue "US-4.2.4: Trajectory writing during simulation" \
"## User Story
> As a simulator, I want trajectory output during simulation.

## Acceptance Criteria
- [ ] Compact binary format (polysim-native)
- [ ] XYZ export (human-readable, debug)
- [ ] DCD or XTC export (VMD/MDAnalysis compatible)
- [ ] Configurable write frequency, optional compression" \
  "Phase 4: Simulation Engines (v2.0+)" \
  "epic:md" "priority:P1" "phase:4"

create_issue "US-4.3.1: Monte Carlo Metropolis with polymer moves" \
"## User Story
> As a simulator, I want MC Metropolis with polymer-adapted moves.

## Acceptance Criteria
- [ ] Moves: translation, pivot, reptation, crankshaft
- [ ] Correct Metropolis criterion
- [ ] Configurable target acceptance rate (auto-adjustment)" \
  "Phase 4: Simulation Engines (v2.0+)" \
  "epic:mc" "priority:P1" "phase:4"

create_issue "US-4.3.2: Configurational-bias Monte Carlo (CBMC)" \
"## User Story
> As a simulator, I want CBMC for chain insertion and regrowth.

## Acceptance Criteria
- [ ] Rosenbluth sampling for chain growth
- [ ] Configurable trial orientations
- [ ] Correct for grand-canonical ensemble" \
  "Phase 4: Simulation Engines (v2.0+)" \
  "epic:mc" "priority:P2" "phase:4"

create_issue "US-4.3.3: Wang-Landau sampling" \
"## User Story
> As a researcher, I want Wang-Landau for density of states and phase transitions.

## Acceptance Criteria
- [ ] 1/t-WL algorithm (proven convergence)
- [ ] Energy histogram
- [ ] g(E) export" \
  "Phase 4: Simulation Engines (v2.0+)" \
  "epic:mc" "priority:P2" "phase:4"

create_issue "US-4.4.1: Dissipative Particle Dynamics (DPD) engine" \
"## User Story
> As a mesoscale simulator, I want DPD for large-scale block copolymer morphology.

## Acceptance Criteria
- [ ] Conservative, dissipative, random forces (Groot-Warren)
- [ ] Built-in DPD thermostat
- [ ] Modified velocity Verlet integrator
- [ ] Repulsion parameters aij mapped from Flory-Huggins chi" \
  "Phase 4: Simulation Engines (v2.0+)" \
  "epic:dpd" "priority:P2" "phase:4"

# =============================================================================
# 3. ISSUES — Phase 5
# =============================================================================
info "Creating Phase 5 issues..."

create_issue "US-5.1.1: Radial distribution function g(r)" \
"## User Story
> As an analyst, I want to compute the RDF from a trajectory.

## Acceptance Criteria
- [ ] Total and per-pair RDF
- [ ] PBC support (minimum image)
- [ ] Frame averaging
- [ ] Correct normalization (ideal gas density)

## New Crate
\`crates/polysim-analysis/\`" \
  "Phase 5: Analysis (v2.1+)" \
  "epic:analysis" "priority:P0" "phase:5"

create_issue "US-5.1.2: Radius of gyration Rg and end-to-end distance R_ee" \
"## User Story
> As an analyst, I want Rg and R_ee for polymer chains.

## Acceptance Criteria
- [ ] Instantaneous and trajectory-averaged Rg
- [ ] R_ee and <R_ee^2>
- [ ] Ratio <R_ee^2> / <Rg^2> (theoretically 6 for ideal chain)
- [ ] Per-chain and ensemble average" \
  "Phase 5: Analysis (v2.1+)" \
  "epic:analysis" "priority:P0" "phase:5"

create_issue "US-5.1.3: Structure factor S(q)" \
"## User Story
> As an analyst, I want S(q) for comparison with SAXS/SANS experiments.

## Acceptance Criteria
- [ ] Debye method (small systems)
- [ ] FFT method (large systems)
- [ ] Configurable q range" \
  "Phase 5: Analysis (v2.1+)" \
  "epic:analysis" "priority:P1" "phase:5"

create_issue "US-5.1.4: Persistence length and Kuhn length" \
"## User Story
> As an analyst, I want to compute persistence length lp and Kuhn length b.

## Acceptance Criteria
- [ ] Bond vector autocorrelation
- [ ] Exponential fit for lp
- [ ] Kuhn length b = 2 * lp" \
  "Phase 5: Analysis (v2.1+)" \
  "epic:analysis" "priority:P1" "phase:5"

create_issue "US-5.2.1: Mean square displacement (MSD) and diffusion coefficient" \
"## User Story
> As an analyst, I want MSD and diffusion coefficient D.

## Acceptance Criteria
- [ ] MSD per atom type and chain center-of-mass
- [ ] Multiple time origin averaging
- [ ] D extraction by linear regression (diffusive regime)
- [ ] Auto-detection of ballistic/Rouse/reptation regimes" \
  "Phase 5: Analysis (v2.1+)" \
  "epic:analysis" "priority:P0" "phase:5"

create_issue "US-5.2.2: Autocorrelation functions (VACF, end-to-end)" \
"## User Story
> As an analyst, I want autocorrelation functions for velocity, dipole, end-to-end vector.

## Acceptance Criteria
- [ ] VACF -> vibrational spectrum
- [ ] End-to-end ACF -> Rouse/reptation relaxation time
- [ ] FFT for efficiency" \
  "Phase 5: Analysis (v2.1+)" \
  "epic:analysis" "priority:P1" "phase:5"

create_issue "US-5.3.1: Tg determination by simulation (cooling ramp)" \
"## User Story
> As a researcher, I want to determine Tg by plotting specific volume vs temperature.

## Acceptance Criteria
- [ ] Automated cooling protocol (NPT at different temperatures)
- [ ] Slope change detection (bi-linear regression)
- [ ] Comparison with group-contribution-predicted Tg" \
  "Phase 5: Analysis (v2.1+)" \
  "epic:analysis" "priority:P2" "phase:5"

create_issue "US-5.3.2: Shear modulus G(t) and viscosity (Green-Kubo)" \
"## User Story
> As a researcher, I want G(t) and viscosity via Green-Kubo method.

## Acceptance Criteria
- [ ] Stress tensor autocorrelation
- [ ] Integration for viscosity eta
- [ ] Long correlation time support (essential for polymers)" \
  "Phase 5: Analysis (v2.1+)" \
  "epic:analysis" "priority:P2" "phase:5"

create_issue "US-5.4.1: CLI \`polysim analyze-trajectory\` command" \
"## User Story
> As a CLI user, I want to analyze a trajectory from the command line.

## Acceptance Criteria
- [ ] \`polysim analyze-trajectory traj.dcd --topology system.top --rdf --msd --rg\`
- [ ] CSV or JSON output
- [ ] Progress bar for long analyses" \
  "Phase 5: Analysis (v2.1+)" \
  "epic:cli" "priority:P1" "phase:5"

create_issue "US-5.5.1: Desktop app — analysis visualization" \
"## User Story
> As a user, I want to visualize analysis curves (RDF, MSD, S(q)) in the desktop app.

## Acceptance Criteria
- [ ] Interactive plots (zoom, pan, range selection)
- [ ] Overlay multiple curves (compare runs)
- [ ] Export PNG/SVG" \
  "Phase 5: Analysis (v2.1+)" \
  "epic:desktop" "priority:P1" "phase:5"

create_issue "US-5.5.2: Desktop app — real-time simulation monitoring" \
"## User Story
> As a simulator, I want to monitor a running simulation with live graphs.

## Acceptance Criteria
- [ ] Real-time temperature, pressure, total energy
- [ ] Alert if simulation diverges
- [ ] Stop/pause button from the interface" \
  "Phase 5: Analysis (v2.1+)" \
  "epic:desktop" "priority:P2" "phase:5"

# =============================================================================
# 3. ISSUES — Phase 6
# =============================================================================
info "Creating Phase 6 issues..."

create_issue "US-6.1.1: Rayon shared-memory parallelism" \
"## User Story
> As a developer, I want Rayon parallelism for force loops and analysis.

## Acceptance Criteria
- [ ] Parallelized force calculation by cell lists
- [ ] Analysis (RDF, MSD) parallelized over frames
- [ ] Linear speedup up to 8-16 cores" \
  "Phase 6: Ecosystem (v2.5+)" \
  "epic:performance" "priority:P0" "phase:6"

create_issue "US-6.1.2: Optimized neighbor lists (cell list + Verlet list)" \
"## User Story
> As a developer, I want optimized neighbor lists for non-bonded interactions.

## Acceptance Criteria
- [ ] Cell list: O(N) construction
- [ ] Verlet list: conditional update (skin distance)
- [ ] Benchmark: > 10^4 atoms at 1 ns/day on desktop" \
  "Phase 6: Ecosystem (v2.5+)" \
  "epic:performance" "priority:P0" "phase:6"

create_issue "US-6.1.3: SIMD vectorization for pair forces" \
"## User Story
> As a developer, I want SIMD vectorization for pair-wise forces.

## Acceptance Criteria
- [ ] SoA (Structure of Arrays) layout for positions
- [ ] \`std::simd\` or \`packed_simd\` intrinsics
- [ ] Benchmark: 2-4x gain on LJ forces" \
  "Phase 6: Ecosystem (v2.5+)" \
  "epic:performance" "priority:P2" "phase:6"

create_issue "US-6.2.1: Python bindings for structures and properties (PyO3)" \
"## User Story
> As a Python chemist, I want to use polysim from Python.

## Acceptance Criteria
- [ ] \`pip install polysim\`
- [ ] \`from polysim import PolymerChain, LinearBuilder, BuildStrategy\`
- [ ] Pythonic API (snake_case, docstrings, type hints)
- [ ] NumPy support for 3D coordinates
- [ ] RDKit integration (SMILES <-> Mol conversion)

## New Crate
\`crates/polysim-python/\`" \
  "Phase 6: Ecosystem (v2.5+)" \
  "epic:python" "priority:P1" "phase:6"

create_issue "US-6.2.2: Python bindings for MD simulation" \
"## User Story
> As a Python researcher, I want to launch MD simulations from Python.

## Acceptance Criteria
- [ ] \`from polysim.md import Simulation\`
- [ ] Python configuration, Rust execution
- [ ] Python callbacks for monitoring
- [ ] MDAnalysis interoperability" \
  "Phase 6: Ecosystem (v2.5+)" \
  "epic:python" "priority:P2" "phase:6"

create_issue "US-6.3.1: Step-by-step tutorials" \
"## User Story
> As a new user, I want tutorials covering main use cases.

## Acceptance Criteria
- [ ] Tutorial 1: PS-b-PEO block copolymer + property prediction
- [ ] Tutorial 2: Polydisperse ensemble + statistical analysis
- [ ] Tutorial 3: PE simulation box + MD run
- [ ] Tutorial 4: Tg by simulation vs prediction" \
  "Phase 6: Ecosystem (v2.5+)" \
  "epic:docs" "priority:P1" "phase:6"

create_issue "US-6.3.2: Exhaustive API documentation with examples" \
"## User Story
> As a Rust developer, I want comprehensive API docs with examples.

## Acceptance Criteria
- [ ] \`#[doc]\` on all public functions/structs/traits
- [ ] Runnable examples (\`/// # Examples\`) for key methods
- [ ] \`cargo doc\` without warnings" \
  "Phase 6: Ecosystem (v2.5+)" \
  "epic:docs" "priority:P1" "phase:6"

create_issue "US-6.3.3: Contributing guide and benchmark suite" \
"## User Story
> As a contributor, I want a contribution guide and benchmarks to validate PRs.

## Acceptance Criteria
- [ ] CONTRIBUTING.md with conventions, architecture, workflow
- [ ] Reproducible benchmark suite (PE melt, PS glass, PS-b-PMMA morphology)
- [ ] Automated comparison with LAMMPS on reference cases" \
  "Phase 6: Ecosystem (v2.5+)" \
  "epic:docs" "priority:P2" "phase:6"

create_issue "US-6.4.1: CLI \`polysim run\` command" \
"## User Story
> As a CLI user, I want to launch a full simulation from a TOML config file.

## Acceptance Criteria
- [ ] \`polysim run simulation.toml\`
- [ ] Progress bar, ETA, real-time stats (T, P, E)
- [ ] Automatic checkpointing (resume after interruption)
- [ ] \`polysim continue checkpoint.bin\`" \
  "Phase 6: Ecosystem (v2.5+)" \
  "epic:cli" "priority:P1" "phase:6"

create_issue "US-6.4.2: CLI \`polysim init\` command" \
"## User Story
> As a CLI user, I want to generate a simulation config template.

## Acceptance Criteria
- [ ] \`polysim init --polymer \"{[]CC[]}\" --ensemble nvt --temperature 300\`
- [ ] Complete TOML with reasonable defaults
- [ ] Explanatory comments in TOML" \
  "Phase 6: Ecosystem (v2.5+)" \
  "epic:cli" "priority:P2" "phase:6"

# =============================================================================
# 4. PROJECT BOARD
# =============================================================================
info "Creating GitHub Project board..."

PROJECT_ID=$(gh project create --owner Peariforme --title "Polysim Roadmap" --format "$(cat <<'BODY'
Roadmap for polysim — the reference polymer simulation tool in Rust.
Organized in 6 phases: Architectures > Properties > 3D+Desktop > Simulation > Analysis > Ecosystem.
BODY
)" 2>&1 | grep -oP '\d+$' || true)

if [[ -n "$PROJECT_ID" ]]; then
  ok "Project created: https://github.com/orgs/Peariforme/projects/$PROJECT_ID"

  # Add custom fields
  info "Adding custom fields to project..."
  gh project field-create "$PROJECT_ID" --owner Peariforme --name "Phase" --data-type "SINGLE_SELECT" \
    --single-select-options "Phase 1,Phase 2,Phase 3,Phase 4,Phase 5,Phase 6" 2>/dev/null || warn "Could not add Phase field"
  gh project field-create "$PROJECT_ID" --owner Peariforme --name "Priority" --data-type "SINGLE_SELECT" \
    --single-select-options "P0 - Critical,P1 - Important,P2 - Nice-to-have" 2>/dev/null || warn "Could not add Priority field"

  # Add all open issues to the project
  info "Adding issues to project board..."
  gh issue list --repo "$REPO" --state open --limit 200 --json number -q '.[].number' | while read -r num; do
    gh project item-add "$PROJECT_ID" --owner Peariforme --url "https://github.com/$REPO/issues/$num" 2>/dev/null && \
      ok "Added #$num to project" || warn "Failed to add #$num"
  done
else
  warn "Could not create project. You may need to create it manually."
  info "To create manually: gh project create --owner Peariforme --title 'Polysim Roadmap'"
  info "Then add issues: gh issue list --repo $REPO --state open --json number -q '.[].number' | xargs -I{} gh project item-add <PROJECT_ID> --owner Peariforme --url 'https://github.com/$REPO/issues/{}'"
fi

echo ""
echo "═══════════════════════════════════════════════════════════"
echo -e "${GREEN}Setup complete!${NC}"
echo "═══════════════════════════════════════════════════════════"
echo ""
echo "Summary:"
gh issue list --repo "$REPO" --state open --json number -q 'length' | xargs -I{} echo "  Issues created: {}"
gh api "repos/$REPO/milestones" -q 'length' | xargs -I{} echo "  Milestones:     {}"
echo ""
echo "View your project board at: https://github.com/Peariforme/polysim/milestones"
echo "View all issues at: https://github.com/Peariforme/polysim/issues"
