# Polysim Roadmap

> Devenir LA reference de la simulation numerique des polymeres en Rust.

## Vision

Polysim ambitionne de combiner en un seul outil Rust :
- La **generation de structures polymeres** a partir de BigSMILES (toutes architectures)
- La **prediction de proprietes physiques** par contributions de groupes (Van Krevelen, Bicerano)
- La **simulation numerique** (MD, MC, DPD, SCFT) avec des force fields polymer-first
- Une **application desktop** pour la visualisation 3D et le pilotage interactif
- Des **bindings Python** (PyO3) pour l'integration dans l'ecosysteme scientifique

L'approche est **progressive** : fondations structurelles d'abord, simulation ensuite.

---

## Etat actuel (v0.2.0)

| Statut | Fonctionnalite |
|--------|----------------|
| :white_check_mark: | Homopolymere lineaire (BigSMILES -> SMILES) |
| :white_check_mark: | Strategies : ByRepeatCount, ByTargetMn, ByExactMass |
| :white_check_mark: | Masse molaire moyenne (IUPAC) + monoisotopique |
| :white_check_mark: | Formule moleculaire (notation Hill) |
| :white_check_mark: | Tg — equation de Fox (melanges/copolymeres) |
| :white_check_mark: | Distributions de longueur de chaine (Flory, Schulz-Zimm, Log-Normal) |
| :white_check_mark: | Ensembles polydisperses (Mn, Mw, PDI) |
| :white_check_mark: | CLI : `analyze` + `generate` |
| :construction: | Copolymeres (random, alterne, bloc) — stubbe |
| :construction: | Architectures ramifiees (peigne, greffe, macromonomere) — stubbe |
| :construction: | Tg Van Krevelen, Tm, solubilite — stubbe |

---

## Architecture workspace cible

```
polysim/
  crates/
    polysim-core/         # Structures, builders, proprietes, contributions de groupes
    polysim-geom/         # Topologie, coordonnees 3D, boite de simulation, PBC
    polysim-ff/           # Force fields et potentiels (OPLS-AA, TraPPE-UA, MARTINI)
    polysim-engine/       # MD, MC, DPD, SCFT, Langevin
    polysim-analysis/     # RDF, MSD, Rg, S(q), autocorrelation
    polysim-io/           # PDB, XYZ, LAMMPS, GROMACS, DCD, XTC
    polysim-cg/           # Coarse-graining : mapping, IBI, force-matching
    polysim-cli/          # CLI : analyze, generate, properties, run, init, convert
    polysim-gui/          # Application desktop (framework TBD)
    polysim-python/       # Bindings PyO3
  docs/tutorials/
  tests/validation/
```

---

## Phase 1 — Fondations : Architectures polymeres completes (v0.3 - v0.5)

> Debloquer toutes les architectures est prerequis a tout le reste.
> Sans copolymeres ni ramifications, impossible de predire correctement Tg, Tm, solubilite.

### Epic 1.1 : Copolymeres lineaires

**US-1.1.1 — Copolymere statistique (random)**
> En tant que chimiste, je veux generer un copolymere statistique a partir d'un BigSMILES
> multi-monomeres avec des fractions molaires, afin de modeliser des copolymeres issus de
> polymerisation radicalaire.

- [ ] `LinearBuilder::random_copolymer(&[f64])` genere un SMILES avec placement aleatoire selon les fractions
- [ ] Les fractions doivent sommer a 1.0 (erreur sinon)
- [ ] Support d'un seed pour reproductibilite
- [ ] Calcul correct de Mn pour la chaine resultante
- [ ] Tests : distribution des monomeres converge vers les fractions cibles (N=1000)

**US-1.1.2 — Copolymere alterne**
> En tant que chimiste, je veux generer un copolymere alterne (A-B-A-B) pour modeliser
> des polymeres comme le SMA (styrene-maleic anhydride).

- [ ] `LinearBuilder::alternating_copolymer()` genere une sequence strictement alternee
- [ ] BigSMILES doit contenir exactement 2 unites repetitives
- [ ] SMILES resultat verifie le pattern ABABAB
- [ ] Tests sur PE-alt-PP, PS-alt-PMA

**US-1.1.3 — Copolymere a blocs**
> En tant que chimiste, je veux generer un copolymere a blocs (AAAA-BBBB) avec des longueurs
> de blocs specifiees, pour modeliser des systemes comme PS-b-PEO.

- [ ] `LinearBuilder::block_copolymer(&[usize])` accepte N blocs avec leurs longueurs
- [ ] Support de 2+ blocs (dibloc, tribloc, multibloc)
- [ ] Gestion correcte des ring closures entre blocs
- [ ] Tests : PS-b-PMMA dibloc, PS-b-PB-b-PS tribloc (SBS)

**US-1.1.4 — Copolymere gradient**
> En tant que chimiste, je veux generer un copolymere gradient ou la composition varie
> progressivement le long de la chaine.

- [ ] Nouvelle methode `LinearBuilder::gradient_copolymer(profile: GradientProfile)`
- [ ] Profils : lineaire, sigmoide
- [ ] La fraction du monomere A passe de `f_start` a `f_end` le long de la chaine
- [ ] Verification statistique de la monotonie du gradient

### Epic 1.2 : Architectures ramifiees

**US-1.2.1 — Polymere en peigne (comb)**
> En tant que chimiste, je veux generer un polymere en peigne avec des branches regulierement
> espacees, pour modeliser des architectures comme le PEG-peigne.

- [ ] `BranchedBuilder::comb_polymer(branch_every)` fonctionne
- [ ] Backbone et branches peuvent etre des monomeres differents
- [ ] SMILES resultat contient les branchements corrects
- [ ] Mn calcule inclut backbone + branches

**US-1.2.2 — Copolymere greffe (graft)**
> En tant que chimiste, je veux generer un copolymere greffe avec des greffons places
> aleatoirement, pour modeliser des systemes comme HIPS.

- [ ] `BranchedBuilder::graft_copolymer(graft_fraction)` place les greffons aleatoirement
- [ ] Seed pour reproductibilite
- [ ] La fraction de greffage est respectee statistiquement

**US-1.2.3 — Polymere etoile (star)**
> En tant que chimiste, je veux generer des polymeres etoiles avec un nombre defini de bras.

- [ ] `BranchedBuilder::star_polymer(arms: usize)`
- [ ] Support de 3 a 12 bras
- [ ] Chaque bras : meme monomere (homostar) ou differents (miktostar)

**US-1.2.4 — Dendrimere**
> En tant que chimiste, je veux generer des dendrimeres de generation specifiee, pour
> modeliser des PAMAM ou poly(propylene imine).

- [ ] `BranchedBuilder::dendrimer(generation: usize, branching_factor: usize)`
- [ ] Generations 1 a 6
- [ ] SMILES correct avec gestion des ring closures profondes
- [ ] Comptage du nombre de monomeres par generation

**US-1.2.5 — Polymere cyclique (ring)**
> En tant que chimiste, je veux generer des polymeres cycliques sans extremites de chaine.

- [ ] `LinearBuilder::cyclic_homopolymer()`
- [ ] SMILES forme un cycle (ring closure premier/dernier atome)
- [ ] Mn ne contient pas de masse de groupements terminaux

### Epic 1.3 : Representation enrichie des chaines

**US-1.3.1 — Composition monomere dans PolymerChain**
> En tant que developpeur, je veux enrichir `PolymerChain` pour stocker la composition
> monomere (types et fractions), afin de permettre les calculs de proprietes sur les copolymeres.

- [ ] Nouveau champ : `composition: Vec<MonomerUnit>` avec `MonomerUnit { smiles, fraction }`
- [ ] Retrocompatible : homopolymeres ont une composition a 1 element
- [ ] Accessible via API publique

**US-1.3.2 — Enum Architecture**
> En tant que developpeur, je veux ajouter un champ architecture a `PolymerChain`
> pour conditionner les calculs de proprietes.

- [ ] Enum `Architecture { Linear, Star { arms }, Comb { branch_spacing }, Dendrimer { generation }, Cyclic }`
- [ ] Chaque builder peuple ce champ automatiquement

**US-1.3.3 — Groupes terminaux explicites**
> En tant que chimiste, je veux specifier des groupes terminaux (end groups) explicitement,
> car ils influencent significativement les proprietes a bas Mn.

- [ ] Support des sections `begin`/`end` du BigSMILES dans les builders
- [ ] Prise en compte dans les calculs de MW, formule, et Tg

### Epic 1.4 : CLI — Nouvelles architectures

**US-1.4.1 — CLI analyze etendu pour copolymeres**
> En tant qu'utilisateur CLI, je veux `polysim analyze` etendu supportant
> `--copolymer random|alternating|block|gradient`.

- [ ] `polysim analyze "{[]CC[]; CC(C)[]}" --copolymer random --fractions 0.7,0.3 --by-repeat 100`
- [ ] `polysim analyze "{[]CC[]; CC(c1ccccc1)[]}" --copolymer block --block-lengths 50,50`
- [ ] Output : composition, architecture, proprietes

**US-1.4.2 — CLI generate pour copolymeres**
> En tant qu'utilisateur CLI, je veux que `polysim generate` supporte les ensembles de
> copolymeres pour etudier la polydispersite compositionnelle.

- [ ] Ensembles avec variation de composition et de longueur de chaine

---

## Phase 2 — Predictions de proprietes physiques (v0.6 - v1.0)

> C'est la raison d'etre principale de polysim a court terme — predire les proprietes
> physiques des polymeres a partir de leur structure, sans simulation couteuse.
> Les methodes de contributions de groupes (Van Krevelen, Bicerano) sont le standard industriel.

### Epic 2.1 : Infrastructure des contributions de groupes

**US-2.1.1 — Systeme de contributions de groupes generique**
> En tant que developpeur, je veux un systeme de contributions de groupes qui decompose
> un SMILES en groupes fonctionnels et somme les contributions.

- [ ] Trait `GroupContributionMethod` avec `fn predict(&self, groups: &[Group]) -> f64`
- [ ] Base de donnees de groupes Van Krevelen (~40 groupes pour les polymeres courants)
- [ ] Decomposition SMILES -> groupes fonctionnels (pattern matching ou SMARTS)
- [ ] Tests : decomposition de PE, PP, PS, PMMA, PET, Nylon-6,6

**US-2.1.2 — Base de donnees extensible via TOML**
> En tant que developpeur, je veux etendre la base de donnees de contributions de groupes
> via un fichier TOML externe.

- [ ] Format TOML pour definir de nouveaux groupes et leurs contributions
- [ ] Merge avec la base de donnees par defaut
- [ ] Validation des contributions (intervalles physiquement raisonnables)

### Epic 2.2 : Proprietes thermiques

**US-2.2.1 — Tg par Van Krevelen**
> En tant que chimiste, je veux predire Tg par la methode de Van Krevelen
> (Yg = somme contributions / masse unite repetitive).

- [ ] `tg_van_krevelen(chain)` retourne Tg en K
- [ ] Precision : erreur < 15 K pour PE, PP, PS, PMMA, PVC, PET
- [ ] Reference : Van Krevelen, chap. 6

**US-2.2.2 — Temperature de fusion Tm**
> En tant que chimiste, je veux predire Tm par contributions de groupes.

- [ ] `tm_van_krevelen(chain) -> Option<f64>`
- [ ] Retourne `None` pour les polymeres amorphes
- [ ] Validation croisee : ratio Tg/Tm ~ 2/3 (regle de Boyer-Beaman)
- [ ] Tests : PE (Tm ~ 411 K), PP iso (Tm ~ 449 K), PET (Tm ~ 538 K)

**US-2.2.3 — Tendance a la cristallisation**
> En tant que chimiste, je veux estimer la tendance a la cristallisation d'un polymere.

- [ ] Implementer `crystallization_tendency()` (actuellement `todo!()`)
- [ ] Basee sur symetrie de chaine et difference Tm - Tg
- [ ] Retourne l'enum `CrystallizationTendency` existant

**US-2.2.4 — Coefficient de dilatation thermique et capacite calorifique**
> En tant que chimiste, je veux predire le CTE et le Cp.

- [ ] `thermal_expansion_coefficient(chain) -> f64` (1/K)
- [ ] `specific_heat_capacity(chain) -> f64` (J/g.K)
- [ ] Contributions de groupes Van Krevelen

### Epic 2.3 : Parametres de solubilite

**US-2.3.1 — Parametre de solubilite de Hildebrand**
> En tant que chimiste, je veux calculer le parametre de Hildebrand (delta)
> pour predire la miscibilite polymere-solvant.

- [ ] `hildebrand_solubility_parameter(chain) -> f64` en (MPa)^0.5
- [ ] Base sur Ecoh et Vw par contributions de groupes
- [ ] Tests : PS (delta ~ 18.5), PE (delta ~ 16.2), PMMA (delta ~ 18.6)

**US-2.3.2 — Parametres de solubilite de Hansen**
> En tant que chimiste, je veux calculer les parametres de Hansen (delta_d, delta_p, delta_h)
> pour une analyse de miscibilite plus fine.

- [ ] `hansen_solubility_parameters(chain) -> HansenParams { d, p, h }`
- [ ] Decomposition en contributions dispersives, polaires, H-bond
- [ ] Distance RED entre polymere et solvant

**US-2.3.3 — Parametre d'interaction de Flory-Huggins chi**
> En tant que chimiste, je veux evaluer la miscibilite entre deux polymeres
> via le parametre de Flory-Huggins chi.

- [ ] `flory_huggins_chi(polymer_a, polymer_b, temperature) -> f64`
- [ ] Base sur chi ~ (delta_A - delta_B)^2
- [ ] chi < 0.5 : miscible, chi > 0.5 : immiscible

### Epic 2.4 : Proprietes mecaniques

**US-2.4.1 — Module de Young et resistance a la traction**
> En tant qu'ingenieur materiaux, je veux predire le module de Young et la
> resistance a la traction d'un polymere.

- [ ] `youngs_modulus(chain) -> f64` (GPa) — via Rao function de Van Krevelen
- [ ] `tensile_strength(chain) -> f64` (MPa)
- [ ] Differenciation amorphe vs semi-cristallin
- [ ] Tests : PS (E ~ 3.0 GPa), PE-HD (E ~ 1.0 GPa), PMMA (E ~ 3.3 GPa)

**US-2.4.2 — Densite**
> En tant qu'ingenieur materiaux, je veux estimer la densite a temperature ambiante.

- [ ] `density(chain) -> f64` (g/cm3) — volume Van der Waals + facteur de packing
- [ ] Tests : PE (0.95), PS (1.05), PMMA (1.18), PVC (1.40)

### Epic 2.5 : Proprietes de transport

**US-2.5.1 — Permeabilite aux gaz**
> En tant qu'ingenieur procedes, je veux predire la permeabilite aux gaz
> (O2, CO2, N2, H2O) d'un film polymere pour le packaging.

- [ ] `gas_permeability(chain, gas: Gas) -> f64` (Barrer)
- [ ] Methode Permachor (Van Krevelen) ou correlations Salame
- [ ] Tests : PE (haute perm O2), PET (bonne barriere), EVOH (excellente barriere)

**US-2.5.2 — Viscosite intrinseque**
> En tant que chimiste, je veux predire la viscosite intrinseque [eta]
> via l'equation de Mark-Houwink.

- [ ] `intrinsic_viscosity(chain, params: MarkHouwinkParams) -> f64` (dL/g)
- [ ] Base de donnees des constantes K, a pour les couples polymere-solvant courants

### Epic 2.6 : Proprietes optiques et electriques

**US-2.6.1 — Indice de refraction**
> En tant qu'ingenieur optique, je veux predire l'indice de refraction d'un polymere.

- [ ] `refractive_index(chain) -> f64` — via refraction molaire de Lorentz-Lorenz
- [ ] Tests : PS (n ~ 1.59), PMMA (n ~ 1.49), PC (n ~ 1.585)

**US-2.6.2 — Constante dielectrique et tension de surface**
> En tant qu'ingenieur, je veux predire la constante dielectrique et la tension de surface.

- [ ] `dielectric_constant(chain) -> f64`
- [ ] `surface_tension(chain) -> f64` (mN/m) — via Parachor
- [ ] Tests : PE (gamma ~ 31 mN/m), PS (gamma ~ 40 mN/m)

### Epic 2.7 : CLI — Rapport de proprietes

**US-2.7.1 — Commande `polysim properties`**
> En tant qu'utilisateur CLI, je veux une commande qui genere un rapport complet
> de toutes les proprietes predites pour un polymere donne.

- [ ] `polysim properties "{[]CC(c1ccccc1)[]}" --by-repeat 100`
- [ ] Sections : Thermique, Mecanique, Solubilite, Transport, Optique
- [ ] Sortie : tableau console, JSON (`--format json`), CSV (`--format csv`)
- [ ] Indicateur de confiance par propriete

**US-2.7.2 — Commande `polysim compare`**
> En tant qu'utilisateur CLI, je veux comparer les proprietes de plusieurs polymeres.

- [ ] `polysim compare "{[]CC[]}" "{[]CC(C)[]}" "{[]CC(c1ccccc1)[]}" --by-repeat 100`
- [ ] Tableau comparatif avec toutes les proprietes

---

## Phase 3 — Generation 3D, topologie et application desktop (v1.1 - v1.5)

> Prerequis indispensable pour la simulation numerique (MD/MC).
> La generation de coordonnees 3D transforme polysim d'un outil QSPR en pre-processeur
> de simulation. C'est aussi le moment ou le CLI atteint ses limites visuelles :
> une application desktop devient necessaire pour la visualisation 3D.

### Epic 3.1 : Representation atomistique et topologie

**US-3.1.1 — Structure de donnees Topology**
> En tant que developpeur, je veux une structure Topology representant atomes, liaisons,
> angles, dihedrales, et impropers d'un systeme polymere.

- [ ] `Atom { id, element, mass, charge, position: [f64; 3] }`
- [ ] `Bond { i, j, bond_type }`, `Angle { i, j, k }`, `Dihedral { i, j, k, l }`
- [ ] Methodes de construction incrementale (add_atom, add_bond, etc.)
- [ ] Detection automatique des angles et dihedrales a partir des liaisons

**US-3.1.2 — Conversion PolymerChain (SMILES) -> Topology**
> En tant que developpeur, je veux convertir un `PolymerChain` en Topology
> avec connectivite complete.

- [ ] Parsing SMILES -> graphe moleculaire (atomes + liaisons)
- [ ] Ajout automatique des hydrogenes implicites
- [ ] Detection liaisons, angles, dihedrales
- [ ] Types d'atomes assignes selon la chimie locale

**US-3.1.3 — Boite de simulation avec PBC**
> En tant que developpeur, je veux un systeme de boite de simulation
> avec conditions aux limites periodiques.

- [ ] `SimulationBox { origin, dimensions, periodicity: [bool; 3] }`
- [ ] Support : orthogonale et triclinique
- [ ] Methodes `wrap(position)` et `minimum_image(r_ij)`

### Epic 3.2 : Generation de coordonnees 3D

**US-3.2.1 — Conformation initiale par random walk**
> En tant que chimiste computationnel, je veux generer une conformation initiale
> par marche aleatoire pour des chaines polymeres.

- [ ] Longueur de liaison fixe (C-C = 1.54 A, etc.)
- [ ] Angles de liaison tetraedriques (109.5 deg)
- [ ] Rotation dihedrales aleatoire (ou distribution de Boltzmann)
- [ ] Pas de chevauchement atomique (verification distances minimales)
- [ ] Seed pour reproductibilite

**US-3.2.2 — Conformation par modele RIS (Rotational Isomeric State)**
> En tant que chimiste computationnel, je veux generer des conformations par le modele
> RIS de Flory pour des conformations plus realistes.

- [ ] Matrices de poids statistiques U pour les paires de dihedrales
- [ ] Parametres RIS pour PE, PP, PS, PMMA (litterature)
- [ ] Generation Monte Carlo avec les poids RIS
- [ ] Calcul du rapport C_inf (rapport caracteristique) comme validation

**US-3.2.3 — Packing d'une boite de simulation**
> En tant que simulateur, je veux peupler une boite avec plusieurs chaines
> a une densite cible.

- [ ] `PackingBuilder::new(box_size, density, chains)`
- [ ] Algorithme : insertion aleatoire avec rejet des chevauchements
- [ ] Alternative : methode "grow and push"
- [ ] Densite resultante dans les 5% de la cible

**US-3.2.4 — Minimisation d'energie rapide**
> En tant que simulateur, je veux une minimisation d'energie (steepest descent / CG)
> pour relaxer les geometries initiales.

- [ ] Steepest descent avec line search
- [ ] Convergence sur gradient max ou energie
- [ ] Support PBC
- [ ] Critere d'arret configurable

### Epic 3.3 : I/O — Formats de fichiers

**US-3.3.1 — Export LAMMPS data file**
> En tant que simulateur, je veux exporter au format LAMMPS data file
> pour lancer des simulations LAMMPS.

- [ ] Header complet (atomes, liaisons, angles, dihedrales)
- [ ] Sections Atoms, Bonds, Angles, Dihedrals, Masses
- [ ] Styles : full, molecular
- [ ] Tests round-trip : ecriture -> relecture -> comparaison

**US-3.3.2 — Export GROMACS (.gro + .top)**
> En tant que simulateur, je veux exporter au format GROMACS.

- [ ] Fichier .gro (coordonnees + velocites)
- [ ] Fichier .top (topologie, force field includes)
- [ ] Support PBC

**US-3.3.3 — Export/Import PDB et XYZ**
> En tant que chimiste, je veux exporter/importer les formats PDB et XYZ.

- [ ] PDB : ATOM/HETATM records, CONECT, CRYST1
- [ ] XYZ : format standard
- [ ] Support lecture et ecriture

**US-3.3.4 — Configuration de simulation en TOML**
> En tant que developpeur, je veux un format de configuration TOML
> pour definir completement un calcul.

- [ ] Specification systeme (polymere, densite, boite)
- [ ] Parametres force field
- [ ] Parametres simulation (integrator, dt, nsteps, T, P)
- [ ] Sorties demandees (trajectoire, energie, proprietes)
- [ ] Schema valide et messages d'erreur clairs

### Epic 3.4 : Application desktop

> A partir de la Phase 3, le rendu CLI est limite pour la visualisation 3D.
> L'application desktop permet de visualiser les structures generees, explorer
> les proprietes, et piloter les simulations de facon interactive.
> Le framework (Tauri, egui, iced, etc.) sera decide au moment de l'implementation.

**US-3.4.1 — Visualisation 3D des chaines polymeres**
> En tant que chimiste, je veux visualiser en 3D les chaines polymeres generees
> pour verifier l'architecture et la conformation.

- [ ] Rendu 3D des atomes et liaisons (ball-and-stick ou space-filling)
- [ ] Coloration par type d'atome, par monomere, ou par charge
- [ ] Rotation, zoom, pan interactifs
- [ ] Affichage de la boite de simulation (wireframe)

**US-3.4.2 — Panneau de proprietes**
> En tant que chimiste, je veux voir les proprietes predites dans un panneau lateral
> a cote de la vue 3D.

- [ ] Toutes les proprietes de Phase 2 affichees dans un panneau structuree
- [ ] Mise a jour en temps reel quand le polymere change
- [ ] Export du rapport en JSON/CSV/PDF

**US-3.4.3 — Constructeur interactif de polymeres**
> En tant que chimiste, je veux construire des polymeres interactivement
> en selectionnant monomeres, architecture, longueur, etc.

- [ ] Selecteur de monomeres (bibliotheque des monomeres courants)
- [ ] Choix de l'architecture (lineaire, etoile, peigne, bloc...)
- [ ] Parametres (longueur, fractions, nombre de bras...)
- [ ] Generation instantanee avec preview 3D

**US-3.4.4 — Visualisation de la boite de simulation**
> En tant que simulateur, je veux visualiser une boite de simulation peuplee
> avec ses conditions aux limites periodiques.

- [ ] Rendu de la boite et des images periodiques
- [ ] Selection de chaines individuelles
- [ ] Affichage des distances, angles, dihedrales par clic
- [ ] Animation des images periodiques (ghost atoms optionnels)

**US-3.4.5 — Editeur de configuration de simulation**
> En tant que simulateur, je veux configurer une simulation (force field, integrator,
> thermostat, etc.) via une interface graphique plutot que d'editer du TOML manuellement.

- [ ] Formulaire pour tous les parametres de simulation
- [ ] Validation en temps reel des parametres
- [ ] Preview du fichier TOML genere
- [ ] Bouton "Lancer la simulation"

**US-3.4.6 — Lecteur de trajectoires**
> En tant que simulateur, je veux visualiser les trajectoires de simulation
> frame par frame avec des controles de lecture.

- [ ] Chargement de trajectoires (DCD, XTC, XYZ)
- [ ] Controles : play, pause, avance/recul frame par frame, vitesse
- [ ] Slider temporel
- [ ] Overlay de proprietes (energie, temperature) sur la timeline

---

## Phase 4 — Moteurs de simulation numerique (v2.0+)

> Transformer polysim en un veritable outil de simulation.
> L'ambition long terme : un LAMMPS-like en Rust, specialise polymeres.

### Epic 4.1 : Force fields et potentiels

**US-4.1.1 — Trait Potential generique**
> En tant que developpeur, je veux un trait Potential pour definir des potentiels d'interaction.

- [ ] `trait Potential { fn energy(&self, r: f64) -> f64; fn force(&self, r: f64) -> f64; }`
- [ ] Implementations : LennardJones, Harmonic (bond/angle), CosineDihedral, FENE
- [ ] Calcul vectorise (SIMD-friendly)

**US-4.1.2 — Force field TraPPE-UA**
> En tant que simulateur, je veux TraPPE-UA pour les simulations de polymeres hydrocarbones.

- [ ] Parametres pour CH4, CH3, CH2, CH groupes
- [ ] Liaisons, angles, dihedrales parametres
- [ ] LJ avec cutoff et tail corrections
- [ ] Tests : energie d'un PE court vs reference LAMMPS

**US-4.1.3 — Force field OPLS-AA**
> En tant que simulateur, je veux OPLS-AA pour les simulations all-atom detaillees.

- [ ] Parametres pour C, H, O, N courants dans les polymeres
- [ ] Regles de combinaison (geometric epsilon, arithmetic sigma)
- [ ] Support des charges partielles

**US-4.1.4 — Force field personnalise via fichier de parametres**
> En tant que simulateur, je veux definir un champ de force personnalise via TOML.

- [ ] Format TOML pour types d'atomes, potentiels, parametres
- [ ] Validation des parametres
- [ ] Mixing rules configurable

### Epic 4.2 : Dynamique moleculaire (MD)

**US-4.2.1 — Integrateur velocity Verlet (NVE)**
> En tant que simulateur, je veux un integrateur velocity Verlet pour MD NVE.

- [ ] Conservation de l'energie totale (drift < 10^-4 kT sur 10^6 pas)
- [ ] Pas de temps configurable (1-2 fs typique)
- [ ] Support PBC
- [ ] Cell lists pour O(N) scaling

**US-4.2.2 — Thermostats (NVT)**
> En tant que simulateur, je veux des thermostats pour les simulations NVT.

- [ ] Nose-Hoover chain (deterministe, correct pour les ensembles)
- [ ] Langevin (stochastique, bon pour la relaxation)
- [ ] Velocity rescaling (Bussi-Donadio-Parrinello, canonique correct)
- [ ] Verification : distribution de Maxwell-Boltzmann des vitesses

**US-4.2.3 — Barostats (NPT)**
> En tant que simulateur, je veux des barostats pour les simulations NPT.

- [ ] Berendsen (equilibrage rapide)
- [ ] Parrinello-Rahman (fluctuations correctes)
- [ ] Support isotrope et anisotrope
- [ ] Verification : densite d'equilibre vs experimentale pour PE

**US-4.2.4 — Ecriture de trajectoires**
> En tant que simulateur, je veux l'ecriture de trajectoires pendant la simulation.

- [ ] Format binaire compact propre a polysim
- [ ] Export XYZ (lisible, debug)
- [ ] Export DCD ou XTC (compatible VMD/MDAnalysis)
- [ ] Frequence configurable, compression optionnelle

### Epic 4.3 : Monte Carlo (MC)

**US-4.3.1 — Monte Carlo Metropolis avec mouvements polymeres**
> En tant que simulateur, je veux un moteur MC avec des mouvements adaptes aux polymeres.

- [ ] Mouvements : translation, pivot, reptation, crankshaft
- [ ] Critere de Metropolis correct
- [ ] Taux d'acceptation cible configurable (ajustement automatique)

**US-4.3.2 — Monte Carlo configurationnel-biais (CBMC)**
> En tant que simulateur, je veux le CBMC pour l'insertion et la regrowth de chaines.

- [ ] Rosenbluth sampling pour la croissance de chaines
- [ ] Nombre de trial orientations configurable
- [ ] Correct pour l'ensemble grand-canonique

**US-4.3.3 — Wang-Landau**
> En tant que chercheur, je veux Wang-Landau pour la densite d'etats et transitions de phase.

- [ ] Algorithme 1/t-WL (convergence prouvee)
- [ ] Histogramme d'energie
- [ ] Export g(E)

### Epic 4.4 : Dissipative Particle Dynamics (DPD)

**US-4.4.1 — Moteur DPD**
> En tant que simulateur mesoscopique, je veux DPD pour la morphologie de copolymeres
> a blocs a grande echelle.

- [ ] Forces conservatives, dissipatives, aleatoires (Groot-Warren)
- [ ] Thermostat DPD integre
- [ ] Integrateur modified velocity Verlet
- [ ] Parametres aij mappes depuis chi (Flory-Huggins)

---

## Phase 5 — Analyse et post-traitement (v2.1+)

> Extraire des proprietes physiques a partir des trajectoires de simulation.

### Epic 5.1 : Analyse structurale

**US-5.1.1 — Fonction de distribution radiale g(r)**
> En tant qu'analyste, je veux calculer la RDF a partir d'une trajectoire.

- [ ] RDF totale et par paires de types d'atomes
- [ ] Support PBC (minimum image)
- [ ] Moyennage sur les frames
- [ ] Normalisation correcte (densite ideale)

**US-5.1.2 — Rayon de giration Rg et distance bout-a-bout R_ee**
> En tant qu'analyste, je veux calculer Rg et R_ee des chaines polymeres.

- [ ] Rg instantane et moyenne sur la trajectoire
- [ ] R_ee et <R_ee^2>
- [ ] Rapport <R_ee^2> / <Rg^2> (theoriquement 6 pour une chaine ideale)
- [ ] Par chaine et moyenne ensemble

**US-5.1.3 — Facteur de structure S(q)**
> En tant qu'analyste, je veux calculer S(q) pour comparer avec SAXS/SANS.

- [ ] Methode Debye (petits systemes)
- [ ] Methode FFT (grands systemes)
- [ ] Plage de q configurable

**US-5.1.4 — Longueur de persistance et longueur de Kuhn**
> En tant qu'analyste, je veux calculer lp et b des chaines.

- [ ] Autocorrelation des vecteurs liaison
- [ ] Fit exponentiel pour lp
- [ ] Kuhn length b = 2 * lp

### Epic 5.2 : Analyse dynamique

**US-5.2.1 — Deplacement carre moyen (MSD) et coefficient de diffusion**
> En tant qu'analyste, je veux calculer le MSD et D.

- [ ] MSD par type d'atome et centre de masse des chaines
- [ ] Multiple time origin averaging
- [ ] Extraction de D par regression lineaire (regime diffusif)
- [ ] Detection automatique regime balistique/Rouse/reptation

**US-5.2.2 — Fonctions d'autocorrelation**
> En tant qu'analyste, je veux calculer des ACFs (vitesses, dipoles, end-to-end).

- [ ] VACF -> spectre vibrationnel
- [ ] End-to-end vector ACF -> temps de relaxation Rouse/reptation
- [ ] FFT pour l'efficacite

### Epic 5.3 : Analyse thermomecanique depuis la simulation

**US-5.3.1 — Determination de Tg par simulation**
> En tant que chercheur, je veux determiner Tg en tracant volume specifique vs temperature.

- [ ] Protocole de refroidissement automatise (NPT multi-T)
- [ ] Detection changement de pente (regression bi-segments)
- [ ] Comparaison avec Tg predit par contributions de groupes

**US-5.3.2 — Module de cisaillement G(t) et viscosite (Green-Kubo)**
> En tant que chercheur, je veux calculer G(t) et la viscosite par Green-Kubo.

- [ ] Autocorrelation du tenseur des contraintes
- [ ] Integration pour la viscosite eta
- [ ] Support des longs temps de correlation

### Epic 5.4 : CLI — Analyse

**US-5.4.1 — Commande `polysim analyze-trajectory`**
> En tant qu'utilisateur CLI, je veux analyser une trajectoire.

- [ ] `polysim analyze-trajectory traj.dcd --topology system.top --rdf --msd --rg`
- [ ] Sortie CSV ou JSON
- [ ] Progress bar pour les longues analyses

### Epic 5.5 : Desktop app — Analyse

**US-5.5.1 — Visualisation des resultats d'analyse**
> En tant qu'utilisateur, je veux visualiser les courbes d'analyse (RDF, MSD, S(q))
> directement dans l'application desktop.

- [ ] Graphiques interactifs (zoom, pan, selection de plage)
- [ ] Overlay de plusieurs courbes (comparaison entre runs)
- [ ] Export PNG/SVG des graphiques

**US-5.5.2 — Monitoring en temps reel de simulation**
> En tant que simulateur, je veux suivre l'evolution d'une simulation en cours
> avec des graphiques live.

- [ ] Temperature, pression, energie totale en temps reel
- [ ] Alerte si la simulation diverge (energie explose)
- [ ] Bouton stop/pause depuis l'interface

---

## Phase 6 — Performance et ecosysteme (v2.5+)

> Rendre polysim competitif en performance et accessible a la communaute scientifique.

### Epic 6.1 : Performance et parallelisme

**US-6.1.1 — Parallelisme shared-memory via Rayon**
> En tant que developpeur, je veux Rayon pour les boucles de force et l'analyse.

- [ ] Calcul des forces parallelise par cell lists
- [ ] Analyse (RDF, MSD) parallelisee sur les frames
- [ ] Speedup lineaire jusqu'a 8-16 coeurs

**US-6.1.2 — Listes de voisins optimisees (cell list + Verlet list)**
> En tant que developpeur, je veux des neighbor lists optimisees.

- [ ] Cell list : O(N) construction
- [ ] Verlet list : mise a jour conditionnelle (skin distance)
- [ ] Benchmark : > 10^4 atomes a 1 ns/jour sur un desktop

**US-6.1.3 — Vectorisation SIMD pour les forces**
> En tant que developpeur, je veux SIMD pour les forces pair-wise.

- [ ] Layout SoA (Structure of Arrays) pour les positions
- [ ] Intrinsics via `std::simd` ou `packed_simd`
- [ ] Benchmark : gain de 2-4x sur les forces LJ

### Epic 6.2 : Bindings Python (PyO3)

**US-6.2.1 — Python bindings pour structures et proprietes**
> En tant que chimiste Python, je veux utiliser polysim depuis Python.

- [ ] `pip install polysim`
- [ ] `from polysim import PolymerChain, LinearBuilder, BuildStrategy`
- [ ] API Pythonique (snake_case, docstrings, type hints)
- [ ] Support NumPy pour les coordonnees 3D
- [ ] Integration RDKit (conversion SMILES <-> Mol)

**US-6.2.2 — Python bindings pour simulation MD**
> En tant que chercheur Python, je veux lancer des simulations MD depuis Python.

- [ ] `from polysim.md import Simulation`
- [ ] Configuration Python, execution Rust
- [ ] Callbacks Python pour le monitoring
- [ ] Interoperabilite MDAnalysis

### Epic 6.3 : Documentation et tutoriels

**US-6.3.1 — Tutoriels pas-a-pas**
> En tant que nouvel utilisateur, je veux des tutoriels couvrant les cas principaux.

- [ ] Tutoriel 1 : Generation d'un copolymere a blocs PS-b-PEO + prediction proprietes
- [ ] Tutoriel 2 : Ensemble polydisperse + analyse statistique
- [ ] Tutoriel 3 : Construction boite PE + run MD
- [ ] Tutoriel 4 : Tg par simulation vs prediction

**US-6.3.2 — Documentation API exhaustive**
> En tant que developpeur Rust, je veux une doc API avec exemples.

- [ ] `#[doc]` sur toutes les fonctions/structs/traits publics
- [ ] Exemples runnable pour chaque methode cle
- [ ] `cargo doc` sans warnings

**US-6.3.3 — Guide de contribution et benchmark suite**
> En tant que contributeur, je veux un guide de contribution et des benchmarks.

- [ ] CONTRIBUTING.md
- [ ] Benchmark suite (PE melt, PS glass, PS-b-PMMA morphology)
- [ ] Comparaison automatisee avec LAMMPS

### Epic 6.4 : CLI avancee

**US-6.4.1 — Commande `polysim run`**
> En tant qu'utilisateur CLI, je veux lancer une simulation depuis un fichier TOML.

- [ ] `polysim run simulation.toml`
- [ ] Progress bar, ETA, statistiques temps reel
- [ ] Checkpointing automatique
- [ ] `polysim continue checkpoint.bin`

**US-6.4.2 — Commande `polysim init`**
> En tant qu'utilisateur CLI, je veux generer un template de configuration.

- [ ] `polysim init --polymer "{[]CC[]}" --ensemble nvt --temperature 300`
- [ ] TOML complet avec valeurs par defaut raisonnables
- [ ] Commentaires explicatifs

---

## Matrice de priorites

### Phase 1 — Fondations (v0.3 - v0.5)

| Priorite | US | Description |
|----------|----|-------------|
| P0 | US-1.3.1 | Composition monomere dans PolymerChain |
| P0 | US-1.3.2 | Enum Architecture |
| P0 | US-1.1.1 | Copolymere random |
| P0 | US-1.1.2 | Copolymere alterne |
| P1 | US-1.1.3 | Copolymere a blocs |
| P1 | US-1.3.3 | Groupes terminaux |
| P1 | US-1.2.1 | Polymere en peigne |
| P1 | US-1.2.3 | Polymere etoile |
| P1 | US-1.4.1 | CLI copolymeres |
| P2 | US-1.1.4 | Copolymere gradient |
| P2 | US-1.2.2 | Copolymere greffe |
| P2 | US-1.2.4 | Dendrimere |
| P2 | US-1.2.5 | Polymere cyclique |
| P2 | US-1.4.2 | CLI generate copolymeres |

### Phase 2 — Proprietes (v0.6 - v1.0)

| Priorite | US | Description |
|----------|----|-------------|
| P0 | US-2.1.1 | Systeme contributions de groupes |
| P0 | US-2.2.1 | Tg Van Krevelen |
| P0 | US-2.2.2 | Tm Van Krevelen |
| P0 | US-2.3.1 | Hildebrand delta |
| P0 | US-2.4.2 | Densite |
| P1 | US-2.2.3 | Cristallisation |
| P1 | US-2.3.2 | Hansen 3D |
| P1 | US-2.4.1 | Module de Young |
| P1 | US-2.6.1 | Indice de refraction |
| P1 | US-2.7.1 | CLI properties |
| P2 | US-2.1.2 | Base extensible TOML |
| P2 | US-2.2.4 | CTE, Cp |
| P2 | US-2.3.3 | Flory-Huggins chi |
| P2 | US-2.5.1 | Permeabilite gaz |
| P2 | US-2.5.2 | Viscosite intrinseque |
| P2 | US-2.6.2 | Dielectrique, tension surface |
| P2 | US-2.7.2 | CLI compare |

### Phase 3 — 3D + Desktop (v1.1 - v1.5)

| Priorite | US | Description |
|----------|----|-------------|
| P0 | US-3.1.1 | Structure Topology |
| P0 | US-3.1.2 | SMILES -> Topology |
| P0 | US-3.1.3 | Boite de simulation PBC |
| P0 | US-3.2.1 | Random walk |
| P0 | US-3.4.1 | Visualisation 3D (desktop) |
| P1 | US-3.2.2 | Modele RIS |
| P1 | US-3.2.3 | Packing boite |
| P1 | US-3.3.1 | Export LAMMPS |
| P1 | US-3.3.4 | Config TOML |
| P1 | US-3.4.2 | Panneau proprietes (desktop) |
| P1 | US-3.4.3 | Constructeur interactif (desktop) |
| P2 | US-3.2.4 | Minimisation energie |
| P2 | US-3.3.2 | Export GROMACS |
| P2 | US-3.3.3 | Export PDB/XYZ |
| P2 | US-3.4.4 | Visualisation boite (desktop) |
| P2 | US-3.4.5 | Editeur config (desktop) |
| P2 | US-3.4.6 | Lecteur trajectoires (desktop) |

### Phase 4 — Simulation (v2.0+)

| Priorite | US | Description |
|----------|----|-------------|
| P0 | US-4.1.1 | Trait Potential |
| P0 | US-4.1.2 | TraPPE-UA |
| P0 | US-4.2.1 | Velocity Verlet NVE |
| P0 | US-4.2.2 | Thermostats NVT |
| P1 | US-4.2.3 | Barostats NPT |
| P1 | US-4.2.4 | Trajectoires |
| P1 | US-4.1.3 | OPLS-AA |
| P1 | US-4.3.1 | MC Metropolis |
| P2 | US-4.1.4 | Force field custom |
| P2 | US-4.3.2 | CBMC |
| P2 | US-4.3.3 | Wang-Landau |
| P2 | US-4.4.1 | DPD |

### Phase 5 — Analyse (v2.1+)

| Priorite | US | Description |
|----------|----|-------------|
| P0 | US-5.1.1 | RDF g(r) |
| P0 | US-5.1.2 | Rg et R_ee |
| P0 | US-5.2.1 | MSD et coefficient de diffusion |
| P1 | US-5.1.3 | Facteur de structure S(q) |
| P1 | US-5.1.4 | Persistance et Kuhn |
| P1 | US-5.2.2 | Autocorrelation |
| P1 | US-5.4.1 | CLI analyze-trajectory |
| P1 | US-5.5.1 | Visualisation analyse (desktop) |
| P2 | US-5.3.1 | Tg par simulation |
| P2 | US-5.3.2 | Green-Kubo viscosite |
| P2 | US-5.5.2 | Monitoring temps reel (desktop) |

### Phase 6 — Ecosysteme (v2.5+)

| Priorite | US | Description |
|----------|----|-------------|
| P0 | US-6.1.1 | Rayon parallelisme |
| P0 | US-6.1.2 | Neighbor lists |
| P1 | US-6.2.1 | Python bindings (structures) |
| P1 | US-6.3.1 | Tutoriels |
| P1 | US-6.3.2 | API docs |
| P1 | US-6.4.1 | CLI run |
| P2 | US-6.1.3 | SIMD |
| P2 | US-6.2.2 | Python bindings (simulation) |
| P2 | US-6.3.3 | Contributing guide |
| P2 | US-6.4.2 | CLI init |

---

## Differenciation vs LAMMPS / GROMACS / HOOMD-blue

1. **Securite Rust** : Pas de segfault dans les neighbor lists, pas de race conditions, `Result<T, E>` au lieu de NaN silencieux
2. **Polymere-first** : BigSMILES natif, architectures et distributions comme citoyens de premiere classe
3. **Prediction integree** : Seul outil combinant contributions de groupes ET simulation atomistique/CG (Materials Studio facture des dizaines de milliers d'euros)
4. **Distribution single-binary** : `cargo install polysim-cli` — pas de conda
5. **Multiscale** : Atomistique -> CG -> DPD -> SCFT avec structures de donnees coherentes
6. **Application desktop integree** : Visualisation 3D, configuration, monitoring — sans logiciel tiers

---

## Validation

### Polymeres de reference
PE, PP, PS, PMMA, PET, Nylon-6,6, PEO, PDMS, PC, PVC

### Precision cible
| Propriete | Precision vs experimental |
|-----------|--------------------------|
| Tg | < 20 K |
| Tm | < 25 K |
| Densite | < 5% |
| Parametre de solubilite | < 0.5 (cal/cm3)^0.5 |
| Module de Young | < 20% |

### Criteres de qualite
- `cargo test --all-features` — tous les tests passent
- `cargo clippy --all-targets -- -D warnings` — zero warning
- `cargo bench` — pas de regression de performance
- Tests de validation numerique contre donnees publiees (Van Krevelen, NIST)