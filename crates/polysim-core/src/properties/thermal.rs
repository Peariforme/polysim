use crate::error::PolySimError;
use crate::polymer::PolymerChain;
use crate::properties::group_contribution::{sum_contribution, GroupDatabase};
use crate::properties::molecular_weight::average_mass;

/// Estimates the glass transition temperature (K) using the Fox equation.
///
/// # Arguments
///
/// `components` — slice of `(weight_fraction, Tg_homopolymer_K)` pairs,
/// one per distinct repeat unit. Weight fractions must sum to 1.0.
///
/// # Reference
///
/// Fox, T. G. (1956). *Bull. Am. Phys. Soc.* **1**, 123.
///
/// # Example
///
/// ```rust
/// use polysim_core::properties::thermal::tg_fox;
///
/// // 50/50 blend of PS (Tg ≈ 373 K) and PMMA (Tg ≈ 378 K)
/// let tg = tg_fox(&[(0.5, 373.0), (0.5, 378.0)]);
/// assert!((tg - 375.4).abs() < 0.2);
/// ```
pub fn tg_fox(components: &[(f64, f64)]) -> f64 {
    let inv_tg: f64 = components.iter().map(|(wi, tgi)| wi / tgi).sum();
    1.0 / inv_tg
}

/// Estimates Tg (K) using the Van Krevelen group-contribution method.
///
/// The glass transition temperature is computed as:
///
/// **Tg = (sum(Ygi) * 1000) / M0**
///
/// where `Ygi` are the molar Tg contributions of each functional group
/// (in K * kg/mol, VK convention) and `M0` is the molar mass of the
/// repeat unit (g/mol).
///
/// # Errors
///
/// Returns `PolySimError::GroupDecomposition` if the SMILES cannot be decomposed.
///
/// # Reference
///
/// Van Krevelen, D. W. & te Nijenhuis, K. (2009).
/// *Properties of Polymers*, 4th ed., Elsevier. Chapter 6.
pub fn tg_van_krevelen(chain: &PolymerChain) -> Result<f64, PolySimError> {
    let groups = GroupDatabase::decompose(chain)?;
    let yg_total = sum_contribution(&groups, |g| g.yg);

    let m0 = repeat_unit_mass(chain)?;

    // Yg values are in K*kg/mol (VK convention), M0 in g/mol.
    let yg_per_repeat = yg_total / chain.repeat_count.max(1) as f64;
    Ok((yg_per_repeat * 1000.0) / m0)
}

/// Estimates Tm (K) using the Van Krevelen group-contribution method.
///
/// The melting temperature is computed as:
///
/// **Tm = (sum(Ymi) * 1000) / M0**
///
/// Returns `None` if the predicted Tm is below 200 K, indicating an
/// amorphous polymer with no meaningful melting point.
///
/// # Errors
///
/// Returns `PolySimError::GroupDecomposition` if the SMILES cannot be decomposed.
///
/// # Reference
///
/// Van Krevelen, D. W. & te Nijenhuis, K. (2009).
/// *Properties of Polymers*, 4th ed., Elsevier. Chapter 7.
pub fn tm_van_krevelen(chain: &PolymerChain) -> Result<Option<f64>, PolySimError> {
    let groups = GroupDatabase::decompose(chain)?;
    let ym_total = sum_contribution(&groups, |g| g.ym);

    let m0 = repeat_unit_mass(chain)?;

    let ym_per_repeat = ym_total / chain.repeat_count.max(1) as f64;
    let tm = (ym_per_repeat * 1000.0) / m0;

    // Below 200 K is considered amorphous (no meaningful Tm).
    if tm < 200.0 {
        Ok(None)
    } else {
        Ok(Some(tm))
    }
}

/// Computes the average molar mass of a single repeat unit (g/mol).
fn repeat_unit_mass(chain: &PolymerChain) -> Result<f64, PolySimError> {
    let m_chain = average_mass(chain);
    let n = chain.repeat_count.max(1) as f64;
    let m0 = m_chain / n;

    if m0 < f64::EPSILON {
        return Err(PolySimError::GroupDecomposition(
            "repeat unit mass is zero".into(),
        ));
    }
    Ok(m0)
}

/// Qualitative tendency of a polymer chain to crystallise.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CrystallizationTendency {
    /// Highly regular chain — expected to crystallise readily (e.g. PE, isotactic PP).
    High,
    /// Moderate regularity — partial crystallisation possible (e.g. syndiotactic PS).
    Medium,
    /// Low regularity — unlikely to crystallise significantly.
    Low,
    /// Fully amorphous — no crystallisation expected (e.g. atactic PS, PMMA).
    Amorphous,
}

/// Estimates the crystallisation tendency of a polymer chain based on
/// the difference between predicted Tm and Tg.
///
/// **Classification rules (Van Krevelen heuristic):**
///
/// | Condition                            | Tendency    |
/// |--------------------------------------|-------------|
/// | Tm is `None` (amorphous)             | `Amorphous` |
/// | `Tm - Tg > 100 K` and high symmetry  | `High`      |
/// | `Tm - Tg > 50 K`                     | `Medium`    |
/// | `Tm - Tg > 0 K`                      | `Low`       |
/// | else                                 | `Amorphous` |
///
/// Symmetry is estimated by counting the number of distinct heavy-atom
/// substituent types on the backbone carbons: fewer substituent types
/// imply a more regular, symmetric chain.
///
/// # Panics
///
/// Returns `Amorphous` if Tg or Tm computation fails (e.g. unrecognised groups).
pub fn crystallization_tendency(chain: &PolymerChain) -> CrystallizationTendency {
    let tg = match tg_van_krevelen(chain) {
        Ok(v) => v,
        Err(_) => return CrystallizationTendency::Amorphous,
    };
    let tm_opt = match tm_van_krevelen(chain) {
        Ok(v) => v,
        Err(_) => return CrystallizationTendency::Amorphous,
    };

    let tm = match tm_opt {
        Some(t) => t,
        None => return CrystallizationTendency::Amorphous,
    };

    let delta = tm - tg;

    if delta > 100.0 && has_high_symmetry(chain) {
        CrystallizationTendency::High
    } else if delta > 50.0 {
        CrystallizationTendency::Medium
    } else if delta > 0.0 {
        CrystallizationTendency::Low
    } else {
        CrystallizationTendency::Amorphous
    }
}

/// Heuristic for backbone symmetry: a chain is considered highly symmetric
/// if the repeat unit contains only C and H atoms (no heteroatom substituents).
fn has_high_symmetry(chain: &PolymerChain) -> bool {
    use crate::properties::group_contribution::GroupDatabase;
    let groups = match GroupDatabase::decompose(chain) {
        Ok(g) => g,
        Err(_) => return false,
    };
    // High symmetry = only aliphatic carbon groups (CH3, CH2, CH, >C<).
    // Any polar/aromatic group breaks high symmetry.
    let aliphatic_names = ["-CH3", "-CH2-", "-CH<", ">C<"];
    groups
        .iter()
        .all(|gm| aliphatic_names.contains(&gm.group.name))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::{linear::LinearBuilder, BuildStrategy};

    fn build_chain(bigsmiles: &str, n: usize) -> PolymerChain {
        let bs = crate::parse(bigsmiles).expect("valid BigSMILES");
        LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(n))
            .homopolymer()
            .expect("build should succeed")
    }

    // --- Tg tests ---

    #[test]
    fn tg_vk_polyethylene() {
        let chain = build_chain("{[]CC[]}", 50);
        let tg = tg_van_krevelen(&chain).unwrap();
        let error_pct = ((tg - 195.0) / 195.0).abs() * 100.0;
        assert!(
            error_pct < 20.0,
            "PE Tg = {tg:.1} K, error = {error_pct:.1}%"
        );
    }

    #[test]
    fn tg_vk_pvc() {
        let chain = build_chain("{[]C(Cl)C[]}", 50);
        let tg = tg_van_krevelen(&chain).unwrap();
        let error_pct = ((tg - 354.0) / 354.0).abs() * 100.0;
        assert!(
            error_pct < 20.0,
            "PVC Tg = {tg:.1} K, error = {error_pct:.1}%"
        );
    }

    #[test]
    fn tg_vk_returns_positive() {
        let chain = build_chain("{[]CC[]}", 10);
        let tg = tg_van_krevelen(&chain).unwrap();
        assert!(tg > 0.0, "Tg should be positive, got {tg}");
    }

    // --- Tm tests ---

    #[test]
    fn tm_vk_polyethylene() {
        // PE: Tm exp ~ 411 K
        let chain = build_chain("{[]CC[]}", 50);
        let tm = tm_van_krevelen(&chain).unwrap();
        assert!(tm.is_some(), "PE should have a Tm");
        let tm = tm.unwrap();
        let error_pct = ((tm - 411.0) / 411.0).abs() * 100.0;
        assert!(
            error_pct < 25.0,
            "PE Tm = {tm:.1} K, error = {error_pct:.1}%"
        );
    }

    #[test]
    fn tm_vk_returns_some_for_crystallizable() {
        // PE is crystallizable
        let chain = build_chain("{[]CC[]}", 50);
        let tm = tm_van_krevelen(&chain).unwrap();
        assert!(tm.is_some(), "PE should return Some(Tm)");
        assert!(tm.unwrap() > 200.0, "PE Tm should be > 200 K");
    }

    #[test]
    fn tm_vk_positive_when_present() {
        let chain = build_chain("{[]C(Cl)C[]}", 50);
        let tm = tm_van_krevelen(&chain).unwrap();
        if let Some(t) = tm {
            assert!(t > 0.0, "Tm should be positive, got {t}");
        }
    }
}
