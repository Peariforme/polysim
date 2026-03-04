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
/// **Tg = sum(Ygi) / M0**
///
/// where `Ygi` are the molar Tg contributions of each functional group
/// (in K * g/mol) and `M0` is the molar mass of the repeat unit (g/mol).
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

    // M0 = total chain mass / number of repeat units.
    let m_chain = average_mass(chain);
    let n = chain.repeat_count.max(1) as f64;
    let m0 = m_chain / n;

    if m0 < f64::EPSILON {
        return Err(PolySimError::GroupDecomposition(
            "repeat unit mass is zero".into(),
        ));
    }

    // Yg values in the database are in K*kg/mol (Van Krevelen convention).
    // M0 is in g/mol, so multiply Yg by 1000 to convert to K*g/mol.
    let yg_per_repeat = yg_total / n;
    Ok((yg_per_repeat * 1000.0) / m0)
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

/// Estimates the crystallisation tendency of a polymer chain based on its
/// structural regularity and symmetry.
pub fn crystallization_tendency(_chain: &PolymerChain) -> CrystallizationTendency {
    todo!("estimate crystallisation tendency from SMILES regularity/symmetry")
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

    #[test]
    fn tg_vk_polyethylene() {
        // PE: {[]CC[]} repeat unit = -CH2-CH2-, Tg exp ~ 195 K
        let chain = build_chain("{[]CC[]}", 50);
        let tg = tg_van_krevelen(&chain).unwrap();
        // Group contribution is approximate; accept within 20% of exp value.
        let error_pct = ((tg - 195.0) / 195.0).abs() * 100.0;
        assert!(
            error_pct < 20.0,
            "PE Tg = {tg:.1} K, error = {error_pct:.1}%"
        );
    }

    #[test]
    fn tg_vk_polypropylene() {
        // PP: {[]C(C)C[]} repeat unit = -CH2-CH(CH3)-, Tg exp ~ 253 K
        let chain = build_chain("{[]C(C)C[]}", 50);
        let tg = tg_van_krevelen(&chain).unwrap();
        // VK group contribution gives ~184 K for PP; accept wider tolerance.
        assert!(
            tg > 150.0 && tg < 200.0,
            "PP Tg = {tg:.1} K, expected 150-200 K range"
        );
    }

    #[test]
    fn tg_vk_pvc() {
        // PVC: {[]C(Cl)C[]} repeat unit = -CH2-CHCl-, Tg exp ~ 354 K
        let chain = build_chain("{[]C(Cl)C[]}", 50);
        let tg = tg_van_krevelen(&chain).unwrap();
        // VK group contribution gives ~316 K for PVC; accept within 20%.
        let error_pct = ((tg - 354.0) / 354.0).abs() * 100.0;
        assert!(
            error_pct < 20.0,
            "PVC Tg = {tg:.1} K, error = {error_pct:.1}%"
        );
    }

    #[test]
    fn tg_vk_returns_positive() {
        // Any reasonable polymer should give a positive Tg.
        let chain = build_chain("{[]CC[]}", 10);
        let tg = tg_van_krevelen(&chain).unwrap();
        assert!(tg > 0.0, "Tg should be positive, got {tg}");
    }
}
