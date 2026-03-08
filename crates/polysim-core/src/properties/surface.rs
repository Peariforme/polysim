//! Surface energy and dielectric property calculations.

use crate::error::PolySimError;
use crate::polymer::PolymerChain;
use crate::properties::group_contribution::{total_parachor, total_pe, total_vw, GroupDatabase};
use crate::properties::mechanical::density;
use crate::properties::molecular_weight::average_mass;

/// Estimates the surface tension (mN/m) using the Parachor method (Sugden).
///
/// The Parachor [P] is a group-additive quantity defined by:
///
/// **[P] = V_molar * gamma^(1/4)**
///
/// Rearranged:
///
/// **gamma = ([P] / V_molar)^4**
///
/// where `V_molar` (cm^3/mol) is the molar volume per repeat unit and
/// `[P]` is the sum of group Parachor contributions.
///
/// # Errors
///
/// Returns `PolySimError::GroupDecomposition` if the SMILES cannot be decomposed.
///
/// # Reference
///
/// Sugden, S. (1924). *J. Chem. Soc.*, **125**, 1177.
///
/// Van Krevelen, D. W. & te Nijenhuis, K. (2009).
/// *Properties of Polymers*, 4th ed., Elsevier. Chapter 8.
pub fn surface_tension(chain: &PolymerChain) -> Result<f64, PolySimError> {
    let groups = GroupDatabase::decompose(chain)?;

    let n = chain.repeat_count.max(1) as f64;

    // Parachor per repeat unit
    let p_total = total_parachor(&groups);
    let p = p_total / n;

    // Molar volume per repeat unit: V = M0 / rho
    let rho = density(chain)?;
    let m0 = average_mass(chain) / n;
    let v = m0 / rho;

    if v < f64::EPSILON {
        return Err(PolySimError::GroupDecomposition(
            "molar volume is zero".into(),
        ));
    }

    // gamma = (P / V)^4  (mN/m since Parachor units are calibrated for this)
    let gamma = (p / v).powi(4);
    Ok(gamma)
}

/// Estimates the static dielectric constant (relative permittivity) using the
/// Clausius-Mossotti relation with group-contribution molar polarization.
///
/// The Clausius-Mossotti equation:
///
/// **(epsilon - 1) / (epsilon + 2) = Pe / V**
///
/// Rearranged to:
///
/// **epsilon = (V + 2 * Pe) / (V - Pe)**
///
/// where `Pe` (cm^3/mol) is the sum of molar electronic polarization
/// contributions and `V` (cm^3/mol) is the molar volume per repeat unit.
///
/// Note: this approximation gives the *electronic* (optical frequency)
/// contribution only (epsilon ≈ n^2 for nonpolar polymers). For polar polymers
/// the static dielectric constant includes orientation polarization and will
/// be systematically underestimated by this method.
///
/// # Errors
///
/// Returns `PolySimError::GroupDecomposition` if the SMILES cannot be
/// decomposed, or if the Pe/V ratio is unphysical (>= 1/3).
///
/// # Reference
///
/// Van Krevelen, D. W. & te Nijenhuis, K. (2009).
/// *Properties of Polymers*, 4th ed., Elsevier. Chapter 11.
pub fn dielectric_constant(chain: &PolymerChain) -> Result<f64, PolySimError> {
    let groups = GroupDatabase::decompose(chain)?;

    let n = chain.repeat_count.max(1) as f64;
    let vw_total = total_vw(&groups);
    let vw_per_repeat = vw_total / n;

    if vw_per_repeat < f64::EPSILON {
        return Err(PolySimError::GroupDecomposition(
            "Van der Waals volume per repeat unit is zero".into(),
        ));
    }

    // Pe per repeat unit
    let pe_total = total_pe(&groups);
    let pe = pe_total / n;

    // Molar volume per repeat unit
    let rho = density(chain)?;
    let m0 = average_mass(chain) / n;
    let v = m0 / rho;

    if v < f64::EPSILON {
        return Err(PolySimError::GroupDecomposition(
            "molar volume is zero".into(),
        ));
    }

    let ratio = pe / v;

    // Guard against unphysical Pe/V >= 1/3 (Clausius-Mossotti diverges at 1/3)
    if ratio >= 1.0 / 3.0 {
        return Err(PolySimError::GroupDecomposition(format!(
            "Pe/V = {ratio:.4} >= 1/3, unphysical dielectric constant"
        )));
    }

    // epsilon = (1 + 2*r) / (1 - r)  where r = Pe/V
    let eps = (1.0 + 2.0 * ratio) / (1.0 - ratio);
    Ok(eps)
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
    fn surface_tension_pe() {
        // PE: gamma exp ~ 31 mN/m. VK Parachor overestimates → ~50 mN/m
        let chain = build_chain("{[]CC[]}", 50);
        let gamma = surface_tension(&chain).unwrap();
        assert!(
            gamma > 20.0 && gamma < 70.0,
            "PE surface tension = {gamma:.1} mN/m, expected in [20, 70]"
        );
    }

    #[test]
    fn surface_tension_ps() {
        // PS: gamma exp ~ 40 mN/m
        let chain = build_chain("{[]CC(c1ccccc1)[]}", 50);
        let gamma = surface_tension(&chain).unwrap();
        assert!(
            gamma > 25.0 && gamma < 60.0,
            "PS surface tension = {gamma:.1} mN/m, expected ~40"
        );
    }

    #[test]
    fn surface_tension_positive() {
        let chain = build_chain("{[]CC[]}", 10);
        let gamma = surface_tension(&chain).unwrap();
        assert!(gamma > 0.0, "surface tension must be positive, got {gamma}");
    }

    #[test]
    fn dielectric_constant_pe() {
        // PE: epsilon exp ~ 2.3 (apolar)
        let chain = build_chain("{[]CC[]}", 50);
        let eps = dielectric_constant(&chain).unwrap();
        assert!(
            eps > 1.5 && eps < 4.0,
            "PE dielectric constant = {eps:.2}, expected ~2.3"
        );
    }

    #[test]
    fn dielectric_constant_greater_than_one() {
        let chain = build_chain("{[]CC[]}", 10);
        let eps = dielectric_constant(&chain).unwrap();
        assert!(eps > 1.0, "dielectric constant must be > 1, got {eps:.3}");
    }

    #[test]
    fn dielectric_constant_physical_range() {
        let polymers = [
            ("{[]CC[]}", "PE"),
            ("{[]CC(C)[]}", "PP"),
            ("{[]CC(c1ccccc1)[]}", "PS"),
            ("{[]CC(C)(C(=O)OC)[]}", "PMMA"),
            ("{[]C(Cl)C[]}", "PVC"),
        ];
        for (bigsmiles, name) in polymers {
            let chain = build_chain(bigsmiles, 50);
            let eps = dielectric_constant(&chain).unwrap();
            assert!(
                eps > 1.5 && eps < 6.0,
                "{name}: dielectric constant = {eps:.2}, out of [1.5, 6.0]"
            );
        }
    }
}
