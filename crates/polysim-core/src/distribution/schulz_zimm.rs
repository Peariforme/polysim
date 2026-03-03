use rand::RngCore;
use rand_distr::{Distribution, Gamma};

use super::ChainLengthDistribution;

/// Schulz-Zimm (gamma) chain length distribution.
///
/// Generalization of the Flory distribution with shape parameter k:
///   k = 1/(PDI − 1),  θ = Xn/k
///
/// Covers a wide range of dispersities:
/// - k → ∞ : monodisperse (PDI → 1.0)
/// - k = 1 : Flory distribution (PDI = 2.0)
/// - k < 1 : broader than Flory (PDI > 2.0)
pub struct SchulzZimm;

impl ChainLengthDistribution for SchulzZimm {
    fn sample(
        &self,
        mn: f64,
        pdi: f64,
        m0: f64,
        num_chains: usize,
        rng: &mut dyn RngCore,
    ) -> Vec<usize> {
        let xn = (mn / m0).max(1.0);
        let pdi_clamped = pdi.max(1.0 + 1e-9);
        let k = (1.0 / (pdi_clamped - 1.0)).min(1e6);
        let theta = xn / k;

        let gamma = Gamma::new(k, theta).expect("valid gamma parameters");
        (0..num_chains)
            .map(|_| {
                let x: f64 = gamma.sample(rng);
                (x.round() as usize).max(1)
            })
            .collect()
    }

    fn name(&self) -> &'static str {
        "schulz_zimm"
    }
}
