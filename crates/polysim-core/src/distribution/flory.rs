use rand::{Rng, RngCore};

use super::ChainLengthDistribution;

/// Flory (most probable) chain length distribution.
///
/// Characteristic of step-growth (polycondensation) polymers.
/// The number distribution is geometric: P(n) = (1−p)·p^(n−1),
/// where p = 1 − 1/Xn is the extent of reaction.
///
/// PDI is inherently ≈ 1+p (approaches 2.0 for high conversion);
/// the `pdi` parameter is ignored.
pub struct Flory;

impl ChainLengthDistribution for Flory {
    fn sample(
        &self,
        mn: f64,
        _pdi: f64,
        m0: f64,
        num_chains: usize,
        rng: &mut dyn RngCore,
    ) -> Vec<usize> {
        let xn = (mn / m0).max(1.0);
        let p = 1.0 - 1.0 / xn;

        (0..num_chains)
            .map(|_| {
                if p <= 0.0 {
                    return 1;
                }
                // Geometric sampling: n = ceil(ln(U) / ln(p))
                let u: f64 = rng.random_range(f64::EPSILON..1.0);
                let n = (u.ln() / p.ln()).ceil() as usize;
                n.max(1)
            })
            .collect()
    }

    fn name(&self) -> &'static str {
        "flory"
    }
}
