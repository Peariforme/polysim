use rand::RngCore;
use rand_distr::{Distribution, LogNormal as LN};

use super::ChainLengthDistribution;

/// Log-normal chain length distribution.
///
/// Flexible distribution parameterized by Mn and PDI:
///   σ² = ln(PDI),  μ = ln(Xn) − σ²/2
///
/// Note: the relation PDI = exp(σ²) is exact for a continuous mass
/// distribution. When end-group mass (m_end) is non-negligible compared
/// to Mn, the actual PDI of the generated ensemble will deviate slightly
/// because MW(n) = n·m0 + m_end introduces a constant offset. The
/// approximation is excellent when m_end << Mn (typical for high-MW polymers).
///
/// Suitable for most synthetic polymers (radical, anionic, etc.).
pub struct LogNormal;

impl ChainLengthDistribution for LogNormal {
    fn sample(
        &self,
        mn: f64,
        pdi: f64,
        m0: f64,
        num_chains: usize,
        rng: &mut dyn RngCore,
    ) -> Vec<usize> {
        let sigma_sq = pdi.max(1.0).ln();
        let sigma = sigma_sq.sqrt();
        let xn = (mn / m0).max(1.0);
        let mu = xn.ln() - sigma_sq / 2.0;

        let dist = LN::new(mu, sigma).expect("valid log-normal parameters");
        (0..num_chains)
            .map(|_| {
                let x: f64 = dist.sample(rng);
                (x.round() as usize).max(1)
            })
            .collect()
    }

    fn name(&self) -> &'static str {
        "log_normal"
    }
}
