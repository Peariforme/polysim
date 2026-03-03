use std::fmt;

use crate::polymer::PolymerEnsemble;

/// Summary statistics for a polymer ensemble.
#[derive(Debug, Clone)]
pub struct EnsembleStats {
    pub num_chains: usize,
    pub mn: f64,
    pub mw: f64,
    pub pdi: f64,
    pub mn_min: f64,
    pub mn_max: f64,
    pub mn_median: f64,
    pub mn_std_dev: f64,
}

impl EnsembleStats {
    /// Computes statistics from an ensemble.
    pub fn from_ensemble(ensemble: &PolymerEnsemble) -> Self {
        let mut masses: Vec<f64> = ensemble.chains().iter().map(|c| c.mn).collect();
        masses.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let n = masses.len();
        let mn = ensemble.mn();
        let mw = ensemble.mw();

        let mn_median = if n % 2 == 1 {
            masses[n / 2]
        } else {
            (masses[n / 2 - 1] + masses[n / 2]) / 2.0
        };

        let variance = masses.iter().map(|m| (m - mn).powi(2)).sum::<f64>() / n as f64;
        let mn_std_dev = variance.sqrt();

        Self {
            num_chains: n,
            mn,
            mw,
            pdi: ensemble.pdi(),
            mn_min: masses[0],
            mn_max: masses[n - 1],
            mn_median,
            mn_std_dev,
        }
    }
}

impl fmt::Display for EnsembleStats {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Ensemble Statistics")?;
        writeln!(f, "  Chains:    {}", self.num_chains)?;
        writeln!(f, "  Mn:        {:.1} g/mol", self.mn)?;
        writeln!(f, "  Mw:        {:.1} g/mol", self.mw)?;
        writeln!(f, "  PDI:       {:.3}", self.pdi)?;
        writeln!(f, "  Median:    {:.1} g/mol", self.mn_median)?;
        writeln!(f, "  Std dev:   {:.1} g/mol", self.mn_std_dev)?;
        write!(
            f,
            "  Mn range:  {:.1} - {:.1} g/mol",
            self.mn_min, self.mn_max
        )
    }
}
