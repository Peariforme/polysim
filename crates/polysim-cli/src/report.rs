/// All data needed to render one analysis report.
pub struct AnalysisResult {
    pub bigsmiles_str: String,
    pub strategy_label: String,
    pub architecture_label: String,
    pub begin_block: Option<String>,
    pub end_block: Option<String>,
    pub smiles: String,
    pub repeat_count: usize,
    pub mn: f64,
    pub mono_mass: f64,
    /// Raw (ASCII) molecular formula, subscript conversion is done at render time.
    pub formula_raw: String,
    pub n_atoms: usize,
    /// Mn − target, present only when `--by-mn` was used.
    pub delta_mn: Option<f64>,
    /// monoisotopic mass − target, present only when `--by-mass` was used.
    pub delta_mass: Option<f64>,
}
