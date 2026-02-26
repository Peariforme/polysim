use bigsmiles::{BigSmiles, BigSmilesSegment};

/// Returns the concatenated SMILES fragments that appear **before** the first
/// stochastic object in the BigSMILES string (initiator / beginning block).
///
/// Returns `None` when no such fragment exists (plain `{…}` without a prefix).
pub fn before_stochastic(bs: &BigSmiles) -> Option<String> {
    let first_stoch = bs
        .segments
        .iter()
        .position(|s| matches!(s, BigSmilesSegment::Stochastic(_)))?;

    let s: String = bs.segments[..first_stoch]
        .iter()
        .filter_map(|seg| match seg {
            BigSmilesSegment::Smiles(mol) => Some(format!("{mol}")),
            BigSmilesSegment::Stochastic(_) => None,
        })
        .collect();

    (!s.is_empty()).then_some(s)
}

/// Returns the concatenated SMILES fragments that appear **after** the last
/// stochastic object in the BigSMILES string (terminator / ending block).
///
/// Returns `None` when no such fragment exists.
pub fn after_stochastic(bs: &BigSmiles) -> Option<String> {
    let last_stoch = bs
        .segments
        .iter()
        .rposition(|s| matches!(s, BigSmilesSegment::Stochastic(_)))?;

    let s: String = bs.segments[last_stoch + 1..]
        .iter()
        .filter_map(|seg| match seg {
            BigSmilesSegment::Smiles(mol) => Some(format!("{mol}")),
            BigSmilesSegment::Stochastic(_) => None,
        })
        .collect();

    (!s.is_empty()).then_some(s)
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use bigsmiles::parse;

    #[test]
    fn no_blocks_for_plain_bigsmiles() {
        let bs = parse("{[]CC[]}").unwrap();
        assert!(before_stochastic(&bs).is_none());
        assert!(after_stochastic(&bs).is_none());
    }

    #[test]
    fn detects_begin_block() {
        let bs = parse("CC{[$]CC[$]}").unwrap();
        assert!(before_stochastic(&bs).is_some());
        assert!(after_stochastic(&bs).is_none());
    }

    #[test]
    fn detects_end_block() {
        let bs = parse("{[$]CC[$]}CC").unwrap();
        assert!(before_stochastic(&bs).is_none());
        assert!(after_stochastic(&bs).is_some());
    }

    #[test]
    fn detects_both_blocks() {
        let bs = parse("CC{[$]CC[$]}CC").unwrap();
        assert_eq!(before_stochastic(&bs).as_deref(), Some("CC"));
        assert_eq!(after_stochastic(&bs).as_deref(), Some("CC"));
    }
}
