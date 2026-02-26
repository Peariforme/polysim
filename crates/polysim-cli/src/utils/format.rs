use comfy_table::Color as TableColor;

/// Returns the sign string and the colour for a delta value.
///
/// The colour is **green** when |Δ / reference| < 0.5 %, **yellow** otherwise.
pub fn delta_style(delta: f64, reference: f64) -> (&'static str, TableColor) {
    let sign = if delta >= 0.0 { "+" } else { "" };
    let relative = if reference.abs() > 1e-9 {
        (delta / reference).abs()
    } else {
        delta.abs()
    };
    let color = if relative < 0.005 {
        TableColor::Green
    } else {
        TableColor::Yellow
    };
    (sign, color)
}

/// Replaces ASCII digits with their Unicode subscript equivalents.
///
/// Example: `"C20H42"` → `"C₂₀H₄₂"`.
pub fn subscript_digits(s: &str) -> String {
    s.chars()
        .map(|c| match c {
            '0' => '₀',
            '1' => '₁',
            '2' => '₂',
            '3' => '₃',
            '4' => '₄',
            '5' => '₅',
            '6' => '₆',
            '7' => '₇',
            '8' => '₈',
            '9' => '₉',
            _ => c,
        })
        .collect()
}

/// Truncates a long string with a mid-string ellipsis `…`.
///
/// If `s.len() <= max_len` the original string is returned unchanged.
pub fn truncate(s: &str, max_len: usize) -> String {
    if s.len() <= max_len {
        s.to_string()
    } else {
        let half = (max_len.saturating_sub(1)) / 2;
        format!("{}…{}", &s[..half], &s[s.len() - half..])
    }
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // subscript_digits --------------------------------------------------------

    #[test]
    fn subscript_digits_converts_all_ten() {
        assert_eq!(subscript_digits("0123456789"), "₀₁₂₃₄₅₆₇₈₉");
    }

    #[test]
    fn subscript_digits_typical_formula() {
        assert_eq!(subscript_digits("C20H42"), "C₂₀H₄₂");
    }

    #[test]
    fn subscript_digits_formula_with_heteroatoms() {
        assert_eq!(subscript_digits("C8H8O2"), "C₈H₈O₂");
    }

    #[test]
    fn subscript_digits_no_digits_unchanged() {
        assert_eq!(subscript_digits("CHONSFClBrI"), "CHONSFClBrI");
    }

    #[test]
    fn subscript_digits_empty_string() {
        assert_eq!(subscript_digits(""), "");
    }

    // truncate ----------------------------------------------------------------

    #[test]
    fn truncate_short_string_returned_unchanged() {
        assert_eq!(truncate("CCCC", 10), "CCCC");
    }

    #[test]
    fn truncate_exact_max_len_returned_unchanged() {
        let s = "A".repeat(20);
        assert_eq!(truncate(&s, 20), s);
    }

    #[test]
    fn truncate_long_string_contains_ellipsis() {
        let long = "A".repeat(100);
        let result = truncate(&long, 20);
        assert!(
            result.contains('…'),
            "truncated string must contain '…', got: '{result}'"
        );
    }

    #[test]
    fn truncate_long_string_is_shorter_than_original() {
        let long = "A".repeat(100);
        let result = truncate(&long, 20);
        assert!(result.len() < 100, "result len={}", result.len());
    }

    #[test]
    fn truncate_preserves_start_and_end() {
        let s = format!("{}{}", "A".repeat(50), "B".repeat(50));
        let result = truncate(&s, 20);
        assert!(result.starts_with('A'), "start should be preserved");
        assert!(result.ends_with('B'), "end should be preserved");
    }

    // delta_style -------------------------------------------------------------

    #[test]
    fn delta_style_small_relative_error_is_green() {
        let (_, color) = delta_style(1.0, 1000.0); // 0.1% < 0.5%
        assert_eq!(color, TableColor::Green);
    }

    #[test]
    fn delta_style_large_relative_error_is_yellow() {
        let (_, color) = delta_style(10.0, 100.0); // 10% > 0.5%
        assert_eq!(color, TableColor::Yellow);
    }

    #[test]
    fn delta_style_positive_gets_plus_sign() {
        let (sign, _) = delta_style(5.0, 100.0);
        assert_eq!(sign, "+");
    }

    #[test]
    fn delta_style_negative_gets_no_extra_sign() {
        let (sign, _) = delta_style(-5.0, 100.0);
        assert_eq!(sign, "");
    }
}
