//! Fast CSV writer reproducing `DataFrame.to_csv(index=False, float_format="%.Ng")`.

use std::fmt::Write;

/// C/Python `%.{p}g` formatting, appended to `out`. `scratch` is reused between
/// calls so that nothing is allocated per value.
pub fn fmt_g_into(out: &mut String, scratch: &mut String, x: f64, p: usize) {
    if x.is_nan() {
        return; // pandas na_rep=""
    }
    if x.is_infinite() {
        out.push_str(if x > 0.0 { "inf" } else { "-inf" });
        return;
    }
    let p = p.max(1);
    if x == 0.0 {
        out.push_str(if x.is_sign_negative() { "-0" } else { "0" });
        return;
    }
    // one exact formatting gives the digits and the exponent after rounding
    scratch.clear();
    write!(scratch, "{:.*e}", p - 1, x).unwrap();
    let bytes = scratch.as_bytes();
    let (neg, rest) = if bytes[0] == b'-' { (true, &bytes[1..]) } else { (false, bytes) };
    let epos = rest.iter().position(|&b| b == b'e').unwrap();
    let exp: i32 = std::str::from_utf8(&rest[epos + 1..]).unwrap().parse().unwrap();
    // significant digits, without the decimal point
    let mut digits = [0u8; 40];
    let mut nd = 0;
    for &b in &rest[..epos] {
        if b != b'.' {
            digits[nd] = b;
            nd += 1;
        }
    }
    while nd > 1 && digits[nd - 1] == b'0' {
        nd -= 1;
    }
    let digits = &digits[..nd];

    if neg {
        out.push('-');
    }
    let push = |out: &mut String, d: &[u8]| out.push_str(std::str::from_utf8(d).unwrap());
    if exp < -4 || exp >= p as i32 {
        out.push(digits[0] as char);
        if nd > 1 {
            out.push('.');
            push(out, &digits[1..]);
        }
        out.push('e');
        out.push(if exp < 0 { '-' } else { '+' });
        let e = exp.unsigned_abs();
        if e < 10 {
            out.push('0');
        }
        write!(out, "{}", e).unwrap();
    } else if exp >= 0 {
        let int_len = exp as usize + 1;
        if nd <= int_len {
            push(out, digits);
            for _ in nd..int_len {
                out.push('0');
            }
        } else {
            push(out, &digits[..int_len]);
            out.push('.');
            push(out, &digits[int_len..]);
        }
    } else {
        out.push_str("0.");
        for _ in 0..(-exp - 1) {
            out.push('0');
        }
        push(out, digits);
    }
}

/// C/Python `%.{p}g` formatting (allocating convenience wrapper, used by tests).
#[cfg(test)]
pub fn fmt_g(x: f64, p: usize) -> String {
    let (mut out, mut scratch) = (String::new(), String::new());
    fmt_g_into(&mut out, &mut scratch, x, p);
    out
}

pub enum Column<'a> {
    Int(&'a [i64]),
    Float(&'a [f64]),
}

impl Column<'_> {
    fn len(&self) -> usize {
        match self {
            Column::Int(v) => v.len(),
            Column::Float(v) => v.len(),
        }
    }
}

/// Header + rows, `\n` terminated. Floats use `%.{precision}g`, ints are verbatim.
pub fn to_csv(names: &[String], cols: &[Column], precision: usize) -> String {
    assert_eq!(names.len(), cols.len());
    let n = cols.first().map_or(0, |c| c.len());
    assert!(cols.iter().all(|c| c.len() == n), "length mismatch");
    let mut out = String::with_capacity(n * cols.len() * 6 + 64);
    let mut scratch = String::with_capacity(48);
    out.push_str(&names.join(","));
    out.push('\n');
    for i in 0..n {
        for (j, c) in cols.iter().enumerate() {
            if j > 0 {
                out.push(',');
            }
            match c {
                Column::Int(v) => write!(out, "{}", v[i]).unwrap(),
                Column::Float(v) => {
                    let before = out.len();
                    fmt_g_into(&mut out, &mut scratch, v[i], precision);
                    // like the csv module, quote an empty field when it is the whole row
                    if cols.len() == 1 && out.len() == before {
                        out.push_str("\"\"");
                    }
                }
            }
        }
        out.push('\n');
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_c_printf() {
        assert_eq!(fmt_g(0.0, 3), "0");
        assert_eq!(fmt_g(1234.0, 3), "1.23e+03");
        assert_eq!(fmt_g(999.5, 3), "1e+03");
        assert_eq!(fmt_g(123.0, 3), "123");
        assert_eq!(fmt_g(12.5, 3), "12.5");
        assert_eq!(fmt_g(0.0001234, 3), "0.000123");
        assert_eq!(fmt_g(0.00001234, 3), "1.23e-05");
        assert_eq!(fmt_g(-2.0, 3), "-2");
        assert_eq!(fmt_g(f64::NAN, 3), "");
    }

    #[test]
    fn csv_layout() {
        let s = to_csv(
            &["pos".into(), "cov".into()],
            &[Column::Int(&[1, 2]), Column::Float(&[10.0, 31.456])],
            3,
        );
        assert_eq!(s, "pos,cov\n1,10\n2,31.5\n");
    }
}
