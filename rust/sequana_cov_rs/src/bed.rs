//! Fast readers for the 3/4-column coverage BED files read by `sequana_coverage`
//! (`chrom <TAB> position <TAB> coverage [<TAB> second coverage]`).
//!
//! Rows are counted like pandas does with `header=None`: blank lines are skipped.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Seek, SeekFrom};
use std::path::Path;

pub struct Scan {
    /// contig names in order of first appearance
    pub names: Vec<String>,
    /// row index (0-based, blank lines excluded) of the first and last row of each contig
    pub first: Vec<u64>,
    pub last: Vec<u64>,
    /// byte offset of the first row of each contig
    pub offsets: Vec<u64>,
    pub total: u64,
}

fn trim_eol(mut line: &[u8]) -> &[u8] {
    while let Some((&last, rest)) = line.split_last() {
        if last == b'\n' || last == b'\r' {
            line = rest;
        } else {
            break;
        }
    }
    line
}

/// Find where each contig starts and ends, without parsing the numbers.
pub fn scan_bed(path: &Path) -> io::Result<Scan> {
    let mut reader = BufReader::with_capacity(1 << 20, File::open(path)?);
    let mut buf = Vec::with_capacity(256);
    let mut scan = Scan { names: vec![], first: vec![], last: vec![], offsets: vec![], total: 0 };
    let mut index: HashMap<Vec<u8>, usize> = HashMap::new();
    let (mut offset, mut row) = (0u64, 0u64);
    let mut current: Option<(usize, Vec<u8>)> = None;

    loop {
        buf.clear();
        let n = reader.read_until(b'\n', &mut buf)?;
        if n == 0 {
            break;
        }
        let line_start = offset;
        offset += n as u64;
        let line = trim_eol(&buf);
        if line.is_empty() {
            continue;
        }
        let name = &line[..line.iter().position(|&b| b == b'\t').unwrap_or(line.len())];
        let id = match &current {
            Some((id, cur)) if cur.as_slice() == name => *id,
            _ => {
                let id = *index.entry(name.to_vec()).or_insert_with(|| {
                    scan.names.push(String::from_utf8_lossy(name).into_owned());
                    scan.first.push(row);
                    scan.last.push(row);
                    scan.offsets.push(line_start);
                    scan.names.len() - 1
                });
                current = Some((id, name.to_vec()));
                id
            }
        };
        scan.last[id] = row;
        row += 1;
    }
    scan.total = row;
    Ok(scan)
}

/// A column of numbers; `is_int` is true when every token was an integer
/// (pandas would then infer int64 instead of float64).
pub struct NumColumn {
    pub values: Vec<f64>,
    pub is_int: bool,
}

pub struct Chunk {
    pub pos: Vec<i64>,
    pub cov: NumColumn,
    pub extra: Option<NumColumn>,
    pub next_offset: u64,
}

fn parse_token(tok: &[u8], is_int: &mut bool) -> Result<f64, String> {
    let s = std::str::from_utf8(tok).map_err(|e| e.to_string())?;
    if *is_int && !tok.iter().enumerate().all(|(i, &b)| b.is_ascii_digit() || (i == 0 && b == b'-' && tok.len() > 1)) {
        *is_int = false;
    }
    s.parse::<f64>().map_err(|_| format!("cannot parse '{}' as a number", s))
}

/// Read up to `nrows` rows starting at byte `offset` (the start of a line).
/// Errors out on anything unusual (other column count, empty or non-numeric
/// field) so that the caller can fall back to pandas.
pub fn read_bed(path: &Path, offset: u64, nrows: usize) -> Result<Chunk, String> {
    let mut file = File::open(path).map_err(|e| e.to_string())?;
    file.seek(SeekFrom::Start(offset)).map_err(|e| e.to_string())?;
    let mut reader = BufReader::with_capacity(1 << 20, file);
    let mut buf = Vec::with_capacity(256);

    let mut pos = Vec::with_capacity(nrows.min(1 << 24));
    let mut cov = NumColumn { values: Vec::with_capacity(nrows.min(1 << 24)), is_int: true };
    let mut extra = NumColumn { values: Vec::new(), is_int: true };
    let mut ncols = 0;
    let mut consumed = 0u64;

    while pos.len() < nrows {
        buf.clear();
        let n = reader.read_until(b'\n', &mut buf).map_err(|e| e.to_string())?;
        if n == 0 {
            break;
        }
        consumed += n as u64;
        let line = trim_eol(&buf);
        if line.is_empty() {
            continue;
        }
        let mut fields = line.split(|&b| b == b'\t');
        fields.next(); // chromosome
        let (p, c) = match (fields.next(), fields.next()) {
            (Some(p), Some(c)) => (p, c),
            _ => return Err("lines must have at least 3 columns".into()),
        };
        let e = fields.next();
        if fields.next().is_some() {
            return Err("more than 4 columns".into());
        }
        let cols = if e.is_some() { 4 } else { 3 };
        if ncols == 0 {
            ncols = cols;
        } else if ncols != cols {
            return Err("inconsistent number of columns".into());
        }
        let mut pos_is_int = true;
        let pv = parse_token(p, &mut pos_is_int)?;
        if !pos_is_int {
            return Err("position is not an integer".into());
        }
        pos.push(pv as i64);
        cov.values.push(parse_token(c, &mut cov.is_int)?);
        if let Some(e) = e {
            extra.values.push(parse_token(e, &mut extra.is_int)?);
        }
    }
    Ok(Chunk {
        pos,
        cov,
        extra: if ncols == 4 { Some(extra) } else { None },
        next_offset: offset + consumed,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn tmp(name: &str, content: &str) -> std::path::PathBuf {
        let p = std::env::temp_dir().join(format!("sequana_cov_rs_{}_{}", std::process::id(), name));
        File::create(&p).unwrap().write_all(content.as_bytes()).unwrap();
        p
    }

    #[test]
    fn scan_finds_contigs() {
        let p = tmp("scan", "a\t1\t5\na\t2\t6\n\nb\t1\t7\nb\t2\t8\nb\t3\t9\n");
        let s = scan_bed(&p).unwrap();
        assert_eq!(s.names, ["a", "b"]);
        assert_eq!((s.first.clone(), s.last.clone()), (vec![0, 2], vec![1, 4]));
        assert_eq!(s.offsets, vec![0, 13]);
        assert_eq!(s.total, 5);
    }

    #[test]
    fn read_from_offset_and_detect_ints() {
        let p = tmp("read", "a\t1\t5\na\t2\t6\nb\t1\t7.5\nb\t2\t8\n");
        let c = read_bed(&p, 12, 10).unwrap();
        assert_eq!(c.pos, vec![1, 2]);
        assert_eq!(c.cov.values, vec![7.5, 8.0]);
        assert!(!c.cov.is_int);
        assert!(c.extra.is_none());
        let c = read_bed(&p, 0, 2).unwrap();
        assert!(c.cov.is_int);
        assert_eq!(c.next_offset, 12);
    }

    #[test]
    fn rejects_unusual_files() {
        let p = tmp("bad", "a\t1\tnan_value\n");
        assert!(read_bed(&p, 0, 1).is_err());
        let p = tmp("bad2", "a\t1\t2\t3\t4\n");
        assert!(read_bed(&p, 0, 1).is_err());
    }
}
