import argparse
import gzip
import io
import re
import sys


def open_maybe_gz(path):
    if path == "-":
        return sys.stdin.buffer
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_embl_records(text):
    # Split records by line that is exactly '//'
    parts = []
    cur = []
    for line in text.splitlines():
        if line.strip() == "//":
            parts.append("\n".join(cur))
            cur = []
        else:
            cur.append(line)
    if cur:
        parts.append("\n".join(cur))
    return parts


def extract_record_info(record_text):
    # Find ID or AC
    id_ = None
    desc_lines = []
    seq_lines = []
    in_seq = False
    for line in record_text.splitlines():
        if line.startswith("ID") and not id_:
            # ID   SCU49845; SV 1; linear; mRNA; STD; PLN; 1859 BP.
            fields = line[2:].strip().split()
            if fields:
                id_ = fields[0].rstrip(";")
        elif line.startswith("AC") and not id_:
            # AC   U49845; U49846;
            ac = line[2:].strip().split()[0].rstrip(";")
            id_ = ac
        elif line.startswith("DE"):
            desc_lines.append(line[2:].strip())
        elif line.startswith("SQ"):
            in_seq = True
        elif in_seq:
            # sequence lines until end
            seq_lines.append(line)
    desc = " ".join([d for d in desc_lines if d]) if desc_lines else ""
    # Clean sequence: remove digits and spaces, keep letters
    seq = "".join(seq_lines)
    # Remove non-letters (digits, spaces, punctuation)
    seq = re.sub(r"[^A-Za-z]", "", seq)
    seq = seq.upper()
    return id_ or "unknown_id", desc, seq


def wrap_sequence(seq, width=60):
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def embl_to_fasta(inpath, outpath):
    # Read entire input (EMBL files are usually manageable); if very large,
    # streaming parsing would be better. This simple script uses read().
    if inpath == "-":
        raw = sys.stdin.read()
    else:
        opener = open_maybe_gz(inpath)
        # opener may be a binary buffer for stdin
        if hasattr(opener, "read") and "b" in getattr(opener, "mode", "t"):
            # binary
            raw = opener.read().decode("utf-8")
        else:
            raw = opener.read()
        opener.close()

    records = parse_embl_records(raw)

    if outpath == "-":
        out = sys.stdout
    else:
        out = open(outpath, "w")

    written = 0
    for rec in records:
        if not rec.strip():
            continue
        identifier, desc, seq = extract_record_info(rec)
        if not seq:
            # skip empty sequences but warn
            print(f"Warning: record {identifier} has no sequence, skipping", file=sys.stderr)
            continue
        header = identifier
        if desc:
            header += " " + desc
        out.write(f">{header}\n")
        out.write(wrap_sequence(seq) + "\n")
        written += 1

    if out is not sys.stdout:
        out.close()
    print(f"Done. Wrote {written} record(s) to {outpath}", file=sys.stderr)


def open_maybe_gz(path):
    if path == "-":
        # stdin is already text
        return sys.stdin
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")
