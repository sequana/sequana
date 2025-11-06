import re


def yield_wig_by_chromosome(filename):
    """

    Given a wig file such as::

        variableStep chrom=1 step=1
        9315 723.0
        9316 723.0

    you can read it using::

        for chrom, entries in yield_wig_by_chromosome("example.wig"):
            print(f"Chromosome: {chrom}, {len(entries)} entries")
            print(entries[:3])

    """
    chrom = None
    current_data = []

    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):  # pragma: no cover
                continue

            if line.startswith("fixedStep"):
                # Yield previous chromosome block if present
                if chrom and current_data:
                    yield chrom, current_data
                    current_data = []

                is_fixed_step = True
                m = re.search(r"chrom=(\S+)", line)
                chrom = m.group(1)

                m = re.search(r"start=(\d+)", line)
                pos = int(m.group(1))

                m = re.search(r"step=(\d+)", line)
                step = int(m.group(1)) if m else 1

            elif line.startswith("variableStep"):
                # If we already have data for a chromosome, yield it
                if chrom and current_data:
                    yield chrom, current_data
                    current_data = []

                is_fixed_step = False
                m = re.search(r"chrom=(\S+)", line)
                chrom = m.group(1)

                pos = None  # variableStep lines give positions per-line

            else:
                if is_fixed_step:
                    try:
                        value = float(line)
                        current_data.append((pos, value))
                        pos += step
                    except ValueError:  # pragma: no cover
                        raise ValueError(f"Invalid fixedStep line: {line}")
                else:
                    try:
                        p, v = line.split()
                        current_data.append((int(p), float(v)))
                    except ValueError:  # pragma: no cover
                        raise ValueError(f"Invalid variableStep line: {line}")

    # Yield the last chromosome
    if chrom and current_data:
        yield chrom, current_data
