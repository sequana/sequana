def CpG(sequence, window=200):
    """
    The Sequence Manipulation Suite: CpG Islands
    Results for 1200 residue sequence "sample sequence" starting "taacatactt".
    CpG islands search using window size of 200.
    Range, value
    32 to 231, the y-value is 1.75 and the %GC content is 50.5
    33 to 232, the y-value is 1.75 and the %GC content is 50.5


    Gardiner-Garden M, Frommer M. J Mol Biol. 1987 Jul 20;196(2):261-82.

    """
    return compute_cpg_content(sequence)


def compute_cpg_content(seq):
    seq = seq.upper()
    c_count = seq.count("C")
    g_count = seq.count("G")
    cg_count = seq.count("CG")
    expected_cg = (c_count * g_count) / len(seq) if len(seq) > 0 else 0
    obs_exp_ratio = cg_count / expected_cg if expected_cg > 0 else 0
    gc_content = (c_count + g_count) / len(seq)
    return {"CpG count": cg_count, "Observed/Expected CpG": obs_exp_ratio, "GC content": gc_content}
