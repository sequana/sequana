from sequana.mpileup import MPileup
from pytest import approx


from . import test_dir



def test_mpileup():
    m = MPileup(f"{test_dir}/data/mpileup/mpileup.txt", f"{test_dir}/data/mpileup/ref.fa")
    df = m.get_mutation()
    assert df['A'].sum() == approx(45.51, 0.01)
    assert df['C'].sum() == approx(16.84, 0.01)
    assert df['G'].sum() == approx(10.10, 0.01)
    assert df['T'].sum() == approx(28.85, 0.01)


    # try with a string as the ref
    ref = "ATGACGGAATATAAGCTGGTGGTGGTGGGCGCCGGCGGTGTGGGCAAGAGTGCGCTGACCATCCAGCTGATCCAGAACCATTTTGTGGACGAATACGACCC"
    m = MPileup(f"{test_dir}/data/mpileup/mpileup.txt", ref)
    df = m.get_mutation()


    # try with a wrong filename
    ref = "ATGACGGAATATAAGCTGGTGGTGGTGGGCGCCGGCGGTGTGGGCAAGAGTGCGCTGACCATCCAGCTGATCCAGAACCATTTTGTGGACGAATACGACCC"
    try:
        m = MPileup(f"wrongmpileup.txt", ref)
        assert False
    except IOError:
        assert True

    m.plot_mutation_matrix()
    m.plot_total_errors(include_deletions=True, include_Ns=True)
    m.plot_stack_bars(include_deletions=True)
