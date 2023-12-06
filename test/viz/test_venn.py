from sequana.viz.venn import plot_venn


def test3():
    A = set([1, 2, 3, 4, 5, 6, 7, 8, 9])
    B = set([6, 7, 8, 9, 10, 11, 12, 13])
    C = set([4, 5, 6, 7, 8, 9])
    plot_venn((A, B, C), labels=("A", "B", "C"))
    plot_venn((A, B, C), labels=("A", "B", "C"), weighted=True)


def test2():
    A = set([1, 2, 3, 4, 5, 6, 7, 8, 9])
    B = set([6, 7, 8, 9, 10, 11, 12, 13])
    plot_venn((A, B), labels=("A", "B"))
    plot_venn((A, B), labels=("A", "B"), weighted=True)
