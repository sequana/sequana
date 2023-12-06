from pylab import normal

from sequana.viz import ANOVA


def test_anova():
    import pandas as pd

    A = normal(0.5, size=10000)
    B = normal(0.25, size=10000)
    C = normal(0, 0.5, size=10000)
    df = pd.DataFrame({"A": A, "B": B, "C": C})
    a = ANOVA(df)
    print(a.anova())
    a.imshow_anova_pairs()
