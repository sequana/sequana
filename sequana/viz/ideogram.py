import matplotlib.patches as patches
import matplotlib.pyplot as plt


class Ideogram:
    """


    from sequana import FastA
    f = FastA("assembly.fa")
    N = len(f)
    L = [ f.get_lengths_as_dict()[str(x)] for x in range(1,N+1)]

    #
    import pandas as pd
    centromeres = pd.read_csv("centromeres.csv", sep=",")
    C = [centromeres.query('chrom==@x')['hic_position'].values[0] for x in range(1,N+1)]

    telomeres = pd.read_csv("telomeres//sequana.telomark.telomark.csv")
    X1 = telomeres['5to3_LHS_position'].values
    X2 = telomeres['3to5_RHS_position'].values
    X2 = L - X2
    """

    def __init__(self, df, gaps={}):

        self.df = df.copy()
        self.gaps = gaps

        self.chrom_width = 0.6
        self.spacing = 1.0
        self.tri_height = 10000  # how far triangles extend

    def plot(self, figsize=(10, 12)):

        N = len(self.df)

        fig, ax = plt.subplots(figsize=figsize)

        C = self.df["centromere"].values
        L = self.df["length"].values
        X1 = self.df["LHS_telomere"].values
        X2 = self.df["RHS_telomere"].values

        for i in range(N):
            xpos = i * self.spacing
            length = float(L[i])
            centromere = float(C[i])
            x1 = float(X1[i])
            x2 = float(X2[i])

            # Bottom chromosome segment (up to centromere)
            chrom_bottom = patches.Rectangle(
                (xpos, 0), self.chrom_width, centromere - self.tri_height, edgecolor="black", facecolor="lightgrey"
            )
            ax.add_patch(chrom_bottom)

            # Top chromosome segment (after centromere)
            chrom_top = patches.Rectangle(
                (xpos, centromere + self.tri_height),
                self.chrom_width,
                length - centromere - self.tri_height,
                edgecolor="black",
                facecolor="lightgrey",
            )
            ax.add_patch(chrom_top)

            # Bottom (inverted) triangle - centromere
            triangle_bottom = patches.Polygon(
                [
                    (xpos, centromere - self.tri_height),
                    (xpos + self.chrom_width / 2, centromere),
                    (xpos + self.chrom_width, centromere - self.tri_height),
                ],
                closed=True,
                color="red",
            )
            ax.add_patch(triangle_bottom)

            # Top triangle - centromere
            triangle_top = patches.Polygon(
                [
                    (xpos, centromere + self.tri_height),
                    (xpos + self.chrom_width / 2, centromere),
                    (xpos + self.chrom_width, centromere + self.tri_height),
                ],
                closed=True,
                color="red",
            )
            ax.add_patch(triangle_top)

            # Telomeres
            tel_left = patches.Rectangle((xpos, 0), self.chrom_width, x1, facecolor="green")
            ax.add_patch(tel_left)

            tel_right = patches.Rectangle((xpos, x2), self.chrom_width, length - x2, facecolor="green")
            ax.add_patch(tel_right)

            if self.gaps and str(i + 1) in self.gaps.keys():
                for position in self.gaps[str(i + 1)][1:]:
                    gap = patches.Rectangle((xpos, position), self.chrom_width, 20000, facecolor="blue", zorder=20)
                    ax.add_patch(gap)

        # Aesthetics
        ax.set_xlim(-0.5, self.spacing * N)
        ax.set_ylim(0, max(L) + 10)
        ax.set_xticks([i * self.spacing + self.chrom_width / 2 for i in range(N)])
        ax.set_xticklabels([str(i + 1) for i in range(N)], rotation="vertical", fontsize=16)
        ax.set_ylabel("Position (bp)", fontsize=16)
        ax.set_xlabel("Chromosome", fontsize=16)
        ax.set_title("Chromosome Ideogram")
