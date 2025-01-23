#  This file is part of Sequana software
#
#  Copyright (c) 2016-2020 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
import gzip
import io
import os
import subprocess
import sys
import tempfile
from pathlib import Path

from tqdm import tqdm

from sequana import logger
from sequana.lazy import numpy as np
from sequana.lazy import pandas as pd
from sequana.lazy import pylab, pysam

__all__ = ["SomyScore"]


class ClusterModels:
    """Cluster GMM to simplify models


    Sometimes, a model may be split in 2 if we provide a large K.
    Instead of using AIC criteria, we can cluster GMM that mean is close
    regarding their variability. This could be useful in some cases
    when we know that models should be clearly separater e.g. in somy studies.


    """

    def __init__(self, mus, sigmas, pis):
        self.mus = mus
        self.sigmas = sigmas
        self.pis = pis

    def __len__(self):
        return len(self.mus)

    def __getitem__(self, index):
        return self.mus[index], self.sigmas[index], self.pis[index]

    def cluster(self):
        # find the most significant
        index = np.argmax(self.pis)
        to_cluster = [index]
        mu, sigma, pi = self[index]
        for i in range(len(self)):
            if i == index:
                pass
            else:
                mu2, sigma2, pi2 = self[i]
                d = abs(mu - mu2)
                if d < 1.2 * sigma + 1.2 * sigma2:
                    to_cluster.append(i)
        print(self.mus)
        print(to_cluster)
        newmu = [self.mus[i] * self.pis[i] for i in to_cluster]
        pis = [self.pis[i] for i in to_cluster]
        newmu = sum(newmu) / sum(pis)

        print(newmu)
        return newmu

    def get_std(self):
        return sum([pi * sigma for pi, sigma in zip(self.pis, self.sigmas)]) / sum(self.pis)


class SomyScore:
    """

    from sequana.somy import SomyScore
    ss = SomyScore("file.bam")
    ss.run()
    # filters
    ss.remove_outliers()
    ss.remove_flanking_regions()
    ss.boxplot()

    # store df for later
    df = ss.df.copy()

    ss2.run(mapq=35)
    # filters
    ss.remove_outliers()
    ss.remove_flanking_regions()
    ss.boxplot()

    ss.df = pd.merge(ss.df, df)
    ss.boxplot()



    """

    def __init__(self, filename=None, window_size=1000):
        self.filename = filename
        self.window_size = window_size

        if self.filename:
            align = pysam.AlignmentFile(self.filename)
            if not align.has_index():
                logger.error(
                    f"Your input file must have an index. Please build one using 'samtools index {filename}'. You can install samtools using 'damona install samtools'"
                )
                sys.exit(1)

        self._df = None
        self.save_em = True
        self.save_em_filename = "sequana_somy_em.png"
        self.info = {}

    def _get_df(self):
        return self._df

    df = property(_get_df)

    def _get_coverage(
        self, use_median=True, chrom=None, flag=1796, threads=4, fast=True, window=1000, mapq=40, tag="default"
    ):

        # TODO: include -F 256 to remove secondary ?
        options = "-x" if fast else ""
        options += f" -t {threads} "
        options += f" --flag {flag}"
        options += " -m " if use_median else ""
        options += f" -c {chrom}" if chrom else ""

        # -n do not output per base depth. faster
        options += " -n "

        with tempfile.TemporaryDirectory() as TEMPDIR:
            cmd = f"mosdepth {options} -b {window}  -Q {mapq} {TEMPDIR}/mos {self.filename}"
            logger.debug(cmd)
            status = subprocess.run(cmd.split())
            assert status.returncode == 0

            with gzip.GzipFile(f"{TEMPDIR}/mos.regions.bed.gz", "r") as zip_ref:
                data = zip_ref.read()
            df = pd.read_csv(io.StringIO(data.decode()), sep="\t", header=None)
            df.columns = ["chr", "start", "stop", "depth"]
            # no need to keep both start/stop
            del df["stop"]

            # add a convenient distance to the border for future filtering
            N = len(df)
            X = list(range(N))
            df["dist"] = [min(x, abs(x - N + 1)) for x in range(len(X))]
            df["tag"] = tag
            return df

    def compute_coverage(
        self,
        use_median=True,
        fast=True,
        mapq=0,
        flag=1796,
        chromosomes=None,
        tag="default",
        threads=4,
        exclude_chromosomes=[],
    ):
        # 3844 also exclude supplementary
        align = pysam.AlignmentFile(self.filename)
        contig_names = [x.contig for x in align.get_index_statistics()]

        dfs = []
        if chromosomes is None:
            chromosomes = contig_names
        else:
            for chrom in chromosomes:
                if chrom not in contig_names:
                    logger.error(f"chromosome/contig {chrom} not found in the list: {contig_names}")
                    sys.exit()

        import multiprocessing

        arguments = []
        for chrom in chromosomes:
            if chrom not in exclude_chromosomes:
                arguments.append(
                    {
                        "chrom": chrom,
                        "use_median": use_median,
                        "fast": fast,
                        "flag": flag,
                        "mapq": mapq,
                        "window": self.window_size,
                    }
                )

        with multiprocessing.Pool(processes=threads) as pool:

            with tqdm(total=len(arguments)) as pbar:

                def update(*_):
                    pbar.update()

                results = []
                for args in arguments:
                    result = pool.apply_async(self._get_coverage, kwds=args, callback=update)
                    results.append(result)

                # collect the results
                results = [result.get() for result in results]

        # build df, or accumulate them with previous runs.
        if self._df is not None:
            self._df = pd.concat([self._df] + results)
        else:
            self._df = pd.concat(results)

    def _estimate_diploid_mean_depth_em(self, k=None, bins=100):

        # we can estimate the mean by brute force
        mu = self.df["depth"].median()
        # however, this may be biased if lots of different somy. We can also try an
        # estimate based on mixture model assuming that diploid is the main population
        from sequana.mixture import EM

        em = EM(self.df["depth"])
        sigma = mu / 3.0

        best_k = 2
        best_LL = 0
        results = {}
        if k is None:
            kmin = 2
            kmax = 6
        else:
            kmin = k
            kmax = k + 1

        for ktry in range(kmin, kmax):
            em.estimate(k=ktry)
            LL = em.results["log_likelihood"]
            results[ktry] = em.results
            if LL > best_LL:
                best_LL = em.results["log_likelihood"]
                best_k = ktry
            print(ktry, LL)

        # identifies main peak (diploid)
        import numpy as np

        best_results = results[best_k]

        print(best_k)
        print(best_results)

        estimate_em_mu = best_results["mus"][np.argmax(best_results["pis"])]

        # need to check somy here. Are they close to 1, 1.5, 2, 3
        # or 0.5,0.75,1,1.5 ? in the second case, it means the estimated mu was not correct
        if self.save_em:
            pylab.clf()
            em.estimate(k=best_k)
            em.plot(bins=bins)
            pylab.savefig(self.save_em_filename)
        self.info["mu"] = self.df["depth"].median()
        self.info["median"] = self.df["depth"].median()
        self.info["estimate_mu_diploid"] = estimate_em_mu
        self.info["em_results"] = best_results

        # Given the mus, we need to figure out what is the best mus for diploid.
        # Taken the minimum is not robust enough since, you may have distribution below.
        # ClusterModels, can help merging.clustering closeby values. the diploid case is suppose
        # to be the most important one but if we test only a subset of chromosomes, this may not work
        # or in extreme case, triploi+tetra are as important.

        # For now, find the main gaussian and consider this is the chromosomes with somy of 2
        best_results["pis"]
        index = np.argmax(best_results["pis"])
        muhat = best_results["mus"][index]
        # if k>=2:
        #    cm = ClusterModels(best_results['mus'], best_results['sigmas'], best_results['pis'])
        #    muhat = cm.cluster()
        # else:
        #    muhat = estimate_em_mu

        # k may be too large and main diploid model split in 2 or 3
        # difficuly and not robust to use AIC
        # we can assume that distance between k is large and so merge models that are closeby
        print(mu, estimate_em_mu, muhat)
        return muhat

    def remove_outliers(self, percentile=0.05):
        def remove_outliers(
            group, col_name="depth", lower_percentile=0.05, upper_percentile=0.95, include_groups=False
        ):
            lower = group[col_name].quantile(lower_percentile)
            upper = group[col_name].quantile(upper_percentile)
            filtered_group = group[(group[col_name] >= lower) & (group[col_name] <= upper)]
            return filtered_group

        # remove outliers
        import pandas as pd

        df = self.df.groupby("chr").apply(remove_outliers, include_groups=True)
        self._df = df.reset_index(drop=True)

    def remove_flanks(self, remove_flanking_regions_kb=10):
        # filter out  data to remove telomeric regions based on distance (number of windows)
        # e.g. 10 kb means 10 regions of 1kb to remove on both sides
        from math import ceil

        N = ceil(remove_flanking_regions_kb * 1000.0 / self.window_size)
        self._df = self._df.query("dist>@N")

    def remove_low_depth(self, threshold=10):
        self._df = self.df.query("depth>@threshold")

    def _estimate_diploid_median_depth(self):
        # iterative estimate
        # first, we take everything
        muhat = self.df["depth"].median()
        # this may be biased by tetraploid chromosomes
        return muhat

    def boxplot(
        self,
        legend=False,
        hlines=[1, 1.5, 2, 2.5],
        normalise=True,
        k=None,
        outfile="sequana_somy.png",
        muhat=None,
        hybrid=False,
        method="em",
        palette=None,
    ):
        """

        :param bool legend: set to 'auto' to add a legend
        """
        import seaborn as sns

        data = self._df.copy()

        if muhat:
            if normalise:
                data["depth"] /= muhat
                muhat = 1
                data["depth"] *= 2
            else:
                muhat = known_diploid_coverage
        elif normalise is True:
            if method == "em":
                try:
                    muhat = self._estimate_diploid_mean_depth_em(k=k)
                    print(f"Estimate mu: {muhat}")
                except Exception as err:
                    print(err)
                    logger.warning("Could not estimate mixture model. Rolling back to simple median estimator")
                    muhat = self.df["depth"].median()
                data["depth"] /= muhat
            elif method == "median":
                muhat = self.df["depth"].median()
                print(f"Estimate mu: {muhat}")
                data["depth"] /= muhat
            elif method == "mean":
                muhat = self.df["depth"].mean()
                print(f"Estimate mu: {muhat}")
                data["depth"] /= muhat
            muhat = 1
            data["depth"] *= 2  # (diploid)
        elif method == "bootstrap":
            pass
        else:
            muhat = self._estimate_diploid_median_depth()
            print(f"Estimate mu: {muhat}")

        if hybrid:
            # data['depth'] /=2
            pass

        # self.df.groupby(["chr", "tag"])['depth'].to_csv("somies.csv")

        estimated_somies = []
        measured_somies = []
        chromosomes = []
        tags = []

        for datum in data.groupby(["chr", "tag"])["depth"]:
            tag = datum[0][1]
            x = datum[1].median()
            if abs(1 - x) < abs(2 - x):
                somy = 1
            elif abs(2 - x) < abs(3 - x):
                somy = 2
            elif abs(x - 3) < abs(4 - x):
                somy = 3
            elif abs(x - 4) < abs(x - 5):
                somy = 4
            else:
                somy = 5
            estimated_somies.append(somy)
            measured_somies.append(round(x, 3))
            chromosomes.append(datum[0][0])
            tags.append(tag)
        somies = pd.DataFrame(
            {
                "chrom": chromosomes,
                "tag": tags,
                "estimated_somies": estimated_somies,
                "measured_somies": measured_somies,
            }
        )
        somies.to_csv("somies.csv", index=None)
        # plot somies
        import seaborn as sns

        pylab.clf()
        ax = sns.regplot(
            x="estimated_somies", y="measured_somies", data=somies, x_jitter=0.0, scatter_kws={"alpha": 0.3}
        )
        ax.set(xlabel="Estimated somies", ylabel="Measured somies")

        m, M = somies["estimated_somies"].min(), somies["estimated_somies"].max()
        pylab.plot([m, M], [m, M], color="r", ls="--", lw=1)
        pylab.savefig("sequana_somies_meas_vs_estim.png")
        self.somies = somies
        error = ((somies["estimated_somies"] - somies["measured_somies"]) ** 2).sum() / len(somies)
        print(f"Error: {error}")

        self.info["error"] = error

        pylab.clf()
        sns.boxplot(data=data, x="chr", y="depth", showfliers=False, hue="tag", legend=legend, palette=palette)
        pylab.xticks(rotation=90)
        for h in hlines:
            if h == 1:
                pylab.axhline(2 * muhat * h, color="r", lw=1, ls="--", label=f"diploid")
            elif h == 1.5:
                pylab.axhline(2 * muhat * h, color="r", lw=1, ls="--", label=f"triploid")
            elif h == 2:
                pylab.axhline(2 * muhat * h, color="r", lw=1, ls="--", label=f"tetraploid")
            elif h == 2.5:
                pylab.axhline(2 * muhat * h, color="r", lw=1, ls="--", label=f"pentaploid")

        if normalise:
            pylab.ylabel("Somy", fontsize=16)

        if hybrid:
            from pylab import fill_betweenx, xlim, ylim

            Ymax = ylim()[1]
            Nchrom = len(self.df["chr"].unique())
            for x in range(int(Nchrom / 2)):
                fill_betweenx([0, Ymax], -0.5 + x * 2, 0.5 + x * 2, color="grey", alpha=0.5)
            ylim([0, Ymax])
            xlim([-0.5, Nchrom - 0.5])

        pylab.savefig(outfile)


def plot_sliding_window_boxplot(data, window_size, step=1000, overlap_percentage=50, facecolor="lightblue"):
    """
    Plot a boxplot with sliding windows from a specified column in a DataFrame.

    Parameters:
        data (pd.DataFrame): The input DataFrame.
        column (str): The column of interest.
        window_size (int): The size of each sliding window.
        overlap_percentage (float): The percentage overlap between windows (default 50).


    Example::

        ss = SomyScore(input_bam)
        ss.compute_coverage()
        chrom_name = 1
        data = ss.df.query("chr==@chrom_name")['depth']
        overlap = 50 #percentage
        Nwin = 20  # number of windows of size (Step) to use per boxplot
        step=ss.window_size # value used in the SomyScore
        plot_sliding_window_boxplot(data, Nwin, step, overlap)

    """
    step_size = int(window_size * (1 - overlap_percentage / 100))  # Calculate step size
    print(f"step is {step_size*step}; step_size={step_size}")
    # Create sliding windows
    windows = [data[i : i + window_size].values for i in range(0, len(data) - window_size + 1, step_size)]
    print(len(data) - window_size + 1, step_size)
    print(f"Number of windows = {len(windows)}")
    if not windows:
        raise ValueError("Window size too large for the given data.")

    # Plot boxplot
    # pylab.clf()
    pylab.boxplot(windows, patch_artist=True, boxprops=dict(facecolor=facecolor))
    pylab.xlabel("Sliding Window")
    pylab.ylabel("Values")
    return windows
