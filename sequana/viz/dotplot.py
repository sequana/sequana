import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def dotplot_from_lastz():

    """
    lastz 5_tofix.fa 5_tofix.fa --notransition --strand=both --step=20 --nogapped --format=rdotplot > temp.txt
    """
    # Read command-line arguments
    infile_name = "temp.txt"
    outfile_name = "temp.pdf"
    text_for_the_title = "TEST"

    # Read the input file (assumes tab-delimited with a header)
    df = pd.read_csv(infile_name, sep="\t")

    # Extract the first two columns
    x = df.iloc[:, 0]
    y = df.iloc[:, 1]

    # Define a function to safely compute the max while ignoring NaNs
    def my_max(series):
        return series.max(skipna=True) if not series.isna().all() else np.nan

    # Compute max values
    maximum_value_of_x = my_max(x)
    maximum_value_of_y = my_max(y)

    # Plot
    plt.figure(figsize=(10, 10))  # inches
    plt.plot(x, y, linestyle="-")
    plt.xlim(0, maximum_value_of_x)
    plt.ylim(0, maximum_value_of_y)
    plt.title(text_for_the_title, fontsize=20)
    plt.xlabel("position (bp)", fontsize=15)
    plt.ylabel("position (bp)", fontsize=15)
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(outfile_name)
