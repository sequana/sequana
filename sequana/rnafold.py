from sequana.lazy import numpy as np


def initialize_matrix(seq):
    """Initialize a matrix filled with zeros."""
    n = len(seq)
    return np.zeros((n, n), dtype=int)


def is_complementary(base1, base2):
    """Check if two bases are complementary."""
    pairs = {"A": "U", "U": "A", "G": "C", "C": "G"}
    return pairs.get(base1) == base2


def calculate_mfe(seq):
    """Calculate the Minimum Free Energy (MFE) score."""
    n = len(seq)
    dp = initialize_matrix(seq)  # Ensure this returns a NumPy array

    # Fill the dp matrix using the dynamic programming approach
    for k in range(1, n):
        for i in range(n - k):
            j = i + k

            if is_complementary(seq[i], seq[j]):

                dp[i, j] = max([dp[i + 1, j - 1] + 1, dp[i + 1, j], dp[i, j - 1]])
            else:

                dp[i, j] = max([dp[i + 1, j], dp[i, j - 1]])

    # Return the MFE score, which is stored in dp[0, n-1]
    return -dp[0, n - 1]
