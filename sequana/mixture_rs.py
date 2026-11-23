#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2026 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Rust-accelerated Gaussian mixture EM (reference: :class:`sequana.mixture.EM`)."""
import numpy as np
from easydev import AttrDict

from sequana import criteria
from sequana.mixture import EM

try:
    import sequana_cov_rs
except ImportError:  # pragma: no cover
    sequana_cov_rs = None


class EMRS(EM):
    """:class:`~sequana.mixture.EM` whose iteration loop runs in Rust."""

    def estimate(self, guess=None, k=2):
        if sequana_cov_rs is None:  # pragma: no cover
            return super().estimate(guess=guess, k=k)

        self.k = k
        if guess is None:
            guess = self.get_guess()

        out = sequana_cov_rs.em_fit(
            np.ascontiguousarray(self.data, dtype=np.float64),
            [float(x) for x in guess[0::3]],
            [float(x) for x in guess[1::3]],
            [float(x) for x in guess[2::3]],
            int(self.max_iter),
        )
        self.status = bool(out["ok"])
        self.gamma = np.asarray(out["gamma"]).reshape(len(out["mu"]), -1)
        mu = np.array(out["mu"])
        self.mus = [mu] * out["counter"]  # the reference code appends the same array each time
        if not self.status:
            print("issue arised at iteration %s" % out["counter"])
            self.debug = {"N": np.array(out["n_k"]), "pis": np.array(out["pi"])}

        self.results = AttrDict(
            x=out["x"], nfev=out["counter"], success=self.status and out["counter"] >= self.max_iter
        )
        self.results.mus = self.results.x[0::3]
        self.results.sigmas = self.results.x[1::3]
        self.results.pis = self.results.x[2::3]

        log_likelihood = sequana_cov_rs.neg_log_likelihood(
            np.ascontiguousarray(self.data, dtype=np.float64), [float(v) for v in self.results.x]
        )
        self.results.log_likelihood = log_likelihood
        self.results.AIC = criteria.AIC(log_likelihood, self.k, logL=True)
        self.results.AICc = criteria.AICc(log_likelihood, self.k, self.data.size, logL=True)
        self.results.BIC = criteria.BIC(log_likelihood, self.k, self.data.size, logL=True)
