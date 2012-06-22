#!/usr/bin/env python3.2

import numpy as np

from scipy.special import beta
from scipy.misc import comb
from scipy.optimize import fmin_cobyla


class BetaBinom(object):

    @staticmethod
    def pmf(n, k, a=None, b=None):
        return comb(n, k) * beta(k + a, n - k + b) / beta(a, b)


class BetaBinomMoMEst(BetaBinom):

    def __init__(self, samples):
        n = len(samples) - 1
        tot = sum(samples)
        m1 = sum(i * samples[i] for i in range(len(samples))) / tot
        m2 = sum((i ** 2) * samples[i] for i in range(len(samples))) / tot
        ah = (n * m1 - m2) / (n * ((m2 / m1) - m1 - 1) + m1)
        bh = (n - m1) * (n - (m2 / m1)) / (n * ((m2 / m1) - m1 - 1) + m1)
        self._ah = ah
        self._bh = bh

    def pmf(self, n, k):
        return super(BetaBinomMoMEst, self).pmf(n, k, self._ah, self._bh)

    @property
    def alpha(self):
        return self._ah

    @property
    def beta(self):
        return self._bh


class BetaBinomMLEst(BetaBinom):

    def __init__(self, samples):
        n = len(samples) - 1
        tot = sum(samples)
        pmf = super(BetaBinomMLEst, self).pmf
        # function to minimize
        def f(x):
            pred = (tot * pmf(n, k, x[0], x[1]) for k in range(len(samples)))
            diff = ((p - samples[k]) ** 2 for k, p in enumerate(pred))
            return sum(diff) # L2 error
        # initial guess with method of moments
        mom = BetaBinomMoMEst(samples)
        x0 = np.array([mom.alpha, mom.beta])
        xn = fmin_cobyla(f, x0, [lambda x: x[0], lambda x: x[1]], rhobeg=1e-1, rhoend=1e-10, disp=0)
        self._ah = xn[0]
        self._bh = xn[1]

    def pmf(self, n, k):
        return super(BetaBinomMLEst, self).pmf(n, k, self._ah, self._bh)

    @property
    def alpha(self):
        return self._ah

    @property
    def beta(self):
        return self._bh


if __name__ == '__main__':
    samples = [5073, 1337, 380, 131, 36, 20, 1] # [3, 24, 104, 286, 670, 1033, 1343, 1112, 829, 478, 181, 45, 7]
    n = len(samples) - 1
    mom = BetaBinomMoMEst(samples)
    mle = BetaBinomMLEst(samples)
    tot = sum(samples)
    momp = [tot * mom.pmf(n, k) for k in range(len(samples))]
    diff = [(p - samples[k]) ** 2 for k, p in enumerate(momp)]
    print(momp)
    print(np.sqrt(sum(diff)))
    mlep = [tot * mle.pmf(n, k) for k in range(len(samples))]
    diff = [(p - samples[k]) ** 2 for k, p in enumerate(mlep)]
    print(mlep)
    print(np.sqrt(sum(diff)))
