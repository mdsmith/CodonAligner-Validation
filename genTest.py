#!/usr/bin/env python3.2

from Bio.Seq import translate

from BioExt import untranslate, randgene
from BioExt.stats import bbinom

from scipy.stats import rv_discrete

samples = [5073, 1337, 380, 131, 36, 20, 1]

n = len(samples) - 1

b = bbinom(*bbinom.fit(samples))

custom = rv_discrete(name='precompbbinom', values=(list(range(n+1)),
            [b.pmf(i) for i in range(n+1)]))

# some multiple of 3 to keep frame
randseq = randgene(9999, custom.ppf)

randseqperm = untranslate(translate(randseq))

print(randseq)
