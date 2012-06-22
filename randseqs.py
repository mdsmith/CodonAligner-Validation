
import numpy as np

from collections import Counter
from operator import itemgetter
from random import randint, uniform

from bbinom import BetaBinomMLEst as MLE

samples = [5073, 1337, 380, 131, 36, 20, 1]

mle = MLE(samples)

tot = sum(samples)

mlep = [mle.pmf(k) for k in range(len(samples))]
# print(mlep)

seq = ''
l = 0
ls = []
char = -1
last_char = -1
avoid = -1
while l < 1000:
    lp = mle.ppf(uniform(0, 1)) + 1
    # avoid making stop codons in frame
    if l % 3 == 1 and seq[-1] == 'T' and lp > 1:
        avoid = 0 # avoid TAA
    if l % 3 == 2:
        if seq[-2:] == 'TA':
            avoid = 2 # avoid TAG
        elif seq[-2:] == 'TG':
            avoid = 0 # avoid TGA
    # don't permit artifically longer homopolymers
    # due to sampling the same letter twice in a row
    while char == last_char or char == avoid:
        char = randint(0, 3)
    # set some vars
    last_char = char
    avoid = -1
    # grow our string
    seq += 'ACGT'[char] * lp
    l += lp
    ls.append(lp)

# props = [v / len(ls) for _, v in sorted(Counter(ls).items(), key=itemgetter(0))]
print(seq)
