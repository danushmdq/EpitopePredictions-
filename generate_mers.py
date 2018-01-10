import skbio
import numpy as np
import pandas as pd
import time
import requests
import io

import HLAPredCache
from dani_helpers import *

import sys

if __name__ == '__main__':
    # seqFn = GIT_PATH + 'transpm/data/Aligned Gag Sequences (N=79) and Subtype References -.fasta'
    # seqFn = '/Users/dmonaco/Desktop/Paper 2 +/Adaptive Evolution +/Gag +/Aligned Gag Sequences (N=79) and Subtype References -.fasta'
    # seqFn = '/home/agartlan/gitrepo/transpm/data/Aligned Gag Sequences (N=79) and Subtype References -.fasta'

    seqFn = sys.argv[1]
    outFn = sys.argv[2]

    seqList = fasta2mers(seqFn)
    lrSeqs = filterLRSeqs(seqList)
    mers = generateMersFromNT(lrSeqs)
    with open(outFn, 'w') as fh:
        for m in mers:
            fh.write('{}\n'.format(m))