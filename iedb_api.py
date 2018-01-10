"""
python iedb_api.py seqs_nt.fasta hlafile.hla
"""
import skbio
import numpy as np
import pandas as pd
import time
import sys

sys.path.append('/home/agartlan/gitrepo')
import HLAPredCache
from dani_helpers import *

if __name__ == '__main__':
    mersFn = sys.argv[1]
    hlaFn = sys.argv[2]
    outFn = sys.argv[3]

    with open(mersFn, 'r') as fh:
        raw = fh.read()
    mers = raw.split('\n')[:-1]

    with open(hlaFn, 'r') as fh:
        raw = fh.read()
    hlas = raw.split('\n')[:-1]

    lengths = sorted(list(set([len(m) for m in mers])))

    results = []
    for l in lengths:
        tmpMers = [p for p in mers if len(p) == l]
        res = HLAPredCache.iedbPepPredict(hlas, tmpMers, method='netmhcpan', timeit=True)
        results.append(res)
    outDf = pd.concat(results, axis=0)

    outDf.to_csv(outFn[['allele', 'peptide', 'ic50']], index=False)





