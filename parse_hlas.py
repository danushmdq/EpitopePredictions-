import skbio
import numpy as np
import pandas as pd
import sys
sys.path.append("/Users/dmonaco/CodeRepository+/")
import HLAPredCache
from dani_helpers import *


if __name__ == '__main__':
    hlaFn = sys.argv[1]
    outFn = sys.argv[2]

    # hlaFn = '/home/agartlan/gitrepo/transpm/data/All HLA Typing Data -.xlsx'
    hlaDf = parseHLAFile(hlaFn)
    uHLAs = extractUniqueHLAs(hlaDf.loc[hlaDf['Site'] != 'R'])
    print(uHLAs)
    print(len(uHLAs))

    '''goodHLAs, badHLAs = HLAPredCache.checkHLAs(['HLA-A*02:01', 'HLA-A*30:14', 'B*58:01', 'HLA-B*58:01'],
                                  lengths=[8],
                                  verbose=True)'''
    
    """goodHLAs, badHLAs = checkHLAs(uHLAs, lengths=[9], verbose=True)"""
    
    goodHLAs = [h for h in uHLAs if not h == 'HLA-A*30:14']

    convertedHLAs = [ hlatoafg  (h) for h in goodHLAs ]
    
    with open(outFn, 'w') as fh:
        for h in convertedHLAs:
            fh.write('{}\n'.format(h))