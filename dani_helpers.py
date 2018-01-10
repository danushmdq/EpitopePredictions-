import skbio
import numpy as np
import pandas as pd

def filterLRSeqs(seqList):
    lrSeqs = []
    for s in seqList:
        if '_' in s.metadata['id'] and s.metadata['id'].split('_')[1] == 'LR':
            lrSeqs.append(s)
    return lrSeqs

def fasta2mers(fastaFilename):
    aln = skbio.io.read(fastaFilename, format='fasta')
    return [seq for seq in aln]

def parseHLAFile(hlaFilename):
    def _makeList(row):
        hlaL = [row[col] for col in hlaCols if not pd.isnull(row[col]) and not row[col] == '']
        newHLAL = []
        for h in hlaL:
            locus, allele = h.split('*')
            digits = allele.split(':')
            if '(' in digits[1]:
                pos = digits[1].find('(')
                lh = digits[0]
                possible = digits[1][pos+1:-1].split(',')
                for p in possible:
                    newHLAL.append('{}*{}:{}'.format(locus, lh, p))
            else:
                newHLAL.append(h)
        return list(set(newHLAL))

    hlaDf = pd.read_excel(hlaFilename)
    hlaCols = ['A1', 'A2', 'B1', 'B2', 'C1', 'C2']
    for col in hlaCols:
        # hlaDf.loc[:, col] = hlaDf[col].map(formatHLA)
        hlaDf.loc[:, col] = hlaDf.loc[:, col].map(formatHLA)
    hlaDf.loc[:, 'hlaL'] = hlaDf.apply(_makeList, axis=1)
    return hlaDf

def extractUniqueHLAs(hlaDf, dropXX=True):
    uHLAs = set()
    for hlaL in hlaDf['hlaL']:
        for h in hlaL:
            uHLAs.add(h)
    uHLAs = sorted([h for h in uHLAs])
    if dropXX:
        uHLAs = [h for h in uHLAs if not 'X' in h]
    return uHLAs

    '''uHLAs = []
    for hlaL in hlaDf['hlaL']:
        uHLAs.extend(hlaL)
    return np.unique(uHLAs)'''

def formatHLA(h):
    if pd.isnull(h):
        return ''
    try:
        locus, allele = h.split('*')
        digits = allele.split(':')
    except AttributeError:
        print(h)
        return str(h)

    """Discard w in Cw"""
    locus = locus[0]

    if len(digits) == 1:
        a = digits[0]
        lh = a[:2]
        if '(' in a:
            parenPos = a.find('(')
            rh = a[parenPos:]
        elif len(digits[0]) == 4:
            rh = digits[0][2:]
        elif len(digits[0]) > 4:
            rh = digits[0][2:4]
        else:
            """Two-digit alleles will have XX"""
            rh = 'XX'
    else:
        lh = digits[0]
        rh = digits[1]
    return 'HLA-{locus}*{lh}:{rh}'.format(locus=locus, lh=lh, rh=rh)

def hlatoafg(h):

    """ HLA-A*02:01      A_0201"""

    h= h.replace("HLA-","")
    h= h.replace(":","")
    h= h.replace("*","_")

    return h
