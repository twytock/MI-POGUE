"""
evec2Similarity.py

This file takes the as input
    --a file containing the labels of the selected eigenvectors:
        'tmp/corr_feat_50-carrera-corr-avg.pkl'
    --matrix mapping genes (rows) to eigengenes (columns), i.e. P 
      of P D P^{-1}=C, where C is the correlation matrix:
        'tmp/evecnz_carrera.pkl'

It (1) converts each selected eigengenes into similarity matrices that are
stored in "expanded_evec_<n>.pkl" where n = 0..8.

Then it (2) exports the similarity files to a (gzipped) R list of dataframes
called "L1.gzip".
"""

import numpy as np
from rpy2.robjects import r
from rpy2.robjects import pandas2ri
import pandas as pd
import cPickle as cp
from glob import glob
import sys

def get_file_name(NAME):
    fn_l = glob('corr_feat_*-%s_feat-avg.pkl' % NAME)
    sfn_l = sorted(fn_l, key = lambda fn: int(fn.split('_')[2].split('-')[0]))
    return sfn_l[-1]

def main(NAME='carrera-corr',Nopt=9):
    feat_fn = get_file_name(NAME)
    feats = pd.read_pickle(feat_fn)
    evecs = pd.read_pickle('ecoli_evecs_%s_df.pkl' % NAME)
    ## only the first 9 are needed for the optimal model
    ccp_l = []
    for ii in range(Nopt):
        X = evecs.loc[:,feats[ii:ii+1]]
        #K = np.array(map(float,featsNoise[ii]))
        V = (X.values)
        CCp = np.dot(V,V.T)
        #CCp[np.diag_indices_from(CCp)]=1
        scale = np.amax(CCp[np.diag_indices_from(CCp)])
        #print scale
        CCp = CCp/scale
        CCp = pd.DataFrame(CCp,index=X.index,columns=X.index)
        ccp_l.append(CCp)
        #CCp.to_pickle('expanded_evec_%d.pkl' % ii)
    
    pandas2ri.activate()
    for ii,data in enumerate(ccp_l):
        rdata = pandas2ri.py2ri(data)
        r.assign('X%d' % ii,rdata)
    ft_list_str = ','.join(['X%d' % i for i in range(Nopt)])
    r("L1=list(%s); save(L1,file='L1.gzip',compress=TRUE)" % ft_list_str)

if __name__ == '__main__':
    if len(sys.argv) > 2:
        NAME=sys.argv[1]
        Nopt = int(sys.argv[2])
    else:
        NAME='carrera-corr'
        Nopt = 9
    main(NAME,Nopt)
