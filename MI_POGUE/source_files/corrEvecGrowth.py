""" 
corrEvecGrowth.py

This file takes in: 
    --the labels of the selected eigengenes in the noisy and noiseless cases:
        'tmp/corr_feat_20-carrera-noise-avg.pkl' (noisy)
        'tmp/corr_feat_50-carrera-corr-avg.pkl' (noiseless)
    --the gene expression file (expressed as eigenvectors):
        'SF1-EcoMAC/ecoli_compendium_corr_df.pkl'
    --the associated annotation file:
        'msb145108-sup-0005-DataSetS2.xlsx'

This file outputs the correlation of these eigengenes with growth:
    --Table-SXX_Rank_eigengene_correlation.xlsx
"""
import pandas as pd
import numpy as np
import sys

def main(NAME='carrera-corr',nopt=20):
    feats50 = pd.read_pickle('corr_feat_%d-%s_feat-avg.pkl' % (nopt,NAME))
    featsNoise = pd.read_pickle('data/corr_feat_20-carrera-noise-avg.pkl')    
    A = pd.read_pickle('ecoli_compendium_gncorr_%s_df.pkl' % NAME)
    ECOLI_GFN = 'data/DataSetS1_CarreraAnnotations.xlsx'
    ecannot = pd.read_excel(ECOLI_GFN,sheetname=u'EcoMAC and EcoPhe',index_col='CEL file name')
    rem_rows = ecannot.index.get_level_values(0).isin(['PU00524','PU00586'])
    ecannot = ecannot[~rem_rows]
    _nz = ecannot['Flag growth']>0
    ecannot_phe = ecannot[_nz]
    ec_gr = ecannot_phe.loc[:,'Growth rate (1/h)']

    Ap = A.loc[ec_gr.index]
    cc_d = {} 
    for ii in Ap.columns:
	cc_d[ii] = Ap.loc[:,ii].corr(ec_gr,method='spearman')
    cc_ser = pd.Series(cc_d)
    DF = pd.DataFrame({'val':cc_ser.loc[np.hstack([feats50[:9],featsNoise[:5]])],
     'rank':cc_ser.abs().rank(ascending=False).loc[np.hstack([feats50[:9],featsNoise[:5]])]})
    print DF
    DF.to_excel('Table-SXX_Rank_eigengene_correlation.xlsx')

if __name__ == '__main__':
    if len(sys.argv)< 2:
        main()
    else:
        main(sys.argv[1],int(sys.argv[2]))
