#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 09:42:35 2018

combine_yeast_gene_expression.py
This takes the gene expression exported from the R script
"exportGrowthData2Python.R"
the gene expression downloaded from Slavov et al.,
"ugrr_data.txt"
and the holstege data downloaded from GEO, compiled in 
"readHolstegeData.py",
and combines them into a single data frame for
processing by "gather_data.py"
@author: thomaspwytock
"""

import pandas as pa
import numpy as np
import sys
def combine_yeast_gene_expression(LO_HUGHES=False):
    """
combine_yeast_gene_expression
combines the five different datasets used for calculating correlations
between yeast genes.
The Holstege dataset needs to be expressed in terms of Yeast ORFs
instead of the in-house probesets before integration with the other
datasets.
The combined dataset is saved, along with the corresponding metadata.
RETURNS:
    0
"""
    allMetaData = pa.read_excel('data/allmetadata.xlsx',index_col=0)
    ## gather gene expression data
    charlesGenExp = pa.read_table('CharlesGenExp.tsv',header=0,index_col=0)
    greshamGenExp = pa.read_table('GreshamGenExp.tsv',header=0,index_col=0)
    if not LO_HUGHES: 
        hughesGenExp = pa.read_table('HughesGenExp.tsv',header=0,index_col=0)
    slavovGenExp = pa.read_table('ugrr_data.txt',header=0,index_col=0)
    holstegeGenExp = pa.read_pickle('holstege_data_wt_responsive.pkl')
    hge_ind = holstegeGenExp.index
    holstegeFormat=pa.read_csv('data/holstege_lab_yeast_microarray_format.csv',index_col=0)
    holstegeGenExp.loc[hge_ind,'SystematicName'] = holstegeFormat.loc[hge_ind,'SYSTEMATIC_NAME']
    holstegeGenExpGB = holstegeGenExp.groupby('SystematicName').mean()
    charlesGenExp.columns = ['Charles_%d' % ii for ii in range(charlesGenExp.shape[1])]
    greshamGenExp.columns = [c.split('.')[0] for c in greshamGenExp.columns]
    slavovGenExp.columns = allMetaData[allMetaData.Group=='Slavov'].index.tolist()
    if not LO_HUGHES: 
        hughesGenExp.columns = allMetaData[allMetaData.Group=='Hughes'].index.tolist()
        common_genes = holstegeGenExpGB.index.intersection(
                           slavovGenExp.index).intersection(
                           hughesGenExp.index).intersection(
                           greshamGenExp.index).intersection(
                           charlesGenExp.index)
        GEXP = pa.concat([charlesGenExp,
                          greshamGenExp,
                          hughesGenExp,
                          slavovGenExp,
                          holstegeGenExpGB],axis=1).loc[common_genes]
        GEXP.to_pickle('all_yeast_data.pkl')
    else:
        common_genes = holstegeGenExpGB.index.intersection(
                           slavovGenExp.index).intersection(
                           greshamGenExp.index).intersection(
                           charlesGenExp.index)
        GEXP = pa.concat([charlesGenExp,
                          greshamGenExp,
                          slavovGenExp,
                          holstegeGenExpGB],axis=1).loc[common_genes]
        GEXP.to_pickle('hughes_removed_yeast_data.pkl')
    z1 = GEXP
    gndata = pa.read_excel('holstege_gene_spreadsheet.xlsx',index_col=1)
    responsive_mutants = gndata[gndata.iloc[:,7]=='responsive mutant']
    holst_group = allMetaData[allMetaData.Group=='Holstege']
    #hughes_group = allMetaData[allMetaData.Group=='Hughes']
    inds2rem = [i for i,row in holst_group.iterrows() if row.loc['Gene1'] not in responsive_mutants.index.tolist()]
    z1p = z1.loc[:,[col for col in z1.columns if col not in inds2rem]]
    metadata_nr_removed = allMetaData.loc[z1p.columns]
    if not LO_HUGHES:    
        z1p.to_pickle('data_nonresponsive_removed.pkl')
        metadata_nr_removed.to_pickle('metadata_nonresponsive_removed.pkl')
    else:
        z1p.to_pickle('data_hughes_removed.pkl')
        metadata_nr_removed.to_pickle('metadata_hughes_removed.pkl')
    return 0

if __name__ == '__main__':
    try:
        LO_HUGHES = sys.argv[1].startswith(('T','t'))
    except ValueError:
        LO_HUGHES=False
    combine_yeast_gene_expression(LO_HUGHES)
