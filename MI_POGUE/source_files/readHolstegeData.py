#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
readHolstegeData.py
reads in Holstege dataset from series matrix files downloaded from GEO

@author: thomaspwytock
"""

import pandas as pd
import gzip as gz
import sys

def count_header_lines(fn):
    ## use for correct handling of gzip files
    op = gz.open if fn.endswith('gz') else open        
    with op(fn,'r') as fh:
        x = fh.readline(); cnt = 0
        while not x.startswith('!series_matrix_table_begin'):
            x = fh.readline()
            cnt+=1
    return cnt
            
def read_series_matrix_files(fn):
    nhead = count_header_lines(fn)
    comp = 'gzip' if fn.endswith('gz') else 'infer'
    X = pd.read_table(fn,compression=comp,sep='\t',header=nhead)
    X = X.iloc[:-1,:]
    X.ID_REF = X.ID_REF.astype(int)
    X = X.set_index("ID_REF")
    X = X.astype(float)
    return X

def main():
    fn_list = ['./GSE42215_series_matrix.txt.gz',
               './GSE42217_series_matrix.txt.gz',
               './GSE42240_series_matrix.txt.gz',
               './GSE42241_series_matrix.txt.gz',
               './GSE42527_series_matrix.txt.gz']
    holstege_data = pd.concat(map(read_series_matrix_files,fn_list),axis=1)
    holstege_data.to_pickle('holstege_data_wt_responsive.pkl')

if __name__ == '__main__':
    main()
