import sys,os,os.path as osp,getopt
import pandas as pa
import numpy as np

help_message = '''
This script takes a tab-delimited file of gene expression data
and a tab-delimited file of growth rate data as input.
The column identifiers of the gene expression file must match
with the row identifiers of the growth rate file.
It applies MI-POGUE to develop a mapping between gene expression
and growth rate.
Possible options:
-h, --help : see this help message
-e, --expression-file=: name of the gene expression file
-N, --name=: Filename extension for the results
-r, --growth-rate=: name of the growth rate file
'''

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

## need to weight by slower-growth states


def get_gsm(fn):
    return osp.split(fn)[1].split('.')[0].split('_')[0]

def annotate_data(FN):
    """
annotate_data takes in gene expression data (represented as a pkl generated from readHolstegeData.py) 
and organizes them into a pandas data frame (saved as a pickle in FN) for analysis in python.
holstege_data_wt_responsive.pkl -- gene expression of yeast strains from the deletion library and wild type
pertdata_targets.txt -- names of the deletions.
all_targets.txt -- which channel was wt and which was the pool
Inputs:
    FN -- the filename of the compressed data frame
Outputs:
    alldat -- the dataframe of all data with rows labelled by genes and columns
    labelled by GSMs and deletions.
"""
    data = pa.read_pickle('../holstege_data_wt_responsive.pkl')
    ma_lay = pa.read_csv('../data/holstege_lab_yeast_microarray_format.csv')
    MI=pa.MultiIndex.from_tuples(zip(ma_lay.ID,ma_lay.PROBE_ID,ma_lay.GENE_SYMBOL))
    pc_targs = pa.read_table('../data/pertdata_targets.txt')
    pc_targs.loc[:,'GSMS'] = map(get_gsm,pc_targs.Names)
    pc_targs = pc_targs.set_index('GSMS')
    pc_gsm = data.columns.isin(pc_targs.index)
    pc_targs = pc_targs.loc[pc_gsm]
    pc_perts = [pc_targs.loc[gsm,'experiment Cy3']
                if not pc_targs.loc[gsm,'experiment Cy3'] == 'ref'
                else pc_targs.loc[gsm,'experiment Cy5'] for gsm in pc_gsm]
    pc_col_tups = zip(pc_gsm.tolist(),pc_perts)
    wtc_targs = pa.read_table('all_targets.txt')
    wtc_targs.loc[:,'GSMS'] = map(get_gsm,wtc_targs.Names)
    wtc_targs = wtc_targs.set_index('GSMS')
    wtc_gsm = data.columns.isin(wtc_targs.index)
    wtc_targs = wtc_targs.loc[wtc_gsm]
    wtc_perts = [wtc_targs.loc[gsm,'experiment Cy3']
                 if not wtc_targs.loc[gsm,'experiment Cy3'] == 'ref'
                 else wtc_targs.loc[gsm,'experiment Cy5'] for gsm in wtc_gsm]
    wtc_col_tups= zip(wtc_gsm,wtc_perts)
    comb_col_tups = wtc_col_tups+pc_col_tups
    comb_col_series = pa.DataFrame(comb_col_tups).set_index(0)
    array2ind = comb_col_series.loc[data.columns].reset_index().values
    col_MI = pa.MultiIndex.from_array(array2ind.T)
    data.index = MI
    data.columns = col_MI
    data.to_pickle(FN)
    return data

def fill_nan(df):
    for ind,row in df.iterrows():
        newrow = row.fillna(row.mean())
        df.loc[ind,:]=newrow
    return df

def load_ypd_growth_data():
    ## Archive_regression was downloaded from:
    ## http://www-deletion.stanford.edu/YDPM/YDPM_index.html
    ## this file is currently not used
    AA = pa.read_table('./Archive_regression/Regression_Tc1_hom.txt')
    AA2 = pa.read_table('./Archive_regression/Regression_Tc2_hom.txt')
    AA = AA.set_index('orf')
    AA2 = AA2.set_index('orf')
    return (AA.YPD + AA2.YPD)/2.

def trim_evals(vals,vecs,TOL=1e-6):
    nzi = np.where(vals>TOL)[0]
    valnz = vals[nzi]
    vecnz = vecs[:,nzi]
    vecnz = vecnz[:,::-1]
    valnz = np.sqrt(valnz[::-1])
    vallbl = ['{0:.6f}'.format(x) for x in valnz]
    return valnz,vallbl,vecnz

def correlated_data_yeast(z1,NAME):
    ## restrict to genes
    if isinstance(z1.index,pa.MultiIndex):
        z1_gene = z1.groupby(level='SystematicName').mean()
    else:
        z1_gene = z1
    CC_gene = np.corrcoef(z1_gene)
    evl_gn,evc_gn = np.linalg.eigh(CC_gene)
    np.save('sacCer_eval_gn_sqr_%s.npy' % NAME,evl_gn)
    np.save('sacCer_evec_gn_%s.npy' % NAME,evc_gn)
    evl_gn_nz,evl_gnlbl_nz,evc_gn_nz = trim_evals(evl_gn,evc_gn)
    evc_gn_df=pa.DataFrame(evc_gn_nz,index=z1_gene.index,columns=evl_gnlbl_nz)
    z1_gene_corr = z1_gene.T.dot(evc_gn_df)
    z1_gene_corr.to_pickle('sacCer_compendium_gncorr_%s_df.pkl' % NAME)
    evc_gn_df.to_pickle('yeast_evecs_%s_df.pkl' % NAME)
    return z1_gene_corr, evc_gn_df

def correlated_data_ecoli(z1_gene,NAME):
    ## restrict to genes
    CC_gene = np.corrcoef(z1_gene)
    evl_gn,evc_gn = np.linalg.eigh(CC_gene)
    np.save('ecoli_eval_gn_sqr_%s.npy' % NAME,evl_gn)
    np.save('ecoli_evec_gn_%s.npy' % NAME,evc_gn)
    evl_gn_nz,evl_gnlbl_nz,evc_gn_nz = trim_evals(evl_gn,evc_gn)
    evc_gn_df=pa.DataFrame(evc_gn_nz,index=z1_gene.index,columns=evl_gnlbl_nz)
    z1_gene_corr = z1_gene.T.dot(evc_gn_df)
    z1_gene_corr.to_pickle('ecoli_compendium_gncorr_%s_df.pkl' % NAME)
    evc_gn_df.to_pickle('ecoli_evecs_%s_df.pkl' % NAME)
    return z1_gene_corr, evc_gn_df


def main():
    try:
        FN = 'SF1-EcoMAC/ecoli_compendium_df.pkl'
        argv = sys.argv
        NAME= 'carrera-corr'
        org = 'ecoli'
        try:
            opts, args = getopt.getopt(argv[1:], "he:N:o:",
                                       ["help","expression-file=","name=",
                                        "organism="])
        except getopt.error, msg:
            raise Usage(msg)
        # option processing
        for option, value in opts:
            ## BOUNDONLEFT True in US, False in England
            WEIGHTS = 'uniform' if option == "-u" else 'distance'
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option in ("-e", "--expression-file"):
                FN = value
            if option in ("-N","--name"):
                NAME=value
            if option in ("-o","--organism"):
                if value.startswith(('E','e')):
                    org = 'ecoli'
                if value.startswith(('Y','y','S','s')):
                    org = 'saccer'
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"


    #alldat = annotate_data(FN.replace('csv','pkl'))
    if FN.endswith('pkl'):
        alldat = pa.read_pickle(FN)
        alldat = fill_nan(alldat)
        #alldat.to_csv(FN.replace('csv','pkl'))
    else:
        alldat = pa.read_table(FN)
        alldat = fill_nan(alldat)
    if FN.endswith('csv'):
        cols= alldat.columns.tolist(); cols[0]='SystematicName'
        alldat.columns = cols; alldat = alldat.set_index('SystematicName')
    if org=='saccer':
        __,__ =correlated_data_yeast(alldat,NAME)
    elif org=='ecoli':
        __,__ =correlated_data_ecoli(alldat,NAME)
    return 0

if __name__ == '__main__':
    main()
