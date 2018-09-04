"""
importCarreraData.py 
converts the text files associated with Carrera et al. Mol. Syst. Biol. 2014
"An integrative multi-scale genome-wide model reveals the phenotypic landscape of E. coli."
into a pandas data frame.

It takes in a text file of the gene expression,
    --SF1-EcoMAC/Ecoli-compendium-assembly-quantilenorm.txt
the column labels
    --SF1-EcoMAC/numbering_arrays-compendium.txt
and the row labels
    --SF1-EcoMAC/numbering_assembly-tfs-enzymes-genes.txt
which are extracted by default into the directory SF1-EcoMAC when the Carrera dataset is downloaded.

This file outputs a pandas data frame of the gene expression data to:
    --SF1-EcoMAC/ecoli_compendium_df.pkl
"""

import pandas as pd

def main():
    X = pd.read_table('SF1-EcoMAC/Ecoli-compendium-assembly-quantilenorm.txt',header=None)
    columns_list = pd.read_table('SF1-EcoMAC/numbering_arrays-compendium.txt',header=None,
                                 sep=' ',index_col=0)
    index_list = pd.read_table('SF1-EcoMAC/numbering_assembly-tfs-enzymes-genes.txt',
                               header=None,index_col=0)
    MI = pd.MultiIndex.from_arrays(index_list.values.T)
    X.index=MI; X.columns=columns_list.values.ravel();
    X.to_pickle('SF1-EcoMAC/ecoli_compendium_df.pkl')

if __name__=='__main__':
    main()
