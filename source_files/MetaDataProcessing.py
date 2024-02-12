# -*- coding: utf-8 -*-
"""
Spyder Editor

MetaDataProcessing.py
gathers gene expression data from several sources
"""
import pandas as pd
import numpy as np

## process the Charles DataFrame
col_lbls = {'Column':{},'Treatment1':{},'Growth_Rate':{},'Time':{},'Units':{},
            'Name':{},'Gene1':{},'Treatment2':{},'Gene2':{},'Group':{}}

kw_l = ['Treatment1','Growth_Rate','Time','Units']
with open("CharlesColMetaData.tsv") as fh:
    #Sample line: Heatshock_0.25_120_min
    ## code "unstress" as zero
    for col,ln in enumerate(fh):
        col_lbls['Column'][col]=col
        col_lbls['Group'][col]='Charles'
        for ii,elt in enumerate(ln.rstrip().split('_')):
            kw = kw_l[ii]
            if elt == 'unstress':
                elt =0
            elt =  elt if kw not in ('Growth_Rate','Time') else float(elt)
            col_lbls[kw][col]=elt
        col_lbls['Gene1'][col]='NA'
        col_lbls['Gene2'][col]='NA'
        col_lbls['Treatment2'][col]='NA'
        col_lbls['Name'][col]='Charles_%d' % col

CharlesMetaData = pd.DataFrame(col_lbls)

## process the Gresham DataFrame
col_lbls = {'Column':{},'Treatment1':{},'Growth_Rate':{},'Time':{},'Units':{},
            'Name':{},'Gene1':{},'Treatment2':{},'Gene2':{},'Group':{}}
kw_l = ['Name','Treatment1','Growth_Rate','Time','Units']
with open("GreshamColMetaData.tsv") as fh:
    #Sample line: CDG5a_Allantoin_0.06
    ## code "unstress" as zero
    for col,ln in enumerate(fh):
        col_lbls['Column'][col]=col
        col_lbls['Group'][col]='Gresham'
        for ii,elt in enumerate(ln.rstrip().split('_')):
            kw = kw_l[ii]
            elt = float(elt) if kw in ('Growth_Rate','Time') else elt            
            col_lbls[kw][col]=elt
        col_lbls['Time'][col]=0
        col_lbls['Units'][col]='h'
        col_lbls['Gene1'][col]='NA'
        col_lbls['Gene2'][col]='NA'
        col_lbls['Treatment2'][col]='NA'
       
GreshamMetaData = pd.DataFrame(col_lbls)

orf2common = pd.read_table("data/orf2common.tsv",sep="\t")

common2orf = orf2common.set_index('GeneSymbol')

## process the Hughes DataFrame
col_lbls = {'Column':{},'Treatment1':{},'Gene1':{},'Treatment2':{},'Gene2':{},
            'Growth_Rate':{},'Time':{},'Units':{},'Name':{},'Group':{}}

with open("HughesGrowth.tsv") as fh:
    gr_l = [np.nan if ln.rstrip()=='NA' else float(ln.rstrip()) for ln in fh]

hughes_gr= np.log(2)/1.5    
with open("HughesColMetaData.tsv") as fh:
    cols = [ln.rstrip() for ln in fh]
    ## idea have a KEY,TREATMENT1,GENE1,TREATMENT2,GENE2
    new_col = []
    for ii,col in enumerate(cols):
        gr = gr_l[ii]
        gr = gr if gr == 'NA' else float(gr)*hughes_gr
        col_lbls['Growth_Rate'][ii]=gr
        col_lbls['Column'][ii]=ii
        col_lbls['Group'][ii]='Hughes'
        col_lbls['Time'][ii]=0
        col_lbls['Units'][ii]='h'
        col_lbls['Name'][ii]=col
        A = '__' in col
        B = col.endswith('_hap')
        C = col.endswith(('_1','_2'))
        D = col.endswith('_a')
        E = col.endswith('_oe')
        if col.upper() in common2orf.index: ## all of A/B/C/D/E are False
            nc = common2orf.loc[col.upper(),"SystematicName"]
            col_lbls['Gene1'][ii]=nc
            col_lbls['Treatment1'][ii]='Diploid Knockout'
            col_lbls['Gene2'][ii]='NA'
            col_lbls['Treatment2'][ii]='NA'
        elif col.upper() in orf2common.SystematicName:
            col_lbls['Gene1'][ii]=col.upper()
            col_lbls['Treatment1'][ii]='Diploid Knockout'
            col_lbls['Gene2'][ii]='NA'
            col_lbls['Treatment2'][ii]='NA'
        elif A:
            ## if '__' in col, its a double knockout
            gn1,gn2 = col.split('__')
            if B:
                colp = gn2.split('_hap')[0].upper()
                if gn1.upper() in common2orf.index:
                    nc1 = common2orf.loc[gn1.upper(),"SystematicName"]
                elif gn1.upper() in orf2common.index:
                    nc1 = gn1.upper()
                nc2 = common2orf.loc[colp,"SystematicName"]
                col_lbls['Gene1'][ii]=nc1
                col_lbls['Gene2'][ii]=nc2
                col_lbls['Treatment1'][ii]='Haploid Knockout'
                col_lbls['Treatment2'][ii]='Haploid Knockout'
            else:
                nc1 = common2orf.loc[gn1.upper(),"SystematicName"]
                nc2 = common2orf.loc[gn2.upper(),"SystematicName"]
                col_lbls['Gene1'][ii]=nc1
                col_lbls['Gene2'][ii]=nc2
                col_lbls['Treatment1'][ii]='Diploid Knockout'
                col_lbls['Treatment2'][ii]='Diploid Knockout'
        elif B:
            ## if endswith('_hap') its a haploid strain
            colp = col.split('_hap')[0].upper()
            if colp in common2orf.index:
                nc = common2orf.loc[colp,"SystematicName"]
                col_lbls['Gene1'][ii]=nc
                col_lbls['Treatment1'][ii]='Haploid Knockout'
                col_lbls['Treatment2'][ii]='NA'
                col_lbls['Gene2'][ii]='NA'
            elif colp in orf2common.index:
                ## need to look up the 
                col_lbls['Gene1'][ii]=colp
                col_lbls['Treatment1'][ii]='Haploid Knockout'
                col_lbls['Treatment2'][ii]='NA'
                col_lbls['Gene2'][ii]='NA'
        elif C:
            ## _1/_2 indicates replicate
            gn1 = col.split('_')[0]
            nc = common2orf.loc[gn1.upper(),"SystematicName"]
            col_lbls['Gene1'][ii]=nc
            col_lbls['Treatment1'][ii]='Diploid Knockout'
            col_lbls['Gene2'][ii]='NA'
            col_lbls['Treatment2'][ii]='Rep %s' % col.split('_')[1]
        elif D:
            ## _a should be -A
            gn1 = col.replace('_','-')
            nc = gn1.upper() #common2orf.loc[gn1.upper(),"SystematicName"]
            col_lbls['Gene1'][ii]=nc
            col_lbls['Treatment1'][ii]='Diploid Knockout'
            col_lbls['Gene2'][ii]='NA'
            col_lbls['Treatment2'][ii]='Rep %s' % col.split('_')[1]
        elif E:
            ## _oe should be coded as an overexpression.
            gn1 = col.split('_')[0]
            nc = common2orf.loc[gn1.upper(),"SystematicName"]
            col_lbls['Gene1'][ii]=nc
            col_lbls['Treatment1'][ii]='Overexpression'
            col_lbls['Gene2'][ii]='NA'
            col_lbls['Treatment2'][ii]='NA'

HughesMetaData = pd.DataFrame(col_lbls)
hughescorrections = pd.read_excel('data/missingHughesIds.xlsx',header=0,index_col=0)
for ii,row in hughescorrections.iterrows():
    if isinstance(row.loc['SystematicName'],float):
        HughesMetaData.Gene1.iloc[ii] = 'NA'
        HughesMetaData.Gene2.iloc[ii] = 'NA'
        HughesMetaData.Treatment1.iloc[ii] =str(HughesMetaData.Name.iloc[ii])
    else:
        HughesMetaData.Gene1.iloc[ii] = str(row.loc['SystematicName'].upper())
        HughesMetaData.Gene2.iloc[ii] = 'NA'

## to combine allMetaData with slavov_combined
##    ExpCondition --- Gene1; replace np.nan with 'NA'
##    CarbonSource ---> Treatment2
##    LimitingNutrient ---> Treatment1
##    GrowthRate ---> Growth_Rate
##    Add Time = 0 Units = h
##    Gene2 = NA
##    Group should go to "Slavov" or "Holstege" depending on whether
##    Column should just be the number of the column
slavov_combined = pd.read_pickle('data/yeast_growth_slavov-holstege_df.pkl')
slavov_combined.ExpCondition = slavov_combined.ExpCondition.replace((np.nan,'WT'),'NA')
slavov_combined.LimitingNutrient = slavov_combined.LimitingNutrient.replace(np.nan,'NA')
col_lbls = {'Column':{},'Treatment1':{},'Gene1':{},'Treatment2':{},'Gene2':{},
            'Growth_Rate':{},'Time':{},'Units':{},'Name':{},'Group':{}}

for ii,(nm,row) in enumerate(slavov_combined.iterrows()):
    grp = 'Holstege' if nm.startswith('GSM') else 'Slavov'
    col_lbls['Group'][ii]=grp
    col_lbls['Name'][ii]=nm
    col_lbls['Column'][ii]=ii
    col_lbls['Treatment1'][ii]=row.loc['LimitingNutrient']
    col_lbls['Gene1'][ii]=row.loc['ExpCondition']
    col_lbls['Treatment2']=row.loc['CarbonSource']
    col_lbls['Gene2'][ii]='NA'
    col_lbls['Time'][ii]=0
    col_lbls['Units'][ii]='h'
    col_lbls['Growth_Rate']=row.loc['GrowthRate']
SlavovHolstegeMetadata = pd.DataFrame(col_lbls)

allMetaData = pd.concat([CharlesMetaData,GreshamMetaData,HughesMetaData,SlavovHolstegeMetadata])
allMetaData = allMetaData.set_index("Name")
#allMetaData.to_excel('allmetadata.xlsx')
allMetaData.to_pickle('data/allmetadata.pkl')
## gather gene expression data
charlesGenExp = pd.read_table('CharlesGenExp.tsv',header=0,index_col=0)
greshamGenExp = pd.read_table('GreshamGenExp.tsv',header=0,index_col=0)
hughesGenExp = pd.read_table('HughesGenExp.tsv',header=0,index_col=0)
slavovGenExp = pd.read_table('ugrr_data.txt',header=0,index_col=0)
holstegeGenExp = pd.read_pickle('data/holstege_data_wt_responsive.pkl')
charlesGenExp.columns = ['Charles_%d' % ii for ii in range(charlesGenExp.shape[1])]
greshamGenExp.columns = [c.split('.')[0] for c in greshamGenExp.columns]
hughesGenExp.columns = HughesMetaData.Name
slavovGenExp.columns = SlavovMetaData.Name ### need to change!
#slavovHolstegeGenExp = pd.read_pickle('corrected_yeast_slavov-holstege_df.pkl')
common_genes = holstegeGenExp.index.intersection(
                   slavovGenExp.index.intersection(
                           hughesGenExp.index).intersection(
                                   greshamGenExp.index).intersection(
                                           charlesGenExp.index))
GEXP = pd.concat([charlesGenExp,greshamGenExp,hughesGenExp,slavovGenExp,holstegeGenExp],axis=1).loc[common_genes]
GEXP.to_pickle('all_yeast_data.pkl')
z1 = GEXP
gndata = pd.read_excel('data/holstege_gene_spreadsheet.xlsx',index_col=1)
responsive_mutants = gndata[gndata.iloc[:,7]=='responsive mutant']
metadata = pd.read_pickle('data/allmetadata.pkl')
holst_group = metadata[metadata.Group=='Holstege']
hughes_group = metadata[metadata.Group=='Hughes']

inds2rem = [i for i,row in holst_group.iterrows() if row.loc['Gene1'] not in responsive_mutants.index.tolist()]
z1p = z1.loc[:,[col for col in z1.columns if col not in inds2rem]]
z1p.to_pickle('data_nonresponsive_removed.pkl')
metadata_nr_removed = metadata.loc[z1p.columns]
metadata_nr_removed.to_pickle('metadata_nonresponsive_removed.pkl')
