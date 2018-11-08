import pandas as pa
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations as combs
from sklearn.neighbors import KNeighborsRegressor as KNR
from sklearn.cross_validation import KFold,LeaveOneOut,StratifiedKFold
from collections import defaultdict
from glob import glob
import getopt
import sys,os,os.path as osp
from scoop import futures
plt.style.use('seaborn-paper')
csrc_marker_d = {'ETH':'^','GLU':'s'}
limnut_color_d = {'ethanol':'#D62728',#,
                  'nitrogen':'#FF7F0E',
                  'phosphate':'#2CA02C',
                  'glucose':'#1F77B4',
                  'leusine':'#9467BD',
                  'uracile':'#8C564B',#,
                  'sulfur':'#E377C2',#''}
                  'heatshock':'#7F7F7F',
                  'allantoin':'#FF7F0E',
                  'arginine':'#FF7F0E',
                  'proline':'#17BECF',
                  'glutamine':'#FF7F0E',
                  'urea':'#FF7F0E',
                  'glutamate':'#FF7F0E'}
grp_shape_d = {'Charles':'^','Gresham':'s','Slavov':'o'}

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
-k, --num-neighbors=: number of neighbors for the KNN algorithm
-l, --lambda=: regularization parameter that balances the
               discretization of the data with the prediction accuracy
-N, --name=: Filename extension for the results
-R, --raw-KNN: regression is NOT weighted by inverse distance;
      If this parameter is not passed, regression will be distance-weighted.
-r, --growth-rate=: name of the growth rate file
'''

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def create_growth_rate_mapping(lim_data,Z,**kwargs):
    """builds a model which maps the gene expression to growth rate
    includes no growth points outside the data"""
    # weights are not calculated here, should be done above this level
    #lim_data = lim_data
    if 'n_neighbors' in kwargs.keys():
        n_neighbors = kwargs['n_neighbors']
    else:
        n_neighbors=8
    allgrdat = lim_data
    allgr = Z
    clf = KNR(n_neighbors,'distance')
    clf.fit(allgrdat.values,allgr.values)
    return clf

def select_most_accurate(d):
    fn_l = glob(osp.join(d,'LOO_results_cv_*.pkl'))
    opd = {}
    for fn in fn_l:
        X = pa.read_pickle(fn)
        if not 'Actual' in X.columns:
            X = X.T
        pct_error = np.abs(X.Actual.sub(X.Predicted)/X.Actual)
        opd[fn] = np.mean(pct_error)
    opdf = pa.Series(opd)
    return opdf.argmin(),opdf.min()

def mse(X):
    return np.mean(np.square(X.Actual-X.Predicted))

def success_rate(X,t=[0.25,0.05]):
    if not 'Actual' in X.columns:
        X = X.T
    pct_error = np.abs(X.Actual.sub(X.Predicted)/X.Actual)
    return dict([(tt,(pct_error<tt).sum()) for tt in t])

def logseo(Z1,gexp,useRND=False,n_neighbors=7):
    gseDF_L = []
    for gse in Z1['GEO accession'].unique():
        test = Z1['GEO accession']==gse
        test_gr = Z1[test]['Growth rate (1/h)']
        train_gr = Z1[~test]['Growth rate (1/h)']
        test_gexp = gexp[test]
        train_gexp = gexp[~test]
        if useRND:
            # update the training growth rate and gene expression
            pass
        clf_logseo= create_growth_rate_mapping(train_gexp,train_gr,n_neighbors=n_neighbors)
        pred = clf_logseo.predict(test_gexp)
        act = test_gr
        gseDF = pa.DataFrame({'Predicted':pa.Series(pred,index=test_gexp.index),
                              'Actual':test_gr})
        gseDF_L.append(gseDF)
    return pa.concat(gseDF_L,axis=0)

def loexpo(Z1,gexp,annot,grp_keys=[],n_neighbors=8):
    DF_L = []
    Z1 = Z1.loc[gexp.index]
    if len(grp_keys)>0:
        if grp_keys[0]!='Growth_Rate':
            GB = annot.groupby(grp_keys)
        else:
            GB = annot.groupby(grp_keys[0])
        for lbl,grp in GB:
            test = gexp.index.isin(grp.index.get_level_values(0))
            test_gexp = gexp[test]
            train_gexp = gexp[~test]
            test_gr = Z1[test]
            train_gr = Z1[~test]
            clf_logrpo = create_growth_rate_mapping(train_gexp,train_gr,n_neighbors=n_neighbors)
            pred = clf_logrpo.predict(test_gexp)
            act = test_gr
            DF = pa.DataFrame({'Predicted':pa.Series(pred,index=test_gexp.index),
                               'Actual':act})
            DF_L.append(DF)
    else:
        for gse in annot.unique():
            test = annot==gse
            test_gr = Z1[test]
            train_gr = Z1[~test]
            test_gexp = gexp[test]
            train_gexp = gexp[~test]
            clf_logseo = create_growth_rate_mapping(train_gexp,train_gr,n_neighbors=n_neighbors)
            pred = clf_logseo.predict(test_gexp)
            act = test_gr
            DF = pa.DataFrame({'Predicted':pa.Series(pred,index=test_gexp.index),
                                  'Actual':test_gr})
            DF_L.append(DF)
    return pa.concat(DF_L,axis=0)

def cv_kfold(Z1,gexp,nfold=5,n_neighbors=8):
    from sklearn.cross_validation import KFold
    kfold_L = []
    for _itrain,_itest in KFold(gexp.shape[0],nfold,shuffle=True):
        train_gexp = gexp.iloc[_itrain]
        train_gr = Z1.iloc[_itrain]
        test_gexp = gexp.iloc[_itest]
        test_gr = Z1.iloc[_itest]
        clf_kfold = create_growth_rate_mapping(train_gexp,train_gr,n_neighbors=n_neighbors)
        pred = clf_kfold.predict(test_gexp)
        act = test_gr
        kfoldDF = pa.DataFrame({'Predicted':pa.Series(pred,index=test_gexp.index),
                              'Actual':test_gr})
        kfold_L.append(kfoldDF)
    return pa.concat(kfold_L,axis=0)

def cv_stratified_kfold(Z1,gexp,nfold=5,nbins=0,n_neighbors=8):
    from sklearn.cross_validation import StratifiedKFold
    if nbins<2: return cv_kfold(Z1,gexp,max_d,all_zeros,nfold)
    qc,bins = pa.cut(Z1,nbins,False,retbins=True)
    stratkfold_L = []
    for _itrain,_itest in StratifiedKFold(qc,nfold,shuffle=True):
        train_gexp = gexp.iloc[_itrain]
        train_gr = Z1.iloc[_itrain]
        test_gexp = gexp.iloc[_itest]
        test_gr = Z1.iloc[_itest]
        clf_stratkfold = create_growth_rate_mapping(train_gexp,train_gr,n_neighbors=n_neighbors)
        pred = clf_stratkfold.predict(test_gexp)
        act = test_gr
        stratkfoldDF = pa.DataFrame({'Predicted':pa.Series(pred,index=test_gexp.index),
                              'Actual':test_gr})
        stratkfold_L.append(stratkfoldDF)
    return pa.concat(stratkfold_L,axis=0)

def cv_loo(Z1,gexp,n_neighbors=8):
    from sklearn.cross_validation import LeaveOneOut
    looDF_L = []
    for _itrain,_itest in LeaveOneOut(gexp.shape[0]):
        train_gexp = gexp.iloc[_itrain]
        train_gr = Z1.iloc[_itrain]
        test_gexp = gexp.iloc[_itest]
        test_gr = Z1.iloc[_itest]
        clf_loo = create_growth_rate_mapping(train_gexp,train_gr,n_neighbors=n_neighbors)
        pred = clf_loo.predict(test_gexp)
        act = test_gr
        looDF = pa.DataFrame({'Predicted':pa.Series(pred,index=test_gexp.index),
                              'Actual':test_gr})
        looDF_L.append(looDF)
    return pa.concat(looDF_L,axis=0)


def select_model(ind,data_l,FEATS,method='STRAT',n_neighbors=8):
    nfeat_pred_d = defaultdict(dict)
    ## restrict to selected features
    annot,gr,ac_gexp = data_l
    if method == 'LO_TO':
        df = loexpo(gr,ac_gexp,annot,grp_keys=['Treatment1'],n_neighbors=n_neighbors)
    elif method == 'LO_GRPO':
        df = loexpo(gr,ac_gexp,annot,['Group'],n_neighbors=n_neighbors)
    elif method == 'LO_GRO':
        df = loexpo(gr,ac_gexp,annot,['Growth_Rate'],n_neighbors=n_neighbors)
    elif method == 'LO_GSEO':
        df = loexpo(gr,ac_gexp,annot,['Growth_Rate'],n_neighbors=n_neighbors)
    elif method =='KFOLD':
        df = cv_kfold(gr,ac_gexp,5,n_neighbors=n_neighbors)
    elif method == 'STRAT':
        df = cv_stratified_kfold(gr,ac_gexp,5,5,n_neighbors=n_neighbors)
    elif method == 'LOO':
        df = cv_loo(gr,ac_gexp,n_neighbors=n_neighbors)
    return df

def bootstrap_R2(df,Nrep=100,nfold=5):
    R2_l = []
    for __ in range(100):
        for _itrain,_itest in KFold(df.shape[0],nfold,shuffle=True):
            sdf = df.iloc[_itrain]
            R2_l.append(do_ols(sdf.Actual,sdf.Predicted)[1])
    return np.mean(R2_l),np.std(R2_l)/np.sqrt(Nrep*nfold)

def do_ols(xvs,yvs):
    import statsmodels.api as sm
    ols_model_lim = sm.OLS(yvs.values,sm.add_constant(xvs.values))
    results_lim = ols_model_lim.fit()
    b,m = results_lim.params
    slope_uncert = np.sqrt(results_lim.cov_params()[1,1])
    r2 = results_lim.rsquared
    RMSE = np.sqrt(np.mean(np.square(xvs-yvs)))
    return RMSE, r2, (m,b), b/slope_uncert

def find_optimum(data_l,FEATS,method,n_neighbors=8):
    from functools import partial
    rmse_l = []; r2_l = []; res_l = []
    for _i in range(1,len(FEATS)+1):
        sel_feat = FEATS[:_i]
        gexp = data_l[-1].loc[:,sel_feat]
        DATA_l = [data_l[0],data_l[1],gexp]
        sm = partial(select_model,data_l=DATA_l,FEATS=sel_feat,method=method,n_neighbors=n_neighbors)
        df_l = list(futures.map(sm,range(20)))
        Act = df_l[0].Actual
        meanPred = np.mean(np.vstack([df.Predicted.loc[Act.index].values for df in df_l]),axis=0)
        Pred = pa.Series(meanPred,df_l[0].Actual.index)
        rmse,r2,__,__ = do_ols(Act,Pred)
        r2_l.append(r2)
        rmse_l.append(rmse)
        res_l.append(pa.DataFrame({'Actual':Act,'Predicted':Pred}))
    ii = np.argmax(r2_l)
    print method,ii,r2_l[ii]
    return r2_l,rmse_l,res_l

def scatterplot_results(df_d):
    fig,(ecax,scax) = plt.subplots(2,1,figsize=(2.25,4))
    name_l = []
    for org,(NAME,df,metadata_phe) in df_d.items():
        name_l.append(NAME)
        if not 'Actual' in df.columns:
            df = df.T
        if org=='ecoli':
            cax = ecax
            color_d = {'ATCC25404': '#7F7F7F', # GRAY
                       'BL21':'#BCBD22', # YELLOW
                       'BW13711': '#8C564B',#, BROWN
                       'BW25113': '#2CA02C',# GREEN
                       'DH5alpha':'#E377C2',# PINK
                       'EDL933':'#1F77B4',# BLUE
                       'MG1655':'#17BECF', #CYAN/LIGHT BLUE
                       'P4X':'#FF7F0E',# ORANGE
                       'P4XB2':'#FF7F0E',# ORANGE
                       'W3110':'#9467BD'} # PURPLE
            grp_shape_d = {'~M. Isalan':'o','M. Isalan':'s'}
        elif org=='saccer':
            color_d = {'ethanol':'#D62728',#, RED
                              'nitrogen':'#FF7F0E',# ORANGE
                              'phosphate':'#2CA02C',# GREEN
                              'glucose':'#1F77B4',# BLUE
                              'leusine':'#9467BD', # PURPLE
                              'uracile':'#8C564B',#, BROWN
                              'sulfur':'#E377C2',# PINK
                              'heatshock':'#7F7F7F', # GRAY
                              'allantoin':'#FF7F0E', # ORANGE
                              'arginine':'#FF7F0E', # ORANGE
                              'proline':'#17BECF', #CYAN/LIGHT BLUE
                              'glutamine':'#FF7F0E', # ORANGE
                              'urea':'#FF7F0E', # ORANGE
                              'glutamate':'#FF7F0E'} # ORANGE
            grp_shape_d = {'Charles':'^','Gresham':'s','Slavov':'o'}
            METANAME = NAME.split('-growth-only')[0]
            metadata = pa.read_pickle('metadata_%s.pkl' % METANAME) # NAME=nonresponsive_removed
            cax = scax
        xvals = np.linspace(0,df.Actual.max()*1.25,20)
        cax.plot(xvals,xvals,color='k',ls='--')
        cax.fill_between(xvals,xvals*.95,xvals*1.05,color='grey',alpha=0.5)
        ## remove the duplicate values before proceeding
        ## want to color the dots by strain
        NCOL=1
        obj_d = {}
        for nm,C in color_d.items():
            if org =='ecoli':
                metadata_sel = metadata_phe[metadata_phe.Strain==nm]
            else:
                metadata_sel = metadata_phe[(metadata_phe.Treatment1==nm)]
            sel_df = df.loc[metadata_sel.index]
            for grp,shp in grp_shape_d.items():
                if grp.startswith('~'):
                    if org=='ecoli':
                        metadata_grp = metadata_sel[metadata_sel.Author!=grp[1:]]
                else:
                    if org=='ecoli':
                        if nm=='MG1655':
                            metadata_grp = metadata_sel[metadata_sel.Author==grp]
                            kw = nm+'-rewired'
                            C = '#D62728'
                        else:
                            metadata_grp = metadata_sel
                    elif org=='saccer':
                        metadata_grp = metadata_sel[metadata_sel.Group==grp]
                        if metadata_grp.shape[0]<1: continue
                    kw = nm
                grp_df = sel_df.loc[metadata_grp.index.get_level_values(0)]
                KW = '%s-%s' % (kw,grp) if org=='saccer' else kw
                obj_d[KW] = cax.scatter(grp_df.Actual,grp_df.Predicted,c=C,marker=shp,alpha=0.7)
        cc = np.corrcoef(df.Predicted,df.Actual)[0,1]**2
        if org=='ecoli':
            cax.set_xlim(0,1.7)
            cax.set_ylim(0,1.7)
            cax.annotate(r'$R^2$ = %.3f' % (cc),(1.1,0),color='k',fontsize=12,family='Times New Roman')
        elif org=='saccer':
            cax.set_xlim(0,0.35)
            cax.set_ylim(0,0.35)
            cax.set_xlabel(r'Measured growth rate (1/h)',family='Times New Roman',size=12)
            cax.annotate(r'$R^2$ = %.3f' % (cc), (0.3,0),color='k',fontsize=12,family='Times New Roman')
        hndls = obj_d.values()
        lbls = obj_d.keys()
        leg = cax.legend(hndls,lbls,loc=0,prop={'size':8,'family':'Times New Roman'},ncol=NCOL)
        cax.set_ylabel(r'Predicted growth rate (1/h)',family='Times New Roman',size=12)
        plt.setp(cax.get_xticklabels(),family='Times New Roman',size=12)
        plt.setp(cax.get_yticklabels(),family='Times New Roman',size=12)
    fig.savefig('figure-3_%s-%s.svg' % tuple(name_l))
    return fig,(ecax,scax)

def figure_3(name_d={'saccer':'nonresponsive_removed','ecoli':'carrera-corr'},MAX_FEAT=20):
    ORG_D = {'saccer':'#2CA02C','ecoli':'#1F77B4'}
    DIR_D = {'saccer':'./','ecoli':'./'}
    EXT_D = name_d
    MAX_D = {'saccer':MAX_FEAT,'ecoli':MAX_FEAT}
    FT_D = {'GENE':'--','CORR':'-'}
    CRIT_D = {'Efficiency':'^','LSQ':'s'}#,'Reg':'o'}
    fig,ax = plt.subplots(1,1,sharex=True,figsize=(3.25,2))
    #ax,ax2 = ax_l
    fig2,ax_l = plt.subplots(2,1,sharex=True,figsize=(3.25,4))
    ax2,ax3 = ax_l #plt.subplots(1,1,sharex=True,figsize=(3.25,2))
    #ax = fig.add_subplot(111)
    L_D = {'saccer':0.001,'ecoli':0.05};
    myd = {}; xmax=np.amax(MAX_D.values())
    org_eff_d = {}
    org_lstsq_d = {}
    org_score_d = {}
    nselft_d = {}
    for O,clr in ORG_D.items():
        d = DIR_D[O]
        ext=EXT_D[O]
        max_feat = MAX_D[O]
        l = L_D[O]
        for F,ls in FT_D.items():
            if F=='GENE': continue
            op_d = {}
            eff_d = {}
            lstsq_d = {}
            # if O=='ecoli':
                # for _i in range(1,max_feat+1):
                    # lstsq = pa.read_pickle('lstsq_%d_%s--r1.pkl' % (_i,ext)))
                    # eff = pa.read_pickle('eff_%d_%s--r1.pkl' % (_i,ext)))
                    # op_d[_i]=(lstsq+l*eff).min()
                    # arg =(lstsq+l*eff).idxmin()
                    # eff_d[_i] = eff.loc[arg]
                    # lstsq_d[_i] = lstsq.loc[arg]
                    # op_ser = pa.Series(op_d)
                # eff_ser = pa.Series(eff_d)
                # lstsq_ser = pa.Series(lstsq_d)
            # elif O=='saccer':
            lstsq = pa.read_pickle('lstsq_results_%s_feat.pkl' % (ext))
            eff = pa.read_pickle('eff_results_%s_feat.pkl' % (ext))
            lstsq_ser = lstsq.min(axis=1)
            new_inds = [ii+1 for ii in lstsq_ser.index.tolist()]
            lstsq_idx = lstsq.idxmin(axis=1)
            lstsq_ser.index = new_inds
            lstsq_idx.index = new_inds
            eff_vals = []
            for ii,ft in lstsq_idx.iteritems():
                eff_vals.append(eff.loc[ii-1,ft])
            eff_ser = pa.Series(eff_vals,index=lstsq_ser.index)
            op_ser =lstsq_ser.loc[range(1,MAX_FEAT+1)] + eff_ser.loc[range(1,MAX_FEAT+1)]*l
            ax.plot(range(1,MAX_FEAT+1),op_ser/op_ser.iloc[0],ls=ls,color=clr)#,marker='o')
            nselft_d[O] = op_ser.idxmin()
        org_eff_d[O] = eff_ser
        org_lstsq_d[O] = lstsq_ser
        org_score_d[O] = op_ser
    org_eff_df =  pa.DataFrame(org_eff_d)
    org_lstsq_df = pa.DataFrame(org_lstsq_d)
    org_score_df = pa.DataFrame(org_score_d)
    ax.set_xlim(0,xmax)
    ax.set_ylim(0,1.25)
    ax.set_xlabel('Number of features',family='Times New Roman',size=12)
    ax.set_ylabel('Scaled objective',family='Times New Roman',size=12)
    plt.setp(ax.get_xticklabels(),family='Times New Roman',size=12)
    plt.setp(ax.get_yticklabels(),family='Times New Roman',size=12)
    fig.savefig('figure-2_%s-%s.svg' % tuple(EXT_D.values()))
    #fig2 = plt.figure(figsize=(3.25,2))
    #ax2 = fig2.add_subplot(111)
    for O,clr in ORG_D.items():
        #df = pl.loc[((O,'CORR')),:,:]
        eff = org_eff_df[O]
        lstsq = org_lstsq_df[O]
        for C,M in CRIT_D.items():
            if C!='LSQ':
                #if O=='ecoli':
                ax3.plot(eff.index,eff,marker=M,markersize=9,ls='',color=clr)
            else:
                ax2.plot(lstsq.index,lstsq.values/lstsq.values[0],marker=M,markersize=9,ls='',color=clr)
    ax3.set_xlabel('Number of features',family='Times new Roman',size=12)
    ax3.set_xlim(0,MAX_FEAT)
    ax2.set_ylabel('Objective',family='Times New Roman',size=12)
    ax3.set_ylabel('Objective',family='Times New Roman',size=12)
    #ax2.set_yscale('log')
    #ax2.set_ylim(0,6)
    ax3.set_ylim(.1,120)
    ax3.set_yscale('log')
    plt.setp(ax3.get_xticklabels(),family='Times New Roman',size=12)
    plt.setp(ax2.get_yticklabels(),family='Times New Roman',size=12)
    plt.setp(ax3.get_yticklabels(),family='Times New Roman',size=12)
    fig2.savefig('figure-2B_%s-%s.svg' % tuple(EXT_D.values()))
    return pa.Series(nselft_d)

def main():
    try:
        N_NEIGHBORS = 8
        FN = 'corrected_yeast_df.csv'
        argv = sys.argv
        YEAST_NAME= ''
        ECOLI_NAME= 'carrera-corr'
        YEAST_GFN = 'metadata_nonresponsive_removed.pkl' 
        ECOLI_GFN = 'data/DataSetS1_CarreraAnnotations.xlsx'
        GROWTH_ONLY = False
        try:
            opts, args = getopt.getopt(argv[1:], "hk:E:Y:e:y:g",
                                       ["help","num-neighbors=","e-coli-name=","yeast-name=",
                                        "e-coli-growth-rate=","yeast-growth-rate=",
                                        "growth-only"])
        except getopt.error, msg:
            raise Usage(msg)
        # option processing
        for option, value in opts:
            WEIGHTS = 'uniform' if option == "-u" else 'distance'
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option in ("-k", "--num-neighbors"):
                N_NEIGHBORS = int(value)
            if option in ("-E","--e-coli-name"):
                ECOLI_NAME=value
            if option in ("-Y","--yeast-name"):
                YEAST_NAME=value
            if option in ("-e","--e-coli-growth-rate"):
                ECOLI_GFN=value ## this is the file name with annotations
            if option in ("-y","--yeast-growth-rate"):
                YEAST_GFN=value ## this is the file name with annotations
            if option in ("-g","--growth-only"):
                GROWTH_ONLY = True
        if GROWTH_ONLY:
            ECOLI_NAME+='-growth-only'
            YEAST_NAME+='-growth-only'

    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
    try:
        ## get the E. coli annotation information
        if ECOLI_GFN.endswith('.pkl'):
            ecannot = pa.read_pickle(GFN)
        elif ECOLI_GFN.endswith('.csv'):
            ecannot = pa.read_table(GFN)
        elif ECOLI_GFN.endswith('.xlsx'):
            ecannot = pa.read_excel(ECOLI_GFN,sheetname=u'EcoMAC and EcoPhe',index_col='CEL file name')
            rem_rows = ecannot.index.get_level_values(0).isin(['PU00524','PU00586'])
            ecannot = ecannot[~rem_rows]
            _nz = ecannot['Flag growth']>0
            ecannot_phe = ecannot[_nz]
            ec_gr = ecannot_phe.loc[:,'Growth rate (1/h)']
        ## get the E. coli expression information
        if ECOLI_NAME.startswith('carrera'):
            #DIR = '../../E-coli_synthetic-rescues/SF1-EcoMAC/'
            ec_gncorr = pa.read_pickle('ecoli_compendium_gncorr_%s_df.pkl' % ECOLI_NAME)
        ## get the yeast annotation information
        if YEAST_GFN.endswith('.pkl'):
            scannot = pa.read_pickle(YEAST_GFN)
        elif YEAST_GFN.endswith('.csv'):
            scannot = pa.read_table(YEAST_GFN)
        ## get the yeast expression information
        if YEAST_NAME == 'nonresponsive_hughes':
            sc_gncorr = pa.read_pickle('sacCer_compendium_gncorr_nonresponsive_removed_df.pkl')
        else:
            sc_gncorr = pa.read_pickle('sacCer_compendium_gncorr_%s_df.pkl' % YEAST_NAME)
        ## for the yeast information, restrict to certain data
        if YEAST_NAME=='holstege':
            sc_gncorr,__ = convert_deletion2orf(sc_gncorr,sc_gr)
        elif YEAST_NAME=='all':
            sc_gr = scannot.Growth_Rate.astype(float)
        elif YEAST_NAME.startswith('nonresponsive_removed'):
            ## remove the holstege and hughes data
            scannot = scannot[[grp not in ('Holstege','Hughes') for grp in scannot.Group]]
            sc_gr = scannot.Growth_Rate.astype(float)
            sc_gncorr = sc_gncorr.loc[sc_gr.index]
        elif YEAST_NAME=='nonresponsive_hughes':
            ## remove the holstege data ONLY
            scannot = scannot[[grp not in ('Holstege',) for grp in scannot.Group]]
            sc_gr = scannot.Growth_Rate.astype(float)
            sc_gncorr = sc_gncorr.loc[sc_gr.index]
        elif YEAST_NAME.startswith('hughes_removed'):
            scannot =scannot[[grp not in ('Holstege', 'Hughes') for grp in scannot.Group]]
            sc_gr = scannot.Growth_Rate.astype(float)
            sc_gncorr = sc_gncorr.loc[sc_gr.index]
        else:
            sc_gr = sc_gr.GrowthRate.astype(float)
        ## eliminate any experiments with growth rate=np.nan
        if np.any(np.isnan(sc_gr)):
            sc_gr = sc_gr[~np.isnan(sc_gr)]
            sc_gncorr = sc_gncorr.loc[sc_gr.index]
            scannot = scannot.loc[sc_gr.index]
    except ValueError,e:
        print "Missing file(s). Run 'gather_yeast_data.py' first."
        raise

    if YEAST_NAME.startswith('nonresponsive_hughes'):
        yeast_bins = np.load('data/growth_bins_nonresponsive_removed.npy')
    else:
        if not GROWTH_ONLY:
            yeast_bins = np.load('data/growth_bins_%s.npy' % YEAST_NAME)
        else:
            fnext = YEAST_NAME.split('-growth-only')[0]
            yeast_bins = np.load('data/growth_bins_%s.npy' % fnext)

    ecoli_bins = np.array([0., 0.215, 0.265, 0.315, 0.365, 0.395, 0.425, 0.455,
                           0.485, 0.515, 0.545, 0.575, 0.605, 0.665, 0.825,
                           1.35, 10. ])
    #np.save('ecoli_growth_bins.npy',ecoli_bins)
    sc_data_l = [scannot,sc_gr,sc_gncorr]
    #for org,d in dir_d.items():
    scfeat_fns_l = glob('corr_feat_*-%s_feat-avg.pkl' %  YEAST_NAME)
    #scbin_fns_l = glob('corr_bins_*-%s_feat-avg.pkl' %  YEAST_NAME)
    scsel_fn =sorted(scfeat_fns_l, key = lambda f: int(f.split('-')[0].split('_')[-1]))[-1]
    scnfeat = int(scsel_fn.split('-')[0].split('_')[-1])
    scbin_fn = 'corr_feat_%d-%s_feat-avg.pkl' % (scnfeat,YEAST_NAME)
    if len(scfeat_fns_l)>1:
        print "Cleaning intermediate saved files."
        #map(os.remove,filter(lambda f: f!=scsel_fn,scfeat_fns_l))
        #map(os.remove,filter(lambda f: f!=scbin_fn,scbin_fns_l))
    YEAST_FEATS = pa.read_pickle(scsel_fn)
    ecfeat_fns_l = glob('corr_feat_*-%s_feat-avg.pkl' %  ECOLI_NAME)
    ecsel_fn= sorted(ecfeat_fns_l, key = lambda f: int(f.split('-')[0].split('_')[-1]))[-1]
    ECOLI_FEATS = pa.read_pickle(ecsel_fn)
    ec_data_l = [ecannot_phe,ec_gr,ec_gncorr.loc[ec_gr.index]]
    org_opt_method_d = {}
    org_fnext_d = {}
    for org in ['ecoli','saccer']:
        opt_method_d = defaultdict(dict)
        FEATS = YEAST_FEATS  if org=='saccer' else ECOLI_FEATS
        data_l = sc_data_l if org=='saccer' else ec_data_l
        N_NEIGHBORS= 8 if org=='saccer' else 7
        for method in ['LO_GRPO','LO_TO','LO_GRO','KFOLD']:
            if org == 'ecoli' and method != 'KFOLD': continue
            r2_l,rmse_l,res_l = find_optimum(data_l,FEATS,method,n_neighbors=N_NEIGHBORS)
            opt_method_d[method]['RMSE']=rmse_l
            opt_method_d[method]['R2']=r2_l
            opt_method_d[method]['RES']=res_l
        opt_method_df =pa.DataFrame(opt_method_d)
        NAME = YEAST_NAME if org=='saccer' else ECOLI_NAME
        opt_method_df.to_pickle('opt_method_d_%s.pkl' % NAME)
        plot_trend(opt_method_df,NAME)
        org_opt_method_d[org] = opt_method_df
        org_fnext_d[org]=NAME
    ## once the optimal parameters are found for lstsq and efficiency
    ## need to get the ec/sc results
    org_optnfeat_d = figure_3(org_fnext_d,MAX_FEAT=scnfeat)
    print org_optnfeat_d
    ec_result = org_opt_method_d['ecoli'].loc['RES','KFOLD'][org_optnfeat_d['ecoli']]
    E,F = bootstrap_R2(ec_result)
    print org_opt_method_d['ecoli'].loc['R2','KFOLD'][org_optnfeat_d['ecoli']],\
          org_opt_method_d['ecoli'].loc['RMSE','KFOLD'][org_optnfeat_d['ecoli']],E,F
    sc_result = org_opt_method_d['saccer'].loc['RES','KFOLD'][org_optnfeat_d['saccer']]
    C,D = bootstrap_R2(sc_result)
    print org_opt_method_d['saccer'].loc['R2','KFOLD'][org_optnfeat_d['saccer']],\
          org_opt_method_d['saccer'].loc['RMSE','KFOLD'][org_optnfeat_d['saccer']], C,D
    print np.sum((ec_result.Actual-ec_result.Predicted)/ec_result.Actual <0.05)
    print np.sum((sc_result.Actual-sc_result.Predicted)/sc_result.Actual <0.05)
    #SEL_ECOLI_FEATS = ECOLI_FEATS[:org_optnfeat_d['ecoli']]
    #SEL_YEAST_FEATS = YEAST_FEATS[:org_optnfeat_d['saccer']]
    #ec_result = select_model(org_optnfeat_d['ecoli'],ec_data_l,SEL_ECOLI_FEATS,'KFOLD')
    #sc_result = select_model(org_optnfeat_d['saccer'],sc_data_l,SEL_YEAST_FEATS,'KFOLD')
    ## need to get the results
    ## we can get the scatter plot
    scatterplot_results({'ecoli':(ECOLI_NAME,ec_result,ecannot_phe),
                         'saccer':(YEAST_NAME,sc_result,scannot)})
    return 0

def plot_trend(opt_method_df,NAME,kw='R2'):
    fig,ax = plt.subplots()
    my_d = {}
    for kk,vv in opt_method_df.loc[kw].iteritems():
        my_d[kk]=vv
    DF = pa.DataFrame(my_d)
    ax = DF.plot(ax=ax)
    ax.set_xlabel('Number of eigengenes')
    if kw=='R2':
        ax.set_ylabel(r'$R^2$')
    elif kw=='RMSE':
        ax.set_ylabel(r'RMSE')
    fig.savefig('figure-1_%s.svg' % NAME)
    return fig,ax

if __name__ == '__main__':
    #figure_3(MAX_FEAT=20)
    main()
