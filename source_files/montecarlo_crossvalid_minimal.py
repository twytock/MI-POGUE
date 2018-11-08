#!/usr/bin/env python
# encoding: utf-8
"""
montecarlo_crossvalid_minimal.py

This script applies MI-POGUE to a gene expression dataset.
"""
import sys,os,os.path as osp, getopt
import numpy as np
import cPickle as cp
import pandas as pa
from itertools import combinations as combs
from functools import partial
from sklearn.neighbors import KNeighborsRegressor as KNR
from sklearn.cross_validation import KFold,LeaveOneOut,StratifiedKFold
from scoop import futures
from glob import glob
from collections import defaultdict
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

## need to weight by slower-growth states
try:
    N_NEIGHBORS = 7
    LAMBDA = 5e-2
    FN = 'SF1-EcoMAC/ecoli_compendium_df.pkl'
    argv = sys.argv
    NAME= 'carrera-corr'
    GFN = 'data/DatasetS1_CarreraAnnotations.xlsx'
    org = 'ecoli'
    WEIGHTS = 'distance'
    BIOMASS_PRECURSORS= False
    GROWTH_ONLY=False
    try:
        opts, args = getopt.getopt(argv[1:], "he:k:l:N:r:o:ubg",
                                   ["help","expression-file=","num-neighbors="
                                    "lambda=","name=","growth-rate=","organism=",
                                    "uniform-weights","biomass-precursors","growth-only"])
    except getopt.error, msg:
        raise Usage(msg)
    # option processing
    for option, value in opts:
        WEIGHTS = 'uniform' if option == "-u" else 'distance'
        if option in ("-h", "--help"):
            raise Usage(help_message)
        if option in ("-e", "--expression-file"):
            FN = value
        if option in ("-k", "--num-neighbors"):
            N_NEIGHBORS = int(value)
        if option in ("-l","--lambda"):
            LAMBDA = float(value)
        if option in ("-N","--name"):
            NAME=value
        if option in ("-r","--growth-rate"):
            GFN=value
        if option in ("-u","--uniform-weights"):
            WEIGHTS='uniform'
        if option in ("-b","--biomass-precursors"):
            BIOMASS_PRECURSORS = True
        if option in ("-g","--growth-only"):
            GROWTH_ONLY = True
        if option in ("-o","--organism"):
            if value.startswith(('E','e')):
                org = 'ecoli'
            if value.startswith(('Y','y','S','s')):
                org = 'saccer'
    if GROWTH_ONLY:
        NAME+='-growth-only'

except Usage, err:
    print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
    print >> sys.stderr, "\t for help use --help"


def SAVE(obj,fn):
    with open(fn,'wb') as fh:
        cp.dump(obj,fh,-1)
    return 0

def group_membership(dat,bd_l):
    if dat.ndim<2:
        dat = dat.to_frame()
    gm = dat.copy()
    for ii,bd in enumerate(bd_l):
        qc = pa.cut(dat.iloc[:,ii],bd)
        col = qc.name
        gm.loc[qc.index,col]=qc.values
    return gm.groupby(gm.columns.tolist())

def convert_deletion2orf(gndat,growth):
    d2orf_df = pa.read_pickle('deletion2orf_df.pkl')
    d2orf_df = d2orf_df.set_index('Entry')
    commonNames = gndat.index.get_level_values(1).unique()
    orfs = d2orf_df.loc[commonNames[1:],'ORF']
    myd = dict(zip(commonNames[1:],orfs))
    myd['BY4742-WT'] = 'BY4742-WT'
    mySer = pa.Series(myd)
    l2 = mySer.lroc[gndat.index.get_level_values(1)]
    l1 = gndat.index.get_level_values(0)
    MI = pa.MultiIndex.from_tuples(zip(l1,l2))
    gndat.index =MI
    gndat_mu = gndat.groupby(level=1).mean()
    growth_mu = growth.groupby(level=0).mean()
    growth_mu = growth_mu.loc[gndat_mu.index]
    growth_mu.loc['BY4742-WT']=1
    growth_mu = growth_mu.dropna()
    gndat_mu = gndat_mu.loc[growth_mu.index]
    return gndat_mu,growth_mu

try:
    if FN.endswith('.pkl'):
        ## use the transpose to put the genes in COLUMNS and the experimental conditions in ROWS
        sc_gndat = pa.read_pickle(FN).T
    elif FN.endswith('.csv'):
        sc_gndat = pa.read_table(FN).T
    if GFN.endswith('.pkl'):
        sc_gr = pa.read_pickle(GFN)
    elif GFN.endswith('.csv'):
        sc_gr = pa.read_table(GFN)
    elif GFN.endswith('.xlsx'):
        if org=='ecoli' and GFN.startswith('data/DataSetS1'):
            sc_gr = pa.read_excel('data/DataSetS1_CarreraAnnotations.xlsx',
                        sheetname='EcoMAC and EcoPhe',index_col=u'CEL file name')
        else:
            sc_gr = pa.read_excel(GFN)
    if org == 'saccer':
        if NAME== 'nonresponsive_hughes' and not BIOMASS_PRECURSORS:
            sc_gncorr = pa.read_pickle('sacCer_compendium_gncorr_nonresponsive_removed_df.pkl')
            evecs = pa.read_pickle('yeast_evecs_nonresponsive_removed_df.pkl')
        elif BIOMASS_PRECURSORS:
            sc_gncorr = sc_gndat
        else:
            sc_gncorr = pa.read_pickle('sacCer_compendium_gncorr_%s_df.pkl' % NAME)
            evecs = pa.read_pickle('yeast_evecs_%s_df.pkl' % NAME)
    elif org=='ecoli':
        if NAME == 'carrera-corr' and not BIOMASS_PRECURSORS:
            sc_gncorr = pa.read_pickle('ecoli_compendium_gncorr_%s_df.pkl' % NAME)
            evecs = pa.read_pickle('ecoli_evecs_%s_df.pkl' % NAME)
        elif BIOMASS_PRECURSORS:
            sc_gncorr = sc_gndat
        else:
            sc_gncorr = pa.read_pickle('ecoli_compendium_gncorr_%s_df.pkl' % NAME)
            evecs = pa.read_pickle('ecoli_evecs_%s_df.pkl' % NAME)

    ## need to put in function to calculate correlations here.
    if org=='saccer':
        if NAME.startswith('holstege'):
            sc_gndat,sc_gr = convert_deletion2orf(sc_gndat,sc_gr)
            sc_gncorr,__ = convert_deletion2orf(sc_gncorr,sc_gr)
        elif NAME.startswith('all'):
            sc_gr = sc_gr.Growth_Rate.astype(float)
        elif NAME.startswith('nonresponsive_removed') and not BIOMASS_PRECURSORS:
            sc_gr = sc_gr[[grp not in ('Holstege','Hughes') for grp in sc_gr.Group]] ## remove the holstege and hughes data
            sc_gr = sc_gr.Growth_Rate.astype(float)
            sc_gndat = sc_gndat.loc[sc_gr.index]
            sc_gncorr = sc_gncorr.loc[sc_gr.index]
        elif NAME.startswith('nonresponsive_removed') and BIOMASS_PRECURSORS:
            sc_gr = sc_gr[[grp not in ('Holstege','Hughes') for grp in sc_gr.Group]] ## remove the holstege and hughes data
            sc_gr = sc_gr.Growth_Rate.astype(float)
            sc_gndat = sc_gndat.loc[sc_gr.index]
            sc_gncorr = sc_gndat
        elif NAME.startswith('nonresponsive_hughes'):
            sc_gr = sc_gr[[grp not in ('Holstege',) for grp in sc_gr.Group]] ## remove the holstege data ONLY
            sc_gr = sc_gr.Growth_Rate.astype(float)
            sc_gndat = sc_gndat.loc[sc_gr.index]
            sc_gncorr = sc_gncorr.loc[sc_gr.index]
        elif NAME.startswith('hughes_removed'):
            sc_gr = sc_gr[[grp not in ('Holstege',) for grp in sc_gr.Group]] ## remove the holstege data ONLY
            sc_gr = sc_gr.Growth_Rate.astype(float)
            sc_gndat = sc_gndat.loc[sc_gr.index]
            sc_gncorr = sc_gncorr.loc[sc_gr.index]
        else:
            sc_gr = sc_gr.Growth_Rate.astype(float)
    elif org=='ecoli':
        if NAME.startswith('carrera-corr') and BIOMASS_PRECURSORS:
            sc_gr = sc_gr[sc_gr['Flag growth']==1]
            sc_gr = sc_gr['Growth rate (1/h)'].astype(float)
            sc_gncorr = sc_gncorr.loc[sc_gr.index]
        elif NAME.startswith('carrera-corr'):
            sc_gr = sc_gr[sc_gr['Flag growth']==1]
            sc_gr = sc_gr['Growth rate (1/h)'].astype(float)
            sc_gncorr = sc_gncorr.loc[sc_gr.index]
    if np.any(np.isnan(sc_gr)):
        sc_gr = sc_gr[~np.isnan(sc_gr)]
        sc_gndat = sc_gndat.loc[sc_gr.index][~np.isnan(sc_gr)]
        sc_gncorr = sc_gncorr.loc[sc_gr.index][~np.isnan(sc_gr)]
    
except ValueError,e:
    print "Missing file(s). Run 'gather_yeast_data.py' first."
    raise

def realistic_noise_data(cols,gr_err,nfold,isCorr=True):
    new_inds = []
    sc_gr_ = sc_gr#.loc[gr_cols,u'Growth rate']
    sc_gndat_ = sc_gndat#.loc[gexp_cols]
    L,M = sc_gndat_.shape
    new_gr = np.empty(L*nfold)
    new_noise=np.empty((L*nfold,M))
    STD = sc_gndat.std()
    for ii in range(nfold):
        new_gr[ii*L:(ii+1)*L] = sc_gr_.loc[sc_gndat_.index] + \
                                gr_err*np.random.randn(L)*sc_gr_.loc[sc_gndat_.index]
        new_inds.extend([gn+'--%d' % ii for gn in sc_gndat_.index.tolist()])
        new_noise[ii*L:(ii+1)*L,:]=STD.values*np.random.randn(L,M)
    rnd_gr = pa.Series(new_gr,index=new_inds)
    if isCorr:
        AV = sc_gndat_.dot(evecs.loc[:,cols])
    else:
        if sum([ft in sc_gndat_.columns.tolist() for ft in cols]) > 0:
            AV = sc_gndat_.loc[:,cols]
        else:
            cols = [col for col in sc_gndat_.columns if col[0] in cols]
            AV = sc_gndat_.loc[:,cols]
    Act_vals = pa.concat([AV for nf in range(nfold)],axis=0);
    Act_vals.index = new_inds
    if isCorr:
        try:
            var = pa.DataFrame(np.dot(new_noise,evecs.loc[:,cols].values),index=new_inds,columns=cols)
        except ValueError:
            print cols
            raise
    else:
        var = pa.DataFrame(new_noise,index=new_inds,columns=sc_gndat_.columns).loc[:,cols]
        Act_vals = Act_vals.loc[:,cols]; var = var.loc[:,cols]
    rnd_gexp = Act_vals + var
    return rnd_gr,rnd_gexp

#weights should be used based on distance away
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

def fixed_bin_divisions(X,nbins=10):
    sp = X.shape[0]//(nbins-1)
    sX = sorted(X)
    candidate_bins = np.hstack([sX[::sp],sX[-1]])
    candidate_bins = np.unique(candidate_bins) # filter duplicates
    c_min = candidate_bins[0] * .9 if candidate_bins[0]>0 \
        else 1.1 * candidate_bins[0]
    c_max = candidate_bins[-1] * 1.1 if candidate_bins[-1]>0 \
        else .9 * candidate_bins[-1]
    candidate_mps = (candidate_bins[:-1] + candidate_bins[1:])/2.
    cw = candidate_mps[1:]-candidate_mps[:-1]
    cm = (candidate_mps[1:]+candidate_mps[:-1])/2.
    inv_bins = cw < 0.05 * cm
    while np.any(inv_bins):
        ri = np.random.randint(0,inv_bins.sum())
        rb = np.where(inv_bins)[0][ri]
        Ns,__=np.histogram(sX,bins=np.hstack([c_min,candidate_mps,c_max]))
        ## rb + 1 is the bin corresponding to rb
        if rb == 0:
            TF = True
        elif rb == inv_bins.shape[0]-1:
            TF = False
        elif Ns[rb]==Ns[rb+2]:
            rn = np.random.uniform()
            TF = True if rn > 0.5 else False
        else:
            TF = np.argmax((Ns[rb],Ns[rb+2]))
        candidate_mps,cw,cm = update_mps(candidate_mps,rb,TF)
        inv_bins = cw < 0.05 * cm
    return pa.cut(X,bins=np.hstack([c_min,candidate_mps,c_max]),retbins=True)

def update_mps(old,ind,TF):
    if TF:
        new = np.hstack([old[:ind],(old[ind]+old[ind+1])/2.,old[ind+2:]])
    else:
        new = np.hstack([old[:ind-1],(old[ind]+old[ind-1])/2.,old[ind+1:]])
    nw = new[1:]-new[:-1]
    nm = (new[1:]+new[:-1])/2.
    return new,nw,nm

def bin_divisions(X,nbins=10):
    if isinstance(nbins,np.ndarray):
        qc,bins = pa.cut(X,nbins,False,retbins=True)
    else:
        #qc,bins = pa.qcut(X,nbins,retbins=True)
        qc,bins = fixed_bin_divisions(X,nbins)
    dat= pa.DataFrame(np.vstack([X,qc]).T,index=X.index,columns=[X.name,'bin_no'])
    return dat,bins

def initialize_model(filt,corr=True,feats=[],**kwargs):
    ## need to update this!
    ##
    ## in both cases need to load this and sort it
    sc_gexp_phe = sc_gncorr#.iloc[gexp_cols,:]
    sc_gr2 = sc_gr#.loc[:,'Growth_Rate']
    if len(feats)>0:
        if sum([ft in sc_gexp_phe.columns.tolist() for ft in feats]) > 0:
            sc_gexp_phe = sc_gexp_phe.loc[:,feats]
        else:
            cols = [col for col in sc_gexp_phe.columns if col[0] in feats]
            sc_gexp_phe = sc_gexp_phe.loc[:,cols]

        #sc_gncorr = sc_gncorr.loc[:,feats]
    zscale = 1 if not 'zscale' in kwargs.keys() else kwargs['zscale']
    n_neighbors=8 if not 'n_neighbors' in kwargs.keys() else kwargs['n_neighbors']
    NUM = -1 if not 'NUM' in kwargs.keys() else kwargs['NUM']
    gr_err = 0.1 if not 'gr_err' in kwargs.keys() else kwargs['gr_err']
    noiseFolds = 4 if not 'noiseFolds' in kwargs.keys() else kwargs['noiseFolds']

    ## if looking for the growth data, then report only the data that
    ## has a corresponding measured growth rate.
    rnd_gr,rnd_phe = realistic_noise_data(feats,gr_err,noiseFolds,corr)
    return rnd_gr, rnd_phe

def mc_data(gr,dat,gn_err,gr_err,nfold):
    new_inds = []
    L = dat.shape[0]
    new_gr = np.empty(L*nfold)
    if dat.ndim>1:
        cols = dat.columns.tolist()
        new_dat=np.empty((L*nfold,dat.shape[1]))
        M = dat.shape[1]
    else:
        M = 1
        cols = dat.index.tolist()
        new_dat=np.empty((L*nfold,1))
    gr_ = gr.loc[dat.index]
    gr_.index = dat.index
    for ii in range(nfold):
        new_gr[ii*L:(ii+1)*L] = gr_ + gr_err*np.random.randn(L)*gr_
        new_dat[ii*L:(ii+1)*L,:] = dat + gn_err*np.random.randn(L,M)*dat
        new_inds.extend([(gsm+'_{0:d}'.format(ii),DEL) for gsm,DEL in dat.index.tolist()])
    rnd_gr = pa.Series(new_gr,index=pa.MultiIndex.from_tuples(new_inds))
    rnd_gr.index.names = ['ORF']
    rnd_gexp =pa.DataFrame(new_dat,index=pa.MultiIndex.from_tuples(new_inds),
                            columns=dat.columns)
    rnd_gexp.index.names = ['ORF']
    return gr_,rnd_gr,rnd_gexp

def do_strat(nfold,**kwargs):
    """
do_strat
estimates the log-state-space occupancy of the dataset and performs stratified k-fold CV
Arguments:
    nfold -- number of randomized profiles for every experiment
    kwargs -- a dictionary of keyword arguments, including:
        z_gr -- pandas series of growth rate, rows=experiments
        rnd_gr -- pandas series of randomized growth rates
        z_phe -- pandas dataframe of gene expression, columns=features, row=experiments
                 the experiments match z_gr
        rnd_phe -- pandas dataframe of randomized gene expression, columns/rows as in z_phe
        all_bins -- numpy array of the growth rate bins used for stratifying the data
Returns:
    pa.DataFrame(OP) -- DataFrame containing the cross-validation results;
                        Rows=experiments; Columns = Actual/Noiseless/Downsampled/Predicted
                        Actual=measured growth rates (h^-1)
                        Noiseless=predictions with the pseudodata excluded from training data
                        Downsampled=predictions with overrepresented growth rates downsampled in the training data
                        Predicted=predictions with no modifications
    fbins            -- list of bins describing the discretization along each feature
    efficiency       -- log state-space occupancy measure
"""
    #print "Begin Stratified KFold"
    z_gr = kwargs['z_gr']
    rnd_gr = kwargs['rnd_gr']
    z_phe = kwargs['z_phe']
    rnd_phe = kwargs['rnd_phe']
    all_bins=kwargs['all_bins']
    N = len(all_bins)
    if N > 1:
        Y = pa.concat([z_gr,z_phe.iloc[:,:N-1]],axis=1)
    else:
        Y = z_gr.to_frame()
    cat,grbins = bin_divisions(z_gr,all_bins[0])
    mem = group_membership(Y,all_bins) ## all_bins contains the bins of all fixed features
    dat,fbins = bin_divisions(z_phe.iloc[:,-1])
    dat.index.name = ''
    fbins[0]=-1e3; fbins[-1]=1e3
    fdivs = len(dat.bin_no.unique())
    if N > 1:
        path_probs = pa.DataFrame(np.zeros((mem.ngroups,fdivs)),
                                  index=pa.MultiIndex.from_tuples(mem.groups.keys()),
                                  columns=dat.bin_no.unique())
    else:
        path_probs = pa.DataFrame(np.zeros((mem.ngroups,fdivs)),
                                  index=mem.groups.keys(),columns=dat.bin_no.unique())
    for gid,gdf in mem:
        if gdf.shape[0]<2:
            try:
                xxx = path_probs.index.get_loc(gid)
                yyy = path_probs.columns.get_loc(dat.loc[gdf.index[0],'bin_no'])
                path_probs.iloc[xxx,yyy]=1
            except KeyError:
                #print '>',z_phe.columns
                #print '>',gid
                #print '>',dat.index
                #print '>',gdf
                #import pdb
                #pdb.set_trace()
                raise
        else:
            ginds = gdf.index
            #=dat.index.get_level_values(level="GSM").isin(ginds.index.get_level_values(level="GSM"))
            sel=dat[dat.index.isin(ginds)]
            datcount = sel.bin_no.value_counts()
            for xi,xv in datcount.iteritems():
                path_probs.loc[gid,xi] = xv
    Nocc = np.sum(path_probs>0).sum()
    NREP=100; REP_L = []
    Nposs = reduce(lambda x,y:x*y,map(lambda x: len(x)-1,all_bins+[fbins]))
    efficiency = np.abs(np.log(Nocc**2/float(Nposs)))
    SC_GROWTH = np.log(2)/1.5
    DS_BINS = np.linspace(0.98,1.02,5)*SC_GROWTH
    for _nrep in range(NREP):
        strat_df_l = []
        for train_ind,test_ind in StratifiedKFold(cat.bin_no,nfold,shuffle=True):
            rnd_ind_l = []
            for _ind in z_phe.iloc[train_ind].index:
                rnd_ind_l.extend([_ind+'--%d' % _ii for _ii in range(nfold)])
            train_ind_ds = downsample_inds(train_ind,z_gr,DS_BINS,max_samp=50)
            test_ind_ds = downsample_inds(test_ind,z_gr,DS_BINS,max_samp=50)
            train_rows = pa.concat([z_phe.iloc[train_ind],rnd_phe.loc[rnd_ind_l]],axis=0)
            test_rows = z_phe.iloc[test_ind]
            train_grs = pa.concat([z_gr.iloc[train_ind],rnd_gr.loc[rnd_ind_l]],axis=0)
            test_grs = z_gr.iloc[test_ind]
            clf = create_growth_rate_mapping(train_rows,train_grs)
            CLF = create_growth_rate_mapping(z_phe.iloc[train_ind],z_gr.iloc[train_ind])
            CLF_DS = create_growth_rate_mapping(z_phe.iloc[train_ind_ds],z_gr.iloc[train_ind_ds])
            pred = clf.predict(test_rows)
            no_noise_pred = CLF.predict(test_rows)
            downsamp_pred = CLF_DS.predict(z_phe.iloc[test_ind_ds])
            arr_id_l = test_rows.index.tolist()
            pred_ser = pa.Series(pred,index=arr_id_l)
            nsls_ser = pa.Series(no_noise_pred,index=arr_id_l)
            ds_ser = pa.Series(downsamp_pred,index=z_phe.iloc[test_ind_ds].index.tolist())
            act = test_grs
            df= pa.DataFrame({'REP':_nrep,'ORF':arr_id_l,
                              'Predicted':pred_ser.loc[arr_id_l],
                              'Actual':act,
                              'Noiseless':nsls_ser.loc[arr_id_l],
                              'Downsampled':ds_ser.loc[arr_id_l]})
            df = df.set_index(['REP','ORF'])
            strat_df_l.append(df)
        DF = pa.concat(strat_df_l,axis=0)
        REP_L.append(DF)
    DDF = pa.concat(REP_L,axis=0)
    OP = {}
    OP['Predicted'] = pa.concat([DDF.xs(_i,level='REP').Predicted for _i in range(NREP)],
                                axis=1,sort=True).mean(axis=1)
    OP['Actual'] = DDF.xs(0,level='REP').Actual
    OP['Noiseless'] = pa.concat([DDF.xs(_i,level='REP').Noiseless for _i in range(NREP)],
                                axis=1,sort=True).mean(axis=1)
    OP['Downsampled']=pa.concat([DDF.xs(_i,level='REP').Downsampled for _i in range(NREP)],
                                  axis=1,sort=True).mean(axis=1)
    return pa.DataFrame(OP),fbins,efficiency


def cv_mc(feats,corr=True,gn_err=0.10,gr_err=0.05,\
              zscale=1,n_neighbors=7,name=[],all_bins=[],recalc=False,\
              useNoise=True,noiseFolds=4,ext=''):
    iinit_d = {'n_neighbors':n_neighbors,'zscale':zscale,\
                   'gn_err':gn_err,'gr_err':gr_err,'useNoise':useNoise,\
                   'noiseFolds':noiseFolds}
    rnd_gr,rnd_phe = initialize_model('growth',corr,feats,**iinit_d)
    if sum([ft in sc_gncorr.columns.tolist() for ft in feats]) > 0:
        sc_gnexp_phe = sc_gncorr.loc[:,feats]
    else:
        cols = [col for col in sc_gncorr.columns if col[0] in feats]
        sc_gnexp_phe = sc_gncorr.loc[:,cols]
    if len(feats)<1:
        print "no feature passed."
        return
    input_d = {'z_gr':sc_gr,'z_phe':sc_gnexp_phe,'all_bins':all_bins,
               'rnd_gr':rnd_gr,'rnd_phe':rnd_phe}
    #subdir = 'cv_stats_rep'
    n1 = len(feats)
    if all_bins == []:
        if org=='ecoli':
            yeast_bins=np.load('data/growth_bins_carrera-corr.npy') 
        else:
            yeast_bins=np.load('data/growth_bins_nonresponsive_removed.npy') 
        yeast_bins[0] = 0; yeast_bins[-1] = 10
        input_d['all_bins']= [yeast_bins]
        df,__,eff = do_strat(noiseFolds,**input_d)
        df.to_pickle('%s_results_%d.pkl' % (name,n1))
        print "%s precursors of biomass accuracy:" % name
        print df.corr().iloc[0,1:]**2
    else:
        df,bins,eff = do_strat(noiseFolds,**input_d)
        sel_feat = feats[-1]
        #df.to_pickle(osp.join(subdir,'%s_%d_%s.pkl' % (sel_feat,n1,ext)))
        D = np.sqrt(np.sum(np.square(df.Actual-df.Predicted)))
        return sel_feat,df,bins,eff


def sel_pairs_predictions(feat_l,bins_l,kind='feat',recalc=False,err=.1,useNoise=True,
                          noiseFolds=4,ext='',filt_thr=1e-6,testAll=True,n_neighbors=7):
    l = evecs.columns.tolist()#pa.read_pickle('sacCer_%s_l.pkl' % kind)
    L = list(set(l) - set(feat_l))
    if not testAll:
        L = filter(lambda x: float(x) > filt_thr,L)
    test_genes = sorted(L,key=lambda l: float(l), reverse=True)
    LOL = [feat_l + [g] for g in sorted(test_genes,key=lambda x: float(x),reverse=True)]
    NN = True if kind.startswith('feat') else False
    p_cvp = partial(cv_mc,corr=NN,all_bins=bins_l,recalc=recalc,gn_err=err,useNoise=True,noiseFolds=4,ext=ext,n_neighbors=n_neighbors)
    L = list(futures.map(p_cvp,LOL))
    return L

def downsample_inds(inds,act,bins,max_samp=50):
    ds_inds = []
    binno_ind_d = defaultdict(list)
    for ind,gr in zip(inds,act.iloc[inds]):
        if gr<np.amin(bins):
            ds_inds.append(ind)
        elif gr>=np.amax(bins):
            ds_inds.append(ind)
        else:
            for no,(start,end) in enumerate(zip(bins[:-1],bins[1:])):
                if (start<=gr<end):
                    binno_ind_d[no].append(ind)
                    break
    for no,inds in dict(binno_ind_d).items():
        if len(inds)> max_samp:
            ds_inds.extend(np.random.permutation(inds)[:max_samp])
        else:
            ds_inds.extend(inds)
    return ds_inds

def compare_columns(ii,df,dnsamp_bins,column):
    sel_inds = downsample_inds(np.arange(df.shape[0]),df.Actual,dnsamp_bins)
    df = df.iloc[sel_inds]
    return np.linalg.norm(df.Actual-df.loc[:,column])

def best_feat(L,ii,name='',p=0.1):
    print "Finding best feature."
    eff_d = {}
    lstsq_d = {}
    dummy = 100
    SC_GROWTH = np.log(2)/1.5
    DS_BINS = np.linspace(0.98,1.02,5)*SC_GROWTH
    for feat,df,bins,eff in L:
        if 'noise' in name:
            col = 'Predicted'
        elif name.startswith('slavov-holstege'):
            col = 'Downsampled'
        else:
            col = 'Noiseless'
        norm = partial(compare_columns,df=df,dnsamp_bins=DS_BINS,column=col)
        D_l = list(futures.map(norm,range(dummy)))
        Dmu = np.mean(D_l)
        #score = D+eff*p
        eff_d[feat]=eff
        lstsq_d[feat]=Dmu
    eff_ser = pa.Series(eff_d)
    eff_ser.to_pickle('eff_%d_%s.pkl' % (ii+1, name))
    lstsq_ser = pa.Series(lstsq_d)
    lstsq_ser.to_pickle('lstsq_%d_%s.pkl' % (ii+1, name))
    score_ser = lstsq_ser + p*eff_ser
    sel = score_ser.idxmin()
    print "The minimum feature at",ii+1,"features is:",sel ,lstsq_ser.loc[sel], lstsq_ser.idxmin(),eff_ser.loc[sel]
    ll = filter(lambda l: l[0]==sel,L)
    return ll[0]#L[sel]

def select_features(ext='',F=16384.,K=7,P=0.1,kind='feat',start=0,nrep=1,err=.1,MAX_FEATS=50):
    all_feats = evecs.columns.tolist()
    if not osp.exists('tmp'):
        os.mkdir('tmp')
    if NAME.startswith('nonresponsive_hughes'):
        yeast_bins = np.load('data/growth_bins_all.npy' )
    else:
        if GROWTH_ONLY:
            fnext = ext.rsplit('_',1)[0].split('-growth-only')[0]
            yeast_bins = np.load('data/growth_bins_%s.npy' % fnext)
        else:
            yeast_bins = np.load('data/growth_bins_%s.npy' % ext.rsplit('_',1)[0])
    yeast_bins[0] = 0; yeast_bins[-1] = 10
    if start > 0:
        all_bins = pa.read_pickle('tmp/corr_bins_%d-%s-avg.pkl' % (start,ext))
        sel_feat= pa.read_pickle('tmp/corr_feat_%d-%s-avg.pkl' % (start,ext))
    else:
        all_bins = [yeast_bins]
        sel_feat = []
    for ii in range(start,min(MAX_FEATS,len(all_feats))):
        #print ii
        for jj in range(nrep):
            #print jj
            testAll = True if jj==0 else False
            if not testAll:
                _fn = './tmp/lstsq_%d_%s--r1.pkl' % (ii+1, ext)
                _X = pa.read_pickle(_fn)
                _Y = pa.read_pickle(_fn.replace('lstsq','eff'))
                _SER = _X + P*_Y
                thr = _SER.min() + (_SER.max()-_SER.min()) * .5
                THR =np.amin([np.amin(np.array(map(float,_SER[_SER<thr].index.tolist())))-1e-6,1])
            else: THR = 1e-6
            L = sel_pairs_predictions(sel_feat,all_bins,kind='%s_%s' % (kind,NAME),
                                      noiseFolds=4,err=err,ext=ext+'-r%d' % (jj+1),
                                      testAll=testAll,filt_thr=THR,n_neighbors=K)
            feat,data,fbins,eff = best_feat(L,ii,ext+'--r%d' % (jj+1),P)
            predictions_df = filter(lambda elt: elt[0]==feat,L)[0][1]
            predictions_df.to_pickle('%s_%d_%s--r%d.pkl' % (feat,ii+1,ext,jj+1))
            #rep_feat2 = sel_feat + [feat]
            #rep_bins2 = all_bins + [fbins]
            #SAVE(rep_feat2,'./tmp/corr_feat_%d-%s-r%d.pkl' % (ii+1,ext,jj+1))
            #SAVE(rep_bins2,'./tmp/corr-bins_%d-%s-r%d.pkl' % (ii+1,ext,jj+1))
            pfn = 'corr_feat_%d-%s-r%d.pkl' % (ii-1,ext,jj-1)
            bfn = 'corr_bins_%d-%s-r%d.pkl' % (ii-1,ext,jj-1)
            if osp.exists(pfn):
                os.remove(pfn)
            if osp.exists(bfn):
                os.remove(bfn)
            #L=glob(osp.join('cv_stats_rep/*_%d_%s-r%d.pkl' % (ii+1,ext,jj+1)))
            #if not osp.exists('cv_stats_rep/%s_%d_r%d_res.pkl' % (ext,ii+1,jj+1)):
            #PLI = pa.Panel(dict([(osp.split(l)[1].split('_')[0],
            #                      pa.read_pickle(l)) for l in L]))
            #PLI.to_pickle('cv_stats_rep/%s_%d_r%d_res.pkl' % (ext,ii+1,jj+1))
            #for fn in L:
            #    if osp.exists(fn):
            #        os.remove(fn)
        #sel_feat.append(feat)
        vc_feat,new_bins = find_feature_and_bins(ii+1,ext,P)
        all_bins = all_bins + [new_bins]
        SAVE(all_bins,'corr_bins_%d-%s-avg.pkl' % (ii+1,ext))
        sel_feat.append(vc_feat)
        SAVE(sel_feat,'corr_feat_%d-%s-avg.pkl' % (ii+1,ext))
    for wd in ['feat','bins']:
        fn_l = glob('corr_%s_*-%s-avg.pkl' % (wd,ext))
        __ = map(os.remove,sorted(fn_l,key=lambda fn: int(fn.split('-')[0].split('_')[-1]))[:-1])
    for wd in ['lstsq','eff','sel_feat']:
        X_d = {}
        for i in range(min(MAX_FEATS,len(all_feats))):
            WD = wd if wd != 'sel_feat' else sel_feat[i]
            fn = '%s_%d_%s--r1.pkl' % (WD,i+1,ext)
            X = pa.read_pickle(fn)
            if WD != 'sel_feat':
                X_d[i]= X
            else:
                X_d[(i,WD)]=X
            os.remove(fn)
        if wd != 'sel_feat':
            df = pa.DataFrame(X_d).T
            df.to_csv('%s_results_%s.csv' % (wd,ext))
        else:
            df = pa.Panel(X_d)
        df.to_pickle('%s_results_%s.pkl' % (wd,ext))
    return 0

def find_feature_and_bins(nfeat,ext,p=0.1):
    data_fn_l = glob('lstsq_%d_%s--r1.pkl' % (nfeat,ext))
    feat_d = {}
    for data_fn in data_fn_l:
        eff_fn = data_fn.replace('lstsq','eff')
        _zz = int(osp.split(data_fn)[1].split('--r')[1].split('.')[0])
        eff = pa.read_pickle(eff_fn)
        delta = pa.read_pickle(data_fn)
        feat_d[_zz] = delta + p*eff
    feat_df = pa.DataFrame(feat_d)
    srt = feat_df.rank(axis=0).mean(axis=1).sort_values()
    feat_VC = srt.index[0]
    dat,fbins = bin_divisions(sc_gncorr.loc[:,feat_VC])
    fbins[0] = -1e3; fbins[-1] = 1e3
    print "Feature %d, Selected: %s" % (nfeat,feat_VC)
    return feat_VC,fbins

def biomass_precursors(kind='yeast7',org='saccer'):
    """
biomass_precursors
runs MI-POGUE with the selected features as biomass precursors
Arguments:
    kind -- identifier of the metabolic model from which the precursors are taken
    org  -- organism on which MI-POGUE is being applied
Returns:
    None (cross validation output is saved by cv_mc)
"""
    if org=='saccer':
        sel_feat = np.load('data/%s_precursors_sysname.npy' % kind)
        cv_mc(sel_feat,False,n_neighbors=8,name='%s' % kind)
    elif org=='ecoli':
        sel_feat = pa.read_pickle('data/%s_precursors_l.pkl' % kind)
        cv_mc(sel_feat,False,n_neighbors=7,name='%s' % kind)


if __name__ == '__main__':
    if BIOMASS_PRECURSORS:
        if org=='saccer':
            KIND='iFF708'
            biomass_precursors(org=org)
            biomass_precursors(KIND)
        else:
            biomass_precursors(kind='iJO1366',org=org)
    else:
        KIND='feat'
        select_features('%s_%s' % (NAME,KIND),K=N_NEIGHBORS,P=LAMBDA,kind=KIND,MAX_FEATS=20)
    sys.exit()
