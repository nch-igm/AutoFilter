#!/usr/bin/env python
# coding: utf-8

# In[2]:


## LOAD MODULES
import sys
import pandas as pd
import numpy as np   
import imp
import vcf 
import math
import pickle 
import xgboost as xgb 
import f  # (local script bank)
import random
try:
    import shap
    featimpflag=1
except:
    print("shap python package not found...will not compute feature importances.")
    featimpflag=0
#loc_cmd='python /igm/home/jbg001/git/germfilt/germline_variant_filtering/kf/autofilt_train.py --traindir /igm/projects/FOR_JEFF/germfilt/train/fulltrain_030920_afternoon_batch_3_combo_1_prepy_0 --mod_paramdict /igm/projects/FOR_JEFF/germfilt/train/fulltrain_030920_afternoon_batch_3_combo_1_prepy_0/mod_params.json --vcf_query_preprocessed /igm/projects/FOR_JEFF/germfilt/tm/03062020_noon_noprepy_batch_3_giab_3/query_preprocessed_ph.vcf.gz --vcf_happy /igm/projects/FOR_JEFF/germfilt/tm/03062020_noon_noprepy_batch_3_giab_3/happyout_ph.vcf.gz'
#loc_cmd='python /igm/home/jbg001/git/germfilt/germline_variant_filtering/kf/autofilt_train.py --traindir /igm/projects/FOR_JEFF/germfilt/train/fulltrain_030920_afternoon_batch_1_combo_2_prepy_0 --mod_paramdict /igm/projects/FOR_JEFF/germfilt/train/fulltrain_030920_afternoon_batch_1_combo_2_prepy_0/mod_params.json --vcf_query_preprocessed /igm/projects/FOR_JEFF/germfilt/tm/03062020_noon_noprepy_batch_1_giab_3/query_preprocessed_ph.vcf.gz --vcf_happy /igm/projects/FOR_JEFF/germfilt/tm/03062020_noon_noprepy_batch_1_giab_3/happyout_ph.vcf.gz'
#loc_cmd='python /igm/home/jbg001/git/germfilt/germline_variant_filtering/kf/autofilt_train.py --traindir /igm/projects/FOR_JEFF/germfilt/train/fulltrain_030920_afternoon_batch_2_combo_1_prepy_0 --mod_paramdict /igm/projects/FOR_JEFF/germfilt/train/fulltrain_030920_afternoon_batch_2_combo_1_prepy_0/mod_params.json --vcf_query_preprocessed /igm/projects/FOR_JEFF/germfilt/tm/03062020_noon_noprepy_batch_2_giab_1/query_preprocessed_ph.vcf.gz --vcf_happy /igm/projects/FOR_JEFF/germfilt/tm/03062020_noon_noprepy_batch_2_giab_1/happyout_ph.vcf.gz'
filestem='run_train_5_1_0'
testargs=f.grab_bashscript_args(filestem)


# In[ ]:





# In[3]:


## READ INPUT PARAMS, DEFINE OUTFILES
if 1==1:
    myargs=sys.argv
    myargs=testargs
    querylocs=[x+1 for x in range(len(myargs)) if myargs[x]=='--vcf_query_preprocessed']
    haplocs=[x+1 for x in range(len(myargs)) if myargs[x]=='--vcf_happy']
    if len(querylocs) != len(haplocs):
        print("Error: # arguments preceded by --vcf_query_preprocessed must equal # arguments preceded by --vcf-happy")
    input_vcf_dict={i:{'query':myargs[querylocs[i]],   'happy':myargs[haplocs[i]] } for i in range(len(querylocs))}
    if '--traindir' in myargs:
        traindir=myargs[myargs.index('--traindir')+1]
    else:
        print("Error - missing training directory parameter:  should follow --traindir in command")
    if '--mod_paramdict' in myargs:
        jsonfile=myargs[myargs.index('--mod_paramdict')+1]
    else:
        print("Warning - No json model param dict provided...using default parameters.")
else:
    input_vcf_dict={0:{'query':'/igm/projects/FOR_JEFF/germfilt/tm/run_02272020/query_preprocessed_jeffsrun.vcf.gz',
                       'happy':'/igm/projects/FOR_JEFF/germfilt/tm/run_02272020/happyout_jeffsrun.vcf.gz'}}
    traindir='/igm/projects/FOR_JEFF/germfilt/tm/run_02272020/training_022720_beforelunch'
    #jsonfile='/igm/projects/FOR_JEFF/germfilt/tm/batch_1_giab_1_02262020_v2/training/mod_params.json'
perffile=traindir+'/allmod_summary.csv'
featimpfile=traindir+'/feature_importances.csv'
joincols=['CHROM','POS','REF']
moddir=traindir+'/mod'
f.quietly_create_directory(traindir)
f.quietly_create_directory(moddir)


# In[ ]:





# In[6]:


## SET MODEL PARAMETERS
imp.reload(f)
jdict={1:'1_5',6:'6_15',16:'16_PLUS'}
vartypes=['SNP']+[x+jdict[j] for x in ['I','D'] for j in [1,6,16]]+['MNP']
req_cols=['modnum','vartype','zyg','moddir_loc']
default_paramvals=f.grab_default_procedural_param_dict()
#default_paramvals['n_estimators']=500
## GRAB PARAMDICT   
if jsonfile is not None:
    mod_paramdict=f.jread(jsonfile)
else:
    mod_paramdict=f.grab_mod_paramdict_generic(traindir)
# VALIDATE, ADD DEFAULT VALUES
for i in mod_paramdict:
    for y in req_cols:
        if y not in mod_paramdict[i]:
            print("Error... model parameter-set "+str(i)+" is missing required parameter "+y)
    for z in default_paramvals:
        if z not in mod_paramdict[i]:
            mod_paramdict[i][z]=default_paramvals[z]


# In[ ]:





# In[7]:


## LOAD VCFS
imp.reload(f)
vcfdict={}
formatted_df_dict={}
sampdict={}
trainflag=1
for i in input_vcf_dict:
    queryfile=input_vcf_dict[i]['query']
    truthfile=input_vcf_dict[i]['happy']
     ## set trainflag=1 to require exactly one sample
    ## READ RAW QUERY...ASSIGN ZYG
    queryfile=input_vcf_dict[0]['query']
    vcfdict['query']=f.grab_query_vcf_new(queryfile,trainflag)
    dfraw=vcfdict['query']['df']
    sampdict[i]=vcfdict['query']['samples'][0]
    samp=sampdict[i]
    vcfdict['query']['keepinds']=f.add_zyg(dfraw,'gtfeat_'+samp+'_GT')
    dfpre=dfraw.loc[vcfdict['query']['keepinds']]
    ## ASSIGN VARTYPE
    f.add_vartype(dfpre)
    ## HC OUTPUT
    vcfdict['happy']=f.grab_happy_vcf_new(truthfile)
    hfraw=vcfdict['happy']['df']
    hfpre=hfraw.loc[hfraw[hfraw['ALT_list'].apply(len)==1].index]
    hfpre['ALT']=hfpre['ALT_list'].apply(lambda x: str(x[0]))
    ## JOIN
    refcols=['CHROM','POS','REF','ALT']
    df=pd.merge(left=dfpre.assign(queryflag=1),
              right=hfpre.assign(truthflag=1).drop(columns=[x for x in hfpre.columns if x in dfpre.columns and x not in refcols]),
              on=refcols,how='inner')
    ## DEFINE TRUTH VALUES
    df.loc[df['BD_TRUTH'].isnull(),'BD_TRUTH']='X'
    zyg_mismatch_inds=df[((df.BD_QUERY=='FP') & (df.BD_TRUTH=='FN')) & (df.BLT_TRUTH != df.BLT_QUERY)].index
    df['y']=0
    df.loc[df.BD_QUERY=='TP','y']=1
    df.loc[zyg_mismatch_inds,'y']=1
    formatted_df_dict[i]=df.copy()

## COMBINE VCF DATA
df=pd.concat([formatted_df_dict[i].rename(columns={x:x.replace(sampdict[i],sampdict[0]) for x in formatted_df_dict[i].columns})  for i in range(len(formatted_df_dict))]).reset_index(drop=True)


# In[ ]:





# In[8]:


## BUILD/FORMAT FEATS
## FORMAT sample-specific df (with zygosity)
imp.reload(f)
df,sampfeats=f.grab_basefeats(df)
cursamp=sampdict[0]
featdict={x:x.replace('gtfeat_'+cursamp,'gtfeat') for x in sampfeats if 'gtfeat_'+cursamp in x or 'infofeat_' in x or 'genfeat_' in x}
gf=df.rename(columns=featdict)
adhocfeats=['gtfeat_total_AD']
## ADD TOTAL AD... we require AD annotation (but this seems reasonable)
gf['gtfeat_total_AD']=gf[[x for x in gf.columns if 'gtfeat_AD' in x]].apply('sum',axis=1)
## ADD SS zygosity
gtcol='gtfeat_'+cursamp+'_GT'
gf['zyg']='het'
gf.loc[gf[gtcol].apply(lambda x: (x.count('1'))>1),'zyg']='hom'
globfeats=[featdict[x] for x in featdict]+adhocfeats
cutfeats=['gtfeat_PS']
globfeats=[x for x in globfeats if x not in cutfeats]
catfeats_hardcoded=['genfeat_REF', 'genfeat_ALT']
for x in [y for y in globfeats if y not in catfeats_hardcoded]:
    gf[x]=gf[x].apply(lambda x: float(x) if x is not None else x)
globfeats=[x for x in globfeats if 'VQS' not in x]


# In[ ]:





# In[9]:


## INDICES
indict={}
for vartype in vartypes:
    indict[vartype]=gf[gf.vartype==vartype].index
for zyg in ['hom','het']:
    indict[zyg]=gf[gf.zyg==zyg].index


# In[ ]:





# In[10]:


## MAIN MODEL LOOP
imp.reload(f)
perfdict={}
perfflag=1
zygperfdict={}
lc=0
shapdict={}
moddir=traindir+'/mod'
f.quietly_create_directory(moddir)
## RESTRICT BY VARTYPE
#for vartype in ['SNP','INDEL']:
#for vartype in [x for x in vartypes if x != 'SNP']:
for vartype in vartypes:
    vardat=gf.loc[indict[vartype]]
    ## FORMAT BY VARTYPE
    if vartype=='SNP':
        locfeats=[x for x in globfeats if x not in ['genfeat_REF_len','genfeat_ALT_len']]
        catfeats=catfeats_hardcoded
    else:
        locfeats=[x for x in globfeats if x not in catfeats_hardcoded]
        catfeats=[]
    ## ONEHOT
    ohfeats=f.onehot_encode_cold(vardat,catfeats)
    modfeats=[x for x in locfeats if x not in catfeats]+ohfeats
    ## RESTRICT BY ZYG
    for zyg in ['het','hom']:
#    for zyg in ['het']:
        locdat_base=vardat.loc[vardat.index.intersection(indict[zyg])]
        locdat_base=locdat_base[(locdat_base.y >=0)]
        modinds_mid=[x for x in mod_paramdict if ((mod_paramdict[x]['zyg']==zyg) and (mod_paramdict[x]['vartype']==vartype))]
        if len(modinds_mid)==0:
            print("No models for " + vartype + " "+zyg+"...skipping")
            continue
        constrained_paramdat=pd.DataFrame({x:mod_paramdict[x] for x in modinds_mid}).transpose()[['nfold','train_prop','training_depth_cutoff','seed']]
        distinct_depth_cutoffs=list(set([mod_paramdict[x]['training_depth_cutoff'] for x in mod_paramdict]))
        distinct_depth_cutoff_inddict={x:locdat_base[locdat_base.gtfeat_total_AD >= x].index for x in distinct_depth_cutoffs}
        for depth in distinct_depth_cutoffs:
            modinds_loc=[x for x in modinds_mid if mod_paramdict[x]['training_depth_cutoff']==depth]
            locdat=locdat_base.loc[distinct_depth_cutoff_inddict[depth]]
            if len(locdat)==0:
                continue
            distinct_trainprops=list(set([mod_paramdict[x]['train_prop'] for x in modinds_loc]))
            for trainprop in distinct_trainprops:
                modinds_prop=[x for x in modinds_loc if mod_paramdict[x]['train_prop']==trainprop]
                fold_dict_glob=f.grab_traininds(locdat,trainprop)
                distinct_max_trainsizes=list(set([mod_paramdict[x]['max_trainsize'] for x in modinds_prop]))
                for max_trainsize in distinct_max_trainsizes:
                    modinds_maxtrainsize=[x for x in modinds_prop if mod_paramdict[x]['max_trainsize']==max_trainsize]
                    nfold_max=max([mod_paramdict[x]['nfold'] for x in modinds_maxtrainsize])
                    if (nfold_max > len(fold_dict_glob)):
                        print("Requested number of folds "+str(nfold_max) + " exceeds that permitted by training proportion "+str(trainprop)+". Only running "+str(len(fold_dict_glob)) +" folds.")
                        nfold_max=len(fold_dict_glob)
                    fold_dict=fold_dict_glob.copy()
                    #print("Handling models with depth-cutoff "+str(depth)+", trainprop "+str(trainprop) + ", max_trainsize "+str(max_trainsize))
                    for fold in range(nfold_max):
                        traindat=locdat.loc[fold_dict[fold]['train']]
                        if len(traindat)>max_trainsize:
                            traindat=traindat.sample(max_trainsize)
                        valdat=locdat.loc[fold_dict[fold]['val']]
                        for j in modinds_maxtrainsize:
                            paramdict=mod_paramdict[j]
                            moddir_loc=paramdict['moddir_loc']
                            f.quietly_create_directory(moddir_loc)
                            specs=[x for x in paramdict.keys() if x in xgb.XGBClassifier()._get_param_names()]
                            xgbdict=xgb.XGBClassifier().get_params()
                            for y in specs:
                                xgbdict[y]=paramdict[y]
                            mymod=xgb.XGBClassifier(**xgbdict)
                            #print("Training")
                            #g.grabtime()
                            mymod.fit(traindat[modfeats],traindat['y'])
                            modfile=moddir_loc+'/mod_'+str(fold)+'.pkl'
                            print("Wrote model to file "+modfile)
                            sys.stdout.flush()
                            #g.grabtime()
                            pickle.dump(mymod, open(modfile, "wb"),protocol=2)
                            valdat['yprob_'+str(j)]=mymod.predict_proba(valdat[modfeats])[:,1]
                        if perfflag==1:
                            perfdat=f.grab_modperf_manymods(valdat,'y',paramdict['PASS_cutoffs'],modinds_loc,vartype,zyg).assign(
                                vartype=vartype,zyg=zyg,modfold=fold,modfile=modfile)
                            perfdict[lc]=perfdat;lc+=1
                            ## FEATURE IMPORTANCES
                        if featimpflag==1:
                            explainer=shap.TreeExplainer(mymod)
                            shapdat_loc=pd.DataFrame(explainer.shap_values(locdat[modfeats]),columns=['shap_'+x for x in modfeats])
                            shapdat=shapdat_loc.apply('mean',0).reset_index().rename(columns={'index':'shapfeat',0:'shapval'})
                            shapdat['featnice']=shapdat['shapfeat'].apply(lambda x: x.replace('shap_genfeat_','').replace('shap_gtfeat_','').replace('shap_infofeat_',''))
                            shapdat['shapval_abs']=shapdat['shapval'].apply(abs)
                            shapdat=shapdat.sort_values('shapval_abs',ascending=False).assign(
                                ).assign(vartype=vartype,zyg=zyg,modnum=paramdict['modnum'],modfold=fold);
                            shapdict[lc]=shapdat[['featnice','shapval','vartype','zyg']].rename(columns={'featnice':'feat'}).assign(modnum=paramdict['modnum'])
## SAVE PERFORMANCE REPORT
if perfflag==1:
    perfdat=pd.concat(perfdict,ignore_index=True)
    perfdat=perfdat.sort_values(['vartype','zyg','modid','thresh','modfold'],ascending=[False,True,True,True,True])
    perfdat.to_csv(perffile,index=False)
    print("Wrote performance file to "+perffile)
## CONSOLIDATE FEATURE IMPORT
if featimpflag==1:
    featimpdat=f.grab_featimp(shapdict)
    featimpdat.to_csv(featimpfile,index=False)
    print("Wrote feature-importance file to "+featimpfile)


# In[ ]:





# In[ ]:





# In[ ]:




