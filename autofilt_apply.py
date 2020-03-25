#!/usr/bin/env python
# coding: utf-8

# In[14]:


## IMPORT
import sys
import pandas as pd
#import numpy as np   
import vcf 
import math  
import pickle  
import imp
#import xgboost as xgb 
import f
myargs=f.grab_bashscript_args('run_apply_3')
## TODO
## chage 'thresh' to 'pass_score'
## mention indel capabilities in MS
## chr21 test files


# In[15]:


## GRAB INPUT PARAMS
jdict={1:'1_5',6:'6_15',16:'16_PLUS'}
vartypes=['SNP','MNP']+[x+jdict[j] for x in ['I','D'] for j in [1,6,16]]
if 1==1:
    #myargs=sys.argv
    infile=myargs[myargs.index('--vcf_in')+1]
    outfile=myargs[myargs.index('--vcf_out')+1]
    moddict={vartype:{zyg:[myargs[i+1] for i,x in enumerate(myargs) if x=='--modfile_'+vartype+'_'+zyg] for zyg in ['hom','het']} for vartype in vartypes}
    threshdict={vartype:{zyg:float(myargs[myargs.index('--thresh_'+vartype+'_'+zyg)+1]) for zyg in ['het','hom']} for vartype in vartypes}
else:
    infile='/igm/home/jbg001/git/germfilt/germline_variant_filtering/kf/tm/temp.vcf.gz'
    outfile='/igm/home/jbg001/git/germfilt/germline_variant_filtering/kf/tm/applied_filter.vcf'
    moddict={x:{y:{} for y in ['hom','het']} for x in ['SNP','INDEL']}
    moddict['SNP']['het']='/igm/projects/FOR_JEFF/germfilt/tm/batch_1_giab_1_02262020_v2/training/mod/mod_SNP_het_0/mod_0.pkl'
    moddict['SNP']['hom']='/igm/projects/FOR_JEFF/germfilt/tm/batch_1_giab_1_02262020_v2/training/mod/mod_SNP_hom_1/mod_0.pkl'
    moddict['INDEL']['het']='/igm/projects/FOR_JEFF/germfilt/tm/batch_1_giab_1_02262020_v2/training/mod/mod_INDEL_het_2/mod_0.pkl'
    moddict['INDEL']['hom']='/igm/projects/FOR_JEFF/germfilt/tm/batch_1_giab_1_02262020_v2/training/mod/mod_INDEL_hom_3/mod_0.pkl'
    threshdict={x:{y:{} for y in ['hom','het']} for x in ['SNP','INDEL']}
    threshdict['SNP']['het']=0.75
    threshdict['SNP']['hom']=0
    threshdict['INDEL']['het']=0.25
    threshdict['INDEL']['hom']=0
 


# In[ ]:





# In[16]:


## LOAD VCFs, SET INITIAL PREDS
imp.reload(f)
## LOAD VCFS
if 1==1:
    vcfdict={}
    vcfdict['query']=f.grab_query_vcf_new(infile,0)
    samps=vcfdict['query']['samples']
    dfraw=vcfdict['query']['df']
## IDENTIFY PROBLEMATIC INDS AT RECORD LEVEL...DEFINE ALT
if 1==1:
    multiallelic_inds=dfraw[dfraw.ALT_list.apply(len)>1].index
    noalt_inds=dfraw[dfraw.ALT_list.apply(len)==0].index
    dftemp=dfraw.loc[dfraw.index.difference(multiallelic_inds.union(noalt_inds))]
    dftemp['ALT']=dftemp.ALT_list.apply(lambda x: x[0])
    dftemp.loc[dftemp.ALT.isnull(),'ALT']='X'
    dftemp['ALT']=dftemp['ALT'].apply(str)
    weirdalt_inds=dftemp[dftemp['ALT'].apply(lambda x:x.replace('A','').replace('C','').replace('G','').replace('T','')).apply(len)>0].index
## LOOK AT GENOTYPES, ASSIGN ZYG
if 1==1:
    zygdict={}
    for samp in samps:
        zygdict[samp]={}
        gtcol='gtfeat_'+samp+'_GT'
        gtlistcol='gtlist_'+samp
        dftemp[gtlistcol]=dftemp[gtcol].apply(lambda x: [int(y) for y in x.replace('|','/').replace('.','0').split('/')])
        nocallinds=dftemp[dftemp[gtlistcol].apply(max)==0].index
        multiallelic_call_inds=dftemp[dftemp[gtlistcol].apply(max)>1].index
        tossinds=nocallinds.union(multiallelic_call_inds)
        hetinds=dftemp[dftemp[gtlistcol].apply(min)==0].index.difference(tossinds)
        hominds=dftemp[dftemp[gtlistcol].apply(min)==1].index.difference(tossinds)
        zygdict[samp]['het']=hetinds.difference(weirdalt_inds)
        zygdict[samp]['hom']=hominds.difference(weirdalt_inds)
        zygdict[samp]['nocall']=nocallinds
    dfpre=dftemp.loc[dftemp.index.difference(weirdalt_inds)]
# ASSIGN VARTYPES
if 1==1:
    for x in ['REF','ALT']:
        dfpre[x+'_len']=dfpre[x].apply(len)
    dfpre['len_diff']=dfpre['ALT_len']-dfpre['REF_len']
    vtdict={}
    vtdict['SNP']=dfpre[(dfpre.len_diff==0) & (dfpre.REF_len==1)].index
    vtdict['MNP']=dfpre[(dfpre.len_diff==0) & (dfpre.REF_len>1)].index
    cutoffdict=[1,6,16,10000000]
    for cutoffind in range(len(cutoffdict)-1):
        jlower=cutoffdict[cutoffind]
        jupper=cutoffdict[cutoffind+1]
        vtdict['I'+jdict[jlower]]=dfpre[(dfpre.len_diff >= jlower) & (dfpre.len_diff < jupper)].index
        vtdict['D'+jdict[jlower]]=dfpre[(dfpre.len_diff <= -jlower) & (dfpre.len_diff > -jupper)].index
    dfpre['vartype']=None
    for x in vtdict:
        dfpre.loc[vtdict[x],'vartype']=x


# In[ ]:





# In[17]:


## FEATS...INITIALIZE PREDS
df = dfpre.where(pd.notnull(dfpre), None)
df,sampfeats=f.grab_basefeats(df)
for samp in samps:
    df['gtfeat_'+samp+'_total_AD']=df[[x for x in df.columns if 'gtfeat_'+samp+'_AD_' in x]].apply('sum',axis=1)
    df['yprob_'+samp]=0.0
    df['ypred_'+samp]=0
cutfeats=['gtfeat_PS']
catfeats_hardcoded=['genfeat_REF', 'genfeat_ALT']


# In[ ]:





# In[18]:


## MAIN LOOP
for vartype in vartypes:
    vardat=df.loc[vtdict[vartype]].copy()
    catfeats=[]
    if vartype=='SNP':
        catfeats=catfeats_hardcoded
    dummyvar=f.onehot_encode_cold_allvals(vardat,catfeats)
    ## COMPUTE 
    for zyg in ['het','hom']:
        for samp in samps:
            print('Applying filter to '+ vartype + ', ' + zyg+', '+samp)
            locdat_base=vardat.loc[vardat.index.intersection(zygdict[samp][zyg])]
            if len(locdat_base)==0:
                continue
            featdict={x:x.replace('gtfeat_'+samp,'gtfeat') for x in locdat_base.columns if 'gtfeat_'+samp in x}
            locdat=locdat_base.rename(columns=featdict)
            print(len(locdat))
            for modcount,modfile in enumerate(moddict[vartype][zyg]):
                with open(modfile, 'rb') as ff:
                    mymod = pickle.load(ff)
                feats=mymod.get_booster().feature_names
                if modcount==0:
                    for x in feats:
                        if x not in locdat.columns:
                            locdat[x]=float('nan')
                        locdat.loc[locdat[x].isnull(), x]=float('nan')
                        locdat[x]=locdat[x].apply(float)
                df.loc[locdat.index,'yprob_'+samp+'_'+str(modcount)]=mymod.predict_proba(locdat[feats])[:,1]
            df.loc[locdat.index,'yprob_'+samp]=df.loc[locdat.index][['yprob_'+samp+'_'+str(modcount) for modcount in range(len(moddict[vartype][zyg]))]].apply('mean',axis=1)
            df.loc[locdat.index,'ypred_'+samp]=1*(df.loc[locdat.index]['yprob_'+samp] >= threshdict[vartype][zyg])
            


# In[ ]:





# In[19]:


## ADD PREDICTIONS TO RAW DATAFRAME
## PASS ALL MULTIALLELIC INDS 
## PASS ALL INDS WITH ALT CONTAINING NON ACGT CHARACTERS 
## (THE NOCALLS WERE MANUALLY FAILED IN THE MODEL-APPLICATION SETUP)
multiallelic_inds,noalt_inds,weirdalt_inds
for samp in samps:
    predcol='ypred_'+samp
    dfraw[predcol]=0
    dfraw.loc[multiallelic_inds,predcol]=1 
    dfraw.loc[weirdalt_inds,predcol]=1    
    dfraw.loc[df.index,predcol]=df['ypred_'+samp]
## PASS ALL POSITIONS WITH A PASSED VARIANT
dfraw['ypred']=dfraw[['ypred_'+x for x in samps]].apply(max,axis=1)


# In[ ]:





# In[20]:


## WRITE OUTPUT
## Write VCF
predcodedict={0:'LowQual',1:'PASS'}
vcf_reader = vcf.Reader(filename=infile,compressed=True)
vcf_reader.filters['LowQual']=vcf_reader.filters[list(vcf_reader.filters.keys())[0]]._make(['LowQual','Flagged as low quality by AutoFilter models'])
if 'PASS' not in vcf_reader.filters:
    vcf_reader.filters['PASS']=vcf_reader.filters[list(vcf_reader.filters.keys())[0]]._make(['PASS','Passed by AutoFilter models'])
vcf_writer= vcf.Writer(open(outfile, 'w'), vcf_reader)
recdict=vcfdict['query']['recdict']
recdict_active=recdict.copy()
for i in range(len(recdict)):
    record=recdict[i]
    record.FILTER=predcodedict[dfraw.loc[i]['ypred']]
    if 'FT' in record.genotype(samp).data._asdict():
        record.genotype(samp).data=record.genotype(samp).data._replace(FT=record.FILTER)
    vcf_writer.write_record(record)
    vcf_writer.flush() 
    i=i+1
print('Wrote annotated vcf to '+outfile)


# In[ ]:





# In[21]:


## HAP
hapfile='/'.join(infile.split('/')[:infile.count('/')])+'/happyout_ph.vcf.gz'#tfile='/igm/projects/FOR_JEFF/germfilt/tm/03062020_noon_noprepy_batch_1_giab_4/applied_granular_train_031720_afternoon_filt_happed.vcf.gz'
hapdict=f.grab_happy_vcf_new(hapfile)
hapdat_raw=hapdict['df']
hapdat=hapdat_raw[~(hapdat_raw['ALT_list'].apply(len)==0)].copy()
hapdat['ALT']=hapdat['ALT_list'].apply(lambda x: ';'.join([str(y) for y in x]))
hapdat.loc[hapdat.BD_TRUTH.isnull(),'BD_TRUTH']='X'
hapdat.loc[hapdat.BD_QUERY.isnull(),'BD_QUERY']='X'
hapdat['y']=0
hapdat.loc[(hapdat.BD_QUERY=='TP'),'y']=1
hapdat.loc[(hapdat.BD_TRUTH=='FN') & (hapdat.BD_QUERY=='FP'),'y']=1


# In[ ]:





# In[22]:


## JOIN TO HAP
filtcols=['ypred_'+samp,'tranchflag_1','tranchflag_2','tranchflag_both','gatkflag','lqflag','glqflag','lqtflag']
filtcols=['ypred_'+samp,'yprob_'+samp]
refcols=['CHROM','POS','REF','ALT']
joindat=pd.merge(left=df[refcols+['vartype']+filtcols],right=hapdat,
                 on=['CHROM','POS','REF','ALT'],how='inner').rename(columns={x+'_'+samp:x for x in ['ypred','yprob']})


# In[23]:


joindat.assign(n=1).groupby(['vartype','ypred','y'])['n'].count().reset_index()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


## EXTRACT FILTS
dfcur=dfraw[dfraw.ALT_list.apply(len)==1].copy()
dfcur['filtstring']=dfcur['FILTER'].apply(lambda x: ';'.join(x))
dfcur['tranchflag_1']=1-1*(dfcur.filtstring.str.contains('VQSRTrancheSNP99.00'))
dfcur['tranchflag_2']=1-1*(dfcur.filtstring.str.contains('VQSRTrancheSNP99.90'))
dfcur['tranchflag_both']=dfcur['tranchflag_1']*dfcur['tranchflag_2']
dfcur['gatkflag']=1-1*(dfcur.filtstring.str.contains('gatkRecommendIndelFilter'))
dfcur['lqflag']=1-1*(dfcur.filtstring.str.contains('LowQual'))
dfcur['glqflag']=dfcur['gatkflag']*dfcur['lqflag']
dfcur['lqtflag']=dfcur['tranchflag_2']*dfcur['lqflag']


# In[ ]:


## JOIN
filtcols=['ypred_'+samp,'tranchflag_1','tranchflag_2','tranchflag_both','gatkflag','lqflag','glqflag','lqtflag']
refcols=['CHROM','POS','REF','ALT']
joindat=pd.merge(left=dfcur[~dfcur.ALT.isnull()][refcols+['vartype']+filtcols],right=hapdat,
                 on=['CHROM','POS','REF','ALT'],how='inner')


# In[ ]:


joindat['ypred_ext']=joindat['lqtflag']
joindat.loc[joindat.vartype == 'INDEL','ypred_ext']=joindat.loc[joindat.vartype == 'INDEL']['glqflag']
joindat['ypred_loc']=joindat['ypred_'+samp]


# In[ ]:


agdat=joindat.assign(n=1).groupby(['vartype','y','ypred_loc','ypred_ext'])['n'].count().reset_index()
for x in ['loc','ext']:
    agdat['TP_'+x]=agdat['ypred_'+x]*agdat['y']*agdat['n']
    agdat['FP_'+x]=agdat['ypred_'+x]*(1-agdat['y'])*agdat['n']
    agdat['TN_'+x]=(1-agdat['ypred_'+x])*(1-agdat['y'])*agdat['n']
    agdat['FN_'+x]=(1-agdat['ypred_'+x])*agdat['y']*agdat['n']
outcols=[x+'_'+y for x in ['TP','FP','FN','TN'] for y in ['loc','ext']]


# In[ ]:


cagdat=agdat.groupby('vartype').agg({x:'sum' for x in outcols}).reset_index()
truthvars=['TP','FP','FN','TN']
tempdict={}
for x in ['loc','ext']:
    tempdict[x]=cagdat[['vartype']+[y+'_'+x for y in truthvars]].rename(columns={y+'_'+x:y for y in truthvars}).assign(filt=x)
#agdat.vartype.unique()
#{x:'sum' for x in outcols}
sdat=pd.concat(tempdict,ignore_index=True)[['vartype','filt']+truthvars]
sdat['recall']=sdat['TP']*1.0/(sdat['TP']+sdat['FN'])
sdat['precision']=sdat['TP']*1.0/(sdat['TP']+sdat['FP'])


# In[ ]:


import matplotlib.pyplot as plt
fig,ax=plt.subplots()
for filt in ['loc','ext']:
    plotdat=sdat[sdat.filt==filt]
    plt.scatter(plotdat.precision,plotdat.recall,label=filt)
ax.legend()


# In[ ]:


sdat


# In[ ]:





# In[ ]:





# In[ ]:




