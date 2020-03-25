#!/usr/bin/env python
# coding: utf-8

# In[1]:


## IMPORT
import sys
import os
import f
import app_paths


# In[ ]:





# In[18]:


## READ COMMAND-LINE ARGUMENTS
if 1==1:
    myargs=sys.argv
    cmdline_params=['vcf_truth','vcf_query','bedfile_truth','bedfile_query',
                          'reffile','scriptfile_preprocess','locdir','runtag']
    argdict={}
    for arg in cmdline_params:
        if '--'+arg in myargs:
            argdict[arg]=myargs[myargs.index('--'+arg)+1]
        else:
            print("Missing required parameter "+arg+", which should be preceded by --"+arg+" in command")
else:
    argdict={'vcf_truth': '/igm/home/pjb004/NIST/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz',
             'vcf_query': '/igm/projects/GATK4.1.2_GIAB/WES/CHURCHILL/Variants/NA12878.HC.vcf.gz',
             'bedfile_truth': '/igm/home/pjb004/NIST/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed',
             'bedfile_query': '/igm/projects/GATK4.1.2_GIAB/WES/CHURCHILL/TEMP/kit_padded.bed',
             'reffile': '/igm/apps/genomes/Homo_sapiens/human_g1k_v37_decoy/human_g1k_v37_decoy.fasta',
             'scriptfile_preprocess': '/igm/home/jbg001/git/germfilt/germline_variant_filtering/kf/script/run_preprocess.sh',
             'locdir': '/igm/projects/FOR_JEFF/germfilt/tm/run_02272020',
             'runtag': 'jeffsrun'}


# In[ ]:





# In[20]:


## LOAD/CREATE SCRIPT VARIABLES


cmd_dict=app_paths.cmd_dict
vcf_truth_raw=argdict['vcf_truth']
vcf_query_raw=argdict['vcf_query']
bedfile_truth=argdict['bedfile_truth']
bedfile_query=argdict['bedfile_query']
reffile=argdict['reffile']
scriptfile_preprocess=argdict['scriptfile_preprocess']
locdir=argdict['locdir']
runtag=argdict['runtag']

bedfile_intersected=locdir+'/giab_intersected_'+runtag+'.bed'

vcf_truth_filtered=locdir+'/truth_filtered_'+runtag+'.vcf.gz'
vcf_query_filtered=locdir+'/query_filtered_'+runtag+'.vcf.gz'
vcf_truth_prepyed=locdir+'/truth_filtered_prepyed_'+runtag+'.vcf.gz'
vcf_query_prepyed=locdir+'/query_filtered_prepyed_'+runtag+'.vcf.gz'
vcf_query_preprocessed=locdir+'/query_preprocessed_'+runtag+'.vcf.gz'
happy_outstring=locdir+'/happyout_'+runtag
vcf_happy=happy_outstring+'.vcf.gz'
happy_summary=happy_outstring+'.summary.csv'


# In[ ]:





# In[22]:


## WRITE PREPROCESSING SCRIPT
f.quietly_create_directory(locdir)
prepflag=1
if bedfile_query is None:
    print("No query bedfile given...assuming this is WGS.")
with open(scriptfile_preprocess,'w') as fx:
    fx.write("#!/bin/bash\n")
    fx.write("#$ -pe smp 8\n")
    fx.write(cmd_dict['bedtools_path']+' intersect -a '+bedfile_truth +' -b '+bedfile_query+' |\\\n')
    fx.write(cmd_dict['bedtools_path']+  '  sort -i -  > ' + bedfile_intersected +'\n')
    fx.write(cmd_dict['bedtools_path']+  ' intersect -a '+ vcf_truth_raw + ' -b ' + bedfile_intersected +' -header '+' | bgzip > '+ vcf_truth_filtered+ '\n')
    fx.write(cmd_dict['bedtools_path']+  ' intersect -a '+ vcf_query_raw + ' -b ' + bedfile_intersected +' -header '+' | bgzip > '+ vcf_query_filtered+ '\n')
    fx.write(cmd_dict['tabix_path'] + ' -f ' + vcf_truth_filtered+'\n')
    fx.write(cmd_dict['tabix_path'] + ' -f ' + vcf_query_filtered+'\n')
    gatkflag=0
    if 'gatk38_path' in cmd_dict:
        if cmd_dict['gatk38_path'] is not None:
            gatkflag=1
    if gatkflag==0:
        print("No GATK 3.8 path given, so TandemRepeat annotations will not be calculated.")
        fx.write('mv '+vcf_query_filtered+' '+vcf_query_preprocessed+'\n')
        fx.write('mv '+vcf_query_filtered+'.tbi '+vcf_query_preprocessed+'.tbi'+'\n')
    else:
        xx=1
        #fx.write('java -jar '+cmd_dict['gatk38_path']+'  --analysis_type VariantAnnotator -V '+
        #        vcf_query_filtered + ' -o ' + vcf_query_preprocessed + ' -R '+reffile + ' -A TandemRepeatAnnotator -L '+ bedfile_intersected +' -ip 50\n')
    fx.write(cmd_dict['happy_path'] + ' ' + vcf_truth_filtered +' ' +vcf_query_preprocessed + ' -r ' + reffile + 
             ' -f '+bedfile_intersected+ ' -o '+happy_outstring + ' --threads 4 --no-roc --no-json -X --output-vtc \n')
os.chmod(scriptfile_preprocess, 0o775)
print(" ")
print("Preprocessing script: "+scriptfile_preprocess)
print("Output file locations")
print("Truth vcf: "+vcf_happy)
print("Query vcf: "+vcf_query_preprocessed)
print("Hap.py summary: "+happy_summary)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




