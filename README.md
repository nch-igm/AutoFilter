
# VCF AutoFilter
* * *

Jeff Gaither <jeffrey.gaither@nationwidechildrens.org>
<br>
<br>

AutoFilter is a pipeline for training VCF germline filters on gold-standard 
[Genome in a Bottle](https://github.com/samtools/htslib) VCFs. 

<br>

<br>

## Machine-learning VCF filters in four easy steps. 

### Step #1 - Record application paths

```python
## app_paths.py 
cmd_dict={
    'bedtools_app':'/btpath/bedtools',
    'happy_app':'/hppath/hap.py',
    'gatk38_app':'/gatkpath/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar', ## optional
    'tabix_app':'/tabixpath/tabix',
    }
```
<br>

<br>


### Step #2 - Generate preprocessing script

```bash
python autofilt_preprocess.py \
--vcf_truth truth.vcf.gz \
--vcf_query query.vcf.gz \
--bedfile_truth truth.bed \
--bedfile_query query.bed \
--reffile reference.fa \
--locdir locdir \
--runtag runid \
--scriptfile_preprocess preprocess.sh
```

<br>

<br>

### Step #3 - Preprocess
```bash
preprocess.sh
```

<br>

<br>

### Step #4 - Train
```bash
python autofilt_train.py \
--vcf_query_preprocessed query_preprocessed.vcf.gz \ 
--vcf_happy hap.py_output.vcf.gz \ 
--traindir traindir  
```

<br>

<br>

### Models are now built...apply them to clinical VCF
```bash
python autofilt_apply.py \
--vcf_in unfiltered.vcf.gz \
--vcf_out filtered.vcf \
--modfile_SNP_het SNP_het_mod.pkl  			    --pass_score_SNP_het .5 \
--modfile_SNP_het SNP_hom_mod.pkl  			    --pass_score_SNP_hom 0 \
--modfile_I1_5_het I1_5_het_mod.pkl   		    --pass_score_I1_5_het .1 \
--modfile_I1_5_hom I1_5_hom_mod.pkl   		    --pass_score_I1_5_hom 0 \
--modfile_I6_15_het I6_15_het_mod.pkl   	    --pass_score_I6_15_het .1 \
--modfile_I6_15_hom I6_15_hom_mod.pkl   	    --pass_score_I6_15_hom 0 \
--modfile_I16_PLUS_het I16_PLUS_het_mod.pkl     --pass_score_I16_PLUS_het .1 \
--modfile_I16_PLUS_hom I16_PLUS_hom_mod.pkl     --pass_score_I16_PLUS_hom 0 \
--modfile_D1_5_het D1_5_het_mod.pkl   		    --pass_score_D1_5_het .1 \
--modfile_D1_5_hom D1_5_hom_mod.pkl    		    --pass_score_D1_5_hom 0 \
--modfile_D6_15_het D6_15_het_mod.pkl   	    --pass_score_D6_15_het .1 \
--modfile_D6_15_hom D6_15_hom_mod.pkl   	    --pass_score_D6_15_hom 0 \
--modfile_D16_PLUS_het D16_PLUS_het_mod.pkl     --pass_score_D16_PLUS_het .1 \
--modfile_D16_PLUS_hom D16_PLUS_hom_mod.pkl     --pass_score_D16_PLUS_hom 0 \
--modfile_MNP_het MNP_het_mod.pkl   		    --pass_score_MNP_het .1 \
--modfile_MNP_hom MNP_hom_mod.pkl   		    --pass_score_MNP_hom 0 
```
<br>
<br>

In the output, low-confidence variants will be marked __"LowQual"__ in the FILTER field. 

```python
...
21      46916516        .       AC      A       2156.73 PASS   ...
21      46916699        .       C       G       20.8    LowQual ...
21      46916737        .       C       G       235.77  PASS ...
...
```
<br>

* * *
<BR>

# AutoFilter Guide and FAQ

## Motivation
AutoFilter is a tool to automatically build false-positive filters for germline VCFs.
AutoFilter achieves the machine-learning-based accuracy of GATK's Variant Quality Score Recalibration tool (the industry standard)
while eliminating most of VQSR's more problematic requirements, including the need for:
1. a large training-set
2. careful selection of which VCF annotations to use 
3. the re-training of fresh models on every VCF to be filtered. 

In short, AutoFilter helps
institutions to build their own accurate VCF filters with minimum fuss.<br><br>

## Required resources
__FILES__<br>
__*Truth VCF from Genome in a Bottle website:*__  A theoretically
100%-accurate set of mutations from one of seven extensively-studied samples (the most
common being NA12878.) You can download it from the Genome in a Bottle ftp site: 
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/ (If 
the link doesn't work, try Internet Explorer.)<br>
__*Query VCF for Genome in a Bottle Sample:*__ A set of variants called 
on a Genome in a Bottle sample using your own institutional pipeline. <br>
__*High-confidence bed file for Genome in a Bottle sample:*__ Available from the
Genome in a Bottle website (link above).<br>
__*Capture bed file for your platform (if WES):*__  The bed file for your 
exon capture kit. <br>
__*Human reference genome file:*__ Hg37/19 or hg38. <br><br>
__APPLICATIONS__<br>
bedtools, hap.py, tabix and bgzip.<br><br>

## Preprocessing step (autofilt_preprocess.py): summary<br>
This step compares your query VCF to the gold-standard truth sets and identifies the true and
false variants. 
1. (WES only) Intersect the Genome in a Bottle bedfile with the bedfile for your capture kit.<br>
2. Filter both VCFs by the bedfile from Step 1.<br>
3. Compare truth and query VCFs using the program hap.py. (Requires Python 2.7)<br><br>

## Training options (autofilt_train.py) <br>
AutoFilter will happily train "vanilla" models if provided only with input and output 
files. However, you can achieve a high degree of model customization by passing in a json
file. <br>

```python
mod_paramdict={'0': {'modnum': 0,
  'vartype': 'SNP',
  'zyg': 'het',
  'moddir_loc': 'traindir/mod/mod_SNP_het_0',            ### opt: if not supplied, model i, fold foldnum is written to path /traindir/mod/mod_<i>_<vartype>_<zyg>/mod_<foldnum>.pkl.
  'PASS_cutoffs': [0,.1,.2,.3,.4,.5,.6,.7,.8,.9],        ### opt: model score cutoffs to use during validation    
  'train_prop': .75,                                     ### opt: proportion of variants to train on (remainder withheld for validation)
  'nfold':1,                                             ### opt: cross-validation folds
  'max_trainsize':200000,                                ### opt: maximimum # variants to be used in training (useful for WGS)    
  'training_depth_cutoff':0,                             ### opt: minimum total depth for a tranining variant (changing not recommended)
  'n_estimators':1000,                                   ### opt: xgboost parameter
  '<any_xgboost_parameter>':.1},                         ### opt: xgboost parameter
              '1': {...}
  ...
      }
```
<br>
<br>

## Reading model output, setting thresholds (autofilt_apply.py)
AutoFilter outputs a __model performance file__ to the path traindir/allmod_summary.csv. Indels
are sub-categorized by insertion/deletion and length. <br>

AutoFilter builds a separate model for each variant-type/zygosity combination. 
Each model scores its applicable variants in the range \[0,1\], with the score representing the level of confidence 
in the variant's being real. During training the model performance is
validated on a held-out set at calling thresholds of 0.1, 0.2, ..., 0.9. <br>

The model performance file shows the number of TP, FP... at different threshholds. These
thresholds are required inputs when applying the models, allowing the user
to intelligently take control over the tradeoff of sensitivity/specificity.<br>

| vartype | zyg | modid | thresh | n | TP | FP | FN | TN | 
| --- | --- | --- | --- | --- | --- | --- | --- | --- |  
| SNP | het | 0 | 0.1 | 456631 | 454229 | 1166 | 80 | 1156 |  
| SNP | hom | 1 | 0.1 | 295853 | 295824 | 29 | 0 | 0 |  
| I1_5 | het | 2 | 0.1 | 27187 | 26409 | 453 | 9 | 316 |  
| I1_5 | hom | 3 | 0.1 | 19533 | 19102 | 353 | 7 | 71 |  
| I6_15 | het | 4 | 0.1 | 2369 | 2286 | 59 | 1 | 23 |  
| I6_15 | hom | 5 | 0.1 | 1420 | 1306 | 73 | 13 | 28 |  
| I16_PLUS | het | 6 | 0.1 | 493 | 461 | 11 | 0 | 21 |  
| I16_PLUS | hom | 7 | 0.1 | 281 | 271 | 5 | 1 | 4 |  
| D1_5 | het | 8 | 0.1 | 31999 | 30306 | 834 | 31 | 828 |  
| D1_5 | hom | 9 | 0.1 | 19461 | 19077 | 323 | 5 | 56 |  
| D6_15 | het | 10 | 0.1 | 2842 | 2661 | 114 | 10 | 57 |  
| D6_15 | hom | 11 | 0.1 | 1544 | 1445 | 77 | 7 | 15 |  
| D16_PLUS | het | 12 | 0.1 | 655 | 604 | 21 | 7 | 23 |  
| D16_PLUS | hom | 13 | 0.1 | 299 | 288 | 7 | 1 | 3 |  
| SNP | het | 0 | 0.2 | 456631 | 454174 | 948 | 135 | 1374 |  
| ... | ... | ... | ... | ... | ... | ... | ... | ... |  
| D16_PLUS | hom | 13 | 0.9 | 299 | 280 | 3 | 9 | 7 |  
<br>

<br>

## FAQ
#### Incorrect removal of a true variant could prevent the successful diagnosis of a patient's disorder. Do you take any responsibility for patient outcomes influenced in part by AutoFilter?
Of course we take no responsibilty for any outcomes that might result, directly or indirectly, from the use of AutoFilter. But in practice,
removing false positives from variant call-set is a necessary evil for which machine-learning algorithms currently represent the best solution.<br>

AutoFilter's summary file will give a great idea on how many true variants AutoFilter is throwing away at any threshold. From here you can 
choose the calling threshold that achieves the sensitivity/specificity tradeoff most suitable to your needs. 
<br>
<br>

#### Can AutoFilter handle multi-sample VCFs?
Sure. A multi-sample variant will be flagged as PASS if it is passed in any sample.
<br>
<br>

#### Does AutoFilter work on all variants?
No. Most importantly, AutoFilter labels all multi-allelic variants as PASS. AutoFilter 
considers a variant to "multi-allelic" if it contains
more than one allele in the ALT field OR its genotype field "X/Y" includes an integer greater than one 
(though the second condition should imply the first). Note that these criteria also apply to 
multi-sample variants, even if the variant has only one alternate allele in any given sample.<br>

Variants with unconventional annotations (such as structural variants enclosed by < >) will always 
be passed, provided they do not have a reference genotype.<br>
<br>

#### During preprocessing, shouldn't you be running pre.py to ensure consistent format for INDELs?
This is indeed the best practice. However, pre.py changes the VCF annotations, meaning 
that a model trained on pre.py-processed VCFs will only run on VCFs that have been subjected to 
the same procedure. Moreover, in our experience the improvements gained from applying pre.py are
minimal. <br>

Since ease of use is one of our primary goals, it seemed better not 
to require our users to apply pre.py to every VCF in their pipeline. That said, users 
wishing to go this route can easily modify the file autofilt_preprocess.py to include 
a pre.py command.
<br>

<br>

#### Can I train on multiple VCFs?
Sure. Just pass multiple "--vcf_query_preprocessed" and  "--vcf_happy" arguments to the training script
autofilt_train.py in 
a consistent order. 
<br>

<br>

#### When building custom models, do I have to include entries for all sixteen types (SNP het, SNP hom, I1_5_het, ..., MNP_hom)?
Nope!
<br>

<br>

#### Why does AutoFilter only support xgboost architecture? I want to use neural networks/logistic regression/random forests/support vector machines!
None of these model-types is robust enough to blindly train on all features in an arbitrary VCF and achieve good results.
The remarkable discerning power of gradient-boosted trees is really what makes AutoFilter possible in the first place. 
One could build an AutoFilter-like package for applying these other architectures to VCFs, but accuracy and ease 
of use would both suffer. 
<br>

<br>

#### You say AutoFilter trains on all the features in a VCF. Won't some of these features just be noise that impairs accuracy?
This is another problem that is largely obviated by the gradient-boosted tree architecture, which focuses on important features and ignores the 
uninformative ones. In our experience, feature-pruning an xgboost model does not lead to significant gains in accuracy. Regardless,
AutoFilter achieves solid performance even when all features are included.
<br>
<br>

#### What version of Python is AutoFilter supposed to run on? And is this the same version required by the preprocessing script?
AutoFilter is designed to work in Python 3.6, though it works in 2.7. The program
hap.py _requires_ python 2.7, meaning the script preprocess.sh in Step 3 must be run in a
python 2.7 environment. 

We recommend the use of a version-flexible distribution of Python such as Anaconda. 
You can seamlessly add your package-loading commands to the pre-processing script-writer
autofilt_preprocess.py as shown below:

```python
## autofilt_preprocess.py
...
with open(scriptfile_preprocess,'w') as fx:
    fx.write("#!/bin/bash\n")
    fx.write("#$ -pe smp 8\n")
    fx.write("conda activate your_python2.7_env")    ##### ADD THIS LINE   ######
    fx.write(cmd_dict['bedtools_path']+' intersect -a '+bedfile_truth +' -b '+bedfile_query+' |\\\n')
    fx.write(cmd_dict['bedtools_path']+  '  sort -i -  > ' + bedfile_intersected +'\n')
```
<br>

#### Can I use AutoFilter to identify somatic variants?
No, this sounds like a terrible idea. Somatic variant calling is much harder
than germline calling owing to the enormously larger genetic diversity of a tumor sample. 
There are many tools
(MuTect, VarScan) specifically designed to address this difficult problem. 


