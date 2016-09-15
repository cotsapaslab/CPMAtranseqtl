#!/bin/bash
##########################################################################################
### Wrapper script for trans-eQTL analysis by Boel Brynedal ##############################
##########################################################################################



# In directory data
cd data
wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2009-01_phaseIII/plink_format/hapmap3_r2_b36_fwd.qc.poly.tar.bz2

bunzip2 -c hapmap3_r2_b36_fwd.qc.poly.tar.bz2 | tar xvf -

mv ./hapmap3_pop/hapmap3_r2_b36_fwd.LWK.qc.poly.map LWK.map
mv ./hapmap3_pop/hapmap3_r2_b36_fwd.YRI.qc.poly.map YRI.map
mv ./hapmap3_pop/hapmap3_r2_b36_fwd.MKK.qc.poly.map MKK.map
mv ./hapmap3_pop/hapmap3_r2_b36_fwd.LWK.qc.poly.ped LWK.ped
mv ./hapmap3_pop/hapmap3_r2_b36_fwd.YRI.qc.poly.ped YRI.ped
mv ./hapmap3_pop/hapmap3_r2_b36_fwd.MKK.qc.poly.ped MKK.ped

cd ..

##########################################################################################
# Genotype data pre-processing

plink --noweb --file ./data/MKK --maf 0.15 --filter-founders --hwe 0.000001 --make-bed ./data/mkk15

plink --noweb --file ./data/YRI --maf 0.15 --filter-founders --hwe 0.000001 --make-bed ./data/yri15

plink --noweb --file ./data/LWK --maf 0.15 --filter-founders --hwe 0.000001 --make-bed ./data/lwk15

# extract common SNPs 

Rscript ./scripts/commonBIM.R ./data/mkk15.bim ./data/yri15.bim ./data/lwk15.bim


plink --noweb --file ./data/mkk15 --extract commonSNPs --make-bed ./data/mkk15_common

plink --noweb --file ./data/yri15 --extract commonSNPs --make-bed ./data/yri15_common

plink --noweb --file ./data/lwk15 --extract commonSNPs --make-bed ./data/lwk15_common


##########################################################################################
### Update SNP annotation to hg38
# How long did this take? 2 days. (3 total)

# Download UCSC SNP data base in hg38, tab delimited
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp142.txt.gz
gunzip snp142.txt.gz

python plinkAnno.py snp142.txt

rm snp142.txt

# Update base pair and chromosome for each SNP in our data sets
plink --bfile lwk15_common --update-map ./data/snp142_bp.txt --make-bed --out lwk15_common_hg38bp
plink --bfile lwk15_common_hg38bp --update-chr ./data/snp142_chr.txt --make-bed --out lwk15_common_hg38

plink --bfile mkk15_common --update-map ./data/snp142_bp.txt --make-bed --out mkk15_common_hg38bp
plink --bfile mkk15_common_hg38bp --update-chr ./data/snp142_chr.txt --make-bed --out mkk15_common_hg38

plink --bfile yri15_common --update-map ./data/snp142_bp.txt --make-bed --out yri15_common_hg38bp
plink --bfile yri15_common_hg38bp --update-chr ./data/snp142_chr.txt --make-bed --out yri15_common_hg38




##########################################################################################
# Gene expression data pre-processing

# Collect data

wget http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-264/E-MTAB-264.processed.1.zip
unzip E-MTAB-264.processed.1.zip

sed '2d' MKK_p3_expression.txt > MKK_expression.txt
sed '2d' YRI_p3_expression.txt > YRI_expression.txt
sed '2d' LWK_p3_expression.txt > LWK_expression.txt

# Filter and normalize

Rscript ./scripts/filtExp3.R

# Creates a plink phenotype file for each population based on the filtered data. 

##########################################################################################
# Gene expression annotation
# The annotation file gencode20.hg38.hgnc.sorted.bed is supplied. The code for how this annotation file was created is included in the Annotation_README.txt file. 


##########################################################################################
# Detection of optimal lambda - number of 'genetic principal components'
# How long did this take? Took me about a week (excluding linear regression, another week) to get to the boxplots. But then I was struggling to get things to work. 

 
# This is done for 100 genes, 50 random and 50 randomly selected among genes high lambda (top 5%).
# Requires flashpca (https://github.com/gabraham/flashpca)

# I did not use seed.. So can not replicate the exact same gene selection as before.

# Linear regression without covariates (needed to calculate lambda)

plink --bfile mkk15_common --linear --pheno MKKphenoQ --all-pheno --out ./MKKassoc/MKK

plink --bfile yri15_common --linear  --pheno YRIphenoQ --all-pheno --out ./YRIassoc/YRI

plink --bfile lwk15_common --linear --pheno LWKphenoQ --all-pheno  --out ./LWKassoc/LWK



# Calculate lambda for all genes in data set.

Rscript ./scripts/calcLambda.R ./LWKfiles/LWK.P .assoc.linear_z 9085 LWKlambda

Rscript ./scripts/calcLambda.R ./MKKfiles/MKK.P .assoc.linear_z 9085 MKKlambda

Rscript ./scripts/calcLambda.R ./YRIfiles/YRI.P .assoc.linear_z 9085 YRIlambda

# Select random set of genes: 50 completely at random, 50 among the highest lambda. 

Rscript ./script/selectRGenes.R LWKlambda.RData lwkset
Rscript ./script/selectRGenes.R MKKlambda.RData mkkset
Rscript ./script/selectRGenes.R YRIlambda.RData yriset

# indecies of random genes in the files lwkset, mkkset and yriset


# Prune genetic data by r2 0.2

for pop in {lwk,yri,mkk}
do
plink --bfile ./data/${pop}15_common --indep-pairwise 100 5 0.2 --out ${pop}15_common
plink --map3  --bfile ${pop}15_common --extract ${pop}15_common.prune.in --make-bed --out ./data/${pop}15_0.2p
done


###### Finding optimal lambda
# This step requires flashpca (https://github.com/gabraham/flashpca)

# manhattan.qq.Boel.R
# run.linear.altpheno.sh
## Usage: sh find.optimal.lambda.sh filename altpheno mpheno
# altpheno is the name of the gene expresion phenotype file
# mpheno is the number of the gene
# Need to have all Mitjas files on github. 

# Below script will execute the find.optimal.lambda step for each gene (30 sets of genome 
# wide linear regressions) consecutively. Recommended to submit these jobs in parallell 
# on a cluster instead. 

while read g; do
./scripts/find.optimal.lambda.linear.sh ./data/mkk_0.2p MKKpheno $g
done < mkkset

while read g; do
./scripts/find.optimal.lambda.linear.sh ./data/lwk_0.2p LWKpheno $g
done < lwkset

while read g; do
./scripts/find.optimal.lambda.linear.sh ./data/yri_0.2p YRIpheno $g
done < yriset

## Investigate results and select number of PCs. 
# Needs lambdas without correction, prefix of lambdas.combined and the gene index file. 
# Rscript optimalLambdaPlot.R initialLambda prefix geneIndex

Rscript ./scripts/optimalLambdaPlot.R MKKlambda.RData mkk_0.2p mkkset

Rscript ./scripts/optimalLambdaPlot.R LWKlambda.RData lwk_0.2p lwkset

Rscript ./scripts/optimalLambdaPlot.R YRIlambda.RData yri_0.2p yriset

# Selected number of PCs after looking at plots. We chose:
# MKK: 20
# LWK: 7
# YRI: 2

##########################################################################################
# linear regression with principal components as covariates.
# Preferably to be run in parallell.  

plink --bfile mkk15_common_Au --linear hide-covar --pheno MKKphenoQ --all-pheno --covar mkk_0.2p.pcs_30PCs.tab.22 --out ./MKKassoc/MKKpc20

plink --bfile yri15_common_Au --linear hide-covar --covar yri_0.2p.pcs_30PCs.tab.4 --pheno YRIphenoQ --all-pheno --out ./YRIassoc/YRIpc2

plink --bfile lwk15_common_Au --linear hide-covar --pheno LWKphenoQ --all-pheno --covar lwk_0.2p.pcs_30PCs.tab.9 --out ./LWKassoc/LWKpc7

##########################################################################################
# Calculate z-score files: Create files with '_z' suffix from plink assoc files. 
# Python 2.7
# Script takes a start and end index of files to process to mediate parallell execution. 

start=1
end=9085

python ./scripts/plink2z.py ./MKKassoc/MKKpc20 $start $end
python ./scripts/plink2z.py ./YRIassoc/YRIpc2 $start $end
python ./scripts/plink2z.py ./LWKassoc/LWKpc7 $start $end


##########################################################################################
# Create binary data sets (eQTL p-values and z-scores)

Ngenes=9085
Nsnps=737867


# Of p-values:
column=4

python ./scripts/merge2binary.py ./YRIassoc/YRIpc2.P .assoc.linear_z $Ngenes YRIp.dat $Nsnps $column

python ./scripts/merge2binary.py ./MKKassoc/MKKpc20.P .assoc.linear_z $Ngenes MKKp.dat $Nsnps $column

python ./scripts/merge2binary.py ./LWKassoc/LWKpc7.P .assoc.linear_z $Ngenes LWKp.dat $Nsnps $column

# Of z-scores:
column=5

python ./scripts/merge2binary.py ./YRIassoc/YRIpc2.P .assoc.linear_z $Ngenes YRIz.dat $Nsnps $column

python ./scripts/merge2binary.py ./MKKassoc/MKKpc20.P .assoc.linear_z $Ngenes MKKz.dat $Nsnps $column

python ./scripts/merge2binary.py ./LWKassoc/LWKpc7.P .assoc.linear_z $Ngenes LWKz.dat $Nsnps $column


# Add the count of number of bytes to the binary files:
./scripts/fixCounts LWKp.dat 26814086780 $Ngenes $Nsnps
./scripts/fixCounts YRIp.dat 26814086780 $Ngenes $Nsnps
./scripts/fixCounts MKKp.dat 26814086780 $Ngenes $Nsnps

./scripts/fixCounts LWKz.dat 26814086780 $Ngenes $Nsnps
./scripts/fixCounts YRIz.dat 26814086780 $Ngenes $Nsnps
./scripts/fixCounts MKKz.dat 26814086780 $Ngenes $Nsnps


# Transform binary file so that data can be accessed by gene as well as by SNP:

./scripts/transpose2 YRIp.dat YRIp_trans.dat 1 0
./scripts/transpose2 LWKp.dat LWKp_trans.dat 1 0
./scripts/transpose2 MKKp.dat MKKp_trans.dat 1 0

./scripts/transpose2 YRIz.dat YRIp_trans.dat 1 0
./scripts/transpose2 LWKz.dat LWKp_trans.dat 1 0
./scripts/transpose2 MKKz.dat MKKp_trans.dat 1 0


##########################################################################################
# Empirical CPMA p-values for each pop
# Calculate CPMA

# Calculate covariance between eQTL vectors of genes.  

Rscript ./scripts/largeCovMatrixCS.R ./scripts/API.R LWKz.dat commonGenes commonSNPs LWK_covCS

Rscript ./scripts/largeCovMatrixCS.R ./scripts/API.R YRIz.dat commonGenes commonSNPs YRI_covCS

/Rscript ./scripts/largeCovMatrixCS.R ./scripts/API.R MKKz.dat commonGenes commonSNPs MKK_covCS
 
 
# Simulate null-distribution of CPMA for each pop
# This takes a long time, should preferably be performed in parallell. 

# We need a vector of mead Z-scores for each gene:

Rscript ./scripts/meanZ.R LWKz.dat commonGenes commonSNPs LWKgeneZ
Rscript ./scripts/meanZ.R YRIz.dat commonGenes commonSNPs YRIgeneZ
Rscript ./scripts/meanZ.R MKKz.dat commonGenes commonSNPs MKKgeneZ

Nsim=1000

Rscript ./scripts/simulateCPMA_mvrnorm.R LWK_covCS.$number.RData $Nsim LWKgeneZ.RData LWKsimCPMA
Rscript ./scripts/simulateCPMA_mvrnorm.R MKK_covCS.$number.RData $Nsim MKKgeneZ.RData MKKsimCPMA
Rscript ./scripts/simulateCPMA_mvrnorm.R YRI_covCS.$number.RData $Nsim YRIgeneZ.RData YRIsimCPMA


# Will create the R objects YRIsimCPMA.RData, LWKsimCPMA.RData, MKKsimCPMA.RData



# Calculate empirical p-values

Rscript ./scripts/empPcalcExact.R YRIcpma.RData YRIsimCPMA.RData YRI

Rscript ./scripts/empPcalcExact.R LWKcpma.RData LWKsimCPMA.RData LWK

Rscript ./scripts/empPcalcExact.R MKKcpma.RData MKKsimCPMA.RData MKK

# Creates the R objects YRIempP.RData, LWKempP.RData, MKKempP.RData

##########################################################################################
# Meta analysis of CPMA stats. Weighted by the size of each study (number of individuals). 

Rscript ./scripts/metaCPMA.R

# Creates the R object METAempP.RData


##########################################################################################
# Clumping

# Merge three data sets:

for pop in {mkk,lwk,yri}
do
echo ${pop}15_common_Au >> merge.txt
done

plink --merge-list merge.txt --make-bed --out common_Au

# Clumping. No p-value threshold, want to get all independent SNPs. 
plink --bfile common_Au --clump METAempP.txt --clump-p1 1 --clump-p2 1 --clump-r2 0.2 --out METAempP

Rscript indSNPs.R
# Creates:
# clumps.RData with the independent SNPs and 
# sign.RData with the significant independent SNPs
# signN.RData with the number of overlapping targets at 0.05 level per independent significant SNP 


##########################################################################################
# Overlapping targets between the different populations

Rscript ./scripts/OLtargets.R MKK YRI sign.RData targetsYvM
Rscript ./scripts/OLtargets.R LWK YRI sign.RData targetsYvL
Rscript ./scripts/OLtargets.R MKK LWK sign.RData targetsLvM

# Creates R objects targetsYvM.RData etc


##########################################################################################
# Pairwise assessment of target overlap

# Calculate the overlap at all 178476 independent SNPs. This is divided into 12 different jobs, 
# meaning 14873 independent SNPs in each job. 


for i in {1..12}
do 
echo "library(mmap)
library(hash)
source('./scripts/API.R')
mkk.api_p=getfuncs('MKKp_trans.dat','commonSNPs','commonGenes')
yri.api_p=getfuncs('YRIp_trans.dat','commonSNPs','commonGenes')
lwk.api_p=getfuncs('LWKp_trans.dat','commonSNPs','commonGenes')

load('sign.RData')
load('clumps.RData')
load('signN.RData')

tOL_YvMn.$i=matrix(data = NA, nrow=length(clumps), ncol=length(sign))
colnames(tOL_YvMn.$i)=new1555
rownames(tOL_YvMn.$i)=clumps[($i*14873-14872):($i*14873)]

tOL_YvLn.$i=matrix(data = NA, nrow=14873, ncol=1555)
colnames(tOL_YvLn.$i)=new1555
rownames(tOL_YvLn.$i)=clumps[($i*14873-14872):($i*14873)]

tOL_LvMn.$i=matrix(data = NA, nrow=14873, ncol=1555)
colnames(tOL_LvMn.$i)=new1555
rownames(tOL_LvMn.$i)=clumps[($i*14873-14872):($i*14873)]


for (i in 1:14873) { 
Mtmp=sort( mkk.api_p\$getRow(clumps[$i*14873-14873+i]),index.return=T)\$ix
Ytmp=sort( yri.api_p\$getRow(clumps[$i*14873-14873+i]),index.return=T)\$ix
Ltmp=sort( lwk.api_p\$getRow(clumps[$i*14873-14873+i]),index.return=T)\$ix

	tOL_YvMn.$i[i,]=apply(sign1555[,c('YRI','MKK')],1, function(x) {
	return( length(intersect( Mtmp[1:x[2]],  Ytmp[1:x[1]])))
	})
	
	tOL_YvLn.$i[i,]=apply(sign1555[,c('YRI','LWK')],1, function(x) {
	return( length(intersect( Ltmp[1:x[2]],  Ytmp[1:x[1]])))
	})
	
	tOL_LvMn.$i[i,]=apply(sign1555[,c('LWK','MKK')],1, function(x) {
	return( length(intersect( Mtmp[1:x[2]],  Ltmp[1:x[1]])))
	})
	}
	
save(tOL_YvMn.$i,file='/scratch/bb447/trans_eQTL/tOL_YvM.$i.RData')
save(tOL_YvLn.$i,file='/scratch/bb447/trans_eQTL/tOL_YvLn.$i.RData')
save(tOL_LvMn.$i,file='/scratch/bb447/trans_eQTL/tOL_LvMn.$i.RData')" > tOL.$i.R
done

# Run the created scripts (These jobs should preferably be executed in parallel, not as below).

for i in {1..12}
do 
Rscript tOL.$i.R
done

# Creates 12 R objects containing the random data for each pairwise population assessment.


########################## Empirical significance of target overlap between two populations:

Rscript tOL.p.R LWK MKK targetsLvM.RData tOL_LvM
Rscript tOL.p.R LWK YRI targetsYvL.RData tOL_YvL
Rscript tOL.p.R MKK YRI targetsYvM.RData tOL_YvM

# creates the objects tOL_YvM.p.RData, tOL_YvL.p.RData and tOL_LvM.p.RData
# Not tested yet. <------- OBS!!!

##########################################################################################
# Pairwise assessment of direction overlap
# Calculate the overlap of the target genes across all 178476 independent SNPs. 
# This is divided into 12 different jobs, meaning 14873 independent SNPs in each job. 

# Create the R scripts:
for number in {1..12}
do
echo "library(mmap)
library(hash)
source('./scripts/API.R')
mkk.api_z=getfuncs('MKKz_trans.dat','commonSNPs','commonGenes')
yri.api_z=getfuncs('YRIz_trans.dat','commonSNPs','commonGenes')
lwk.api_z=getfuncs('LWKz_trans.dat','commonSNPs','commonGenes')
load('/scratch/bb447/trans_eQTL/signN.RData')
load('/scratch/bb447/trans_eQTL/clumps.RData')
load('/scratch/bb447/trans_eQTL/sign.RData')
load('targetsYvM.RData')
load('targetsYvL.RData')
load('targetsLvM.RData')

OLdir_YvMn.$number=matrix(NA,nrow=14873, ncol=length(sign))
OLdir_YvLn.$number=matrix(NA,nrow=14873, ncol=length(sign))
OLdir_LvMn.$number=matrix(NA,nrow=14873, ncol=length(sign))

for (j in 1:14873){
Mtmp=mkk.api_z\$getRow(clumps[($number*14873-14873+j)])
Ytmp=yri.api_z\$getRow(clumps[($number*14873-14873+j)])
Ltmp=lwk.api_z\$getRow(clumps[($number*14873-14873+j)])
OLdir_YvMn.$number[j,]=sapply(targetsYvM, function(x) {
if (length(x)==0) { 
return(0)
} else {
if( length( which( (Ytmp[x]<0 & Mtmp[x]<0) |  (Ytmp[x]>0 & Mtmp[x]>0) )) > length( which( (Ytmp[x]<0 & Mtmp[x]>0) | (Ytmp[x]<0 & Mtmp[x]>0) ))) { #if they have more of same then opposite direction
return(length( which( (Ytmp[x]<0 & Mtmp[x]<0) | (Ytmp[x]>0 & Mtmp[x]>0) )))
} else {
return(length( which( (Ytmp[x]<0 & Mtmp[x]>0) | (Ytmp[x]<0 & Mtmp[x]>0) )))
}}})
OLdir_YvLn.$number[j,]=sapply(targetsYvL, function(x) {
if (length(x)==0) { 
return(0)
} else {
if( length( which( (Ytmp[x]<0 & Ltmp[x]<0) | (Ytmp[x]>0 & Ltmp[x]>0) )) > length( which( (Ytmp[x]<0 & Ltmp[x]>0) | (Ytmp[x]<0 & Ltmp[x]>0) ))) { #if they have more of same then opposite direction
return(length( which( (Ytmp[x]<0 & Ltmp[x]<0) | (Ytmp[x]>0 & Ltmp[x]>0) )))
} else {
return(length( which( (Ytmp[x]<0 & Ltmp[x]>0) | (Ytmp[x]<0 & Ltmp[x]>0) )))
}}})
OLdir_LvMn.$number[j,]=sapply(targetsLvM, function(x) {
if (length(x)==0) { 
return(0)
} else {
if( length( which( (Ltmp[x]<0 & Mtmp[x]<0) | (Ltmp[x]>0 & Mtmp[x]>0) )) > length( which( (Ltmp[x]<0 & Mtmp[x]>0) | (Ltmp[x]<0 & Mtmp[x]>0) ))) { #if they have more of same then opposite direction
return(length( which( (Ltmp[x]<0 & Mtmp[x]<0) | (Ltmp[x]>0 & Mtmp[x]>0) )))
} else {
return(length( which( (Ltmp[x]<0 & Mtmp[x]>0) | (Ltmp[x]<0 & Mtmp[x]>0) )))
}}})
}

save(OLdir_LvMn.$number,file='OLdir_LvMn.$number.RData')
save(OLdir_YvMn.$number,file='OLdir_YvMn.$number.RData')
save(OLdir_YvLn.$number,file='OLdir_YvLn.$number.RData')" > OLdirn.$number.R
done		
		

# Run the R scripts:
for i in {1..12}
do
Rscript OLdirn.$i.R 
done


# Creates 12 R objects containing the random data for each pairwise population assessment.


########################## Directionality overlap between pops:

Rscript tOL.p.R LWK MKK targetsLvM.RData sign.RData OLdir_LvMn
Rscript tOL.p.R LWK YRI targetsYvL.RData sign.RData OLdir_YvLn
Rscript tOL.p.R MKK YRI targetsYvM.RData sign.RData OLdir_YvMn

# creates the R objects OLdir_LvMn.p.RData, OLdir_YvLn.p.RData, and OLdir_YvMn.p.RData containing a vector of p-values for the SNP in sign.RData


##########################################################################################
# Selection of SNPs that fulfill our requirements of an trans-eQTL <--- OBS, make into R script!

#sign=commonSNPs[which(METAempP<=0.05)][which(commonSNPs[which(METAempP<=0.05)] %in% clumps)]

load(tOL_*v*.p)
load(OLdir_*v*.p)

theFinal=sign[which(tOL_YvL.p<=0.05 & tOL_YvM.p<=0.05 & tOL_LvM.p<=0.05 & OLdir_YvL.p<=0.05 & OLdir_YvM.p<=0.05 & OLdir_LvM.p<=0.05)]

save(theFinal,file='theFinal.RData')

##########################################################################################
# Gene ontology enrichment. For each SNP run:

Rscript GObpStats.R geneAnno.RData assoc.matrix.RData snp 

# Tested, works. 


##########################################################################################
# Transcription factor binding enrichment

### get the TF data:
wget http://encodenets.gersteinlab.org/enets8.GM_proximal_filtered_network.txt


### Mean expression of each gene:

Rscript meanE.R
# creates the object mE.RData with mean intensity for each gene across the three populations.


### calculate random TF overlap. For each snp run:

Rscript TFrOverlap.R snp geneAnno.RData mE.RData assoc.matrix.RData
# generates R object TF.rOverlap.$SNP.RData with the 1000 null-values of the overlaps. 

### Calculate empirical p-values for final set of SNPs

Rscript TFoverlapP.R
# Creates R object empP.TF.RData


#### PINTS analysis not included
