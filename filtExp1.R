# R script to:
# 1. filters gene expression data set according to mean interquantile range and mean intensity 
# 2. quantile normalize
# 3. outputs plink phenotype-file. 
# Rscript ./scripts/filtExp1.R geneExpressionFile plinkFamFile phenoFileName
# by Boel Brynedal

args<-commandArgs(TRUE)
source('./scripts/filtFun.R')
library(mclust)
library(limma)

Exp=read.table(args[1],as.is=T,header=T,sep='\t',row.names=1)


###### filtering 

IQR=filtIQR(Exp,0.8)
E=filtE(Exp,0.8)

filt=reduceFilt(IQR,E)
names(filt)=rownames(Exp)

Exp=Exp[filt,]


######### Quantile normalize the expression data post filtering

Exp=t(normalize.quantiles(Exp)


######### Output plink-pheno files

fam=read.table(args[2],as.is=T)


write.table(cbind(fam[which(fam[,2] %in% colnames(Exp)),1:2],t(Exp[,fam[which(fam[,2] %in% colnames(Exp)),2]])),quote=F,row.names=F,col.names=F,file=args[3])

