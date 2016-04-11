# R script to:
# 3. filters according to mean interquantile range and mean intensity in each pop separately
# 4. takes the intersect
# 5. outputs a wenn-diagram of filtering and overlap 
# 6. quantile normalize
# 7. outputs plink phenotype-files. 
# Rscript ./scripts/filtExp3.R
# by Boel Brynedal

#args<-commandArgs(TRUE)
source('./scripts/filtFun.R')
library(mclust)
library(limma)

MKKexp=read.table('./data/MKK_expression.txt',as.is=T,header=T,sep='\t',row.names=1)
YRIexp=read.table('./data/YRI_expression.txt',as.is=T,header=T,sep='\t',row.names=1)
LWKexp=read.table('./data/LWK_expression.txt',as.is=T,header=T,sep='\t',row.names=1)


###### filtering pop1 (YRI)

yIQR=filtIQR(YRIexp,0.8)
yE=filtE(YRIexp,0.8)

filtYRI=reduceFilt(yIQR,yE)
names(filtYRI)=rownames(YRIexp)

###### filtering pop2 (LWK)

lIQR=filtIQR(LWKexp,0.8)
lE=filtE(LWKexp,0.8)

filtLWK=reduceFilt(lIQR,lE)
names(filtLWK)=rownames(LWKexp)


###### filtering pop3 (MKK)

mIQR=filtIQR(MKKexp,0.8)
mE=filtE(MKKexp,0.8)

filtMKK=reduceFilt(mIQR,mE)
names(filtMKK)=rownames(MKKexp)


######## Produce Venn diagram of filtering overlaps

vc <- vennCounts(cbind(filtMKK[names(filtYRI)],filtYRI,filtLWK[names(filtYRI)]))
vennDiagram(vc)

today <- Sys.Date()

jpeg(paste('VennDiagram_filter',format(today, format="%Y%m%d"),'.jpg',sep='')
vennDiagram(vc)
dev.off()

# filter all three pops according to above filters
MKKexp=MKKexp[names(filtYRI)[which(filtMKK[names(filtYRI)] & filtYRI & filtLWK[names(filtYRI)])],]
YRIexp=YRIexp[which(filtMKK[names(filtYRI)] & filtYRI & filtLWK[names(filtYRI)]),]
LWKexp=LWKexp[names(filtYRI)[which(filtMKK[names(filtYRI)] & filtYRI & filtLWK[names(filtYRI)])],]

######### Quantile normalize the expression data post filtering

MKKexp=t(normalize.quantiles(MKKexp)
LWKexp=t(normalize.quantiles(LWKexp)
YRIexp=t(normalize.quantiles(YRIexp)



######### Output plink-pheno files

lwkfam=read.table('lwk15_common.fam',as.is=T)
yrifam=read.table('yri15_common.fam',as.is=T)
mkkfam=read.table('mkk15_common.fam',as.is=T)


write.table(cbind(yrifam[which(yrifam[,2] %in% colnames(YRIexp)),1:2],t(YRIexp[,yrifam[which(yrifam[,2] %in% colnames(YRIexp)),2]])),quote=F,row.names=F,col.names=F,file='YRIpheno')


write.table(cbind(lwkfam[which(lwkfam[,2] %in% colnames(LWKexp)),1:2],t(LWKexp[,lwkfam[which(lwkfam[,2] %in% colnames(LWKexp)),2]])),quote=F,row.names=F,col.names=F,file='LWKpheno')


write.table(cbind(mkkfam[which(mkkfam[,2] %in% colnames(MKKexp)),1:2],t(MKKexp[,mkkfam[which(mkkfam[,2] %in% colnames(MKKexp)),2]])),quote=F,row.names=F,col.names=F,file='MKKpheno')

