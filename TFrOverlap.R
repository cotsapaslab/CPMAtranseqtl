# Script to permute overlap betweens TF gene sets and 1000 random gene sets with the same mean intensity distribution as the index SNPs target.
# Rscript TFrOverlap.R SNP geneAnno.RData mE.RData assoc.matrix.RData

source('./scripts/functions_cpma.eQTL.R')

args <- commandArgs(TRUE)
SNP=args[1]
geneAnno=loadRData(args[2])
mE=loadRData(args[3])
assoc.matrix=loadRData(args[4])

load('geneAnno.RData') 
load('mE.RData') 
load('assoc.matrix.RData')

p=getPermProb(mE,which(assoc.matrix[$SNP,]==1))
p[which(is.na(p))]=min(p,na.rm=T)
GM.filtered=read.table('enets8.GM_proximal_filtered_network.txt',sep='\t',stringsAsFactors=F,header=F,strip.white=T)
TFgeneset=unique(intersect(geneAnno[,'HGNC'],union(GM.filtered[,1],GM.filtered[,2])))
TFs=unique(GM.filtered[,1])

TF.rOverlap=matrix(NA,nrow=1000,ncol=length(TFs))
colnames(TF.rOverlap.$SNP)=TFs
for (i in 1:1000) { # for each iteration I want the overlap between each TF and random gene set -> 50*1000 matrix
idx=sample(1:9085,length(which(assoc.matrix[$SNP,]==1)),prob=p)
TF.rOverlap[i,]=apply(as.matrix(TFs,ncol=1),1, function(x) length(intersect(intersect(geneAnno[idx,2],TFgeneset),GM.filtered[which(GM.filtered[,1]==x),2])))
}

assign(paste('TF.rOverlap.',SNP,sep=''),TF.rOverlap)

save(list=paste('TF.rOverlap.',SNP,sep=''),file=paste('TF.rOverlap.',SNP,'.RData',sep=''))

