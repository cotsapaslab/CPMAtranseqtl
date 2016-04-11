# Calculate empirical p-values for overlap with TF gene sets. 
# Rscript TFoverlapP.R

source('./scripts/functions_cpma.eQTL.R')
load('theFinal.RData')

GM.filtered=read.table('enets8.GM_proximal_filtered_network.txt',sep='\t',stringsAsFactors=F,header=F,strip.white=T)
TFgeneset=unique(intersect(geneAnno[,'HGNC'],union(GM.filtered[,1],GM.filtered[,2])))
TFs=unique(GM.filtered[,1])

empP.TF=matrix(NA,ncol=length(TFs),nrow=length(theFinal))
rownames(empP.TF)=theFinal
colnames(empP.TF)=TFs

for (i in 1:length(theFinal)){
r.matrix=getTFoverlapRD(theFinal[i])
observed=apply(as.matrix(TFs,ncol=1),1, function(x) length(intersect(intersect(geneAnno[which(assoc.matrix[snp,]==1),2],TFgeneset),GM.filtered[which(GM.filtered[,1]==x),2])))
for (i in 1:ncol(r.matrix)){
empP.TF[snp,i]=(length(which(r.matrix[,i]>=observed[i]))+1)/(nrow(r.matrix)+1) }}

save(empP.TF,file='empP.TF.RData')
