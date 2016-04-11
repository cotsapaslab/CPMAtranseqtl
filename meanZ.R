args <- commandArgs(TRUE)
zScores=args[1]
Genes=args[2]
SNPs=args[3]
outPrefix=args[4]

library(mmap)
library(hash)
source('./scripts/API.R')

api_gene=getfuncs(zScores,Genes,SNPs)

assign(outPrefix,sapply(1:length(Genes), function(x) {
return(mean(api_gene$getRow(Genes[x]))) }))

save(list=outPrefix,file=paste(outPrefix,'.RData',sep=''))
