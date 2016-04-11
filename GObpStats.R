# Script that calculates gene ontology biological pathways enrichments within gene sets. 
# Rscript GObpStats.R geneAnno.RData assoc.matrix.RData snp 

library(GOstats)
library(org.Hs.eg.db)
source('./scripts/functions_cpma.eQTL.R')

args <- commandArgs(TRUE)
geneAnno=loadRData(args[1])
assoc.matrix=loadRData(args[2])
snp=args[3]

universe=unique(as.character(geneAnno[,'HGNC']))
univmap <- select(org.Hs.eg.db, universe, "ENTREZID", "SYMBOL")

genset=unique(as.character(geneAnno[,'HGNC'][which(assoc.matrix[snp,]==1)]))
genemap <- select(org.Hs.eg.db, genset, "ENTREZID", "SYMBOL")
  
params <- new("GOHyperGParams", geneIds=unique(genemap[which(!is.na(genemap[,2])),2]), universeGeneIds=unique(univmap[which(!is.na(univmap[,2])),2]), ontology="BP", pvalueCutoff=0.05, conditional=TRUE, testDirection="over", annotation='org.Hs.eg.db')

hgOver <- hyperGTest(params)

assign(paste(snp,'_GO',sep=''), summary(hgOver,categorySize=10))
save(list=paste(snp,'_GO',sep=''),file=paste(snp,'_GO.RData',sep=''))
