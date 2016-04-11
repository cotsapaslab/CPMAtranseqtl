# Rscript OLtargets.R pop1 pop2 sign.RData out
# Script that identifies the targets that overlap at one trans-eQTL between two populations, and creates a list of those. 

library(mmap)
library(hash)
source('./scripts/API.R')
source('./scripts/functions_cpma.eQTL.R')

args <- commandArgs(TRUE)
pop1=args[1]
pop2=args[2]
sign=loadRData(args[3])
out=args[4]


pop1.api_p=getfuncs(paste(pop1,'p_trans.dat',sep=''),'commonSNPs','commonGenes')
pop2.api_p=getfuncs(paste(pop2,'p_trans.dat',sep=''),'commonSNPs','commonGenes')

targets=setTargets(sign,pop1.api_p,pop2.api_p)

assign(out,targets)
save(list=out,file=paste(out,'RData',sep='.'))
