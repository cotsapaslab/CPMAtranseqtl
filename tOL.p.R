# Script that calculates the significance of target overlap between two populations
# Rscript tOL.p.R pop1 pop2 targets.RData randomPrefix
# Where randomPrefix is the prefix used when creating the random overlaps at all independent SNPs

args <- commandArgs(TRUE)
pop1=args[1]
pop2=args[2]
randomPrefix=args[4]

library(mmap)
library(hash)
source('./scripts/API.R')
source('./scripts/functions_cpma.eQTL.R')


pop1.api_z=getfuncs(paste(pop1,'z_trans.dat',sep=''),'commonSNPs','commonGenes')
pop2.api_z=getfuncs(paste(pop2,'z_trans.dat',sep=''),'commonSNPs','commonGenes')

targets=loadRData(args[3])

bg=loadRData(paste(randomPrefix,'.1.RData',sep=''))

for (i in 2:12){
file=paste(randomPrefix,i,'RData',sep='.')
new=loadRData(file)
bg=rbind(bg,new)
}

p=vector('numeric',length=length(targets))
for (i in 1:length(targets)){
p[i]=length(which(bg[,i]>=length(targets[[i]])))/(length(clump)+1) }

assign(paste(randomPrefix,'.p',sep=''), p)

save(list=paste(randomPrefix,'.p',sep=''),file=paste(randomPrefix,'.p.RData',sep=''))
