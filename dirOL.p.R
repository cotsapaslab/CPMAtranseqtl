# Script that calculates the significance of target overlap between two populations
# Rscript tOL.p.R pop1 pop2 targets.RData sign.RData randomPrefix
# Where randomPrefix is the prefix used when creating the random overlaps at all independent SNPs

args <- commandArgs(TRUE)
pop1=args[1]
pop2=args[2]
randomPrefix=args[5]

library(mmap)
library(hash)
source('./scripts/API.R')
source('./scripts/functions_cpma.eQTL.R')


pop1.api_z=getfuncs(paste(pop1,'z_trans.dat',sep=''),'commonSNPs','commonGenes')
pop2.api_z=getfuncs(paste(pop2,'z_trans.dat',sep=''),'commonSNPs','commonGenes')

targets=loadRData(args[3])
sign=loadRData(args[4])


bg=loadRData(paste(randomPrefix,'.1.RData',sep=''))

for (i in 2:12){
file=paste(randomPrefix,i,'RData',sep='.')
new=loadRData(file)
bg=rbind(bg,new)
}

p=sapply(1:length(sign), function(i) {
getDirOLp(sign[i],pop1.api_z,pop2.api_z,targets[[i]],bg[,i]) })

assign(paste(randomPrefix,'.p',sep=''), p)

save(list=paste(randomPrefix,'.p',sep=''),file=paste(randomPrefix,'.p.RData',sep=''))
