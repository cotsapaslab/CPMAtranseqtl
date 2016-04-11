# Script that outputs two things: 
# 1. R object with the independent index SNPs from clumping 
# 2. R object with the independent index SNPs from clumping that have a METAempP <= 0.05

clumps=read.table('METAempP.clumped',as.is=T,header=T)[,3]
save(clumps,file='clumps.RData')

commonSNPs=read.table('commonSNPs',as.is=T)[,1]
load('METAempP.RData')

sign=commonSNPs[which(METAempP<=0.05)][which(commonSNPs[which(METAempP<=0.05)] %in% clumps)]
save(sign,file='sign.RData'

library(mmap)
library(hash)
source('./scripts/API.R')

mkk.api_p=getfuncs('MKKp_trans.dat','commonSNPs','commonGenes')
yri.api_p=getfuncs('YRIp_trans.dat','commonSNPs','commonGenes')
lwk.api_p=getfuncs('LWKp_trans.dat','commonSNPs','commonGenes')


signN=sapply(sign,function(x) {
Ytmp=yri.api_p$getRow(x)
Mtmp=mkk.api_p$getRow(x)
Ltmp=lwk.api_p$getRow(x)
return(c(length(which(Ytmp<=0.05)),length(which(Mtmp<=0.05)),length(which(Ltmp<=0.05)))) })

signN=t(signN)
colnames(signN)=c('YRI','MKK','LWK')

save(signN,file='signN.RData')
