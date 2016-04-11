# Script to calculate the mean expression across multiple data sets. 
# Rscript meanE.R

MKKphenoQ=read.table('MKKpheno',as.is=T)
LWKphenoQ=read.table('LWKpheno',as.is=T)
YRIphenoQ=read.table('YRIpheno',as.is=T)


mE.Y=apply(YRIpheno[,3:ncol(YRIpheno)],2,mean) # first two columns are FID and IID. 
mE.L=apply(LWKpheno[,3:ncol(LWKpheno)],2,mean)
mE.M=apply(MKKpheno[,3:ncol(MKKpheno)],2,mean)

mE=apply(cbind(mE.Y,mE.L,mE.M),1,mean)

save(mE,file='mE.RData')
