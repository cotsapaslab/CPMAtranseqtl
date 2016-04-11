# Script to collect common SNPs from three (3) plink bim-files, creates output file 'commonSNPs'
# by Boel Brynedal

args<-commandArgs(TRUE)

pop1=read.table(args[1],as.is=T,header=F)
pop2=read.table(args[2],as.is=T,header=F)
pop3=read.table(args[3],as.is=T,header=F)

write.table(intersect(intersect(pop1[,2],pop2[,2]),pop3[,2]), file='commonSNPs',quote=F,col.names=F,row.names=F)
