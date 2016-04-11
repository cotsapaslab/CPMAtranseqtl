# Script to select 100 random genes for a population:
# 1. 50 random among all genes
# 2. 50 random among those within the highest 5% of lambda
# Rscript selectRGenes.R lambda.RData outset

args<-commandArgs(TRUE)
lambda=get(load(args[1]))
 
set=unique(c(sample(1:length(lambda),50),sample(which(lambda>quantile(lambda,0.95)),100)))[1:100]

# write files of indecies

write.table(set,quote=F,row.names=F,col.names=F,file=args[2])
