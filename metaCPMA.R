# Script that calculates meta p-values from three sets of p-values, weighted by study size.
# Needs p-value objects and plink fam files

############# 	0. Load Data

load('YRIempP.RData')
load('LWKempP.RData')
load('MKKempP.RData')


############# 	1. Convert empirical p-values to z-scores for all pops.

YRIempZ=sapply(YRIempP, function(x) return(qnorm(1 - x)))
MKKempZ=sapply(MKKempP, function(x) return(qnorm(1 - x)))
LWKempZ=sapply(LWKempP, function(x) return(qnorm(1 - x)))

#############    2. Calculate meta-z

Nyri=nrow(read.table('yri15_common_Au.fam'))
Nmkk=nrow(read.table('mkk15_common_Au.fam'))
Nlwk=nrow(read.table('lwk15_common_Au.fam'))
Ntot=Nyri+Nmkk+Nlwk

 
w_i=c(sqrt(Nyri/Ntot),sqrt(Nmkk/Ntot),sqrt(Nlwk/Ntot))


METAempZ=apply(cbind(YRIempZ,MKKempZ,LWKempZ),1, function(x) {
metaz = (x[1]*w_i[1] + x[2]*w_i[2] + x[3]*w_i[3])
return(metaz)})


############# 	3. Convert to p-value for the meta analysis.

SNPs=read.table('commonSNPs',stringsAsFactors=F)[,1]

METAempP=pnorm(-METAempZ)
names(METAempP)=SNPs

save(METAempP,file='METAempP.RData')
