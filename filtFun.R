# R-functions used when filtering normalized expression data
# by Boel Brynedal

filtIQR=function(Exp,prob1) {
iqr=apply(Exp,1, function(x) IQR(x, na.rm = TRUE))
mmIQR <- Mclust(iqr, G=2)
filt=vector('logical',length=nrow(Exp))
filt[which(mmIQR$z[,2]>=prob1)]=T
return(filt) }

filtE=function(Exp,prob2) {
mE=apply(Exp,1, mean)
mmE <- Mclust(mE, G=2)
filt=vector('logical',length=nrow(Exp))
filt[which(mmE$z[,2]>=prob2)]=T
return(filt) }

reduceFilt=function(filt1,filt2) {
return(Reduce('|',list(filt1,filt2))) }







