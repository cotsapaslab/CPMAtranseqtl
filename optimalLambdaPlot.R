# Rscript optimalLambdaPlot.R initialLambda prefix geneIndex
# Creates a box plot of optimal lambda results.
# by Boel Brynedal 

### Define functions

lget=function(x) {
tmp=load(x)
return(get(tmp)) }

llorder=function(prefix,idx) {
lambdas=read.table(paste(prefix,'.',idx,'.lambdas.combined',sep=''),as.is=T)
lambdas=cbind(lambdas,sapply(lambdas[,3],function(x) as.numeric(strsplit(x,paste(idx,'.',sep=''))[[1]][2])))
lambdas=lambdas[order(lambdas[,4]),]
return(lambdas[,1]) }

args<-commandArgs(TRUE)

genes=read.table(args[3])[,1]

Lambdas=llorder(args[2],genes[1])

for (i in 2:100){
Lambdas=cbind(Lambdas,llorder(args[2],genes[i]))}

rownames(Lambdas)=1:30
colnames(Lambdas)=genes

# Include original lambda without PC correction
start=lget(args[1])
start=start[genes]

Lambdas=rbind(start,Lambdas)

# Create box plot of lambdas for all genes among random or 95%-random sets when using 
# 0-30 PCs.

today <- Sys.Date()

print(paste('creating figure: ',paste(args[2],'optLambda',format(today, format="%Y%m%d"),'.jpg',sep=''),sep=''))

jpeg(paste(args[2],'optLambda',format(today, format="%Y%m%d"),'.jpg',sep=''),height=500, width=1000)
par(mfrow=c(1,2))
boxplot(t(Lambdas[,1:50]),main='50 random genes')
boxplot(t(Lambdas[,51:100]),main='50 random genes among 95% lambda')
dev.off()
