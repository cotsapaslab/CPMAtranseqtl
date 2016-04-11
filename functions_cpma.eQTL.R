# Functions for trans-eQTL CPMA african project

loadRData <- function(fileName){
#loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}


setTargets=function(snps,pop1.api_p,pop2.api_p) {
targets=lapply(snps, function(x) {
	tmp1=pop1.api_p$getRow(x)
	tmp2=pop2.api_p$getRow(x)
	return(which(tmp1<=0.05 & tmp2<=0.05))})
return(targets)}

getDirOLp=function(snp,pop1.api_z,pop2.api_z,Tvector, BGvector) {
	tmp1=pop1.api_z$getRow(snp)
	tmp2=pop2.api_z$getRow(snp)
	if (length( which( (tmp1[Tvector]<0 & tmp2[Tvector]<0) |  (tmp1[Tvector]>0 & tmp2[Tvector]>0) )) > length( which( (tmp1[Tvector]<0 & tmp2[Tvector]>0) | (tmp1[Tvector]<0 & tmp2[Tvector]>0) ))) { #if they have more of same then opposite direction
		ol=length( which( (tmp1[Tvector]<0 & tmp2[Tvector]<0) | (tmp1[Tvector]>0 & tmp2[Tvector]>0) ))
	} else {
		ol=length( which( (tmp1[Tvector]<0 & tmp2[Tvector]>0) | (tmp1[Tvector]<0 & tmp2[Tvector]>0) )) }
	return(length(which(BGvector>=ol))/(length(BGvector)+1)) } 


require(vioplot)
require(devtools)
require(digest)

vioplot2 <- function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL, 
                      horizontal = FALSE, col = "magenta", border = "black", lty = 1, 
                      lwd = 1, rectCol = "black", colMed = "white", pchMed = 19, 
                      at, add = FALSE, wex = 1, drawRect = TRUE, side="both") 
{
  datas <- list(x, ...)
  n <- length(datas)
  if (missing(at)) 
    at <- 1:n
  upper <- vector(mode = "numeric", length = n)
  lower <- vector(mode = "numeric", length = n)
  q1 <- vector(mode = "numeric", length = n)
  q2 <- vector(mode = "numeric", length = n)
  q3 <- vector(mode = "numeric", length = n)
  med <- vector(mode = "numeric", length = n)
  base <- vector(mode = "list", length = n)
  height <- vector(mode = "list", length = n)
  baserange <- c(Inf, -Inf)
  args <- list(display = "none")
  radj <- ifelse(side == "right", 0, 1)
  ladj <- ifelse(side == "left", 0, 1)
  if (!(is.null(h))) 
    args <- c(args, h = h)
  med.dens <- rep(NA, n)
  for (i in 1:n) {
    data <- datas[[i]]
    data.min <- min(data)
    data.max <- max(data)
    q1[i] <- quantile(data, 0.25)
    q2[i] <- quantile(data, 0.5)
    q3[i] <- quantile(data, 0.75)
    med[i] <- median(data)
    iqd <- q3[i] - q1[i]
    upper[i] <- min(q3[i] + range * iqd, data.max)
    lower[i] <- max(q1[i] - range * iqd, data.min)
    est.xlim <- c(min(lower[i], data.min), max(upper[i], 
                                               data.max))
    smout <- do.call("sm.density", c(list(data, xlim = est.xlim), 
                                     args))
    med.dat <- do.call("sm.density", 
                           c(list(data, xlim=est.xlim,
                              eval.points=med[i], display = "none")))
    med.dens[i] <- med.dat$estimate
    hscale <- 0.4/max(smout$estimate) * wex
    base[[i]] <- smout$eval.points
    height[[i]] <- smout$estimate * hscale
    med.dens[i] <- med.dens[i] * hscale
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1], t[1])
    baserange[2] <- max(baserange[2], t[2])
  }
  if (!add) {
    xlim <- if (n == 1) 
      at + c(-0.5, 0.5)
    else range(at) + min(diff(at))/2 * c(-1, 1)
    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  if (is.null(names)) {
    label <- 1:n
  }
  else {
    label <- names
  }
  boxwidth <- 0.05 * wex
  if (!add) 
    plot.new()
  if (!horizontal) {
    if (!add) {
      plot.window(xlim = xlim, ylim = ylim)
      axis(2)
      axis(1, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(x = c(at[i] - radj*height[[i]], rev(at[i] + ladj*height[[i]])), 
              y = c(base[[i]], rev(base[[i]])), 
              col = col, border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
              lty = lty)
        rect(at[i] - radj*boxwidth/2, 
             q1[i], 
             at[i] + ladj*boxwidth/2, 
             q3[i], col = rectCol)
        # median line segment
        lines(x = c(at[i] - radj*med.dens[i], 
                  at[i], 
                  at[i] + ladj*med.dens[i]),
              y = rep(med[i],3))
      }
    }
  }
  else {
    if (!add) {
      plot.window(xlim = ylim, ylim = xlim)
      axis(1)
      axis(2, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(c(base[[i]], rev(base[[i]])), 
              c(at[i] - radj*height[[i]], rev(at[i] + ladj*height[[i]])), 
              col = col, border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, 
              lty = lty)
        rect(q1[i], at[i] - radj*boxwidth/2, q3[i], at[i] + 
               ladj*boxwidth/2, col = rectCol)
        lines(y = c(at[i] - radj*med.dens[i], 
                    at[i], 
                    at[i] + ladj*med.dens[i]),
              x = rep(med[i],3))
      }
    }
  }
  invisible(list(upper = upper, lower = lower, median = med, 
                 q1 = q1, q3 = q3))
}



getPermProb = function(meanExp, obs.set) {
dObs=density(meanExp[obs.set])
Prob <- sapply(meanExp,function(x) approx(dObs$x,dObs$y,x)$y)
dBG <- density(meanExp)
Prob <- apply(as.matrix(cbind(meanExp,Prob),ncol=2),1,function(x){ x[2] / approx(dBG$x,dBG$y,x[1])$y } )
return(Prob) }

getTFoverlapRD <- function(x) {
file=paste('TF.rOverlap.',as.character(x),'.RData',sep='')
load(file) 
return(get(paste('TF.rOverlap.',as.character(x),sep='')))
}




