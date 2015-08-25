#####  helper.R
#####  Functions Implementing SBM and MMSBM Fitting in C:
#####    sbm - Fits SBM Model
#####    mmsbm.C - Fits MMSBM Model
#####    mmsbm - Fits MMSBM choosing intelligent starting values
#####
##############################################################

max.acf <- function(flatTable,kk.max=100,make.plot=TRUE){
  cols <- ncol(flatTable)
  rows <- nrow(flatTable)
  out <- array(NA,c(kk.max+1,cols))

  for(ii in 1:cols){
    out[,ii] <- as.vector(acf(flatTable[,ii],lag.max=kk.max,plot=F)$acf)
  }

  out2 <- apply(out,1,max)
  if(make.plot){
    plot(out2,type="h")
    abline(h=.2,col="blue",lwd=2)
  }

  return(out2)

}


plot.ll <- function(obj,add=FALSE){
  if(add){
    lines(obj$chain$logLik,type="l")
  }else{
    plot(obj$chain$logLik,type="l")
  }
}


rdirichlet <- function(n,alpha){
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  x/as.vector(sm)
}

