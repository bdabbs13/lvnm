#####  helper.R
#####  Functions Implementing SBM and MMSBM Fitting in C:
#####    sbm - Fits SBM Model
#####    mmsbm.C - Fits MMSBM Model
#####    mmsbm - Fits MMSBM choosing intelligent starting values
#####
##############################################################

get.weighted.adj <- function(ig){
  elm <- igraph::get.edgelist(ig)
  elm <- cbind(elm,igraph::E(ig)$count)
  aa <- as.matrix(igraph::get.adjacency(ig))
  for(ii in 1:nrow(elm)){
    aa[elm[ii,1],elm[ii,2]] <- elm[ii,3]
  }
  return(aa)
}

import.weighted.gml <- function(fn){
    if(!grepl(".gml",fn)){
        stop("Filename does not have extension .gml")
    }
    ig <- igraph::read_graph(fn,format="gml")
    return(get.weighted.adj(ig))
}

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

