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


###########################################################################
####################     TIME TABLE CONVERSIONS     #######################
###########################################################################



convert.time.table <- function(tab,cutoffs,nn,labs,sorted=FALSE){

    t0 <- proc.time()

    TT <- length(cutoffs) - 1

    if(missing(nn)){
        if(missing(labs)){
            labs <- unique(tab[,1])
        }
        nn <- length(labs)
    }
    labs <- sort(labs)

    inv <- function(x){
        return(which(labs == x))
    }
    cutoffs <- c(cutoffs,Inf)

    AA.mat <- array(0,c(nn,nn,TT))
    if(!sorted){
        if(all(labs == 1:nn)){
            for(ii in 1:nrow(tab)){
                tt <- 0
                while(tab[ii,3] > cutoffs[tt+1]) tt <- tt + 1
                if(tt > 0 && tt <= TT){
                    AA.mat[tab[ii,1],tab[ii,2],tt] <- AA.mat[tab[ii,1],tab[ii,2],tt] + 1
                }
            }
        }else{
            for(ii in 1:nrow(tab)){
                if(ii %% 1e5 == 0){
                    message("Iteration: ",ii,", Time Elapsed: ",(proc.time() - t0)[3])
                }
                tt <- 0
                while(tab[ii,3] > cutoffs[tt+1]) tt <- tt + 1
                if(tt > 0 && tt <= TT){
                    ss <- inv(tab[ii,1])
                    rr <- inv(tab[ii,2])

                    AA.mat[ss,rr,tt] <- AA.mat[ss,rr,tt] + 1
                }
            }
        }
    }else{
        tt <- 0
        for(ii in 1:nrow(tab)){
            if(ii %% 1e5 == 0){
                message("Iteration: ",ii,", Time Elapsed: ",(proc.time() - t0)[3])
            }

            if(tab[ii,3] > cutoffs[tt+1]) {
                tt <- tt + 1
                if(tt > TT) break
            }
            ss <- inv(tab[ii,1])
            rr <- inv(tab[ii,2])

            AA.mat[ss,rr,tt] <- AA.mat[ss,rr,tt] + 1
        }
    }
    dimnames(AA.mat) <- list(sender=labs,receiver=labs,time=1:TT)
    return(AA.mat)
}







#####  Assumes day cutoffs are integers (so with base funcitonality 0,5,7 is Monday through Friday, then Sunday)  monday is zero indexed, so if the first monday is on the 3rd, set monday to 2.
day.time.cutoffs <- function(days,time.cutoffs,week.cutoffs,
                             monday=0){

    day.secs <- (0:days) * 24 * 3600
    cutoffs <- rowSums(expand.grid(day.secs,time.cutoffs))
    min.time <- 0
    max.time <- days*24*3600

    valid.co <- (cutoffs >= min.time) & (cutoffs <= max.time)
    final.cutoffs <- sort(cutoffs[valid.co])

    TT <- length(final.cutoffs) - 1
    if(missing(week.cutoffs)) week.cutoffs <- 0
    ll <- length(time.cutoffs)
    ee <- ll * length(week.cutoffs)
    time.map <- rep(NA,TT)
    day.class <- 0
    dow.vec <- sod.vec <- rep(NA,TT)

    for(ii in 1:TT){
        sod <- final.cutoffs[ii] %% 86400
        day <- final.cutoffs[ii] %/% 86400
        dow <- (day - monday) %% 7

        sod.vec[ii] <- sod
        dow.vec[ii] <- dow

        if(any(dow == week.cutoffs)){
            day.class <- which(dow == week.cutoffs) %% length(week.cutoffs)
        }

        time.map[ii] <- which(time.cutoffs == sod) + ll * day.class

    }
    hours <- as.numeric(diff(c(time.cutoffs,time.cutoffs[1]+86400)))/3600


    return(list(cutoffs=final.cutoffs,time.map=c(0,time.map),
                dow=c(0,dow.vec),sod=c(time.cutoffs[1],sod.vec),
                ee=ee,TT=TT,hours.vec=hours))
}



read.and.convert.time.table <- function(fn,start,duration,
                                        time.cutoffs,
                                        week.cutoffs,
                                        sender.col="pickup_cd",
                                        receiver.col="dropoff_cd",
                                        time.col="pickup_datetime",
                                        time.format="%Y-%m-%d %H:%M:%S",
                                        max.rows=5e7,monday=0,
                                        labs,sorted=TRUE){
    t0 <- proc.time()

    start.time <- strptime(start,time.format,"UTC")
    time.from.start <- function(str){
        str.time <- strptime(str,time.format,"UTC")
        return(difftime(str.time,start.time,units="secs"))
    }

    message((proc.time() - t0)[3],":  ","Compiling Cutoffs...")

    cutoff.info <- day.time.cutoffs(days=ceiling(duration / 86400),
                                    time.cutoffs=time.cutoffs,
                                    week.cutoffs=week.cutoffs,
                                    monday=monday)

    message((proc.time() - t0)[3],":  ","Reading Data File...")
    tab <- raw.data <- read.csv(gzfile(fn),nrows = 10,skip=0,
                                stringsAsFactors = FALSE)
    classes <- sapply(raw.data,class)
    classes[] <- "NULL"
    classes[sender.col] <- "numeric"
    classes[receiver.col] <- "numeric"
    classes[time.col] <- "character"

    tab <- read.csv(gzfile(fn),nrows = max.rows,skip=0,
                    colClasses=classes)

    message("Data File has ",nrow(tab)," rows")

    message((proc.time() - t0)[3],":  ","Converting to Time from Start...")
    tab <- tab[,c(sender.col,receiver.col,time.col)]
    tab[,3] <- sapply(tab[,3],time.from.start)
    if(missing(labs)) labs <- unique(tab[,1])

    if(sorted){
        message((proc.time() - t0)[3],":  ","Sorting data.frame...")
        tab <- tab[order(tab[,3]),]
    }

    message((proc.time() - t0)[3],":  ","Converting to Multimat Format...")
    AA.total <- convert.time.table(tab=tab,labs=labs,
                                   cutoff.info$cutoffs,sorted=sorted)
    message((proc.time() - t0)[3],":  ","Complete.")
    return(list(AA=AA.total,
                time.map=cutoff.info$time.map[-1],
                hours.vec=cutoff.info$hours.vec,
                cutoff.info=cutoff.info))

}



convert.to.edgelist <- function(AA){
    nn <- nrow(AA)
    el <- array(NA,c(nn*nn,3))
    rn <- 0
    for(ii in 1:nn){
        for(jj in 1:nn){
            rn <- rn + 1
            el[rn,] <- c(ii,jj,AA[ii,jj])
        }
    }
    return(el)
}
