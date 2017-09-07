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


symmetrize <- function(mat,lower=TRUE){

    nn <- min(ncol(mat),nrow(mat))
    if(lower){
        for(ii in 1:(nn-1)){
            for(jj in (ii+1):nn){
                mat[ii,jj] <- mat[jj,ii]
            }
        }
    }else{
        for(ii in 1:(nn-1)){
            for(jj in (ii+1):nn){
                mat[jj,ii] <- mat[ii,jj]
            }
        }
    }

    return(mat)
}

weighted.mle <- function(mmb,net,kk){

    nn <- length(mmb)
    BB <- array(0.0,c(kk,kk))
    BB.count <- array(0.0,c(kk,kk))
    for(ii in 1:nn){
        for(jj in 1:nn){
            if(!is.na(net[ii,jj])){
                BB[mmb[ii],mmb[jj]] <- BB[mmb[ii],mmb[jj]] + net[ii,jj]
                BB.count[mmb[ii],mmb[jj]] <- BB.count[mmb[ii],mmb[jj]] + 1
            }
        }
    }
    return(BB/BB.count)
}



ll.pmat <- function(net,pmat,self.ties=TRUE){
    if(!self.ties){
        diag(net) <- NA
    }
    return(sum(log(pmat)*net,na.rm=TRUE) - sum(pmat[!is.na(net)]))
}

ll.pmat.pois <- function(net,pmat,self.ties=TRUE){
    if(!self.ties){
        diag(net) <- NA
    }
    return(sum(dpois(net,pmat,log=TRUE),na.rm=TRUE))
}





##############################################################################
###########################  Plot Device Helpers  ############################
##############################################################################


size.device <- function(nn,scale=1.5){
    if(dev.cur() == 1){
        mfrow <- get.grid.dim(sqrt(nn),sqrt(nn),
                              nn)
        dev.new(height=mfrow[1]*scale,width=mfrow[2]*scale)
        return(mfrow)
    }else{
        d.size <- dev.size()
        mfrow <- get.grid.dim(d.size[1],d.size[2],nn)
    }
    return(mfrow)
}

get.grid.dim <- function(width,height,nn){

    pp <- width*height
    w2 <- width / sqrt(pp/nn); h2 <- height / sqrt(pp/nn)

    opt.set <- rbind(c(ceiling(h2),ceiling(w2)),
                     c(floor(h2), ceiling(w2)),
                     c(ceiling(h2),floor(w2)))

    opt.set <- cbind(opt.set,apply(opt.set,1,prod))
    opt.set <- opt.set[opt.set[,3] >= nn,,drop=FALSE]
    opt.set <- opt.set[opt.set[,3] == min(opt.set[,3]),,drop=FALSE]
    if(nrow(opt.set) > 1){
        if(height <= width){
            opt.set <- opt.set[opt.set[,1] <= opt.set[,2],1:2,drop=FALSE]
        }else{
            opt.set <- opt.set[opt.set[,1] > opt.set[,2],1:2,drop=FALSE]
        }
    }
    return(opt.set[1,1:2])

}


min.threshold <- function(vec,threshold=1e-5){
    return(sapply(vec,FUN=function(x) return(max(x,threshold))))
}





parse.mcmc.control <- function(mcmc.control){
    if(is.null(mcmc.control$total))
        mcmc.control$total <- 1000
    if(is.null(mcmc.control$burn.in))
        mcmc.control$burn.in <- 0
    if(is.null(mcmc.control$thin))
        mcmc.control$thin <- 1

    if(is.null(mcmc.control$extend.alpha))
        mcmc.control$extend.alpha <- 0.01
    if(is.null(mcmc.control$extend.max))
        mcmc.control$extend.max <- 10

    if(is.null(mcmc.control$label.switch.mode))
        mcmc.control$label.switch.mode <- "ad-hoc"
    if(is.null(mcmc.control$label.switch.max))
        mcmc.control$label.switch.max <- 200

    return(mcmc.control)
}

parse.init.control <- function(init.control){
    if(is.null(init.control$spectral.start))
        init.control$spectral.start <- FALSE

    if(is.null(init.control$multistart))
        init.control$multistart <- 0

    if(is.null(init.control$multistart.total))
        init.control$multistart.total <- 20

    return(init.control)
}

format.vector <- function(vec,digits=3,max.width=78,...){

    form <- format(vec,digits=digits,...)
    len <- nchar(form[1])
    width <- ceiling(log10(length(vec)+1))
    ncol <- floor((max.width - width - 3)/(len+1))

    npad <- ncol - (length(form) %% ncol)
    if(npad < ncol) form <- c(form,rep(" ",npad))


    form.mat <- matrix(form,ncol=ncol,byrow=TRUE)
    form.vec <- apply(form.mat,1,paste,collapse=" ")
    form.vec <- paste("[",format(seq(1,length(vec),by=ncol),width=width),
                      "] ",form.vec,sep="")
    form.out <- paste(form.vec,collapse="\n")
    return(form.out)

}

#' @export
network.plot.matrix <- function(x,pal=grey((50:1)/50),
                                node.order=NULL, ...){

    if(is.null(node.order)) node.order <- 1:dim(x)[1]

    adj.image.plot(mat=x,ord=node.order,pal=pal,...)

}


network.plot.array <- function(x,pal=grey((50:1)/50),
                               node.order=NULL, ...){

    if(length(dim(x)) == 2){
        message("This is a matrix.")
        network.plot.matrix(x,pal=pal,node.order=node.order,...)
    }else if(length(dim(x)) == 3){
        nn <- dim(x)[1]
        if(dim(x)[2] != nn) stop("dimensions 1 and 2 of x must be equal")
        TT <- dim(x)[3]

        mfrow <- size.device(TT,scale=5)
        par(mfrow=mfrow,mar=c(3,1,1,1))
        for(tt in 1:TT){
            network.plot(x[,,tt],pal=pal,node.order=node.order,...)
        }
    }

}


adj.image.plot <- function(mat, ord, pal,
                           xaxt="n", yaxt="n", xlab="", ylab="",
                           ...){
    nn <- nrow(mat)
    image(1:nn,1:nn,mat[ord,ord[nn:1]],col=pal,
          xaxt=xaxt,yaxt=yaxt,xlab=xlab,ylab=ylab, ...)
}



ll.chain.plot <- function(x,main="Log-Likelihood Chain",
                          xlab="Iteration",ylab="Log-Likelihood", ...){

    old.par <- par(cex.main=2,cex.lab=1.5,cex.axis=1.5)
    plot(x$chain$logLik,type="l",
         xlab=xlab,ylab=ylab,main=main,...)
    par(old.par)

}



node.diag.plot <- function(mat,
                           xlim=NULL,ylim=NULL,
                           main.prefix=NULL,main.suffix=NULL,
                           scale.ylim=FALSE,
                           nodes=NULL, ...){

    nn <- dim(mat)[1]
    total <- dim(mat)[2]

    if(is.null(ylim) && scale.ylim) ylim <- range(mat)
    if(is.null(xlim)) xlim <- c(0,total)
    if(is.null(main.prefix)) main.prefix <- "["
    if(is.null(main.suffix)) main.suffix <- "]"
    if(is.null(nodes)) nodes <- 1:nn



    mfrow <- size.device(length(nodes))

    old.par <- par(mfrow=mfrow,mar=c(3,1,1,1),
                   cex.main=1,oma=c(0,1,5,1))
    for(ii in nodes){
        plot(mat[ii,],ylim=ylim,xlim=xlim,type="l",
             xlab="",ylab="",
             main=paste0(main.prefix,ii,main.suffix),
             ...)
    }
    par(old.par)

}


gamma.diagnostic <- function(mat,
                                ...){

    total <- nrow(mat)
    alpha <- 1.5 - .5^(1/sqrt(total))

    plot(t(mat),xlab="Alpha",ylab="Beta",
         col=heat.colors(total,alpha=alpha)[total:1],...)

}




node.gamma.plot <- function(mat,
                            xlim=NULL,ylim=NULL,
                            main.prefix=NULL,main.suffix=NULL,
                            scale.xlim=FALSE, scale.ylim=FALSE,
                            nodes=NULL, ...){

    nn <- dim(mat)[1]
    total <- dim(mat)[2]

    if(is.null(ylim) && scale.ylim) ylim <- range(mat[,2,])
    if(is.null(xlim) && scale.xlim) xlim <- range(mat[,1,])
    if(is.null(main.prefix)) main.prefix <- "["
    if(is.null(main.suffix)) main.suffix <- "]"
    if(is.null(nodes)) nodes <- 1:nn



    mfrow <- size.device(length(nodes))

    old.par <- par(mfrow=mfrow,mar=c(3,1,1,1),
                   cex.main=1,oma=c(0,1,5,1))
    for(ii in nodes){
        gamma.diagnostic(mat[ii,,],
                         ylim=ylim,xlim=xlim,
                         main=paste0(main.prefix,ii,main.suffix),
                         ...)
    }
    par(old.par)

}



#####  Function for drawing from posterior predictive Gamma Distributions
#' @export
gamma.post.pred <- function(samp,iter.draws=100){
    total <- nrow(samp)
    out <- rep(NA,iter.draws*total)
    for(ii in 1:total){
        draws <- rgamma(iter.draws,samp[ii,1],samp[ii,2])
        out[((ii-1)*iter.draws + 1):(ii*iter.draws)] <- draws
    }
    return(out)
}


#####  Function for plotting posterior Gamma Distributions
#' @export
density.mat.plot <- function(mat,col="steelblue2",trans=0.3,xlim,
                             add=FALSE, ...){
    TT <- ncol(mat)
    first <- !add
    upper.d <- min(mean(mat)*10,max(mat))
    if(missing(xlim)){
        xlim <- c(0,upper.d)
    }
    for(ii in 1:TT){
        if(first){
            plot(density(mat[,ii],from=0,to=upper.d,n=1e4),
                 col=scales::alpha(col,trans),
                 xlim=xlim, ...)
            first <- FALSE
        }else{
            lines(density(mat[,ii],from=0,to=upper.d,n=1e4),
                  col=scales::alpha(col,trans), ...)
        }
    }
}





#' Generate Adjacency Matrix
#'
#' Generates an adjacency matrix given a tie intensity or tie probability
#' matrix
#'
#' @param x object of class matrix
#' @param self.ties if FALSE, the diagonals are structurally missing
#' @param d distribution for sampling
#'
#' @return Returns a matrix representing the adjacency matrix for a network
#' drawn with the given tie intensities/probabilities
#'
#' @examples
#' nn <- 50
#' pmat <- array(5,c(nn,nn))
#' net.weighted <- net.gen(pmat)
#'
#' pmat <- array(.25,c(nn,nn))
#' net.binary <- net.gen(pmat,d="b")
#'
#' @seealso \code{\link{net.gen.wsbm}}, \code{\link{net.gen.dynsbm}}
#' @export
net.gen.matrix <- function(x,self.ties=TRUE,
                           d=c("poisson","bernoulli")){
    if(dim(x)[1] != dim(x)[2])
        stop("Network intensity matrices must be symmetric")
    nn <- dim(x)[1]

    d <- match.arg(d)
    if(d == "poisson"){
        net <- array(rpois(nn*nn,x),c(nn,nn))
    }else if(d == "bernoulli"){
        net <- array(rbinom(nn*nn,1,x),c(nn,nn))
    }
    if(!(self.ties)){
        diag(net) <- NA
    }
    return(net)
}
