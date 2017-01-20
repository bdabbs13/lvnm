 #####  sbm.R
#####  Functions Implementing SBM and MMSBM Fitting in C:
#####    sbm - Fits SBM Model
#####    mmsbm.C - Fits MMSBM Model
#####    mmsbm - Fits MMSBM choosing intelligent starting values
#####
##############################################################
                                        #library(svd)
#####  Function to generate data from WSBM model

####  Create some good examples for data generation
rgamma.prior <- function(n=100,p=0.99,q=1.01,r=1,s=1,
                         thin=100,burn.in=100,
                         alpha.init = 1, beta.init = 1, prop.sd=1.0){

    ##  Allocating Memory for C++
    alpha.vec <- beta.vec <- rep(0,n)
##    print(prop.sd)
    out <- .C("RGammaPrior",
              as.integer(n),as.integer(thin),as.integer(burn.in),
              as.double(alpha.init),as.double(beta.init),
              as.double(p),as.double(q),as.double(r),as.double(s),
              as.double(prop.sd),as.double(alpha.vec),as.double(beta.vec))

    return(cbind(out[[11]],out[[12]]))
}




dgamma.prior <- function(alpha,beta,p=1,q=2,r=2,s=2){
    l.num <- (alpha-1)*log(p) - beta*q + alpha*s*log(beta)
    l.den <- r*lgamma(alpha)

    return(exp(l.num - l.den))
}



data.gen.dynsbm <- function(nn=40,TT=4,tmap=rep(1:2,2),hours.vec=c(1,1),
                            BB.prior,SS.prior,RR.prior,
                            mmb.prior=c(.5,.5),mmb,
                            self.ties=TRUE){


    ## Loading and checking mmb
    if(missing(mmb)){
        mmb <- sort(sample(kk,nn,replace=TRUE,prob=mmb.prior))
        kk <- length(mmb.prior)
    }else{
        if(length(mmb) != nn) stop("mmb must have length n")
        kk <- max(mmb)
    }

    ##  Checking dimension of Priors
    ee <- max(tmap)
    if(any(dim(BB.prior) != c(kk,kk,2,ee))){
        stop("BB.prior must have dimension (kk,kk,2,ee)")
    }
    if(any(dim(SS.prior) != c(nn,2,ee))){
        stop("SS.prior must have dimension (nn,2,ee)")
    }
    if(any(dim(RR.prior) != c(nn,2,ee))){
        stop("RR.prior must have dimension (nn,2,ee)")
    }

    ##  Creating containers for network and parameters
    AA.mat <- pmat.mat <- array(NA,c(nn,nn,TT))
    BB.mat <- array(NA,c(kk,kk,TT))
    SS.mat <- array(NA,c(nn,TT))
    RR.mat <- array(NA,c(nn,TT))

    for(tt in 1:TT){
        BB.mat[,,tt] <- rgamma(kk^2,BB.prior[,,1,tmap[tt]],
                               BB.prior[,,2,tmap[tt]])
        SS.mat[,tt] <- rgamma(nn,SS.prior[,1,tmap[tt]],SS.prior[,2,tmap[tt]])
        RR.mat[,tt] <- rgamma(nn,RR.prior[,1,tmap[tt]],RR.prior[,2,tmap[tt]])
        hours <- hours.vec[tmap[tt]]
        gen.single <- data.gen.wsbm(nn=nn,BB=BB.mat[,,tt],mmb=mmb,
                                    SS=SS.mat[,tt],RR=RR.mat[,tt],
                                    self.ties=self.ties,hours=hours)
        AA.mat[,,tt] <- gen.single$net
        pmat.mat[,,tt] <- gen.single$pmat
    }

    return(list(net.mat=AA.mat,hours=hours.vec,
                mmb=mmb,BB=BB.mat,SS=SS.mat,RR=RR.mat,
                BB.prior=BB.prior,SS.prior=SS.prior,RR.prior=RR.prior,
                pmat.mat=pmat.mat,pmat=apply(pmat.mat,c(1,2),mean)))
}



convert.time.table <- function(tab,cutoffs,nn,labs){
    TT <- length(cutoffs) - 1
    if(missing(nn)){
        if(missing(labs)){
            labs <- unique(tab[,1])
        }
        nn <- length(labs)
    }

    AA.mat <- array(NA,c(nn,nn,TT))
    for(ii in 1:nrow(tab)){
        tt <- 0
        while(tab[ii,3] > cutoffs[tt+1]) tt <- tt + 1
        if(tt > 0 && tt <= TT){
            AA.mat[tab[ii,1],tab[ii,2],tt] <-AA.mat[tab[ii,1],tab[ii,2],tt] + 1
        }
    }
    return(AA.mat)
}



dynsbm.priors <- function(BB.hyperprior=c(.9,1.1,1,1),
                          SS.hyperprior=c(.9,1.1,1,1),
                          RR.hyperprior=c(.9,1.1,1,1)){

    if(any(BB.hyperprior <= 0)){
        stop("Hyperprior for BB must have positive componenets")
    }
    if(any(SS.hyperprior <= 0)){
        stop("Hyperprior for SS must have positive componenets")
    }
    if(any(RR.hyperprior <= 0)){
        stop("Hyperprior for RR must have positive componenets")
    }


    return(list(BB=BB.hyperprior,
                SS=BB.hyperprior,
                RR=BB.hyperprior))
}


min.threshold <- function(vec,threshold=1e-5){
    return(sapply(vec,FUN=function(x) return(max(x,threshold))))
}


###  Wrapper for C Implementation of DYNSBM
dynsbm <- function(net.mat, kk=3, tmap, hours.vec, self.ties=TRUE,
                   total=1000, burn.in=0, thin=1,
                   hyperpriors=dynsbm.priors(), eta=rep(1/kk,kk),
                   init.vals=NULL, #spectral.start=FALSE,
                   clean.out=TRUE, verbose=0, max.runs=200,
                   label.switch.mode = c("adhoc","kl-loss"),
                   autoconverge=wsbm.convergence(),
                   multiImpute=FALSE){

    if(length(dim(net.mat)) < 3){
        stop("This is a single time point network. Consider using wsbm or sbm.")
    }else if(length(dim(net.mat)) > 3){
        stop("net.mat has too many dimensions")
    }

    nn <- dim(net.mat)[1]
    if(dim(net.mat)[2] != nn) stop("net.mat must be square in first two dims")
    TT <- dim(net.mat)[3]
    if(TT != length(tmap)){
        stop("tmap must have same length as number of adjacency matrices")
    }

    ee <- length(unique(tmap))

    if(missing(hours.vec)){
        hours.vec <- rep(1,ee)
    }else{
        if(length(hours.vec) != ee) stop("hours.vec has improper length")
        if(any(hours.vec <= 0 )) stop("hours.vec must have positive elements")
    }

    ##  Formatting Adjacency Matrix for C
    label.switch.mode <- match.arg(label.switch.mode)
    multi.int <- as.integer(ifelse(multiImpute,1,0))
    net.mat.clean <- net.mat
    if(!self.ties){
        for(tt in 1:TT) diag(net.mat.clean[,,tt]) <- NA
    }
    net.mat.clean[is.na(net.mat.clean)] <- -1


    short.total <- total ##  Length of chain to store
    total <- short.total*thin + burn.in ## Actual Number of Iterations

    ##  Storage for Priors
    flatBB.prior = double(short.total * 2 * kk*kk * ee)
    flatSS.prior = double(short.total * 2 * nn * ee)
    flatRR.prior = double(short.total * 2 * nn * ee)

    ##  Storage for Parameters
    flatSS = double(short.total * nn * TT)
    flatRR = double(short.total * nn * TT)
    flatBB = double(short.total * kk*kk * TT)
    flatMMB = integer(short.total * nn)

    ##  Storage for Performance Indicators
    ll.vec <- double(short.total)
    flatHH <- double(short.total*kk*nn)


    ##  Spectral Clustering Initialization Not Currently Supported
    ## if(is.null(init.vals) & spectral.start){
    ##     if(verbose > 0) message("Initializing with Spectral Clustering...")
    ##     init.vals <- sbm.spectral(net=net,kk=kk,weighted=TRUE)
    ##     if(verbose > 0) message("Initialization Complete.")
    ## }

    ##  Loading Initial Values
    if(is.null(init.vals$mmb)){
        init.vals$mmb <- sample(kk,nn,replace=TRUE,prob=eta)
    }
    if(is.null(init.vals$BB)){
        init.vals$BB.mat <- array(rbeta(kk^2*TT,1,1),c(kk,kk,TT))
    }
    if(is.null(init.vals$SS.mat) || is.null(init.vals$RR.mat) ||
       (dim(init.vals$SS.mat) != c(nn,TT)) || dim(init.vals$RR.mat != c(nn,TT))){
        ss.init <- apply(net.mat,3,rowMeans,na.rm=TRUE)
        rr.init <- apply(net.mat,3,colMeans,na.rm=TRUE)

        if(any(ss.init == 0) || any(rr.init == 0)){
            ss.init <- ss.init + 0.001
            rr.init <- rr.init + 0.001
        }

        normalize <- function(x) return(x/mean(x))

        ss.init <- apply(ss.init,2,normalize)
        rr.init <- apply(rr.init,2,normalize)

        init.vals$SS.mat <- ss.init
        init.vals$RR.mat <- rr.init
    }
    if(is.null(init.vals$HH)){
        init.vals$HH <- mmb.to.PI(init.vals$mmb)
    }
    if(is.null(init.vals$logLik)){
        init.vals$logLik <- dynsbm.log.like.net(net.mat,mmb=init.vals$mmb,
                                                BB.mat=init.vals$BB.mat,
                                                SS.mat=init.vals$SS.mat,
                                                RR.mat=init.vals$RR.mat,
                                                self.ties=self.ties,
                                                hours.vec=hours.vec[tmap])
    }

    start = 1
    for(tt in 1:TT){
        flatBB[1:(kk*kk) + (tt-1)*short.total*kk^2] <- init.vals$BB[,,tt]
        flatSS[1:nn + (tt-1)*short.total*nn] <- init.vals$SS[,tt]
        flatRR[1:nn + (tt-1)*short.total*nn] <- init.vals$RR[,tt]
    }
    for(tt in 1:ee){
        BB.add <- (tt-1)*short.total*2*kk*kk
        flatBB.prior[1:(kk^2) + BB.add] <- min.threshold(apply(init.vals$BB[,,tmap==tt,drop=FALSE],
                                                               c(1,2),mean),
                                                         1e-5)
        flatBB.prior[1:(kk^2) + (kk^2) + BB.add] <- 1
        SR.add <- (tt-1)*short.total*2*nn
        flatSS.prior[1:nn + SR.add] <- min.threshold(rowMeans(init.vals$SS[,tmap==tt]),
                                                     1e-5)
        flatSS.prior[1:nn + nn + SR.add] <- 1
        flatRR.prior[1:nn + SR.add] <- min.threshold(rowMeans(init.vals$RR[,tmap==tt]),
                                                     1e-5)
        flatRR.prior[1:nn + nn + SR.add] <- 1
    }
##    browser()
    flatMMB[1:nn] <- init.vals$mmb
    flatHH[1:(kk*nn)] <- init.vals$HH
    ll.vec[1] <- init.vals$logLik


    extend.max <- autoconverge$extend.max
    shift.size <- autoconverge$shift.size
    qq <- as.double(qnorm(1 - autoconverge$alpha/2))


    ##  Calling C Implementation of MCMC
    out <- .C("dynsbm",
              as.integer(total),                    # 1: NUMBER of MCMC DRAWS
              as.integer(burn.in), as.integer(thin),# 2-3: Burnin, Thin Control
              as.integer(start),                    # 4: Starting Iteration
              extend.max, shift.size, qq,           # 5-7: Autoconverge Control
              as.integer(nn),                       # 8: Nodes in Networks
              as.integer(net.mat.clean),            # 9: Networks
              as.integer(kk),                       # 10: Number of Blocks
              multi.int,                            # 11: Multi Impute Indicator
              as.integer(TT), as.integer(ee),       # 12-13: Time / Equiv Pts
              as.integer(tmap),                     # 14: Time/Equiv Map
              as.double(hours.vec),                 # 15: Hours Per Equiv
              as.double(hyperpriors$SS),            # 16: Sender Hyperprior
              as.double(hyperpriors$RR),            # 17: Receiver Hyperprior
              as.double(hyperpriors$BB),            # 18: BlockMatrix Hyperprior
              as.double(eta),                       # 19: BlockMemb Prior
              flatSS.prior,                         # 20: Sender Priors
              flatRR.prior,                         # 21: Receiver Priors
              flatBB.prior,                         # 22: BlockMatrix Priors
              as.integer(flatMMB),                  # 23: BlockMemb Params
              flatSS,                               # 24: Sender Params
              flatRR,                               # 25: Receiver Params
              flatBB,                               # 26: BlockMatrix Params
              ll.vec,                               # 27: Log-Likelihood
              flatHH,                               # 28: PosteriorMemb Matrix
              as.integer(verbose))                  # 29: Verbose Indicator

    SS.prior <- array(out[[20]],c(nn,2,short.total,ee))
    RR.prior <- array(out[[21]],c(nn,2,short.total,ee))
    BB.prior <- array(out[[22]],c(kk,kk,2,short.total,ee))
    MMB <- array(out[[23]],c(nn,short.total))

    SS.mat <- array(out[[24]],c(nn,short.total,TT))
    RR.mat <- array(out[[25]],c(nn,short.total,TT))
    BB.mat <- array(out[[26]],c(kk,kk,short.total,TT))

    return(list(init.vals=init.vals,
                SS.prior=SS.prior,RR.prior=RR.prior,
                BB.prior=BB.prior,MMB=MMB,logLik=out[[27]],
                SS.mat=SS.mat,RR.mat=RR.mat,BB.mat=BB.mat))

    ##                                     #  browser()
    ## ##  Pulling the flat matrices from the C output
    ## BB.flat <- matrix(out[[9]],nrow=kk*kk)
    ## BB.flat <- apply(BB.flat,2,transpose.vector,nrow=kk)
    ## BB <- array(BB.flat,c(kk,kk,short.total))

    ## mmb <- array(out[[10]],c(nn,short.total))
    ## SS <- array(out[[11]],c(nn,short.total))
    ## RR <- array(out[[12]],c(nn,short.total))


    ## ll.vec <- as.vector(out[[17]])
    ## HH.flat <- out[[21]]
    ## HH <- array(HH.flat,c(nn,kk,short.total))

    ## pmat.mat <- NULL  ## I might want to compute this later...
    ## diag(net.clean) <- -1


    ## dynsbm.out <- structure(list(BB=BB,mmb=mmb,
    ##                            SS=SS,RR=RR,
    ##                            net.mat=net.mat,logLik=ll.vec,
    ##                            burn.in=burn.in,thin=thin,self.ties=self.ties,
    ##                            pmat=pmat.mat,HH=HH),class="dynsbm")

    ## if(label.switch.mode == "kl-loss"){
    ##     dynsbm.out <- switch.labels(dynsbm.out,max.runs=max.runs,verbose=verbose)
    ## }else{
    ##     ##  Do nothing
    ## }
    ##                                     #  browser()
    ## ## Summarizing MCMC Chain
    ## dynsbm.summ <- summary(dynsbm.out)
    ## dynsbm.summ$priors = priors


    ## if(clean.out){
    ##     dynsbm.summ$chain <- list(logLik=dynsbm.out$logLik)
    ##     dynsbm.summ$clean <- TRUE
    ## }else{
    ##     dynsbm.summ$chain <- dynsbm.out
    ##     dynsbm.summ$chain$net <- NULL
    ##     dynsbm.summ$clean <- FALSE
    ## }

    ## return(dynsbm.summ)
}



#################################################################
######################  dynsbm class functions  ####################
#################################################################


summary.dynsbm <- function(object,...){
###  browser()
    total <- dim(object$BB)[3]; kk <- dim(object$BB)[1]

    BB.hat <- apply(object$BB,c(1,2),mean)
    ##  PI.mean <- apply(object$PI[,,my.sub],c(1,2),mean)
    PI.mean <- t(apply(object$mmb,1,tabulate,nbins=kk) / total)
    mmb <- apply(PI.mean,1,which.max)

    SS.hat <- rowMeans(object$SS)
    RR.hat <- rowMeans(object$RR)

    summ.obj <- structure(list(mmb=mmb,BB=BB.hat,
                                SS=SS.hat,RR=RR.hat,PI.mean=PI.mean),class="dynsbm")
    summ.obj$self.ties <- object$self.ties

    summ.obj$pmat <- predict(summ.obj)
    summ.obj$net <- object$net
    ## Calculating DIC

    #######################  FUNCTION OBSOLETE
    summ.obj$logLik <- with(summ.obj,dynsbm.ll.pmat(net,pmat,
                                                    self.ties=self.ties))
    #######################
    summ.obj$DIC <- 2*summ.obj$logLik - 4 * mean(object$logLik)
    summ.obj$burn.in <- object$burn.in; summ.obj$thin <- object$thin

    return(summ.obj)

}

predict.dynsbm <- function(object,...){
    PI <- mmb.to.PI(object$mmb,kk=ncol(object$BB))
    pmat <- PI %*% object$BB %*% t(PI)
    if(!object$self.ties) diag(pmat) <- NA

    sr.mat <- as.matrix(object$SS) %*% t(as.matrix(object$RR))
    pmat <- pmat * sr.mat
    return(pmat)
}


dynsbm.metric <- function(graph,kk=2,total=1500,
                       thin=1,burn.in=500,verbose=0,...){

    ##    mode <- match.arg(mode)
    ##    if(mode == "mcmc"){
    dynsbm.fit <- dynsbm(total=total,net=graph,kk=kk,verbose=verbose,
                   thin=thin,burn.in=burn.in)
    ##    }
    return(dynsbm.fit$pmat)
}



plot.net.dynsbm <- function(dynsbm.obj,
                         ord = order(dynsbm.obj$mmb),...){
    nn <- nrow(dynsbm.obj$net)
    image(1:nn,1:nn,dynsbm.obj$net[ord[nn:1],ord],col=grey((50:1)/50),
          ylab="",xlab="",yaxt="n",xaxt="n",...)
}

plot.fit.dynsbm <- function(dynsbm.obj,
                          ord = order(dynsbm.obj$mmb),...){
    nn <- nrow(dynsbm.obj$net)
    image(1:nn,1:nn,dynsbm.obj$pmat[ord[nn:1],ord],col=grey((50:1)/50),
          ylab="",xlab="",yaxt="n",xaxt="n",...)
}





##########################################################
##################  ROTATION FUNCTIONS  ##################
##########################################################

dynsbm.log.like.net <- function(net.mat,mmb,BB.mat,SS.mat,RR.mat,
                                self.ties=TRUE,hours.vec){

    kk <- dim(BB.mat)[1]
    TT <- dim(net.mat)[3]
    if(missing(hours.vec)){
        hours.vec <- rep(1,TT)
    }
    ll.total <- 0
    for(tt in 1:TT){
        ll.single <- wsbm.log.like.net(net=net.mat[,,tt],
                                                 mmb=mmb,BB=BB.mat[,,tt],
                                                 SS=SS.mat[,tt],RR=RR.mat[,tt],
                                                 self.ties=self.ties,
                                                 hours=hours.vec[tt])
        ll.total <- ll.total + ll.single
    }

    return(ll.total)
}

## dynsbm.ll.pmat <- function(net,pmat,self.ties=TRUE){
##     if(!self.ties){
##         diag(net) <- NA
##     }
##     return(sum(log(pmat)*net,na.rm=TRUE) - sum(pmat[!is.na(net)]))
## }

dynsbm.marginal.ll.single <- function(net,BB,theta=rep(1/ncol(BB),ncol(BB))){
    nn <- ncol(net)
    kk <- ncol(BB)

    PI.t <- rmultinom(nn,1,theta)
    PP <- t(PI.t) %*% BB %*% PI.t
    diag(PP) <- NA
    ll <- sum(log(PP * net + (1-PP) * (1 - net)),na.rm=TRUE)
    return(ll)
}

dynsbm.marginal.log.like.net <- function(net,BB,theta=rep(1/ncol(BB),ncol(BB)),
                                      iter.max=1e4){
    ll.vec <- replicate(iter.max,marginal.ll.single(net,BB,theta))
    ll.max <- max(ll.vec)
    ll.vec <- ll.vec - ll.max
    return(log(mean(exp(ll.vec))) + ll.max)
}



dynsbm.load.init.vals <- function(init.vals,nn,kk){

    if(is.null(init.vals$BB) | is.null(init.vals$PI)){
        stop("init.vals must be a list containing BB and PI")
    }
    ## Checking BB
    if(any(dim(init.vals$BB) != kk)) stop("BB must be a kk by kk matrix")
    if(any(init.vals$BB < 0 | init.vals$BB > 1))
        stop("BB must have values between 0 and 1")

    ## Checking PI
    if(any(dim(init.vals$PI)!=c(nn,kk))) stop("PI must be an nn by kk matrix")

    BB.init <- double(kk^2)
    PI.init <- double(kk*nn)
    for(jj in 1:kk){
        BB.init[((jj-1) * kk + 1):(jj*kk)] <- init.vals$BB[jj,]
    }
    for(jj in 1:nn){
        PI.init[((jj-1)*kk + 1):(jj*kk)] <- init.vals$PI[jj,]
    }
    flatTable <- array(c(BB.init,PI.init),c(1,kk*(kk+nn)))
    return(flatTable)
}
