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
                         alpha.init = 1, beta.init = 1, prop.sd=1.0,
                         method=c("mv","alpha-beta","norm")){

    method <- match.arg(method)
    ##  Allocating Memory for C++
    alpha.vec <- beta.vec <- rep(0,n)
    ##    print(prop.sd)
    if(method == "alpha-beta"){
        out <- .C("RGammaPrior",
                  as.integer(n),as.integer(thin),as.integer(burn.in),
                  as.double(alpha.init),as.double(beta.init),
                  as.double(p),as.double(q),as.double(r),as.double(s),
                  as.double(prop.sd),as.double(alpha.vec),as.double(beta.vec))
    }else if(method == "mv"){
        out <- .C("RGammaPriorMV",
                  as.integer(n),as.integer(thin),as.integer(burn.in),
                  as.double(alpha.init),as.double(beta.init),
                  as.double(p),as.double(q),as.double(r),as.double(s),
                  as.double(prop.sd),as.double(alpha.vec),as.double(beta.vec))
    }else if(method == "norm"){
        out <- .C("RGammaPriorNorm",
                  as.integer(n),as.integer(thin),as.integer(burn.in),
                  as.double(alpha.init),as.double(beta.init),
                  as.double(p),as.double(q),as.double(r),as.double(s),
                  as.double(prop.sd),as.double(alpha.vec),as.double(beta.vec))
    }
    return(cbind(out[[11]],out[[12]]))
}




dgamma.prior <- function(alpha,beta,p=1,q=2,r=2,s=2){
    l.num <- (alpha-1)*log(p) - beta*q + alpha*s*log(beta)
    l.den <- r*lgamma(alpha)

    return(exp(l.num - l.den))
}

ldgamma.prior <- function(alpha,beta,p=1,q=2,r=2,s=2){
    out <- .C("ldGammaPrior",
              as.double(alpha),as.double(beta),
              as.double(p),as.double(q),as.double(r),as.double(s),
              double(1))

    return(out[[7]])
}



data.gen.dynsbm <- function(nn=40,TT=4,tmap=rep(1:2,2),hours.vec=c(1,1),
                            BB.prior,SS.prior,RR.prior,
                            mmb.prior=c(.5,.5),mmb,
                            self.ties=TRUE,debug=FALSE,normalize=FALSE){
    if(debug) browser()

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
        if(normalize){
            mmb.count <- ss.norm <- rr.norm <- rep(NA,kk)

            for(ll in 1:kk){
                ind.ll <- which(mmb == ll)
                mmb.count[ll] <- length(ind.ll)
                ss.norm[ll] <- sum(SS.mat[ind.ll,tt])
                rr.norm[ll] <- sum(RR.mat[ind.ll,tt])
                SS.mat[ind.ll,tt] <- SS.mat[ind.ll,tt] * (mmb.count[ll] / ss.norm[ll])
                RR.mat[ind.ll,tt] <- RR.mat[ind.ll,tt] * (mmb.count[ll] / rr.norm[ll])
            }
            ## for(ii in 1:kk){
            ##     for(jj in 1:kk){
            ##         BB.mat[ii,jj,tt] <- BB.mat[ii,jj,tt] / (mmb.count[ii] / ss.norm[ii]) / (mmb.count[jj] / rr.norm[jj])
            ##     }
            ## }
        }
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



## convert.time.table <- function(tab,cutoffs,nn,labs){
##     TT <- length(cutoffs) - 1
##     if(missing(nn)){
##         if(missing(labs)){
##             labs <- unique(tab[,1])
##         }
##         nn <- length(labs)
##     }

##     AA.mat <- array(0,c(nn,nn,TT))
##     for(ii in 1:nrow(tab)){
##         tt <- 0
##         while(tab[ii,3] > cutoffs[tt+1]) tt <- tt + 1
##         if(tt > 0 && tt <= TT){
##             AA.mat[tab[ii,1],tab[ii,2],tt] <-AA.mat[tab[ii,1],tab[ii,2],tt] + 1
##         }
##     }
##     return(AA.mat)
## }



dynsbm.priors <- function(BB.hyperprior=c(.99,1.01,1,1),
                          SS.hyperprior=c(.99,1.01,1,1),
                          RR.hyperprior=c(.99,1.01,1,1)){

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
                SS=SS.hyperprior,
                RR=RR.hyperprior))
}


min.threshold <- function(vec,threshold=1e-5){
    return(sapply(vec,FUN=function(x) return(max(x,threshold))))
}


###  Wrapper for C Implementation of DYNSBM
dynsbm <- function(net.mat, kk=3, tmap, hours.vec, self.ties=TRUE,
                   total=1000, burn.in=0, thin=1,
                   hyperpriors=dynsbm.priors(), eta=rep(1/kk,kk),
                   init.vals=NULL, #spectral.start=FALSE,
                   clean.out=FALSE, verbose=0, max.runs=200,
                   label.switch.mode = c("adhoc","kl-loss"),
                   autoconverge=wsbm.convergence(),
                   multiImpute=FALSE,update.mmb=TRUE){

##################################################
################# ERROR CHECKING #################
##################################################


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


##################################################
######## Creating Containers for Output ##########
##################################################


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
    ll.vec <- double(short.total*TT)
    flatHH <- double(short.total*kk*nn)


    ##  Spectral Clustering Initialization Not Currently Supported
    ## if(is.null(init.vals) & spectral.start){
    ##     if(verbose > 0) message("Initializing with Spectral Clustering...")
    ##     init.vals <- sbm.spectral(net=net,kk=kk,weighted=TRUE)
    ##     if(verbose > 0) message("Initialization Complete.")
    ## }


##################################################
############# Initializing Values  ###############
##################################################

    ##  Initializing
    init.vals <- dynsbm.load.init.vals(init.vals=init.vals,net.mat=net.mat,kk=kk,tmap=tmap,
                                       eta=eta,self.ties=self.ties,hours.vec=hours.vec)

    ##  Flattening Initial Values for Pass to C++
    start = 1
    for(tt in 1:TT){
        flatBB[1:(kk*kk) + (tt-1)*short.total*kk^2] <- init.vals$BB.mat[,,tt]
        flatSS[1:nn + (tt-1)*short.total*nn] <- init.vals$SS.mat[,tt]
        flatRR[1:nn + (tt-1)*short.total*nn] <- init.vals$RR.mat[,tt]
    }


    for(tt in 1:ee){
        BB.add <- (tt-1)*short.total*2*kk*kk
        flatBB.prior[1:(kk^2) + BB.add] <- init.vals$BB.prior[,,1,tt]
        flatBB.prior[1:(kk^2) + (kk^2) + BB.add] <- init.vals$BB.prior[,,2,tt]

        SR.add <- (tt-1)*short.total*2*nn
        flatSS.prior[1:nn + SR.add] <- init.vals$SS.prior[,1,tt]
        flatSS.prior[1:nn + nn + SR.add] <- init.vals$SS.prior[,2,tt]

        flatRR.prior[1:nn + SR.add] <- init.vals$RR.prior[,1,tt]
        flatRR.prior[1:nn + nn + SR.add] <- init.vals$RR.prior[,2,tt]
    }
    flatMMB[1:nn] <- init.vals$mmb
    flatHH[1:(kk*nn)] <- init.vals$HH

    ##  MCMC Convergence Checking Parameters
    extend.max <- autoconverge$extend.max
    shift.size <- autoconverge$shift.size
    qq <- as.double(qnorm(1 - autoconverge$alpha/2))


##################################################
############  Call to C++ Function  ##############
##################################################
    out <- .C("dynsbm",
              as.integer(nn),                       # 1: Nodes in Networks
              as.integer(kk),                       # 2: Number of Blocks
              as.integer(TT), as.integer(ee),       # 3-4: Time / Equiv Pts
              multi.int,                            # 5: Multi Impute Indicator

              as.integer(net.mat.clean),            # 6: Networks
              as.integer(tmap),                     # 7: Time/Equiv Map
              as.double(hours.vec),                 # 8: Hours Per Equiv

              as.double(hyperpriors$SS),            # 9: Sender Hyperprior
              as.double(hyperpriors$RR),            # 10: Receiver Hyperprior
              as.double(hyperpriors$BB),            # 11: BlockMatrix Hyperprior
              as.double(eta),                       # 12: BlockMemb Prior

              flatSS.prior,                         # 13: Sender Priors
              flatRR.prior,                         # 14: Receiver Priors
              flatBB.prior,                         # 15: BlockMatrix Priors
              as.integer(flatMMB),                  # 16: BlockMemb Params

              flatSS,                               # 17: Sender Params
              flatRR,                               # 18: Receiver Params
              flatBB,                               # 19: BlockMatrix Params

              ll.vec,                               # 20: Log-Likelihood
              flatHH,                               # 21: PosteriorMemb Matrix

              as.integer(verbose),                  # 22: Verbose Indicator
              as.integer(update.mmb),                # 23: Update Membership Vector

              as.integer(total),                    # 24: NUMBER of MCMC DRAWS
              as.integer(burn.in), as.integer(thin),# 25-6: Burnin, Thin Control
              as.integer(start),                    # 27: Starting Iteration
              extend.max, shift.size, qq)           # 28-30: Autoconverge Control


##################################################
###########  Processing C++ Output  ##############
##################################################
    SS.prior <- array(out[[13]],c(nn,2,short.total,ee))
    RR.prior <- array(out[[14]],c(nn,2,short.total,ee))
    BB.prior <- array(out[[15]],c(kk,kk,2,short.total,ee))
    MMB <- array(out[[16]],c(nn,short.total))

    SS.mat <- array(out[[17]],c(nn,short.total,TT))
    RR.mat <- array(out[[18]],c(nn,short.total,TT))
    BB.mat <- array(out[[19]],c(kk,kk,short.total,TT))
    HH.mat <- array(out[[21]],c(nn,kk,short.total))

    ll.mat <- array(out[[20]],c(TT,short.total))
    ll.vec <- colSums(ll.mat)

    ##  Combining output into single structure
    dynsbm.obj <- structure(list(BB.mat=BB.mat,mmb=MMB,HH=HH.mat,
                                 SS.mat=SS.mat,RR.mat=RR.mat,
                                 BB.prior=BB.prior,
                                 SS.prior=SS.prior,RR.prior=RR.prior,
                                 ll.mat=ll.mat,ll.vec=ll.vec,
                                 net.mat=net.mat,self.ties=self.ties,
                                 init.vals=init.vals,
                                 tmap=tmap,hours.vec=hours.vec,
                                 thin=thin,burn.in=burn.in),
                            class="dynsbm")

    ##  Cleaning up output from C++ Call
    rm(out)


######  Label Switching Handler
    if(label.switch.mode == "kl-loss"){
        ##        dynsbm.out <- switch.labels(dynsbm.out,max.runs=max.runs,verbose=verbose)
        message("kl-loss mode is currently disabled for dynsbm")
    }else{
        ##  Do nothing
    }

##################################################
########  Generating Summary Statistics  #########
##################################################

    dynsbm.summ <- summary(dynsbm.obj)
    dynsbm.summ$hyperpriors = hyperpriors


    if(clean.out){
        dynsbm.summ$chain <- list(ll.vec=dynsbm.obj$ll.vec)
        dynsbm.summ$clean <- TRUE
    }else{
        dynsbm.summ$chain <- dynsbm.obj
        dynsbm.summ$chain$net.mat <- NULL
        dynsbm.summ$clean <- FALSE
    }

    return(dynsbm.summ)
}


get.iter.list.dynsbm <- function(object,iter){
    if(is.null(object$chain)){
        stop(paste0("Object does not contain MCMC chain.  Call dynsbm with clean.out=FALSE"))
    }
    if(dim(object$chain$mmb)[2] < iter){
        stop(paste0("Chain has less than ",iter," iterations"))
    }

    BB.prior <- array(NA,dim(object$chain$BB.prior)[-4])
    SS.prior <- RR.prior <- array(NA,dim(object$chain$SS.prior)[-3])
    ee <- dim(BB.prior)[4]
    for(tt in 1:ee){
        BB.prior[,,,tt] <- object$chain$BB.prior[,,,iter,tt]
        SS.prior[,,tt] <- object$chain$SS.prior[,,iter,tt]
        RR.prior[,,tt] <- object$chain$RR.prior[,,iter,tt]
    }

    iter.list <- list(mmb=object$chain$mmb[,iter],
                      BB.mat=object$chain$BB.mat[,,iter,],
                      SS.mat=object$chain$SS.mat[,iter,],
                      RR.mat=object$chain$RR.mat[,iter,],
                      BB.prior=BB.prior,
                      SS.prior=SS.prior,RR.prior=RR.prior)

    return(iter.list)
}

#################################################################
######################  dynsbm class functions  ####################
#################################################################


summary.dynsbm <- function(object,...){
###    browser()
    total <- dim(object$BB.mat)[3]; kk <- dim(object$BB.mat)[1]

###  Summary of Parameters
    BB.hat <- apply(object$BB.mat,c(1,2,4),mean)
    SS.hat <- apply(object$SS.mat,c(1,3),mean)
    RR.hat <- apply(object$RR.mat,c(1,3),mean)
    ##  PI.mean <- apply(object$PI[,,my.sub],c(1,2),mean)
    PI.mean <- t(apply(object$mmb,1,tabulate,nbins=kk) / total)
    mmb <- apply(PI.mean,1,which.max)

###  Summary of Priors
    BB.prior.hat <- apply(object$BB.prior,c(1,2,3,5),mean)
    BB.prior.mean <- apply(object$BB.prior[,,1,,,drop=FALSE]/object$BB.prior[,,2,,,drop=FALSE],c(1,2,5),mean)
    BB.prior.var <- apply(object$BB.prior[,,1,,,drop=FALSE]/(object$BB.prior[,,2,,,drop=FALSE]^2),c(1,2,5),mean)

    SS.prior.hat <- apply(object$SS.prior,c(1,2,4),mean)
    SS.prior.mean <- apply(object$SS.prior[,1,,,drop=FALSE]/object$SS.prior[,2,,,drop=FALSE],c(1,4),mean)
    SS.prior.var <- apply(object$SS.prior[,1,,,drop=FALSE]/(object$SS.prior[,2,,,drop=FALSE]^2),c(1,4),mean)

    RR.prior.hat <- apply(object$RR.prior,c(1,2,4),mean)
    RR.prior.mean <- apply(object$RR.prior[,1,,,drop=FALSE]/object$RR.prior[,2,,,drop=FALSE],c(1,4),mean)
    RR.prior.var <- apply(object$RR.prior[,1,,,drop=FALSE]/(object$RR.prior[,2,,,drop=FALSE]^2),c(1,4),mean)


    summ.obj <- structure(list(mmb=mmb,BB=BB.hat,SS=SS.hat,RR=RR.hat,
                               BB.prior=BB.prior.hat,
                               BB.prior.mean=BB.prior.mean,
                               BB.prior.var=BB.prior.var,
                               SS.prior=SS.prior.hat,
                               SS.prior.mean=SS.prior.mean,
                               SS.prior.var=SS.prior.var,
                               RR.prior=RR.prior.hat,
                               RR.prior.mean=RR.prior.mean,
                               RR.prior.var=RR.prior.var,
                               PI.mean=PI.mean),
                          class="dynsbm")
    summ.obj$self.ties <- object$self.ties
    summ.obj$net.mat <- object$net.mat
    summ.obj$tmap <- object$tmap
    summ.obj$hours.vec <- object$hours.vec


    ## Calculating DIC
    ## summ.obj$logLik <- with(summ.obj,dynsbm.ll.pmat(net,pmat,
    ##                                                 self.ties=self.ties))
    ## summ.obj$pmat <- predict(summ.obj)
    ## summ.obj$DIC <- 2*summ.obj$logLik - 4 * mean(object$logLik)

    summ.obj$burn.in <- object$burn.in; summ.obj$thin <- object$thin

    return(summ.obj)

}


## predict.dynsbm <- function(object,...){
##     PI <- mmb.to.PI(object$mmb,kk=ncol(object$BB))
##     pmat <- PI %*% object$BB %*% t(PI)
##     if(!object$self.ties) diag(pmat) <- NA

##     sr.mat <- as.matrix(object$SS) %*% t(as.matrix(object$RR))
##     pmat <- pmat * sr.mat
##     return(pmat)
## }


## dynsbm.metric <- function(graph,kk=2,total=1500,
##                        thin=1,burn.in=500,verbose=0,...){

##     ##    mode <- match.arg(mode)
##     ##    if(mode == "mcmc"){
##     dynsbm.fit <- dynsbm(total=total,net=graph,kk=kk,verbose=verbose,
##                    thin=thin,burn.in=burn.in)
##     ##    }
##     return(dynsbm.fit$pmat)
## }




##########################################################
##################  HELPER FUNCTIONS  ####################
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

## dynsbm.marginal.ll.single <- function(net,BB,theta=rep(1/ncol(BB),ncol(BB))){
##     nn <- ncol(net)
##     kk <- ncol(BB)

##     PI.t <- rmultinom(nn,1,theta)
##     PP <- t(PI.t) %*% BB %*% PI.t
##     diag(PP) <- NA
##     ll <- sum(log(PP * net + (1-PP) * (1 - net)),na.rm=TRUE)
##     return(ll)
## }

## dynsbm.marginal.log.like.net <- function(net,BB,theta=rep(1/ncol(BB),ncol(BB)),
##                                       iter.max=1e4){
##     ll.vec <- replicate(iter.max,marginal.ll.single(net,BB,theta))
##     ll.max <- max(ll.vec)
##     ll.vec <- ll.vec - ll.max
##     return(log(mean(exp(ll.vec))) + ll.max)
## }



dynsbm.load.init.vals <- function(init.vals,net.mat,kk,tmap,
                                  eta=rep(1/kk,kk),self.ties=TRUE,hours.vec){
    nn <- dim(net.mat)[1]
    TT <- dim(net.mat)[3]
    ee <- length(unique(tmap))

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
        init.vals$HH <- mmb.to.PI(init.vals$mmb,kk)
    }
    if(is.null(init.vals$logLik)){
        init.vals$logLik <- dynsbm.log.like.net(net.mat,mmb=init.vals$mmb,
                                                BB.mat=init.vals$BB.mat,
                                                SS.mat=init.vals$SS.mat,
                                                RR.mat=init.vals$RR.mat,
                                                self.ties=self.ties,
                                                hours.vec=hours.vec[tmap])
    }

    ##  Loading Hyperparameters
    if(is.null(init.vals$BB.prior)){
        init.vals$BB.prior <- array(NA,c(kk,kk,2,ee))
        for(tt in 1:ee){
            init.vals$BB.prior[,,1,tt] <- min.threshold(apply(init.vals$BB.mat[,,tmap==tt,drop=FALSE],
                                                              c(1,2),mean),
                                                        1e-5)
            init.vals$BB.prior[,,2,tt] <- 1
        }
    }
    if(is.null(init.vals$SS.prior) || any(dim(init.vals$SS.prior) != c(nn,2,ee))){
        if(!is.null(init.vals$SS.prior)){
            message("SS.prior has dimension: (",paste0(dim(init.vals$SS.prior),collapse=","),")")
            message("Expected dimension was: (",paste(nn,2,ee,sep=","),")")
            message("Reinitializing SS.prior...")
        }
        init.vals$SS.prior <- array(NA,c(nn,2,ee))
        for(tt in 1:ee){
            init.vals$SS.prior[,1,tt] <- min.threshold(rowMeans(init.vals$SS.mat[,tmap==tt,drop=FALSE]),
                                                       1e-5)
            init.vals$SS.prior[,2,tt] <- 1
        }
    }
    if(is.null(init.vals$RR.prior) || any(dim(init.vals$RR.prior) != c(nn,2,ee))){
        if(!is.null(init.vals$RR.prior)){
            message("RR.prior has dimension: (",paste0(dim(init.vals$RR.prior),collapse=","),")")
            message("Expected dimension was: (",paste(nn,2,ee,sep=","),")")
            message("Reinitializing RR.prior...")
        }
        init.vals$RR.prior <- array(NA,c(nn,2,ee))
        for(tt in 1:ee){
            init.vals$RR.prior[,1,tt] <- min.threshold(rowMeans(init.vals$RR.mat[,tmap==tt,drop=FALSE]),
                                                       1e-5)
            init.vals$RR.prior[,2,tt] <- 1
        }
    }

    return(init.vals)
}



################################################################################
#####################  Data Visualization Functions  ###########################
################################################################################

BlockMat.plot <- function(object, tclass, remap, true.mat,
                          ss.ind, rr.ind, add.mle=FALSE,...){

    kk <- dim(object$BB)[1]
    nn <- dim(object$SS)[1]
    TT <- dim(object$net.mat)[3]

    if(missing(ss.ind)){
        ss.ind <- 1:kk
    }
    if(missing(rr.ind)){
        rr.ind <- 1:kk
    }

    if(missing(remap)){
        remap <- 1:kk
    }else{
        if(length(remap) != kk){
            warning("remap has incorrect length.  Using default of 1:kk")
            remap <- 1:kk
        }
    }

    leg.txt=c("Prior Mean")
    leg.col = c("blue")
    if(add.mle){
        mle.mat <- apply(object$net.mat, 3, weighted.mle,
                         mmb=object$mmb, kk=kk)
        mle.mat <- array(mle.mat,c(kk,kk,TT))/object$hours
        leg.txt <- c(leg.txt,"Naive MLE")
        leg.col <-  c(leg.col,"red")
    }

    if(!missing(true.mat)){
        leg.txt <- c(leg.txt,"True Value")
        leg.col <- c(leg.col,"green")
    }

    tt.vec <- which(object$tmap == tclass)
    par(mfrow=c(length(ss.ind),length(rr.ind)))

    for(ii in ss.ind){#1:kk){
        for(jj in rr.ind){#1:kk){
            boxplot(object$chain$BB.mat[remap[ii],remap[jj],,tt.vec], ...)

            if(add.mle){
                segments(x0=1:length(tt.vec) - .3,
                         x1=1:length(tt.vec) + .3,
                         y0=mle.mat[remap[ii],remap[jj],tt.vec],lwd=3,col="red")
            }

            if(!missing(true.mat)){
                segments(x0=1:length(tt.vec) - .3,
                         x1=1:length(tt.vec) + .3,
                         y0=true.mat[ii,jj,tt.vec],col="green",lwd=3)

            }
            abline(h=object$BB.prior.mean[remap[ii],remap[jj],tclass],lwd=3,col="blue")
        }
    }

    legend("topleft",legend=leg.txt,col=leg.col,lwd=3)
    par(mfrow=c(1,1))
}


SenderMat.plot <- function(object, tclass, true.mat, horizontal=FALSE, nodes, block, ...){

    kk <- dim(object$BB)[1]
    nn <- dim(object$SS)[1]
    TT <- dim(object$net.mat)[3]

    if(missing(nodes)){
        if(missing(block)){
            nodes <- 1:nn
        }else{
            if(block > kk) stop("block must be between 1 and ",kk)
            nodes <- which(object$mmb == block)
        }
    }else{
        if(any(nodes > nn) | any(nodes < 1)){
            stop("nodes must contain integers between 1 and ",nn)
        }
        if(!missing(block)){
            warning("The argument for block is ignored when nodes is used.")
        }
    }


    tt.vec <- which(object$tmap == tclass)
    if(horizontal){
        par(mfrow=c(1,length(tt.vec)))
    }else{
        par(mfrow=c(length(tt.vec),1))
    }
    for(ii in 1:length(tt.vec)){
        boxplot(t(object$chain$SS.mat[nodes,,tt.vec[ii]]),names=nodes, horizontal=horizontal, ...)

        if(!missing(true.mat)){
            if(horizontal){
                segments(y0=1:length(nodes) - .3,
                         y1=1:length(nodes) + .3,
                         x0=true.mat[nodes,tt.vec[ii]],col="green",lwd=3)
            }else{
                segments(x0=1:length(nodes) - .3,
                         x1=1:length(nodes) + .3,
                         y0=true.mat[nodes,tt.vec[ii]],col="green",lwd=3)
            }
        }
        if(horizontal){
            segments(y0=1:length(nodes) - .3,
                     y1=1:length(nodes) + .3,
                     x0=object$SS.prior.mean[nodes,tclass],col="blue",lwd=3)
        }else{
            segments(x0=1:length(nodes) - .3,
                     x1=1:length(nodes) + .3,
                     y0=object$SS.prior.mean[nodes,tclass],col="blue",lwd=3)
        }
    }
}



ReceiverMat.plot <- function(object, tclass, true.mat, horizontal=FALSE, nodes, block, ...){

    kk <- dim(object$BB)[1]
    nn <- dim(object$RR)[1]
    TT <- dim(object$net.mat)[3]

    if(missing(nodes)){
        if(missing(block)){
            nodes <- 1:nn
        }else{
            if(block > kk) stop("block must be between 1 and ",kk)
            nodes <- which(object$mmb == block)
        }
    }else{
        if(any(nodes > nn) | any(nodes < 1)){
            stop("nodes must contain integers between 1 and ",nn)
        }
        if(!missing(block)){
            warning("The argument for block is ignored when nodes is used.")
        }
    }

    tt.vec <- which(object$tmap == tclass)
    if(horizontal){
        par(mfrow=c(1,length(tt.vec)))
    }else{
        par(mfrow=c(length(tt.vec),1))
    }
    for(ii in 1:length(tt.vec)){
        boxplot(t(object$chain$RR.mat[nodes,,tt.vec[ii]]), names=nodes,horizontal=horizontal,
                xaxt='n', ...)

        if(!missing(true.mat)){
            if(horizontal){
                segments(y0=1:length(nodes) - .3,
                         y1=1:length(nodes) + .3,
                         x0=true.mat[nodes,tt.vec[ii]],col="green",lwd=3)
            }else{
                segments(x0=1:length(nodes) - .3,
                         x1=1:length(nodes) + .3,
                         y0=true.mat[nodes,tt.vec[ii]],col="green",lwd=3)
            }
        }
        if(horizontal){
            segments(y0=1:length(nodes) - .3,
                     y1=1:length(nodes) + .3,
                     x0=object$RR.prior.mean[nodes,tclass],col="blue",lwd=3)
        }else{
            segments(x0=1:length(nodes) - .3,
                     x1=1:length(nodes) + .3,
                     y0=object$RR.prior.mean[nodes,tclass],col="blue",lwd=3)
        }
    }
}


BlockPrior.plot <- function(object, remap, true.mat ,
                            senders, receivers, ...){

    kk <- dim(object$BB)[1]
    nn <- dim(object$SS)[1]
    ee <- dim(object$SS.prior)[3]

    if(missing(senders)){
        senders <- 1:kk
    }
    if(missing(receivers)){
        receivers <- 1:kk
    }

    if(missing(remap)){
        remap <- 1:kk
    }else{
        if(length(remap) != kk){
            warning("remap has incorrect length.  Using default of 1:kk")
            remap <- 1:kk
        }
    }

    old.par <- par(mfrow=c(length(senders),length(receivers)),oma=c(0,0,2,0),
                   cex.axis=1.5,cex.main=2)
    for(ii in senders){
        for(jj in receivers){
            boxplot(object$chain$BB.prior[remap[ii],remap[jj],1,,] /
                    object$chain$BB.prior[remap[ii],remap[jj],2,,] ,
                    main=paste0("B[",remap[ii],",",remap[jj],"]"), ...)

            if(!missing(true.mat)){
                segments(x0=1:ee - .3,
                         x1=1:ee + .3,
                         y0=true.mat[ii,jj,1,]/true.mat[ii,jj,2,],col="green",lwd=3)
            }
        }
    }
    title("Block Intensity Posteriors",outer=TRUE)
    par(old.par)
}



#####  Function for plotting posterior Gamma Distributions
density.mat.plot <- function(mat,col="steelblue2",trans=0.3,xlim, add=FALSE, ...){
    TT <- ncol(mat)
    first <- !add
    upper.d <- min(mean(mat)*10,max(mat))
    if(missing(xlim)){
        xlim <- c(0,upper.d)
    }
    for(ii in 1:TT){
        if(first){
            plot(density(mat[,ii],from=0,to=upper.d,n=1e4),col=scales::alpha(col,trans),
                 xlim=xlim, ...)
            first <- FALSE
        }else{
            lines(density(mat[,ii],from=0,to=upper.d,n=1e4),col=scales::alpha(col,trans), ...)
        }
    }
}

#####  Function for drawing from posterior predictive Gamma Distributions
gamma.post.pred <- function(samp,iter.draws=100){
    total <- nrow(samp)
    out <- rep(NA,iter.draws*total)
    for(ii in 1:total){
        draws <- rgamma(iter.draws,samp[ii,1],samp[ii,2])
        out[((ii-1)*iter.draws + 1):(ii*iter.draws)] <- draws
    }
    return(out)
}


#####  Wrapper Functions to perform the above functions for specific classes
#####  and indices within the block probability matrix.
post.pred.block <- function(object,tclass,ss.ind,rr.ind,
                            remap,iter.draws=100){

    kk <- dim(object$BB)[1]
    total <- dim(object$chain$BB.prior)[4]

    if(missing(ss.ind)){
        ss.ind=1:kk
    }
    if(missing(rr.ind)){
        rr.ind = 1:kk
    }
    if(missing(remap)){
        remap=1:kk
    }

    samp <- array(NA,c(length(ss.ind), length(rr.ind), total*iter.draws))
    BB.prior <- object$chain$BB.prior[,,,,tclass]

    for(ii in 1:length(ss.ind)){
        ss <- ss.ind[ii]
        for(jj in 1:length(rr.ind)){
            rr <- rr.ind[jj]
            samp[ii,jj,] <- gamma.post.pred(t(BB.prior[remap[ss],remap[rr],,]),iter.draws=iter.draws)
        }
    }

    return(samp)

}



#####  Wrapper Functions to perform the above functions for specific classes
#####  and indices within the block probability matrix.
post.pred.receiver <- function(object,tclass,nodes,iter.draws=100){

    nn <- dim(object$RR)[1]
    total <- dim(object$chain$RR.prior)[3]

    if(missing(nodes)){
        nodes <- 1:nn
    }

    samp <- array(NA,c(length(nodes), total*iter.draws))
    prior <- object$chain$RR.prdior[,,,tclass]

    for(ii in 1:length(nodes)){
        node <- nodes[ii]
        samp[ii,] <- gamma.post.pred(t(prior[node,,]),iter.draws=iter.draws)
    }

    return(samp)
}


post.pred.sender <- function(object,tclass,nodes,iter.draws=100){

    nn <- dim(object$SS)[1]
    total <- dim(object$chain$SS.prior)[3]

    if(missing(nodes)){
        nodes <- 1:nn
    }

    samp <- array(NA,c(length(nodes), total*iter.draws))
    prior <- object$chain$SS.prior[,,,tclass]

    for(ii in 1:length(nodes)){
        node <- nodes[ii]
        samp[ii,] <- gamma.post.pred(t(prior[node,,]),iter.draws=iter.draws)
    }

    return(samp)
}

BlockMat.Post.Density.Plot <- function(object,tclass,ss.ind=1,rr.ind=1,
                                       remap,col="steelblue2",trans, add=FALSE,
                                       include.pred=FALSE, iter.draws=100, pred.lwd=2,
                                       ...){

    if(length(ss.ind) > 1){
        stop("ss.ind must be a single index")
    }
    if(length(rr.ind) > 1){
        stop("rr.ind must be a single index")
    }
    if(missing(remap)){
        kk <- dim(object$BB)[1]
        remap=1:kk
    }

    tt.vec <- which(object$tmap == tclass)

    if(missing(trans)){
        trans <- 0.5 / sqrt(length(tt.vec))
    }

    density.mat.plot(object$chain$BB.mat[remap[ss.ind],remap[rr.ind],,tt.vec],
                     col=col,trans=trans,add=add, ...)
    if(include.pred){
        samp <- post.pred.block(object=object,tclass=tclass,
                                ss.ind=ss.ind, rr.ind=rr.ind, remap=remap,
                                iter.draws=iter.draws)
        abline(v=quantile(samp,c(0.025,0.975)),col=col,lty=2,lwd=pred.lwd)
        lines(density(samp,from=0,n=1e4),col=col,lwd=pred.lwd)
    }

}


Sender.Post.Density.Plot <- function(object,tclass,node,
                                     col="steelblue2",trans, add=FALSE,
                                     include.pred=FALSE, iter.draws=100,
                                     ...){

    if(length(node) > 1){
        stop("node must be a single index")
    }

    tt.vec <- which(object$tmap == tclass)

    if(missing(trans)){
        trans <- 0.5 / sqrt(length(tt.vec))
    }

    density.mat.plot(object$chain$SS.mat[node,,tt.vec],
                     col=col,trans=trans,add=add, ...)
    if(include.pred){
        samp <- post.pred.sender(object=object,tclass=tclass,
                                nodes=node,
                                iter.draws=iter.draws)
        abline(v=quantile(samp,c(0.025,0.975)),col=col,lty=2)
        lines(density(samp,from=0,n=1e4),col=col,lwd=2)
    }

}




Receiver.Post.Density.Plot <- function(object,tclass,node,
                                       col="steelblue2",trans, add=FALSE,
                                       include.pred=FALSE, iter.draws=100,
                                       ...){

    if(length(node) > 1){
        stop("node must be a single index")
    }

    tt.vec <- which(object$tmap == tclass)

    if(missing(trans)){
        trans <- 0.5 / sqrt(length(tt.vec))
    }

    density.mat.plot(object$chain$RR.mat[node,,tt.vec],
                     col=col,trans=trans,add=add, ...)
    if(include.pred){
        samp <- post.pred.receiver(object=object,tclass=tclass,
                                nodes=node,
                                iter.draws=iter.draws)
        abline(v=quantile(samp,c(0.025,0.975)),col=col,lty=2)
        lines(density(samp,from=0,n=1e4),col=col,lwd=2)
    }

}



#####  Wrapper Functions to perform the above functions for specific classes
#####  and indices within the block probability matrix.
post.pred.receiver <- function(object,tclass,nodes,iter.draws=100){

    nn <- dim(object$RR)[1]
    total <- dim(object$chain$RR.prior)[3]

    if(missing(nodes)){
        nodes <- 1:nn
    }

    samp <- array(NA,c(length(nodes), total*iter.draws))
    prior <- object$chain$RR.prior[,,,tclass]

    for(ii in 1:length(nodes)){
        node <- nodes[ii]
        samp[ii,] <- gamma.post.pred(t(prior[node,,]),iter.draws=iter.draws)
    }

    return(samp)
}


post.pred.sender <- function(object,tclass,nodes,iter.draws=100){

    nn <- dim(object$SS)[1]
    total <- dim(object$chain$SS.prior)[3]

    if(missing(nodes)){
        nodes <- 1:nn
    }

    samp <- array(NA,c(length(nodes), total*iter.draws))
    prior <- object$chain$SS.prior[,,,tclass]

    for(ii in 1:length(nodes)){
        node <- nodes[ii]
        samp[ii,] <- gamma.post.pred(t(prior[node,,]),iter.draws=iter.draws)
    }

    return(samp)
}


