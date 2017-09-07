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


#' Create Dynamic Stochastic Block Model Object
#'
#' Generates a Model Object for a Dynamic Stochastic Block Model
#'
#' @param nn number of nodes in the network
#' @param TT number of networks to generate
#' @param tmap mapping from time points 1:TT to equivalence classes
#' @param hours.vec number of hours for each equivalence class
#' @param BB.prior block effect prior parameters for each equivalence class
#' @param SS.prior sender effect prior parameters for each equivalence class
#' @param RR.prior receiver effect prior parameters for each equivalence class
#' @param mmb.prior prior distribution for multinomial distribution over
#' block memberships
#' @param mmb block memberhsip vector
#' @param self.ties if TRUE self ties are allowed
#' @param normalize if TRUE, sender/receiver effects are normalized to have a
#' mean of 1 within each block
#' @param gen if TRUE also samples a network from the model and includes it as a
#' component of the model object
#'
#' @return Returns a "dynsbm" model object containing all the parameters needed
#' for sampling networks from the model.  If gen = TRUE, the object also
#' contains a network sampled from this parameterization.
#'
#' @examples
#' kk <- 3
#' nn <- kk^2 * 10
#' TT <- 10; ee <- 2
#' tmap <- rep(1:ee,TT/ee)
#' hours.vec <- c(12,12)
#' mmb <-  rep(1:kk,each=nn/kk)
#' BB.prior <- array(1,c(kk,kk,2,ee))
#' BB.prior[cbind(1:3,c(2,3,1),1,1)] <- 2
#' BB.prior[cbind(1:3,c(3,1,2),1,2)] <- 2
#'
#' SS.prior <- array(rgamma(nn*2*ee,10,10),c(nn,2,ee))
#' RR.prior <- array(rgamma(nn*2*ee,10,10),c(nn,2,ee))
#' dat <- dynsbm(nn=nn,TT=TT,tmap=tmap,hours.vec=hours.vec,
#'               BB.prior=BB.prior,SS.prior=SS.prior,RR.prior=RR.prior
#'               mmb=mmb,self.ties=TRUE,normalize=TRUE,gen=TRUE)
#'
#' @seealso \code{\link{dynsbm}}
#' @export
dynsbm <- function(nn=40,TT=4,tmap=rep(1:2,2),hours.vec=c(1,1),
                   BB.prior=NULL,SS.prior=NULL,RR.prior=NULL,
                   BB.mat,SS.mat,RR.mat,
                   mmb.prior=c(.5,.5),mmb,
                   self.ties=TRUE,normalize=FALSE,gen=FALSE){


    ## browser()
    ##  Getting Number of Equivalence Classes
    ee <- max(tmap)


    ##  Loading BB.mat
    if(missing(BB.mat)){
        if(is.null(BB.prior))
            stop("At least one of BB.prior or BB.mat must be provided")
        kk <- dim(BB.prior)[1]
        if(any(dim(BB.prior) != c(kk,kk,2,ee))){
            stop("BB.prior must have dimension (kk,kk,2,ee)")
        }

        BB.mat <- array(NA,c(kk,kk,TT))
        for(tt in 1:TT){
            BB.mat[,,tt] <- rgamma(kk^2,BB.prior[,,1,tmap[tt]],
                                   BB.prior[,,2,tmap[tt]])
        }
    }else{
        kk <- dim(BB.mat)[1]
        if(any(dim(BB.mat) != c(kk,kk,TT))){
            stop("BB.mat must have dimension (kk,kk,TT)")
        }
    }

    ## Loading and checking mmb
    if(missing(mmb)){
        mmb <- sort(sample(kk,nn,replace=TRUE,prob=mmb.prior))
    }else{
        if(length(mmb) != nn) stop("mmb must have length n")
        if(max(mmb) > kk) stop("mmb must contain integers between 1 and kk")
    }

    ## Loading SS
    if(missing(SS.mat)){
        if(is.null(SS.prior))
            stop("At least one of SS.prior or SS.mat must be provided")
        if(any(dim(SS.prior) != c(nn,2,ee))){
            stop("SS.prior must have dimension (nn,2,ee)")
        }

        SS.mat <- array(NA,c(nn,TT))
        for(tt in 1:TT){
            SS.mat[,tt] <- rgamma(nn,SS.prior[,1,tmap[tt]],
                                  SS.prior[,2,tmap[tt]])
        }
    }else{
        if(any(dim(SS.mat) != c(nn,TT))){
            stop("SS.mat must have dimension (nn,TT)")
        }
    }

    ## Loading RR
    if(missing(RR.mat)){
        if(is.null(RR.prior))
            stop("At least one of RR.prior or RR.mat must be provided")

        if(any(dim(RR.prior) != c(nn,2,ee))){
            stop("RR.prior must have dimension (nn,2,ee)")
        }

        RR.mat <- array(NA,c(nn,TT))
        for(tt in 1:TT){
            RR.mat[,tt] <- rgamma(nn,RR.prior[,1,tmap[tt]],
                                  RR.prior[,2,tmap[tt]])
        }
    }else{
        if(any(dim(RR.mat) != c(nn,TT))){
            stop("RR.mat must have dimension (nn,TT)")
        }
    }


    if(normalize){
        for(tt in 1:TT){
            mmb.count <- ss.norm <- rr.norm <- rep(NA,kk)
            for(ll in 1:kk){
                ind.ll <- which(mmb == ll)
                mmb.count[ll] <- length(ind.ll)
                ss.norm[ll] <- sum(SS.mat[ind.ll,tt])
                rr.norm[ll] <- sum(RR.mat[ind.ll,tt])

                ssnorm <- (mmb.count[ll] / ss.norm[ll])
                SS.mat[ind.ll,tt] <- SS.mat[ind.ll,tt] * ssnorm
                rrnorm <- (mmb.count[ll] / rr.norm[ll])
                RR.mat[ind.ll,tt] <- RR.mat[ind.ll,tt] * rrnorm
            }
        }
    }


    obj <- structure(list(mmb=mmb,
                          BB.prior=BB.prior,
                          SS.prior=SS.prior,RR.prior=RR.prior,
                          BB.mat=BB.mat,SS.mat=SS.mat,RR.mat=RR.mat,
                          hours.vec=hours.vec,tmap=tmap,
                          self.ties=self.ties,
                          nn=nn,kk=kk,TT=TT,ee=ee),
                     class="dynsbm")

    if(gen){
        obj$net.mat <- net.gen(obj)
    }

    return(obj)

}


#' @export
is.dynsbm <- function(x) inherits(x,"dynsbm")


#' Generate Set of Dynamic Stochastic Block Model Networks
#'
#' Generates networks from the weighted stochastic block model
#' with or without degree effects.
#'
#' @param nn number of nodes in the network
#' @param TT number of networks to generate
#' @param tmap mapping from time points 1:TT to equivalence classes
#' @param hours.vec number of hours for each equivalence class
#' @param BB.prior block effect prior parameters for each equivalence class
#' @param SS.prior sender effect prior parameters for each equivalence class
#' @param RR.prior receiver effect prior parameters for each equivalence class
#' @param mmb.prior prior distribution for multinomial distribution over
#' block memberships
#' @param mmb block memberhsip vector
#' @param self.ties if TRUE self ties are allowed
#' @param normalize if TRUE, sender/receiver effects are normalized to have a
#' mean of 1 within each block
#'
#' @return Returns a list containing the simulated matrix of networks as the
#' first component 'net.mat'. The list also contains all of the parameters
#' used to simulate the networks.
#'
#' @examples
#' kk <- 3
#' nn <- kk^2 * 10
#' TT <- 10; ee <- 2
#' tmap <- rep(1:ee,TT/ee)
#' hours.vec <- c(12,12)
#' mmb <-  rep(1:kk,each=nn/kk)
#' BB.prior <- array(1,c(kk,kk,2,ee))
#' BB.prior[cbind(1:3,c(2,3,1),1,1)] <- 2
#' BB.prior[cbind(1:3,c(3,1,2),1,2)] <- 2
#'
#' SS.prior <- array(rgamma(nn*2*ee,10,10),c(nn,2,ee))
#' RR.prior <- array(rgamma(nn*2*ee,10,10),c(nn,2,ee))
#' dat <- data.gen.dynsbm(nn=nn,TT=TT,tmap=tmap,hours.vec=hours.vec,
#'                        BB.prior=BB.prior,SS.prior=SS.prior,RR.prior=RR.prior
#'                        mmb=mmb,self.ties=TRUE,normalize=TRUE)
#'
#' @seealso \code{\link{dynsbm}}
#' @export
net.gen.dynsbm <- function(x,marginal=FALSE){

    nn <- x$nn; TT <- x$TT; tmap <- x$tmap; kk <- x$kk

    if(marginal){
        BB.mat <- array(NA,c(kk,kk,TT))
        for(tt in 1:TT){
            BB.mat[,,tt] <- rgamma(kk^2,x$BB.prior[,,1,tmap[tt]],
                                   x$BB.prior[,,2,tmap[tt]])
        }

        SS.mat <- array(NA,c(nn,TT))
        for(tt in 1:TT){
            SS.mat[,tt] <- rgamma(nn,x$SS.prior[,1,tmap[tt]],
                                  x$SS.prior[,2,tmap[tt]])
        }

        RR.mat <- array(NA,c(nn,TT))
        for(tt in 1:TT){
            RR.mat[,tt] <- rgamma(nn,x$RR.prior[,1,tmap[tt]],
                                  x$RR.prior[,2,tmap[tt]])
        }

        x$BB.mat <- BB.mat
        x$SS.mat <- SS.mat
        x$RR.mat <- RR.mat
    }

    wsbm.list <- create.wsbm.list(x)

    ##  Creating containers for network and parameters
    net.mat <- array(NA,c(x$nn,x$nn,x$TT))

    for(tt in 1:x$TT){
        net.mat[,,tt] <- net.gen(wsbm.list[[tt]])
    }
    return(net.mat)
}


#' Dynamic Weighted Stochastic Block Model Estimator
#'
#' Estimates parameters for the dynamic weighted stochastic block model using
#' an MCMC algorithm.
#'
#' @param net.mat a 3-dimensional array containing a sequence of adjacency
#' matrices.
#' @param kk number of blocks
#' @param tmap mapping between networks and equivalence classes
#' @param hours.vec vector of scaling parameters for each equivalence class
#' @param self.ties if true assumes self ties are possible
#' @param priors priors for posterior inference.  See dynsbm.priors for more
#' details
#' @param init.control control parameters for initialization routine.
#' If init.control contains a member named init.vals it is interpreted as a
#' list of parameters from which to start the MCMC chain.
#' @param mcmc.control control parameters for the MCMC algorithm
#' @param clean.out if true, removes MCMC chain from output and only
#' returns posterior means
#' @param verbose higher values correspond to more informative output as the
#' sampler runs.  Setting verbose = 0 suppresses all output.
#'
#' @return Returns a dynsbm object that has posterior means for each parameter
#' in the dynamic, weighted, directed,  degree corrected stochastic block model.
#' If clean.out is FALSE, this object also contains an element 'chain' which
#' contains the thinned and burned in draws from the MCMC sampler.  The prior
#' parameters from the call to dynsbm are also included.
#'
#' @seealso \code{\link{dynsbm.priors}}
#' @export
dynsbm.fit <- function(net.mat, kk=3, tmap, hours.vec, self.ties=TRUE,
                   priors=dynsbm.priors(eta=rep(1/kk,kk)),
                   init.control=list(spectral.start=FALSE,
                                     multistart=0,
                                     multistart.total=10),
                   mcmc.control=list(total=1000,burn.in=1000,thin=10,
                                     extend.alpha=.001,extend.max=10,
                                     label.switch.mode="adhoc",
                                     label.switch.max=200,
                                     multi.impute=FALSE),
                   update.mmb=TRUE,clean.out=FALSE, verbose=1){


    ##  Parsing Control Parameters and Filling with Defaults
    mcmc.control <- parse.mcmc.control(mcmc.control)
    init.control <- parse.init.control(init.control)

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
    multi.int <- as.integer(ifelse(mcmc.control$multiImpute,1,0))
    net.mat.clean <- net.mat
    if(!self.ties){
        for(tt in 1:TT) diag(net.mat.clean[,,tt]) <- NA
    }
    net.mat.clean[is.na(net.mat.clean)] <- -1


##################################################
######## Creating Containers for Output ##########
##################################################


    total <- mcmc.control$total
    thin <- mcmc.control$thin
    burn.in <- mcmc.control$burn.in

    ## total <- short.total*thin + burn.in ## Actual Number of Iterations

    ##  Storage for Priors
    flatBB.prior = double(total * 2 * kk*kk * ee)
    flatSS.prior = double(total * 2 * nn * ee)
    flatRR.prior = double(total * 2 * nn * ee)

    ##  Storage for Parameters
    flatSS = double(total * nn * TT)
    flatRR = double(total * nn * TT)
    flatBB = double(total * kk*kk * TT)
    flatMMB = integer(total * nn)

    ##  Storage for Performance Indicators
    ll.vec <- double(total*TT)
    flatHH <- double(total*kk*nn)



##################################################
############# Initializing Values  ###############
##################################################

    ##  Initializing
    dynsbm.init <- dynsbm.load.init.vals(init.control=init.control,
                                         priors=priors,
                                         net.mat=net.mat,kk=kk,
                                         hours.vec=hours.vec,tmap=tmap,
                                         self.ties=self.ties)
    dynsbm.init$net.mat <- net.mat


    ##  Flattening Initial Values for Pass to C++
    start = 1
    for(tt in 1:TT){
        flatBB[1:(kk*kk) + (tt-1)*total*kk^2] <- dynsbm.init$BB.mat[,,tt]
        flatSS[1:nn + (tt-1)*total*nn] <- dynsbm.init$SS.mat[,tt]
        flatRR[1:nn + (tt-1)*total*nn] <- dynsbm.init$RR.mat[,tt]
    }

    for(tt in 1:ee){
        BB.add <- (tt-1)*total*2*kk*kk
        flatBB.prior[1:(kk^2) + BB.add] <- dynsbm.init$BB.prior[,,1,tt]
        flatBB.prior[1:(kk^2) + (kk^2) + BB.add] <- dynsbm.init$BB.prior[,,2,tt]

        SR.add <- (tt-1)*total*2*nn
        flatSS.prior[1:nn + SR.add] <- dynsbm.init$SS.prior[,1,tt]
        flatSS.prior[1:nn + nn + SR.add] <- dynsbm.init$SS.prior[,2,tt]

        flatRR.prior[1:nn + SR.add] <- dynsbm.init$RR.prior[,1,tt]
        flatRR.prior[1:nn + nn + SR.add] <- dynsbm.init$RR.prior[,2,tt]
    }
    flatMMB[1:nn] <- dynsbm.init$mmb
    ll.vec[1] <- logLik(dynsbm.init)
    flatHH[1:(kk*nn)] <- mmb.to.PI(dynsbm.init$mmb,kk)

    ##  MCMC Convergence Checking Parameters
    extend.max <- as.integer(mcmc.control$extend.max)
    ## shift.size <- 0
    extend.qq <- as.double(qnorm(1 - mcmc.control$extend.alpha/2))

##################################################
############  Call to C++ Function  ##############
##################################################
    #' @useDynLib lvnm dynsbm_R
    out <- .C(dynsbm_R,
              as.integer(nn),                       # 1: Nodes in Networks
              as.integer(kk),                       # 2: Number of Blocks
              as.integer(TT), as.integer(ee),       # 3-4: Time / Equiv Pts
              multi.int,                            # 5: Multi Impute Indicator

              as.integer(net.mat.clean),            # 6: Networks
              as.integer(tmap),                     # 7: Time/Equiv Map
              as.double(hours.vec),                 # 8: Hours Per Equiv

              as.double(priors$SS),                 # 9: Sender Hyperprior
              as.double(priors$RR),                 # 10: Receiver Hyperprior
              as.double(priors$BB),                 # 11: BlockMatrix Hyperprior
              as.double(priors$eta),                # 12: BlockMemb Prior

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
              extend.max, extend.qq)                # 28-29: Autoconverge Control


##################################################
###########  Processing C++ Output  ##############
##################################################
    SS.prior <- array(out[[13]],c(nn,2,total,ee))
    RR.prior <- array(out[[14]],c(nn,2,total,ee))
    BB.prior <- array(out[[15]],c(kk,kk,2,total,ee))
    MMB <- array(out[[16]],c(nn,total))

    SS.mat <- array(out[[17]],c(nn,total,TT))
    RR.mat <- array(out[[18]],c(nn,total,TT))
    BB.mat <- array(out[[19]],c(kk,kk,total,TT))
    HH.mat <- array(out[[21]],c(nn,kk,total))

    ll.mat <- array(out[[20]],c(TT,total))
    ll.vec <- colSums(ll.mat) - sum(lfactorial(net.mat))

    ##  Combining output into single structure
    dynsbm.chain.obj <- structure(list(mmb=MMB,BB.mat=BB.mat,
                                       SS.mat=SS.mat,RR.mat=RR.mat,HH=HH.mat,
                                       BB.prior=BB.prior,
                                       SS.prior=SS.prior,RR.prior=RR.prior,
                                       logLik=ll.vec,ll.mat=ll.mat,
                                       tmap=tmap,hours.vec=hours.vec,
                                       self.ties=self.ties,
                                       nn=nn,kk=kk,TT=TT,ee=ee,
                                       init.control=init.control,
                                       mcmc.control=mcmc.control,priors=priors,
                                       clean=FALSE),
                                  class="dynsbm.chain")

    ##  Cleaning up output from C++ Call
    rm(out)

######  Label Switching Handler
    if(mcmc.control$label.switch.mode == "kl-loss"){
        ## dynsbm.out <- switch.labels(dynsbm.out,
        ##                             max.runs=mcmc.control$max.runs,
        ##                             verbose=verbose)
        message("kl-loss mode is currently disabled for dynsbm")
    }else{
        ##  Do nothing
    }

##################################################
########  Generating Summary Statistics  #########
##################################################

    dynsbm.obj <- summary(dynsbm.chain.obj)
    dynsbm.obj$net.mat <- net.mat

    dynsbm.fit.obj <- dynsbm.obj
    dynsbm.fit.obj$chain <- dynsbm.chain.obj
    class(dynsbm.fit.obj) <- c("dynsbm.fit","dynsbm.mcmc","dynsbm")

    if(clean.out){
        dynsbm.fit.obj$chain$mmb <- NULL
        dynsbm.fit.obj$chain$BB.mat <- NULL
        dynsbm.fit.obj$chain$SS.mat <- NULL
        dynsbm.fit.obj$chain$RR.mat <- NULL
        dynsbm.fit.obj$chain$HH.mat <- NULL

        dynsbm.fit.obj$chain$BB.prior <- NULL
        dynsbm.fit.obj$chain$SS.prior <- NULL
        dynsbm.fit.obj$chain$RR.prior <- NULL
        dynsbm.fit.obj$chain$clean <- TRUE
    }

    return(dynsbm.fit.obj)
}


#' Prior Parameters for dynsbm
#'
#' Generating prior parameters for pass to dynsbm
#'
#' @param eta prior on block memberships
#' @param block.alpha shape parameter for Block Effects
#' @param block.beta rate parameter for Block Effects
#' @param sender.alpha shape parameter for Sender Effects
#' @param sender.beta rate parameter for Sender Effects
#' @param receiver.alpha shape parameter for Receiver Effects
#' @param receiver.beta rate parameter for Receiver Effects
#'
#' @return Returns a list containing the parameters with useful defaults.  eta
#' is the only required parameter.
#'
#' @seealso \code{\link{dynsbm}}
#' @export
dynsbm.priors <- function(eta,
                          BB.hyperprior=c(.99,1.01,1,1),
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
                RR=RR.hyperprior,
                eta=eta))
}


######################################################################
#####################  PREDICTION FUNCTIONS  #########################
######################################################################


#' Predict Method for DynSBM Fits
#'
#' Generates a matrix of predicted intensities for each node pair
#'
#' @param object a fitted object of the class "dynsbm"
#'
#' @return Returns an array of predicted intensities.  These intensities are
#' the predicted values given the posterior mean parameter estimates.
#'
#' @examples
#' data(commuter30)
#' fit.wsbm <- wsbm(commuter30,kk=3,hours=1)
#' pred.mat <- predict(fit.wsbm)
#'
#' @seealso \code{\link{dynsbm}}
#' @export
predict.dynsbm <- function(x, marginal=TRUE, reps = 1000, ...){

    if(marginal){
        pmm <- array(NA,c(x$nn,x$nn,x$ee,reps))
        x$TT <- x$ee
        x$tmap <- 1:x$ee
        for(ii in 1:reps){
            pmm[,,,ii] <- net.gen(x,marginal=marginal)
        }
        pmat <- apply(pmm,c(1,2,3),mean)
        return(pmat)
    }else{
        wsbm.list <-  create.wsbm.list(x)
        TT <- x$TT
        pmat <- array(NA,dim(x$net.mat))

        for(tt in 1:TT){
            pmat[,,tt] <- predict(wsbm.list[[tt]])
        }
        return(pmat)
    }
}



#' Posterior Predictive Distribution Method for DynSBM Fits
#'
#' Generates a posterior predictive distribution for each tie in a network.
#'
#' @param x a fitted object of the class "dynsbm".  Requires clean.out = FALSE
#' in the call to dynsbm.
#' @param marginal if TRUE, values block, sender, and receiver effects are
#' marginalized over using the estimated prior parameters.
#'
#' @return Returns a 3 dimensional array of dimension n x n x total, where n
#' is the number of nodes in the network and total is the number of draws in
#' the MCMC chain generated using wsbm.
#'
#' @examples
#' data(commuter30)
#' fit.wsbm <- wsbm(commuter30,kk=3,hours=1)
#'
#' pp.dist <- post.pred(fit.wsbm)
#' density.ppd <- apply(post.pred,3,mean)
#' summary(density.ppd)
#'
#' ##  Initializing with a specific block membership
#' init.vals <- list(mmb=c(rep(1,15),rep(2,10),rep(3,5)))
#' fit.wsbm.2 <- wsbm(commuter30,kk=3,hours=1,
#'                    init.control=list(init.vals=init.vals))
#' fit.wsbm.2
#'
#' @seealso \code{\link{dynsbm}} \code{\link{predict.dynsbm}}
#' @export
post.predict.dynsbm.mcmc <- function(x, marginal=TRUE, ...){
    if(x$chain$clean) stop("Use clean.out = TRUE to keep MCMC chains")
    nn <- x$nn; TT <- x$TT; ee <- x$ee
    total <- x$chain$mcmc.control$total


    if(marginal) pmat.post <- array(NA,c(nn,nn,ee,total))
    else pmat.post <- array(NA,c(nn,nn,TT,total))

    for(ii in 1:total){
        tmp.obj <- get.iter(x,ii)
        if(marginal){
            tmp.obj$TT <- ee
            tmp.obj$tmap <- 1:ee
        }

        pmat.post[,,,ii] <- net.gen(tmp.obj,marginal=marginal)
    }

    return(pmat.post)
}




#################################################################
###################  dynsbm class functions  ####################
#################################################################


summary.dynsbm.chain <- function(x,...){
###    browser()
    TT <- dim(x$net.mat)[3]
    total <- dim(x$BB.mat)[3];
    nn <- x$nn; TT <- x$TT; ee <- x$ee; kk <- x$kk

###  Summary of Parameters
    BB.hat <- apply(x$BB.mat,c(1,2,4),mean)
    SS.hat <- apply(x$SS.mat,c(1,3),mean)
    RR.hat <- apply(x$RR.mat,c(1,3),mean)

    PI.mean <- t(apply(x$mmb,1,tabulate,nbins=kk) / total)
    mmb <- apply(PI.mean,1,which.max)

###  Summary of Priors
    BB.prior.hat <- apply(x$BB.prior,c(1,2,3,5),mean)
    SS.prior.hat <- apply(x$SS.prior,c(1,2,4),mean)
    RR.prior.hat <- apply(x$RR.prior,c(1,2,4),mean)


    dynsbm.obj <- dynsbm(nn=nn,TT=TT,mmb=mmb,
                         BB.prior=BB.prior.hat,BB.mat=BB.hat,
                         SS.prior=SS.prior.hat,SS.mat=SS.hat,
                         RR.prior=RR.prior.hat,RR.mat=RR.hat,
                         tmap=x$tmap,hours.vec=x$hours.vec,
                         self.ties=x$self.ties)

    return(dynsbm.obj)

    ##  Calculating Posterior Mean/Variance
    BB.prior.mean <- apply(x$BB.prior[,,1,,,drop=FALSE]/
                           x$BB.prior[,,2,,,drop=FALSE],c(1,2,5),mean)
    BB.prior.var <- apply(x$BB.prior[,,1,,,drop=FALSE]/
                          (x$BB.prior[,,2,,,drop=FALSE]^2),c(1,2,5),mean)

    SS.prior.mean <- apply(x$SS.prior[,1,,,drop=FALSE]/
                           x$SS.prior[,2,,,drop=FALSE],c(1,4),mean)
    SS.prior.var <- apply(x$SS.prior[,1,,,drop=FALSE]/
                          (x$SS.prior[,2,,,drop=FALSE]^2),c(1,4),mean)

    RR.prior.mean <- apply(x$RR.prior[,1,,,drop=FALSE]/
                           x$RR.prior[,2,,,drop=FALSE],c(1,4),mean)
    RR.prior.var <- apply(x$RR.prior[,1,,,drop=FALSE]/
                          (x$RR.prior[,2,,,drop=FALSE]^2),c(1,4),mean)
}


#' Log-Likelihood method for DynSBM Fits
#'
#' Returns log-likelihood for posterior mean estimates.
#'
#' @param x an object of class "dynsbm"
#'
#' @return returns the log-likelihood of the data given the posterior mean
#' parameter estimates.
#'
#' @seealso \code{\link{dynsbm}}
#' @export
logLik.dynsbm <- function(x){
    return(sum(sapply(create.wsbm.list(x),logLik)))
}



#' AIC method for DynSBM Fits
#'
#' Returns Akaike Information Criterion for DynSBM Fit
#'
#' @param x an object of class "dynsbm"
#'
#' @return returns the Akaike Information Criterion using the posterior mean
#' parameters as approximations to the MLE.  The number of parameters is taken
#' to be T*(2*(n-k) + k*k) + n.
#'
#' @seealso \code{\link{dynsbm}}
#' @export
AIC.dynsbm <- function(x){
    TT <- x$TT # dim(x$net.mat)[3]
    nn <- x$nn # dim(x$net.mat)[1]
    kk <- x$kk # dim(x$BB)[1]
    edf <- TT*(2*(nn-kk) + kk^2) +nn
    return(-2*logLik(x) + 2 * edf)
}



#' BIC method for DynSBM Fits
#'
#' Returns Bayesian Information Criterion for DynBM Fit
#'
#' @param x an object of class "dynsbm"
#'
#' @return returns the Bayesian Information Criterion using the posterior mean
#' parameters as approximations to the MLE.  The number of parameters is taken
#' to be T*(2*(n-k) + k*k) + n and the sample size is taken to be T*n*n.
#'
#' @seealso \code{\link{dynsbm}}
#' @export
BIC.dynsbm <- function(x){
    TT <- x$TT # dim(x$net.mat)[3]
    nn <- x$nn # dim(x$net.mat)[1]
    kk <- x$kk # dim(x$BB)[1]
    edf <- TT*(2*(nn-kk) + kk^2) +nn
    return(-2*logLik(x) + (2*log(nn) + log(TT)) * edf)
}



#' DIC method for DynSBM Fits
#'
#' Returns Deviance Information Criterion for DynSBM Fit
#'
#' @param x an object of class "dynsbm"
#'
#' @return returns the Deviance Information Criterion.
#'
#' @seealso \code{\link{dynsbm}}
#' @export
DIC.dynsbm <- function(x){
    if(!is.null(x$chain$logLik)){
        DIC.df <- (2*logLik(x) - 2*mean(x$chain$logLik))
        return(-2*logLik(x) + 2 * DIC.df)
    }else{
        return(NULL)
    }
}


## dynsbm.metric <- function(graph,kk=2,total=1500,
##                        thin=1,burn.in=500,verbose=0,...){

##     ##    mode <- match.arg(mode)
##     ##    if(mode == "mcmc"){
##     dynsbm.fit <- dynsbm(total=total,net=graph,kk=kk,verbose=verbose,
##                    thin=thin,burn.in=burn.in)
##     ##    }
##     return(dynsbm.fit$pmat)
## }



############################################################
#################  FORMATTING FUNCTIONS  ###################
############################################################

#' @export
print.dynsbm <- function(x, ...){
    cat(format(x, ...), "\n")
}

#' @export
format.dynsbm <- function(x, digits=4, max.width=78,...){
    head.str <- "Posterior Sample for DynSBM"
    time.str <- paste("Network Count:",dim(x$net.mat)[3])
    map.str <- paste0("Equivalence Classes:\n",
                      format.vector(x$tmap,digits=digits,max.width=max.width,...),
                      "\n")

    node.str <- paste("Nodes:",dim(x$net.mat)[1])
    block.str <- paste("Blocks:",dim(x$BB)[1])
    mmb.str <- paste0("Block Membership:\n",
                     format.vector(x$mmb,digits=digits,max.width=max.width,...))

    ll.str <- paste0("\n","Log-Likelihood:  ",format(logLik(x),digits=7))
    ic.str <- paste("BIC: ",format(BIC(x),digits=7),
                     "   AIC: ",format(AIC(x), digits=7),
                     "   DIC: ",format(DIC(x),digits=7))

    return(paste(head.str,time.str,map.str,
                 node.str,block.str,mmb.str,
                 ll.str,ic.str,sep="\n"))

}


############################################################
###################  HELPER FUNCTIONS  #####################
############################################################

## create.wsbm.list <- function(net.mat,mmb,BB.mat,SS.mat,RR.mat,
##                              hours.vec,tmap,self.ties){
##     TT <- dim(net.mat)[3]
##     wsbm.list <- NULL
##     for(tt in 1:TT){
##         obj <- structure(list(mmb=mmb,
##                               BB=BB.mat[,,tt],
##                               SS=SS.mat[,tt],
##                               RR=RR.mat[,tt],
##                               hours=hours.vec[tmap[tt]],
##                               self.ties=self.ties,
##                               net=net.mat[,,tt]),
##                          class="wsbm")
##         wsbm.list[[tt]] <- obj
##     }
##     return(wsbm.list)

## }

#' @export
create.wsbm.list <- function(x){
    TT <- x$TT
    wsbm.list <- NULL
    for(tt in 1:TT){
        obj <- with(x,wsbm(nn=nn,mmb=mmb,BB=BB.mat[,,tt],
                           SS=SS.mat[,tt],RR=RR.mat[,tt],
                           self.ties=self.ties,hours=hours.vec[tmap[tt]]))
        if(!is.null(x$net.mat)){
            obj$net <- x$net.mat[,,tt]
        }
        wsbm.list[[tt]] <- obj
    }
    return(wsbm.list)
}


dynsbm.load.init.vals <- function(init.control,priors,net.mat,kk,
                                  hours.vec,tmap,self.ties){

    init.vals <- init.control$init.vals

    ##  Spectral Clustering Initialization Not Currently Supported
    if(is.null(init.vals) & init.control$spectral.start){
        ## if(verbose > 0) message("Initializing with Spectral Clustering...")
        ## init.vals <- sbm.spectral(net=net,kk=kk,weighted=TRUE)
        ## if(verbose > 0) message("Initialization Complete.")
    }

    nn <- dim(net.mat)[1]
    TT <- dim(net.mat)[3]
    ee <- length(unique(tmap))


    if(init.control$multistart > 0){
        ll.best <- -Inf
        multi.mcmc <- list(total=init.control$multistart.total,
                           burn.in=0,thin=1,extend.max=0,
                           label.switch.mode="adhoc")
        multi.init <- list(spectral.start=FALSE,multistart=0)
        for(ii in 1:init.control$multistart){


            ##  Update to new call specification
            fit.tmp <- dynsbm.fit(net.mat=net.mat,kk=kk,
                                  hours.vec=hours.vec,tmap=tmap,
                                  self.ties=self.ties,
                                  priors=priors,
                                  mcmc.control=multi.mcmc,
                                  init.control=multi.init,
                                  clean.out=FALSE,verbose=0)

            obj.tmp <- get.iter.dynsbm.mcmc(fit.tmp,
                                            init.control$multistart.total)

            ll.tmp <- logLik(obj.tmp)
            if(ll.tmp > ll.best){
                ll.best <- ll.tmp
                obj.best <- obj.tmp
            }
        }
        return(obj.best)
    }else{

        ##  Loading Initial Values
        if(is.null(init.vals$mmb)){
            init.vals$mmb <- sample(kk,nn,replace=TRUE,prob=priors$eta)
        }

        if(is.null(init.vals$BB.mat)){
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

        ##  Loading Hyperparameters
        if(is.null(init.vals$BB.prior)){
            init.vals$BB.prior <- array(NA,c(kk,kk,2,ee))
            for(tt in 1:ee){
                bp <- apply(init.vals$BB.mat[,,tmap==tt,drop=FALSE],
                            c(1,2),mean)
                init.vals$BB.prior[,,1,tt] <- min.threshold(bp,1e-5)
                init.vals$BB.prior[,,2,tt] <- 1
            }
        }
        if(is.null(init.vals$SS.prior) ||
           any(dim(init.vals$SS.prior) != c(nn,2,ee))){
            if(!is.null(init.vals$SS.prior)){
                message("SS.prior has dimension: (",
                        paste0(dim(init.vals$SS.prior),collapse=","),")")
                message("Expected dimension was: (",paste(nn,2,ee,sep=","),")")
                message("Reinitializing SS.prior...")
            }
            init.vals$SS.prior <- array(NA,c(nn,2,ee))
            for(tt in 1:ee){
                sp <- rowMeans(init.vals$SS.mat[,tmap==tt,drop=FALSE])
                init.vals$SS.prior[,1,tt] <- min.threshold(sp,1e-5)
                init.vals$SS.prior[,2,tt] <- 1
            }
        }
        if(is.null(init.vals$RR.prior) ||
           any(dim(init.vals$RR.prior) != c(nn,2,ee))){
            if(!is.null(init.vals$RR.prior)){
                message("RR.prior has dimension: (",
                        paste0(dim(init.vals$RR.prior),collapse=","),")")
                message("Expected dimension was: (",paste(nn,2,ee,sep=","),")")
                message("Reinitializing RR.prior...")
            }
            init.vals$RR.prior <- array(NA,c(nn,2,ee))
            for(tt in 1:ee){
                rp <- rowMeans(init.vals$RR.mat[,tmap==tt,drop=FALSE])
                init.vals$RR.prior[,1,tt] <- min.threshold(rp,1e-5)
                init.vals$RR.prior[,2,tt] <- 1
            }
        }
    }
    ## browser()
    init.dynsbm <- dynsbm(nn=nn,TT=TT,mmb=init.vals$mmb,
                          BB.prior=init.vals$BB.prior, BB.mat=init.vals$BB.mat,
                          SS.prior=init.vals$SS.prior, SS.mat=init.vals$SS.mat,
                          RR.prior=init.vals$RR.prior, RR.mat=init.vals$RR.mat,
                          tmap=tmap,hours.vec=hours.vec,self.ties=self.ties)

    return(init.dynsbm)
}



################################################################################
#####################  Data Visualization Functions  ###########################
################################################################################



#' Plotting Method for DynSBM Fits
#'
#' Plots the expected mean tie probabilities for a given dynsbm object.
#'
#' @param x object of class dynsbm
#' @param node.order determines ordering of nodes.  "default" uses the original
#' order of the adjacency matrix.  "membership" first orders the data using the
#' estimated block membership vector.
#' @param xaxt a character which specifies the axis type
#' @param yaxt a character which specifies the axis type
#' @param pal pallete of colors to use ordered from low to high intensity
#'
#' @details Plots a greyscale image of the adjacency matrix.  Darker colors
#' indicate larger expected intensities.
#'
#' @examples
#'
#' @seealso \code{\link{dynsbm}}
#' @export
plot.dynsbm <- function(x,marginal=TRUE,
                      node.order=c("default","membership"),
                      pal=grey((50:1)/50), ...){
    nn <- nrow(x$net)

    node.order <- match.arg(node.order)
    if(marginal){

        mfrow <- size.device(x$ee,scale=5)
        old.par <- par(mfrow=mfrow,mar=c(1,1,2,1),cex.main=2)

        pmat <- predict(x,marginal=marginal)
        if(node.order == "default"){
            ord <- 1:nn
        }else if(node.order == "membership"){
            ord <- order(x$mmb)
        }

        for(tt in 1:x$ee){
            adj.image.plot(pmat[,,tt],ord=ord,pal=pal,
                           main=paste0("Class ",tt," Average Intensity"),
                           ...)
        }

        par(old.par)
    }else{
        mfrow <- size.device(x$TT,scale=5)
        old.par <- par(mfrow=mfrow,mar=c(1,1,2,1),cex.main=2)

        wsbm.list <- create.wsbm.list(x)
        for(tt in 1:x$TT){
            plot(wsbm.list[[tt]],node.order=node.order,
                 pal=pal,main=paste0("Network ",tt," Intensity"),
                 ...)
        }
        par(old.par)
    }
}


#' @export
network.plot.dynsbm <- function(x, pal=grey((50:1)/50), node.order=NULL,
                                ...){
    net.mat <- x$net.mat
    if(is.null(net.mat))
        stop("x must have a net component.  Try calling wsbm with gen=TRUE")

    mfrow <- size.device(x$TT,scale=5)
    old.par <- par(mfrow=mfrow,mar=c(1,1,2,1),cex.main=2)

    for(tt in 1:x$TT){
        network.plot(net.mat[,,tt], pal=pal,node.order=node.order,
                     main=paste0("Network ",tt), ...)
    }

    par(old.par)
}





#' Diagnostic Plotting Method for DynSBM fits
#'
#' Plots MCMC chains for various parameters in WSBM model.
#'
#' @param x object of class "dynsbm"
#' @param t network or class number to plot diagnostics for
#' @param marginal if TRUE, plots diagnostics for class level parameters,
#' otherwise plots diagnostics for network level parameters
#' @param param which parameters diagnostics to plot.  The default value "block"
#' plots the block intensity parameters.  Passing "sender" and "receiver" plots
#' the chains for the sender and receiver effects.  The last option, "ll",
#' plots the Log-Likelihood for each step of the chain.
#' @param scale.ylim logical indicating whether each chain should have the same
#' range for the y-axis.  This can be suppressed by passing ylim.
#' @param ... other graphical parameters to pass to the plotting functions
#'
#' @details Plots MCMC chains for the specified set of parameters.
#'
#' @examples
#' data(commuter30)
#' fit.wsbm <- wsbm(commuter30,kk=3,hours=1)
#'
#' diagnostic.plot(fit.wsbm,param="block")
#' diagnostic.plot(fit.wsbm,param="sender")
#' diagnostic.plot(fit.wsbm,param="receiver")
#' diagnostic.plot(fit.wsbm,param="ll")
#'
#' @seealso \code{\link{wsbm}}
#' @export
diagnostic.plot.dynsbm.mcmc <- function(x,t,
                                        marginal=TRUE,
                                        param=c("block",
                                               "sender","receiver","ll"),
                                        ...){
    if(x$chain$clean) stop("clean.out must be FALSE to produce diagnostics")

    param <- match.arg(param)

    if(marginal){
        if(param == "block"){
            block.diagnostic.mar.dynsbm.mcmc(x,t, ...)
        }else if(param == "sender"){
            sender.diagnostic.mar.dynsbm.mcmc(x,t, ...)
        }else if(param == "receiver"){
            receiver.diagnostic.mar.dynsbm.mcmc(x,t, ...)
        }else if(param == "ll"){
            ll.diagnostic.mar.dynsbm.mcmc(x, ...)
        }

    }else{
        if(param == "block"){
            block.diagnostic.cond.dynsbm.mcmc(x,t, ...)
        }else if(param == "sender"){
            sender.diagnostic.cond.dynsbm.mcmc(x,t, ...)
        }else if(param == "receiver"){
            receiver.diagnostic.cond.dynsbm.mcmc(x,t, ...)
        }else if(param == "ll"){
            ll.diagnostic.cond.dynsbm.mcmc(x, ...)
        }
    }
}

block.diagnostic.cond.dynsbm.mcmc <- function(x,t,
                                             xlim=NULL,ylim=NULL,
                                             scale.ylim=FALSE,
                                             blocks=NULL,...){

    if(t > x$TT) stop(paste0("t is greater than the number of networks ",x$TT))

    kk <- x$kk
    if(is.null(ylim) && scale.ylim) ylim <- range(x$chain$BB.mat)
    if(is.null(blocks)) blocks <- 1:kk

    if(dev.cur() == 1){
        len <- max(length(blocks)*2,7)
        dev.new(height=len,width=len)
    }

    old.par <- par(mfrow = c(length(blocks),length(blocks)),
                   cex.main=2,cex.lab=1.5,cex.axis=1.5,
                   mar=c(5,4,2,1),oma=c(0,1,3,1))

    for(ii in blocks){
        for(jj in blocks){
            plot(x$chain$BB.mat[ii,jj,,t],ylim=ylim,xlim=xlim,type="l",
                 xlab="Iteration",ylab="",
                 main=paste0("B[",ii,",",jj,"]"),
                 ...)
        }
    }
    title(main=paste("Block Effect Chains for Network",t),outer=TRUE,cex.main=3)
    par(old.par)
}


sender.diagnostic.cond.dynsbm.mcmc <- function(x,t,
                                              xlim=NULL,ylim=NULL,
                                              scale.ylim=FALSE,
                                              nodes=NULL,...){
    if(t > x$TT) stop(paste0("t is greater than the number of networks ",x$TT))

    node.diag.plot(x$chain$SS.mat[,,t],
                   xlim=xlim,ylim=ylim,scale.ylim=scale.ylim,
                   main.prefix="S[",main.suffix="]",nodes=nodes)
    title(main=paste("Sender Effect Chains for Network",t),outer=TRUE,
          cex.main=2.5,line=-3)
}



receiver.diagnostic.cond.dynsbm.mcmc <- function(x,t,
                                              xlim=NULL,ylim=NULL,
                                              scale.ylim=FALSE,
                                              nodes=NULL,...){
    if(t > x$TT) stop(paste0("t is greater than the number of networks ",x$TT))

    node.diag.plot(x$chain$RR.mat[,,t],
                   xlim=xlim,ylim=ylim,scale.ylim=scale.ylim,
                   main.prefix="R[",main.suffix="]",nodes=nodes)
    title(main=paste("Receiver Effect Chains for Network",t),outer=TRUE,
          cex.main=2.5,line=-3)
}



ll.diagnostic.cond.dynsbm.mcmc <- function(x,...){
    ll.chain.plot(x,...)
}


#####  Marginal Diagnostic Plots


block.diagnostic.mar.dynsbm.mcmc <- function(x,t,
                                             xlim=NULL,ylim=NULL,
                                             scale.xlim=FALSE,scale.ylim=FALSE,
                                             blocks=NULL,...){

    if(t > x$ee) stop(paste0("t is greater than the number of classes ",x$ee))
    total <- x$chain$mcmc.control$total

    kk <- x$kk
    if(is.null(ylim) && scale.ylim) ylim <- range(x$chain$BB.mat[,,2,,t])
    if(is.null(xlim) && scale.xlim) ylim <- range(x$chain$BB.mat[,,1,,t])
    if(is.null(blocks)) blocks <- 1:kk

    if(dev.cur() == 1){
        len <- max(length(blocks)*2,7)
        dev.new(height=len,width=len)
    }

    old.par <- par(mfrow = c(length(blocks),length(blocks)),
                   cex.main=2,cex.lab=1.5,cex.axis=1.5,
                   mar=c(5,5,2,1),oma=c(0,1,3,1))

    for(ii in blocks){
        for(jj in blocks){
            gamma.diagnostic(x$chain$BB.prior[ii,jj,,,t],
                                ylim=ylim,xlim=xlim,
                                main=paste0("B[",ii,",",jj,"]"),
                                ...)
        }
    }
    title(main=paste("Block Prior Chains for Class",t),outer=TRUE,cex.main=3)
    par(old.par)
}

sender.diagnostic.mar.dynsbm.mcmc <- function(x,t,
                                              xlim=NULL,ylim=NULL,
                                              scale.ylim=FALSE,
                                              nodes=NULL,...){
    if(t > x$TT) stop(paste0("t is greater than the number of networks ",x$TT))

    node.gamma.plot(x$chain$SS.prior[,,,t],
                    xlim=xlim,ylim=ylim,scale.ylim=scale.ylim,
                    main.prefix="S[",main.suffix="]",nodes=nodes)
    title(main=paste("Sender Prior Chains for Class",t),outer=TRUE,
          cex.main=2.5,line=-3)
}


receiver.diagnostic.mar.dynsbm.mcmc <- function(x,t,
                                                xlim=NULL,ylim=NULL,
                                                scale.ylim=FALSE,
                                                nodes=NULL,...){
    if(t > x$TT) stop(paste0("t is greater than the number of networks ",x$TT))

    node.gamma.plot(x$chain$RR.prior[,,,t],
                    xlim=xlim,ylim=ylim,scale.ylim=scale.ylim,
                    main.prefix="R[",main.suffix="]",nodes=nodes)
    title(main=paste("Receiver Prior Chains for Class",t),outer=TRUE,
          cex.main=2.5,line=-3)
}


ll.diagnostic.mar.dynsbm.mcmc <- function(x,...){
    ll.chain.plot(x,...)
}



#####  Plotting Parameter Summaries

#' Parameter Plotting Method for DynSBM fits
#'
#' Plots boxplots of posterior samples for parameters in DYNSBM model.
#'
#' @param x object of class "dynsbm"
#' @param param which parameters diagnostics to plot.  The default value "block"
#' plots the block intensity parameters.  Passing "sender" and "receiver" plots
#' the chains for the sender and receiver effects.  The last option, "ll",
#' plots the Log-Likelihood for each step of the chain.
#' @param ... other graphical parameters to pass to the plotting functions
#'
#' @details Plots MCMC chains for the specified set of parameters.
#'
#' @examples
#' data(commuter30)
#' fit.dynsbm <- dynsbm.fit(commuter30,kk=3,hours=1)
#'
#' param.plot(fit.dynsbm,param="block")
#' param.plot(fit.dynsbm,param="sender")
#' param.plot(fit.dynsbm,param="receiver")
#' param.plot(fit.dynsbm,param="ll")
#'
#' @seealso \code{\link{dynsbm.fit}}
#' @export
param.plot.dynsbm.mcmc <- function(x,tclass,
                                   param=c("block","sender","receiver"),
                                   ...){
    if(x$chain$clean) stop("clean.out must be FALSE to produce diagnostics")

    param <- match.arg(param)
    if(param == "block"){
        block.param.plot.dynsbm.mcmc(x, tclass, ...)
    }else if(param == "sender"){
        sender.param.plot.dynsbm.mcmc(x, tclass,  ...)
    }else if(param == "receiver"){
        receiver.param.plot.dynsbm.mcmc(x, tclass, ...)
    }

}


block.param.plot.dynsbm.mcmc <- function(x, tclass, remap=NULL,
                                         truth=NULL, blocks=NULL,...){


    kk <- x$kk; x$nn; x$TT

    if(is.null(blocks)) blocks <- 1:kk
    if(is.null(remap)){
        remap <- 1:kk
    }else{
        if(length(remap) != kk){
            warning("remap has incorrect length.  Using default of 1:kk")
            remap <- 1:kk
        }
    }

    leg.txt=c("Prior Mean")
    leg.col = c("blue")

    if(!is.null(truth)){
        leg.txt <- c(leg.txt,"True Value")
        leg.col <- c(leg.col,"green")
    }

    tt.vec <- which(x$tmap == tclass)

    if(dev.cur() == 1){
        len <- max(length(blocks)*2,7)
        dev.new(height=len,width=len)
    }

    old.par <- par(mfrow = c(length(blocks),length(blocks)),
                   cex.main=1.5,cex.lab=1.2,cex.axis=1.5,
                   mar=c(5,4,2,1),oma=c(0,1,3,1))

    for(ii in blocks){
        for(jj in blocks){
            boxplot(x$chain$BB.mat[remap[ii],remap[jj],,tt.vec],
                    names=tt.vec,xlab="Network",
                    main=paste0("B[",ii,",",jj,"]"),ylab="Effect",...)

            if(!is.null(truth)){
                segments(x0=1:length(tt.vec) - .3,
                         x1=1:length(tt.vec) + .3,
                         y0=truth[ii,jj,tt.vec],col="green",lwd=3)

            }
            ab <- x$BB.prior[remap[ii],remap[jj],,tclass]
            abline(h=ab[1]/ab[2],lwd=3,col="blue")
        }
    }

    legend("topleft",legend=leg.txt,col=leg.col,lwd=3)
    title(main=paste("Block Effect Posteriors for Class",tclass),
          outer=TRUE,cex.main=3)

    par(old.par)
}



sender.param.plot.dynsbm.mcmc <- function(x, tclass, truth,
                                          nodes, block,
                                          scale.ylim=TRUE,ylim=NULL,
                                          ...){
    kk <- x$kk; nn <- x$nn; TT <- x$TT

    if(missing(nodes)){
        if(missing(block)){
            nodes <- 1:nn
        }else{
            if(block > kk) stop("block must be between 1 and ",kk)
            nodes <- which(x$mmb == block)
        }
    }else{
        if(any(nodes > nn) | any(nodes < 1)){
            stop("nodes must contain integers between 1 and ",nn)
        }
        if(!missing(block)){
            warning("The argument for block is ignored when nodes is used.")
        }
    }

    if(is.null(ylim) && scale.ylim) ylim <- range(x$chain$SS.mat[nodes,,])

    mfrow <- size.device(length(nodes))
    old.par <- par(mfrow=mfrow,cex.main=2,
                   mar=c(2,1.5,2,1),oma=c(0,1,3,1))

    tt.vec <- which(x$tmap == tclass)
    for(ii in nodes){
        boxplot(x$chain$SS.mat[ii,,tt.vec],
                names=tt.vec,xlab="",ylim=ylim,
                main=paste0("S[",ii,"]"),ylab="Effect",...)

        if(!missing(truth)){
            segments(x0=1:length(nodes) - .3,
                     x1=1:length(nodes) + .3,
                     y0=truth[nodes,tt.vec[ii]],col="green",lwd=3)
        }
        ab <- x$SS.prior[ii,,tclass]
        mm <- ab[1] / ab[2]
        abline(h=mm,col="blue",lwd=3)
    }
    title(main=paste("Sender Effect Posteriors for Class",tclass),
          outer=TRUE,cex.main=3)
    par(old.par)
}


receiver.param.plot.dynsbm.mcmc <- function(x, tclass, truth,
                                            nodes, block,
                                            scale.ylim=TRUE,ylim=NULL,
                                            ...){
    kk <- x$kk; nn <- x$nn; TT <- x$TT

    if(missing(nodes)){
        if(missing(block)){
            nodes <- 1:nn
        }else{
            if(block > kk) stop("block must be between 1 and ",kk)
            nodes <- which(x$mmb == block)
        }
    }else{
        if(any(nodes > nn) | any(nodes < 1)){
            stop("nodes must contain integers between 1 and ",nn)
        }
        if(!missing(block)){
            warning("The argument for block is ignored when nodes is used.")
        }
    }


    mfrow <- size.device(length(nodes))
    old.par <- par(mfrow=mfrow,cex.main=2,
                   mar=c(2,1.5,2,1),oma=c(0,1,3,1))

    tt.vec <- which(x$tmap == tclass)

    if(is.null(ylim) && scale.ylim) ylim <- range(x$chain$RR.mat[nodes,,])
    for(ii in nodes){
        boxplot(x$chain$RR.mat[ii,,tt.vec],
                names=tt.vec,xlab="",ylim=ylim,
                main=paste0("R[",ii,"]"),ylab="Effect",...)

        if(!missing(truth)){
            segments(x0=1:length(nodes) - .3,
                     x1=1:length(nodes) + .3,
                     y0=truth[nodes,tt.vec[ii]],col="green",lwd=3)
        }
        ab <- x$RR.prior[ii,,tclass]
        mm <- ab[1] / ab[2]
        abline(h=mm,col="blue",lwd=3)
    }
    title(main=paste("Receiver Effect Posteriors for Class",tclass),
          outer=TRUE,cex.main=3)
    par(old.par)
}



## SenderMat.plot <- function(x, tclass, truth,
##                                           horizontal=FALSE,
##                                           nodes, block, ...){

##     kk <- x$kk; nn <- x$nn; TT <- x$TT

##     if(missing(nodes)){
##         if(missing(block)){
##             nodes <- 1:nn
##         }else{
##             if(block > kk) stop("block must be between 1 and ",kk)
##             nodes <- which(x$mmb == block)
##         }
##     }else{
##         if(any(nodes > nn) | any(nodes < 1)){
##             stop("nodes must contain integers between 1 and ",nn)
##         }
##         if(!missing(block)){
##             warning("The argument for block is ignored when nodes is used.")
##         }
##     }


##     tt.vec <- which(x$tmap == tclass)
##     if(horizontal){
##         par(mfrow=c(1,length(tt.vec)))
##     }else{
##         par(mfrow=c(length(tt.vec),1))
##     }

##     for(ii in 1:length(tt.vec)){
##         boxplot(t(x$chain$SS.mat[nodes,,tt.vec[ii]]),
##                 names=nodes, horizontal=horizontal, ...)

##         if(!missing(truth)){
##             if(horizontal){
##                 segments(y0=1:length(nodes) - .3,
##                          y1=1:length(nodes) + .3,
##                          x0=truth[nodes,tt.vec[ii]],col="green",lwd=3)
##             }else{
##                 segments(x0=1:length(nodes) - .3,
##                          x1=1:length(nodes) + .3,
##                          y0=truth[nodes,tt.vec[ii]],col="green",lwd=3)
##             }
##         }
##         ab <- x$SS.prior[nodes,,tclass]
##         mm <- ab[,1] / ab[,2]
##         if(horizontal){
##             segments(y0=1:length(nodes) - .3,
##                      y1=1:length(nodes) + .3,
##                      x0=mm,col="blue",lwd=3)
##         }else{
##             segments(x0=1:length(nodes) - .3,
##                      x1=1:length(nodes) + .3,
##                      y0=mm,col="blue",lwd=3)
##         }
##     }
## }
##
## ReceiverMat.plot <- function(x, tclass, true.mat,
##                              horizontal=FALSE, nodes, block, ...){

##     kk <- dim(x$BB)[1]
##     nn <- dim(x$RR)[1]
##     TT <- dim(x$net.mat)[3]

##     if(missing(nodes)){
##         if(missing(block)){
##             nodes <- 1:nn
##         }else{
##             if(block > kk) stop("block must be between 1 and ",kk)
##             nodes <- which(x$mmb == block)
##         }
##     }else{
##         if(any(nodes > nn) | any(nodes < 1)){
##             stop("nodes must contain integers between 1 and ",nn)
##         }
##         if(!missing(block)){
##             warning("The argument for block is ignored when nodes is used.")
##         }
##     }

##     tt.vec <- which(x$tmap == tclass)
##     if(horizontal){
##         par(mfrow=c(1,length(tt.vec)))
##     }else{
##         par(mfrow=c(length(tt.vec),1))
##     }
##     for(ii in 1:length(tt.vec)){
##         boxplot(t(x$chain$RR.mat[nodes,,tt.vec[ii]]),
##                 names=nodes,horizontal=horizontal,
##                 xaxt='n', ...)

##         if(!missing(true.mat)){
##             if(horizontal){
##                 segments(y0=1:length(nodes) - .3,
##                          y1=1:length(nodes) + .3,
##                          x0=true.mat[nodes,tt.vec[ii]],col="green",lwd=3)
##             }else{
##                 segments(x0=1:length(nodes) - .3,
##                          x1=1:length(nodes) + .3,
##                          y0=true.mat[nodes,tt.vec[ii]],col="green",lwd=3)
##             }
##         }
##         if(horizontal){
##             segments(y0=1:length(nodes) - .3,
##                      y1=1:length(nodes) + .3,
##                      x0=x$RR.prior.mean[nodes,tclass],col="blue",lwd=3)
##         }else{
##             segments(x0=1:length(nodes) - .3,
##                      x1=1:length(nodes) + .3,
##                      y0=x$RR.prior.mean[nodes,tclass],col="blue",lwd=3)
##         }
##     }
## }


#' @export
block.prior.plot.dynsbm.mcmc <- function(x, remap=NULL, truth=NULL,
                                         blocks=NULL,...){

    kk <- x$kk; nn <- x$nn; ee <- x$ee

    if(is.null(blocks)){
        blocks <- 1:kk
    }

    if(missing(remap)){
        remap <- 1:kk
    }else{
        if(length(remap) != kk){
            warning("remap has incorrect length.  Using default of 1:kk")
            remap <- 1:kk
        }
    }


    if(dev.cur() == 1){
        len <- max(length(blocks)*2,7)
        dev.new(height=len,width=len)
    }

    old.par <- par(mfrow = c(length(blocks),length(blocks)),
                   cex.main=2,cex.lab=1.5,cex.axis=1.5,
                   mar=c(5,4,2,1),oma=c(0,1,3,1))

    for(ii in blocks){
        for(jj in blocks){
            boxplot(x$chain$BB.prior[remap[ii],remap[jj],1,,] /
                    x$chain$BB.prior[remap[ii],remap[jj],2,,] ,
                    main=paste0("B[",remap[ii],",",remap[jj],"]"), ...)

            if(!missing(truth)){
                segments(x0=1:ee - .3,
                         x1=1:ee + .3,
                         y0=truth[ii,jj,1,]/truth[ii,jj,2,],
                         col="green",lwd=3)
            }
        }
    }
    title("Block Prior Mean Intensity Posteriors",outer=TRUE,cex.main=3)
    par(old.par)
}



#####  Posterior Prediction for Parameters

#' Posterior Prediction of Future Parameter Values
#'
#' Generates a posterior predictive distribution for parameters in the
#' DynSBM Model
#'
#' @param x a fitted object of the class "dynsbm".  Requires clean.out = FALSE
#' in the call to dynsbm.
#' @param tclass which equivalence class to compute posterior predictive
#' distribution for.
#' @param param which parameter to compute posterior predictive distribution
#' for.
#'
#' @return Returns an array containing posterior predictive samples for the
#' given parameters.
#'
#' @examples
#'
#' @seealso \code{\link{dynsbm}} \code{\link{predict.dynsbm}}
#' @export
param.post.predict.dynsbm.mcmc <- function(x,tclass,
                                           param=c("block",
                                                   "sender","receiver"),
                                           ...){
    if(x$chain$clean) stop("clean.out must be FALSE to produce diagnostics")

    param <- match.arg(param)

    if(param == "block"){
        post.predict.block.dynsbm.mcmc(x,tclass,...)
    }else if(param == "sender"){
        post.predict.sender.dynsbm.mcmc(x,tclass,...)
    }else if(param == "receiver"){
        post.predict.receiver.dynsbm.mcmc(x,tclass,...)
    }

}


post.predict.block.dynsbm.mcmc <- function(x,tclass,blocks,
                                           remap,iter.draws=1,...){

    kk <- x$kk
    total <- x$chain$mcmc.control$total

    if(missing(blocks)) blocks <- 1:kk
    if(missing(remap)) remap <-1:kk

    samp <- array(NA,c(length(blocks), length(blocks), total*iter.draws))
    BB.prior <- x$chain$BB.prior[,,,,tclass]

    for(ii in 1:length(ss.ind)){
        ss <- ss.ind[ii]
        for(jj in 1:length(rr.ind)){
            rr <- rr.ind[jj]
            mat <- t(BB.prior[remap[ss],remap[rr],,])
            samp[ii,jj,] <- gamma.post.pred(mat,iter.draws=iter.draws)
        }
    }
    return(samp)
}



post.predict.sender.dynsbm.mcmc <- function(x,tclass,blocks,
                                            iter.draws=1,...){

    nn <- x$nn
    total <- x$chain$mcmc.control$total

    if(missing(nodes)) nodes <- 1:nn

    samp <- array(NA,c(length(nodes), total*iter.draws))
    prior <- x$chain$SS.prior[,,,tclass]

    for(ii in 1:length(nodes)){
        mat <- t(prior[nodes[ii],,])
        samp[ii,] <- gamma.post.pred(mat,iter.draws=iter.draws)
    }

    return(samp)
}


post.predict.receiver.dynsbm.mcmc <- function(x,tclass,blocks,
                                              iter.draws=1,...){

    nn <- x$nn
    total <- x$chain$mcmc.control$total

    if(missing(nodes)) nodes <- 1:nn

    samp <- array(NA,c(length(nodes), total*iter.draws))
    prior <- x$chain$RR.prior[,,,tclass]

    for(ii in 1:length(nodes)){
        mat <- t(prior[nodes[ii],,])
        samp[ii,] <- gamma.post.pred(mat,iter.draws=iter.draws)
    }

    return(samp)
}


#' @export
BlockMat.Post.Density.Plot <- function(x,tclass,
                                       ss.ind=1,rr.ind=1,
                                       remap,col="steelblue2",
                                       trans, add=FALSE,
                                       include.pred=TRUE,
                                       iter.draws=1, pred.lwd=2,
                                       ...){

    if(length(ss.ind) > 1){
        stop("ss.ind must be a single index")
    }
    if(length(rr.ind) > 1){
        stop("rr.ind must be a single index")
    }
    if(missing(remap)){
        kk <- x$kk
        remap=1:kk
    }

    tt.vec <- which(x$tmap == tclass)

    if(missing(trans)){
        trans <- 0.5 / sqrt(length(tt.vec))
    }

    density.mat.plot(x$chain$BB.mat[remap[ss.ind],remap[rr.ind],,tt.vec],
                     col=col,trans=trans,add=add, ...)
    if(include.pred){
        samp <- post.pred.block(x=x,tclass=tclass,
                                ss.ind=ss.ind, rr.ind=rr.ind, remap=remap,
                                iter.draws=iter.draws)
        abline(v=quantile(samp,c(0.025,0.975)),col=col,lty=2,lwd=pred.lwd)
        lines(density(samp,from=0,n=1e4),col=col,lwd=pred.lwd)
    }

}


#' @export
Sender.Post.Density.Plot <- function(x,tclass,node,
                                     col="steelblue2",trans, add=FALSE,
                                     include.pred=FALSE, iter.draws=100,
                                     ...){

    if(length(node) > 1){
        stop("node must be a single index")
    }

    tt.vec <- which(x$tmap == tclass)

    if(missing(trans)){
        trans <- 0.5 / sqrt(length(tt.vec))
    }

    density.mat.plot(x$chain$SS.mat[node,,tt.vec],
                     col=col,trans=trans,add=add, ...)
    if(include.pred){
        samp <- post.pred.sender(x=x,tclass=tclass,
                                nodes=node,
                                iter.draws=iter.draws)
        abline(v=quantile(samp,c(0.025,0.975)),col=col,lty=2)
        lines(density(samp,from=0,n=1e4),col=col,lwd=2)
    }

}




#' @export
Receiver.Post.Density.Plot <- function(x,tclass,node,
                                       col="steelblue2",trans, add=FALSE,
                                       include.pred=FALSE, iter.draws=100,
                                       ...){

    if(length(node) > 1){
        stop("node must be a single index")
    }

    tt.vec <- which(x$tmap == tclass)

    if(missing(trans)){
        trans <- 0.5 / sqrt(length(tt.vec))
    }

    density.mat.plot(x$chain$RR.mat[node,,tt.vec],
                     col=col,trans=trans,add=add, ...)
    if(include.pred){
        samp <- post.pred.receiver(x=x,tclass=tclass,
                                nodes=node,
                                iter.draws=iter.draws)
        abline(v=quantile(samp,c(0.025,0.975)),col=col,lty=2)
        lines(density(samp,from=0,n=1e4),col=col,lwd=2)
    }

}



#####  Wrapper Functions to perform the above functions for specific classes
#####  and indices within the block probability matrix.
#' @export
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


#' @export
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



#' @export
get.iter.dynsbm.mcmc <- function(object,iter){
    if(is.null(object$chain)){
        stop(paste("Object does not contain MCMC chain.",
                   "Call dynsbm with clean.out=FALSE"))
    }
    chain <- object$chain

    BB.prior <- array(NA,dim(chain$BB.prior)[-4])
    SS.prior <- RR.prior <- array(NA,dim(chain$SS.prior)[-3])
    ee <- dim(BB.prior)[4]
    for(tt in 1:ee){
        BB.prior[,,,tt] <- chain$BB.prior[,,,iter,tt]
        SS.prior[,,tt] <- chain$SS.prior[,,iter,tt]
        RR.prior[,,tt] <- chain$RR.prior[,,iter,tt]
    }

    dynsbm.obj <- dynsbm(nn=object$nn,TT=object$TT,mmb=chain$mmb[,iter],
                         BB.prior=BB.prior, BB.mat=chain$BB.mat[,,iter,],
                         SS.prior=SS.prior, SS.mat=chain$SS.mat[,iter,],
                         RR.prior=RR.prior, RR.mat=chain$RR.mat[,iter,],
                         tmap=object$tmap,hours.vec=object$hours.vec,
                         self.ties=object$self.ties)
    dynsbm.obj$net.mat <- object$net.mat

    return(dynsbm.obj)

    ## if(dim(chain$mmb)[2] < iter){
    ##     stop(paste0("Chain has less than ",iter," iterations"))
    ## }


    ## iter.list <- list(mmb=chain$mmb[,iter],
    ##                   BB.mat=chain$BB.mat[,,iter,],
    ##                   SS.mat=chain$SS.mat[,iter,],
    ##                   RR.mat=chain$RR.mat[,iter,],
    ##                   BB.prior=BB.prior,
    ##                   SS.prior=SS.prior,RR.prior=RR.prior)

    ## return(iter.list)
}
