#####  dynsbm.R
#####  Functions Implementing SBM and MMSBM Fitting in C:
#####    sbm - Fits SBM Model
#####    mmsbm.C - Fits MMSBM Model
#####    mmsbm - Fits MMSBM choosing intelligent starting values
#####
##############################################################

#####  Function to generate data from WSBM model

####  Create some good examples for data generation


#' @include generics.R helper.R wsbm.R
NULL

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
#' @param BB.mat block effect parameters for each network.  Sampled from
#'   BB.prior if not provided.
#' @param SS.mat sender effect parameters for each network.  Sampled from
#'   SS.prior if not provided.
#' @param RR.mat receiver effect parameters for each network.  Sampled from
#'   RR.prior if not provided.
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
#'               BB.prior=BB.prior,SS.prior=SS.prior,RR.prior=RR.prior,
#'               mmb=mmb,self.ties=TRUE,normalize=TRUE,gen=TRUE)
#'
#' plot(dat)
#' plot(dat,marginal=FALSE)
#' network.plot(dat)
#'
#' @seealso \code{\link{dynsbm.fit}}
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


#' @describeIn dynsbm Test whether an object is of class dynsbm
#' @param object object that could be of class "dynsbm"
#' @export
is.dynsbm <- function(object) inherits(object,"dynsbm")


#' Generate Set of Dynamic Stochastic Block Model Networks
#'
#' Generates networks from the weighted stochastic block model
#' with or without degree effects.
#'
#' @param object object of class "dynsbm"
#' @param marginal if TRUE, parameters are sampled from prior distributions
#' before simulating a new network
#' @param ... additional parameters
#'
#' @return Returns an n x n x T array representing a collection of T
#' networks with n nodes corresponding to the dynsbm object.
#'
#' @examples
#' kk <- 3
#' nn <- kk^2 * 10
#' TT <- 10; ee <- 2
#' tmap <- rep(1:ee,TT/ee)
#' hours.vec <- c(12,12)
#' mmb <-  rep(1:kk,each=nn/kk)
#'
#' ###  Setting up Reasonable Prior Distributions
#' BB.prior <- array(1,c(kk,kk,2,ee))
#' BB.prior[cbind(1:3,c(2,3,1),1,1)] <- 2
#' BB.prior[cbind(1:3,c(3,1,2),1,2)] <- 2
#'
#' SS.prior <- array(rgamma(nn*2*ee,10,10),c(nn,2,ee))
#' RR.prior <- array(rgamma(nn*2*ee,10,10),c(nn,2,ee))
#'
#' ###  Setting up Object with Dynsbm Parameters
#' dat <- dynsbm(nn=nn,TT=TT,tmap=tmap,hours.vec=hours.vec,
#'               BB.prior=BB.prior,SS.prior=SS.prior,RR.prior=RR.prior,
#'               mmb=mmb,self.ties=TRUE,normalize=TRUE)
#'
#' ###  Sampling New Edges Only
#' net.mat.1 <- net.gen(dat)
#' network.plot(net.mat.1)
#'
#' net.mat.2 <- net.gen(dat)
#' network.plot(net.mat.2)
#'
#' net.mat.3 <- net.gen(dat)
#' network.plot(net.mat.3)
#'
#' ###  Resampling Parameters from Priors and Edges
#' net.mat.marg.1 <- net.gen(dat,marginal=TRUE)
#' network.plot(net.mat.marg.1)
#'
#' net.mat.marg.2 <- net.gen(dat,marginal=TRUE)
#' network.plot(net.mat.marg.2)
#'
#' net.mat.marg.3 <- net.gen(dat,marginal=TRUE)
#' network.plot(net.mat.marg.3)
#'
#' @seealso \code{\link{dynsbm}}, \code{\link{dynsbm.fit}}
#' @export
net.gen.dynsbm <- function(object,marginal=FALSE, ...){

    nn <- object$nn; TT <- object$TT; tmap <- object$tmap; kk <- object$kk

    if(marginal){
        BB.mat <- array(NA,c(kk,kk,TT))
        for(tt in 1:TT){
            BB.mat[,,tt] <- rgamma(kk^2,object$BB.prior[,,1,tmap[tt]],
                                   object$BB.prior[,,2,tmap[tt]])
        }

        SS.mat <- array(NA,c(nn,TT))
        for(tt in 1:TT){
            SS.mat[,tt] <- rgamma(nn,object$SS.prior[,1,tmap[tt]],
                                  object$SS.prior[,2,tmap[tt]])
        }

        RR.mat <- array(NA,c(nn,TT))
        for(tt in 1:TT){
            RR.mat[,tt] <- rgamma(nn,object$RR.prior[,1,tmap[tt]],
                                  object$RR.prior[,2,tmap[tt]])
        }

        object$BB.mat <- BB.mat
        object$SS.mat <- SS.mat
        object$RR.mat <- RR.mat
    }

    wsbm.list <- create.wsbm.list(object)

    ##  Creating containers for network and parameters
    net.mat <- array(NA,c(object$nn,object$nn,object$TT))

    for(tt in 1:object$TT){
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
#' @param update.mmb if FALSE, the block membership vector is not updated
#'   during the MCMC routine.  This can be used if the block memberships are
#'   known in advance.
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
#' @seealso \code{\link{dynsbm}}, \code{\link{dynsbm.priors}}
#' @export
dynsbm.fit <- function(net.mat, kk=3, tmap, hours.vec, self.ties=TRUE,
                   priors=dynsbm.priors(),
                   init.control=list(spectral.start=FALSE,
                                     multistart=0,
                                     multistart.total=10),
                   mcmc.control=list(total=1000,burn.in=1000,thin=10,
                                     extend.alpha=.001,extend.max=10,
                                     label.switch.mode="adhoc",
                                     label.switch.max=200,
                                     multi.impute=FALSE),
                   update.mmb=TRUE,clean.out=FALSE, verbose=1){


    start.time <- Sys.time()

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

    ## Checking for Block Membership Prior
    if(is.null(priors$eta)) priors$eta <- rep(1/kk,kk)

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

    end.time <- Sys.time()
    dynsbm.fit.obj$chain$time <- (end.time - start.time)

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
#' @param eta prior on block memberships.
#' @param BB.hyperprior hyperprior parameters for conjugate Gamma
#' distribution on block effect priors
#' @param SS.hyperprior hyperprior parameters for conjugate Gamma
#' distribution on sender effect priors
#' @param RR.hyperprior hyperprior parameters for conjugate Gamma
#' distribution on receiver effect priors
#'
#' @return Returns a list containing the parameters with useful defaults.  eta
#' is the only required parameter.
#'
#' @seealso \code{\link{dynsbm}}
#' @export
dynsbm.priors <- function(eta=NULL,
                          BB.hyperprior=c(.9999,.0011,.001,.001),
                          SS.hyperprior=c(.999,.011,.01,.01),
                          RR.hyperprior=c(.999,.011,.01,.01)){

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
#' @param marginal if TRUE, prediction is made using the marginal model.
#' Otherwise, the model conditional on parameter estimates is used
#' @param reps Number of samples from the marginal model to use when marginal
#' = TRUE
#' @param ... additional parameters
#'
#' @return Returns an array of predicted intensities.  These intensities are
#' the predicted values given the posterior mean parameter estimates.
#'
#' @seealso \code{\link{dynsbm.fit}}
#' @export
predict.dynsbm <- function(object, marginal=TRUE, reps = 1000, ...){

    if(marginal){
        pmm <- array(NA,c(object$nn,object$nn,object$ee,reps))
        object$TT <- object$ee
        object$tmap <- 1:object$ee
        for(ii in 1:reps){
            pmm[,,,ii] <- net.gen(object,marginal=marginal)
        }
        pmat <- apply(pmm,c(1,2,3),mean)
        return(pmat)
    }else{
        wsbm.list <-  create.wsbm.list(object)
        TT <- object$TT
        pmat <- array(NA,dim(object$net.mat))

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
#' @param object a fitted object of the class "dynsbm".  Requires
#' clean.out = FALSE in the call to dynsbm.
#' @param marginal if TRUE, values block, sender, and receiver effects are
#' marginalized over using the estimated prior parameters.
#' @param ... additional parameters
#'
#' @return Returns a 3 dimensional array of dimension n x n x total, where n
#' is the number of nodes in the network and total is the number of draws in
#' the MCMC chain generated using wsbm.
#'
#' @seealso \code{\link{dynsbm}} \code{\link{predict.dynsbm}}
#' @export
post.predict.dynsbm.mcmc <- function(object, marginal=TRUE, ...){
    if(object$chain$clean) stop("Use clean.out = TRUE to keep MCMC chains")
    nn <- object$nn; TT <- object$TT; ee <- object$ee
    total <- object$chain$mcmc.control$total


    if(marginal) pmat.post <- array(NA,c(nn,nn,ee,total))
    else pmat.post <- array(NA,c(nn,nn,TT,total))

    for(ii in 1:total){
        tmp.obj <- get.iter(object,ii)
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


summary.dynsbm.chain <- function(object,...){
###    browser()
    TT <- dim(object$net.mat)[3]
    total <- dim(object$BB.mat)[3];
    nn <- object$nn; TT <- object$TT; ee <- object$ee; kk <- object$kk

###  Summary of Parameters
    BB.hat <- apply(object$BB.mat,c(1,2,4),mean)
    SS.hat <- apply(object$SS.mat,c(1,3),mean)
    RR.hat <- apply(object$RR.mat,c(1,3),mean)

    PI.mean <- t(apply(object$mmb,1,tabulate,nbins=kk) / total)
    mmb <- apply(PI.mean,1,which.max)

###  Summary of Priors
    BB.prior.hat <- apply(object$BB.prior,c(1,2,3,5),mean)
    SS.prior.hat <- apply(object$SS.prior,c(1,2,4),mean)
    RR.prior.hat <- apply(object$RR.prior,c(1,2,4),mean)


    dynsbm.obj <- dynsbm(nn=nn,TT=TT,mmb=mmb,
                         BB.prior=BB.prior.hat,BB.mat=BB.hat,
                         SS.prior=SS.prior.hat,SS.mat=SS.hat,
                         RR.prior=RR.prior.hat,RR.mat=RR.hat,
                         tmap=object$tmap,hours.vec=object$hours.vec,
                         self.ties=object$self.ties)

    return(dynsbm.obj)

    ##  Calculating Posterior Mean/Variance
    BB.prior.mean <- apply(object$BB.prior[,,1,,,drop=FALSE]/
                           object$BB.prior[,,2,,,drop=FALSE],c(1,2,5),mean)
    BB.prior.var <- apply(object$BB.prior[,,1,,,drop=FALSE]/
                          (object$BB.prior[,,2,,,drop=FALSE]^2),c(1,2,5),mean)

    SS.prior.mean <- apply(object$SS.prior[,1,,,drop=FALSE]/
                           object$SS.prior[,2,,,drop=FALSE],c(1,4),mean)
    SS.prior.var <- apply(object$SS.prior[,1,,,drop=FALSE]/
                          (object$SS.prior[,2,,,drop=FALSE]^2),c(1,4),mean)

    RR.prior.mean <- apply(object$RR.prior[,1,,,drop=FALSE]/
                           object$RR.prior[,2,,,drop=FALSE],c(1,4),mean)
    RR.prior.var <- apply(object$RR.prior[,1,,,drop=FALSE]/
                          (object$RR.prior[,2,,,drop=FALSE]^2),c(1,4),mean)
}


#' Log-Likelihood method for DynSBM Fits
#'
#' Returns log-likelihood for posterior mean estimates.
#'
#' @param object an object of class "dynsbm"
#' @param ... additional arguments to logLik
#'
#' @return returns the log-likelihood of the data given the posterior mean
#' parameter estimates.
#'
#' @seealso \code{\link{dynsbm}}
#' @export
logLik.dynsbm <- function(object, ...){
    return(sum(sapply(create.wsbm.list(object),logLik)))
}



#' AIC method for DynSBM Fits
#'
#' Returns Akaike Information Criterion for DynSBM Fit
#'
#' @param object an object of class "dynsbm"
#' @param ... additional arguments to AIC
#' @param k numeric, the penalty per parameter to be used; the default k = 2
#' is the classical AIC
#'
#' @return returns the Akaike Information Criterion using the posterior mean
#' parameters as approximations to the MLE.  The number of parameters is taken
#' to be T*(2*(n-k) + k*k) + n.
#'
#' @seealso \code{\link{dynsbm}}
#' @export
AIC.dynsbm <- function(object, ..., k=2){
    TT <- object$TT
    nn <- object$nn
    kk <- object$kk
    edf <- TT*(2*(nn-kk) + kk^2) +nn
    return(-2*logLik(object) + k * edf)
}



#' BIC method for DynSBM Fits
#'
#' Returns Bayesian Information Criterion for DynBM Fit
#'
#' @param object an object of class "dynsbm"
#' @param ... additional arguments to BIC
#'
#' @return returns the Bayesian Information Criterion using the posterior mean
#' parameters as approximations to the MLE.  The number of parameters is taken
#' to be T*(2*(n-k) + k*k) + n and the sample size is taken to be T*n*n.
#'
#' @seealso \code{\link{dynsbm}}
#' @import stats
#' @export
BIC.dynsbm <- function(object, ...){
    TT <- object$TT
    nn <- object$nn
    kk <- object$kk
    edf <- TT*(2*(nn-kk) + kk^2) +nn
    return(-2*logLik(object) + (2*log(nn) + log(TT)) * edf)
}



#' DIC method for DynSBM Fits
#'
#' Returns Deviance Information Criterion for DynSBM Fit
#'
#' @param object an object of class "dynsbm"
#' @param ... additional parameters
#'
#' @return returns the Deviance Information Criterion.
#'
#' @seealso \code{\link{dynsbm}}
#' @export
DIC.dynsbm <- function(object, ...){
    if(!is.null(object$chain$logLik)){
        DIC.df <- (2*logLik(object) - 2*mean(object$chain$logLik))
        return(-2*logLik(object) + 2 * DIC.df)
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
                      format.vector(x$tmap,digits=digits,
                                    max.width=max.width,...),
                      "\n")

    node.str <- paste("Nodes:",dim(x$net.mat)[1])
    block.str <- paste("Blocks:",dim(x$BB)[1])
    mmb.str <- paste0("Block Membership:\n",
                      format.vector(x$mmb,digits=digits,
                                    max.width=max.width,...))

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

#' Convert DynSBM Object to List of WSBM Objects
#'
#' @param object object of class "dynsbm"
#'
#' @return Returns a list of length TT where each object in the list is of
#'   class "wsbm".
#' @export
create.wsbm.list <- function(object){
    TT <- object$TT
    wsbm.list <- NULL
    for(tt in 1:TT){
        obj <- with(object,wsbm(nn=nn,mmb=mmb,BB=BB.mat[,,tt],
                           SS=SS.mat[,tt],RR=RR.mat[,tt],
                           self.ties=self.ties,hours=hours.vec[tmap[tt]]))
        if(!is.null(object$net.mat)){
            obj$net <- object$net.mat[,,tt]
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
#' @param marginal if TRUE average tie probabilities for the equivalence
#'   classes are plotted.  If marginal is FALSE, the average tie probability
#'   for each individual network is plotted instead.
#' @param node.order determines ordering of nodes.  "default" uses the original
#' order of the adjacency matrix.  "membership" first orders the data using the
#' estimated block membership vector.
#' @param pal color palette used to indicate tie intensity
#' @param ... additional graphical parameters
#'
#' @details Plots a greyscale image of the mean tie intensity matrix.  The
#'   individual network parameters are marginalized over each equivalence
#'   class when marginal = TRUE.
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

#' @describeIn network.plot Plots heatmap for each network in net.mat
#' @export
network.plot.dynsbm <- function(x, pal=grey((50:1)/50), node.order=NULL,
                                ...){
    net.mat <- x$net.mat
    if(is.null(net.mat))
        stop("x must have a net component.  Call wsbm with gen=TRUE")

    mfrow <- size.device(x$TT,scale=5)
    old.par <- par(mfrow=mfrow,mar=c(1,1,2,1),cex.main=2)

    for(tt in 1:x$TT){
        network.plot(net.mat[,,tt], pal=pal,node.order=node.order,
                     main=paste0("Network ",tt), ...)
    }

    par(old.par)
}





#' @describeIn diagnostic.plot Plots MCMC chains for various
#'   parameters in DynSBM model.
#'
#' @param marginal if TRUE, plots diagnostics for class level parameters,
#' otherwise plots diagnostics for network level parameters
#' @param t if marginal = TRUE, t is the class number for which to plot
#'   diagnostics.  If marginal = FALSE, t is the network number.
#'
#' @examples
#' ### Dynamic SBM Example
#' data(dc60)
#' fit.dc <- dynsbm.fit(dc60$net.mat,kk=3,
#'                      tmap=dc60$tmap,hours.vec=dc60$hours.vec)
#'
#' diagnostic.plot(fit.dc,t=1,param="block")
#' diagnostic.plot(fit.dc,t=2,param="block")
#'
#' diagnostic.plot(fit.dc,t=1,param="sender")
#' diagnostic.plot(fit.dc,t=1,param="receiver")
#' diagnostic.plot(fit.dc,t=1,param="ll")
#'
#' diagnostic.plot(fit.dc,t=1,param="block",marginal=FALSE)
#' @export
diagnostic.plot.dynsbm.mcmc <- function(object,marginal=TRUE,t,
                                        param=c("block",
                                                "sender","receiver","ll"),
                                        ...){
    if(object$chain$clean) stop("clean.out must be FALSE to use diagnostics")

    param <- match.arg(param)

    if(marginal){
        if(param == "block"){
            block.diagnostic.mar.dynsbm.mcmc(object,t, ...)
        }else if(param == "sender"){
            sender.diagnostic.mar.dynsbm.mcmc(object,t, ...)
        }else if(param == "receiver"){
            receiver.diagnostic.mar.dynsbm.mcmc(object,t, ...)
        }else if(param == "ll"){
            ll.diagnostic.mar.dynsbm.mcmc(object, ...)
        }

    }else{
        if(param == "block"){
            block.diagnostic.cond.dynsbm.mcmc(object,t, ...)
        }else if(param == "sender"){
            sender.diagnostic.cond.dynsbm.mcmc(object,t, ...)
        }else if(param == "receiver"){
            receiver.diagnostic.cond.dynsbm.mcmc(object,t, ...)
        }else if(param == "ll"){
            ll.diagnostic.cond.dynsbm.mcmc(object, ...)
        }
    }
}

block.diagnostic.cond.dynsbm.mcmc <- function(object,t,
                                             xlim=NULL,ylim=NULL,
                                             scale.ylim=FALSE,
                                             blocks=NULL,...){

    if(t > object$TT) stop(paste0("t is greater than the number of networks ",
                                  object$TT))

    kk <- object$kk
    if(is.null(ylim) && scale.ylim) ylim <- range(object$chain$BB.mat)
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
            plot(object$chain$BB.mat[ii,jj,,t],ylim=ylim,xlim=xlim,type="l",
                 xlab="Iteration",ylab="",
                 main=paste0("B[",ii,",",jj,"]"),
                 ...)
        }
    }
    title(main=paste("Block Effect Chains for Network",t),outer=TRUE,cex.main=3)
    par(old.par)
}


sender.diagnostic.cond.dynsbm.mcmc <- function(object,t,
                                              xlim=NULL,ylim=NULL,
                                              scale.ylim=FALSE,
                                              nodes=NULL,...){
    if(t > object$TT) stop(paste0("t is greater than the number of networks ",
                                  object$TT))

    node.diag.plot(object$chain$SS.mat[,,t],
                   xlim=xlim,ylim=ylim,scale.ylim=scale.ylim,
                   main.prefix="S[",main.suffix="]",nodes=nodes)
    title(main=paste("Sender Effect Chains for Network",t),outer=TRUE,
          cex.main=2.5,line=0)
}



receiver.diagnostic.cond.dynsbm.mcmc <- function(object,t,
                                              xlim=NULL,ylim=NULL,
                                              scale.ylim=FALSE,
                                              nodes=NULL,...){
    if(t > object$TT) stop(paste0("t is greater than the number of networks ",
                                  object$TT))

    node.diag.plot(object$chain$RR.mat[,,t],
                   xlim=xlim,ylim=ylim,scale.ylim=scale.ylim,
                   main.prefix="R[",main.suffix="]",nodes=nodes)
    title(main=paste("Receiver Effect Chains for Network",t),outer=TRUE,
          cex.main=2.5,line=0)
}



ll.diagnostic.cond.dynsbm.mcmc <- function(object,...){
    ll.chain.plot(object,...)
}


#####  Marginal Diagnostic Plots


block.diagnostic.mar.dynsbm.mcmc <- function(object,t,
                                             xlim=NULL,ylim=NULL,
                                             scale.xlim=FALSE,scale.ylim=FALSE,
                                             blocks=NULL,...){

    if(t > object$ee) stop(paste0("t is greater than the number of classes ",
                                  object$ee))
    total <- object$chain$mcmc.control$total

    kk <- object$kk
    if(is.null(ylim) && scale.ylim) ylim <- range(object$chain$BB.mat[,,2,,t])
    if(is.null(xlim) && scale.xlim) ylim <- range(object$chain$BB.mat[,,1,,t])
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
            gamma.diagnostic(object$chain$BB.prior[ii,jj,,,t],
                                ylim=ylim,xlim=xlim,
                                main=paste0("B[",ii,",",jj,"]"),
                                ...)
        }
    }
    title(main=paste("Block Prior Chains for Class",t),outer=TRUE,cex.main=3)
    par(old.par)
}

sender.diagnostic.mar.dynsbm.mcmc <- function(object,t,
                                              xlim=NULL,ylim=NULL,
                                              scale.ylim=FALSE,
                                              nodes=NULL,...){
    if(t > object$TT) stop(paste0("t is greater than the number of networks ",
                                  object$TT))

    node.gamma.plot(object$chain$SS.prior[,,,t],
                    xlim=xlim,ylim=ylim,scale.ylim=scale.ylim,
                    main.prefix="S[",main.suffix="]",nodes=nodes)
    title(main=paste("Sender Prior Chains for Class",t),outer=TRUE,
          cex.main=2.5,line=0)
}


receiver.diagnostic.mar.dynsbm.mcmc <- function(object,t,
                                                xlim=NULL,ylim=NULL,
                                                scale.ylim=FALSE,
                                                nodes=NULL,...){
    if(t > object$TT) stop(paste0("t is greater than the number of networks ",
                                  object$TT))

    node.gamma.plot(object$chain$RR.prior[,,,t],
                    xlim=xlim,ylim=ylim,scale.ylim=scale.ylim,
                    main.prefix="R[",main.suffix="]",nodes=nodes)
    title(main=paste("Receiver Prior Chains for Class",t),outer=TRUE,
          cex.main=2.5,line=0)
}


ll.diagnostic.mar.dynsbm.mcmc <- function(object,...){
    ll.chain.plot(object,...)
}



#####  Plotting Parameter Summaries

# Parameter Plotting Method for DynSBM fits
#
# Plots boxplots of posterior samples for parameters in DYNSBM model.
#' @describeIn param.plot Displays boxplots of the posterior samples for
#'   parameters within one of the equivalence classes specified by tclass.
#' @param tclass class for which to plot parameters.
#'
#' @export
param.plot.dynsbm.mcmc <- function(object,tclass,
                                   param=c("block","sender","receiver"),
                                   ...){
    if(object$chain$clean) stop("clean.out must be FALSE to use diagnostics")

    param <- match.arg(param)
    if(param == "block"){
        block.param.plot.dynsbm.mcmc(object, tclass, ...)
    }else if(param == "sender"){
        sender.param.plot.dynsbm.mcmc(object, tclass,  ...)
    }else if(param == "receiver"){
        receiver.param.plot.dynsbm.mcmc(object, tclass, ...)
    }

}


block.param.plot.dynsbm.mcmc <- function(object, tclass, remap=NULL,
                                         truth=NULL, blocks=NULL,...){


    kk <- object$kk; object$nn; object$TT

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

    tt.vec <- which(object$tmap == tclass)

    if(dev.cur() == 1){
        len <- max(length(blocks)*2,7)
        dev.new(height=len,width=len)
    }

    old.par <- par(mfrow = c(length(blocks),length(blocks)),
                   cex.main=1.5,cex.lab=1.2,cex.axis=1.5,
                   mar=c(5,4,2,1),oma=c(0,1,3,1))

    for(ii in blocks){
        for(jj in blocks){
            boxplot(object$chain$BB.mat[remap[ii],remap[jj],,tt.vec],
                    names=tt.vec,xlab="Network",
                    main=paste0("B[",ii,",",jj,"]"),ylab="Effect",...)

            if(!is.null(truth)){
                segments(x0=1:length(tt.vec) - .3,
                         x1=1:length(tt.vec) + .3,
                         y0=truth[ii,jj,tt.vec],col="green",lwd=3)

            }
            ab <- object$BB.prior[remap[ii],remap[jj],,tclass]
            abline(h=ab[1]/ab[2],lwd=3,col="blue")
        }
    }

    legend("topleft",legend=leg.txt,col=leg.col,lwd=3)
    title(main=paste("Block Effect Posteriors for Class",tclass),
          outer=TRUE,cex.main=3)

    par(old.par)
}



sender.param.plot.dynsbm.mcmc <- function(object, tclass, truth,
                                          nodes, block,
                                          scale.ylim=TRUE,ylim=NULL,
                                          ...){
    kk <- object$kk; nn <- object$nn; TT <- object$TT

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

    if(is.null(ylim) && scale.ylim) ylim <- range(object$chain$SS.mat[nodes,,])

    mfrow <- size.device(length(nodes))
    old.par <- par(mfrow=mfrow,cex.main=2,
                   mar=c(2,1.5,2,1),oma=c(0,1,3,1))

    tt.vec <- which(object$tmap == tclass)
    for(ii in nodes){
        boxplot(object$chain$SS.mat[ii,,tt.vec],
                names=tt.vec,xlab="",ylim=ylim,
                main=paste0("S[",ii,"]"),ylab="Effect",...)

        if(!missing(truth)){
            segments(x0=1:length(nodes) - .3,
                     x1=1:length(nodes) + .3,
                     y0=truth[nodes,tt.vec[ii]],col="green",lwd=3)
        }
        ab <- object$SS.prior[ii,,tclass]
        mm <- ab[1] / ab[2]
        abline(h=mm,col="blue",lwd=3)
    }
    title(main=paste("Sender Effect Posteriors for Class",tclass),
          outer=TRUE,cex.main=3)
    par(old.par)
}


receiver.param.plot.dynsbm.mcmc <- function(object, tclass, truth,
                                            nodes, block,
                                            scale.ylim=TRUE,ylim=NULL,
                                            ...){
    kk <- object$kk; nn <- object$nn; TT <- object$TT

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


    mfrow <- size.device(length(nodes))
    old.par <- par(mfrow=mfrow,cex.main=2,
                   mar=c(2,1.5,2,1),oma=c(0,1,3,1))

    tt.vec <- which(object$tmap == tclass)

    if(is.null(ylim) && scale.ylim) ylim <- range(object$chain$RR.mat[nodes,,])
    for(ii in nodes){
        boxplot(object$chain$RR.mat[ii,,tt.vec],
                names=tt.vec,xlab="",ylim=ylim,
                main=paste0("R[",ii,"]"),ylab="Effect",...)

        if(!missing(truth)){
            segments(x0=1:length(nodes) - .3,
                     x1=1:length(nodes) + .3,
                     y0=truth[nodes,tt.vec[ii]],col="green",lwd=3)
        }
        ab <- object$RR.prior[ii,,tclass]
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


#' Block Prior Plotting Method for DynSBM Class
#'
#' Method for plotting posterior distribution of block prior means across
#'   DynSBM classes
#'
#' @param object object of class "dynsbm"
#' @param remap relabeling of blocks that is a permuation of 1:kk
#' @param truth true values for comparison if known
#' @param blocks subset of the blocks 1:kk for which to plot posteriors
#' @param ... additional parameters to pass to plotting functions
#'
#' @export
block.prior.plot.dynsbm.mcmc <- function(object, remap=NULL, truth=NULL,
                                         blocks=NULL,...){

    kk <- object$kk; nn <- object$nn; ee <- object$ee

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
            boxplot(object$chain$BB.prior[remap[ii],remap[jj],1,,] /
                    object$chain$BB.prior[remap[ii],remap[jj],2,,] ,
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
#' @param object a fitted object of the class "dynsbm".  Requires
#' clean.out = FALSE
#' in the call to dynsbm.
#' @param tclass which equivalence class to compute posterior predictive
#' distribution for.
#' @param param which parameter to compute posterior predictive distribution
#' for.
#' @param ... additional parameters
#'
#' @return Returns an array containing posterior predictive samples for the
#' given parameters.
#'
#'
#' @seealso \code{\link{dynsbm}} \code{\link{predict.dynsbm}}
#' @export
param.post.predict.dynsbm.mcmc <- function(object,tclass,
                                           param=c("block",
                                                   "sender","receiver"),
                                           ...){
    if(object$chain$clean) stop("clean.out must be FALSE to use diagnostics")

    param <- match.arg(param)

    if(param == "block"){
        post.predict.block.dynsbm.mcmc(object,tclass,...)
    }else if(param == "sender"){
        post.predict.sender.dynsbm.mcmc(object,tclass,...)
    }else if(param == "receiver"){
        post.predict.receiver.dynsbm.mcmc(object,tclass,...)
    }

}


post.predict.block.dynsbm.mcmc <- function(object,tclass,blocks=NULL,
                                           sblocks=NULL,rblocks=NULL,
                                           remap=NULL,iter.draws=1,...){

    kk <- object$kk
    total <- object$chain$mcmc.control$total

    if(is.null(blocks)) blocks <- 1:kk
    if(is.null(sblocks)) sblocks <- blocks
    if(is.null(rblocks)) rblocks <- blocks
    if(is.null(remap)) remap <-1:kk

    samp <- array(NA,c(length(sblocks), length(rblocks), total*iter.draws))
    BB.prior <- object$chain$BB.prior[,,,,tclass]

    for(ii in 1:length(sblocks)){
        ss <- sblocks[ii]
        for(jj in 1:length(rblocks)){
            rr <- rblocks[jj]
            mat <- t(BB.prior[remap[ss],remap[rr],,])
            samp[ii,jj,] <- post.pred.gamma(mat,iter.draws=iter.draws)
        }
    }
    return(samp)
}



post.predict.sender.dynsbm.mcmc <- function(object,tclass,
                                            nodes=NULL,block=NULL,
                                            iter.draws=1,...){
    nn <- object$nn; kk <- object$kk
    total <- object$chain$mcmc.control$total

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

    samp <- array(NA,c(length(nodes), total*iter.draws))
    prior <- object$chain$SS.prior[,,,tclass]

    for(ii in 1:length(nodes)){
        mat <- t(prior[nodes[ii],,])
        samp[ii,] <- post.pred.gamma(mat,iter.draws=iter.draws, ...)
    }

    return(samp)
}


post.predict.receiver.dynsbm.mcmc <- function(object,tclass,
                                              nodes=NULL,block=NULL,
                                              iter.draws=1,...){

    nn <- object$nn; kk <- object$kk
    total <- object$chain$mcmc.control$total

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

    samp <- array(NA,c(length(nodes), total*iter.draws))
    prior <- object$chain$RR.prior[,,,tclass]

    for(ii in 1:length(nodes)){
        mat <- t(prior[nodes[ii],,])
        samp[ii,] <- post.pred.gamma(mat,iter.draws=iter.draws)
    }

    return(samp)
}

#' Posterior Predictive Plotting for DynSBM Parameters
#'
#' @param object object of class "dynsbm"
#' @param tclass equivalence class for which to plot posterior predictive
#'   distribution
#' @inheritParams param.plot
#' @export
post.pred.plot <- function(object, tclass,
                           param=c("block","sender","receiver"),
                           ...){

    if(object$chain$clean) stop("clean.out must be FALSE to produce diagnostics")
    param <- match.arg(param)

    if(param == "block"){
        block.post.pred.plot(object, tclass, ...)
    }else if(param == "sender"){
        sender.post.pred.plot(object, tclass, ...)
    }else if(param == "receiver"){
        receiver.post.pred.plot(object, tclass, ...)
    }
}



block.post.pred.plot <- function(object,tclass,
                                 sblock,rblock,
                                 remap=NULL, col="steelblue2",add=FALSE,
                                 include.post=TRUE,
                                 iter.draws=50, lwd=2,
                                 alpha=.05,xlim=NULL,
                                 xlab=NULL,ylab=NULL,main=NULL,
                                 ...){

    if(length(sblock) > 1 || sblock > object$kk || sblock < 1)
        stop("sblock must be a single integer between 1 and kk")
    if(length(rblock) > 1 || rblock > object$kk || rblock < 1)
        stop("rblock must be a single integer between 1 and kk")
    if(is.null(remap)) remap <- 1:object$kk

    tt.vec <- which(object$tmap == tclass)


    samp <- param.post.predict(object=object,tclass=tclass,param="block",
                               sblocks=sblock,rblocks=rblock,
                               iter.draws=iter.draws,...)

    mat <- NULL
    if(include.post){
        mat <- object$chain$BB.mat[remap[sblock],
                                   remap[rblock],,]
        mat <- mat[,tt.vec,drop=FALSE]
    }

    if(is.null(xlab)) xlab <- paste0("B[",sblock,",",rblock,"]")

    pp.plot(samp=samp,mat=mat,include.post=include.post,add=add,
            lwd=lwd,alpha=alpha,xlim=xlim,col=col,
            xlab=xlab,ylab=ylab,main=main,...)

}



sender.post.pred.plot <- function(object,tclass,node,
                                  remap=NULL, col="steelblue2",add=FALSE,
                                  include.post=TRUE,
                                  iter.draws=50, lwd=2,
                                  alpha=.05,xlim=NULL,
                                  xlab=NULL,ylab=NULL,main=NULL,
                                  ...){

    if(length(node) > 1 || node < 1 || node > object$nn)
        stop("node must be a single integer between 1 and nn")

    tt.vec <- which(object$tmap == tclass)


    samp <- param.post.predict(object=object,tclass=tclass,param="sender",
                               nodes=node,iter.draws=iter.draws,...)
    mat <- NULL
    if(include.post){
        mat <- object$chain$SS.mat[node,,]
        mat <- mat[,tt.vec,drop=FALSE]
    }

    if(is.null(xlab)) xlab <- paste0("S[",node,"]")
    pp.plot(samp=samp,mat=mat,include.post=include.post,add=add,
            lwd=lwd,alpha=alpha,xlim=xlim,col=col,
            xlab=xlab,ylab=ylab,main=main,...)
}



receiver.post.pred.plot <- function(object,tclass,node,
                                    remap=NULL, col="steelblue2",add=FALSE,
                                    include.post=TRUE,
                                    iter.draws=50, lwd=2,
                                    alpha=.05,xlim=NULL,
                                    xlab=NULL,ylab=NULL,main=NULL,
                                    ...){

    if(length(node) > 1 || node < 1 || node > object$nn)
        stop("node must be a single integer between 1 and nn")

    tt.vec <- which(object$tmap == tclass)


    samp <- param.post.predict(object=object,tclass=tclass,param="receiver",
                               nodes=node,iter.draws=iter.draws,...)

    mat <- NULL
    if(include.post){
        mat <- object$chain$RR.mat[node,,]
        mat <- mat[,tt.vec,drop=FALSE]
    }

    if(is.null(xlab)) xlab <- paste0("R[",node,"]")
    pp.plot(samp=samp,mat=mat,include.post=include.post,add=add,
            lwd=lwd,alpha=alpha,xlim=xlim,col=col,
            xlab=xlab,ylab=ylab,main=main,...)
}


######################################################################
############## Performing Posterior Predictive Test  #################
######################################################################


#' Posterior Predictive Plotting for DynSBM Parameters
#'
#' @param object.1 base object of class "dynsbm.mcmc"
#' @param object.2 test object of class "dynsbm.mcmc"
#' @inheritParams post.pred.plot
#'
#' @return Returns result of a posterior predictive test on the networks in
#'   object.2 for the given parameters.  The test is to determine if the
#'   parameters for each network in object.2 came from the same distribution
#'   as the networks in object.1.  It is assumed that all of the networks in
#'   object.1 come from the same underlying distribution for each equivalence
#'   class.
#' @export
post.pred.test <- function(object.1,object.2,
                           param=c("block","sender","receiver","all"),
                           ...){

    if(object.1$chain$clean) stop("clean.out must be FALSE to use posterior")
    if(object.2$chain$clean) stop("clean.out must be FALSE to use posterior")

    param <- match.arg(param)

    if(param == "block"){
        block.post.pred.test(object.1,object.2, ...)
    }else if(param == "sender"){
        node.post.pred.test(object.1,object.2, param="sender", ...)
    }else if(param == "receiver"){
        node.post.pred.test(object.1,object.2, param="receiver", ...)
    }else if(param == "all"){
        out.df <- block.post.pred.test(object.1,object.2,...)
        out.df <- rbind(out.df,
                        node.post.pred.test(object.1,object.2,
                                            param="sender",...))
        out.df <- rbind(out.df,
                        node.post.pred.test(object.1,object.2,
                                            param="receiver",...))
        return(out.df)
    }
}





single.block.test <- function(object.1, object.2, tclass,
                           alpha=0.05,iter.draws=10){
    ## browser()
    if(missing(object.1)) stop("object.1 must be provided")
    if(missing(object.2)) object.2 <- object.1
    if(tclass > object.1$ee) stop("tclass is greater than ee: ",
                                  object.1$ee)

    kk <- object.1$kk; nn <- object.1$nn; TT <- object.1$TT

    samp <- param.post.predict(object=object.1,tclass=tclass,param="block",
                               iter.draws=iter.draws)

    anom.blocks <- data.frame(Parameter=character(),
                              Sender=numeric(),Receiver=numeric(),
                              Class=numeric(),
                              Intensity=numeric(),Level=numeric(),
                              stringsAsFactors=FALSE)
    for(ss in 1:kk){
        for(rr in 1:kk){

            qq <- quantile(samp[ss,rr,], c(alpha/2, 1 - (alpha/2)))
            tt.vec <- which(object.2$tmap == tclass)

            subsamp <- samp[ss,rr,]
            mat <- object.2$chain$BB.mat[ss,rr,,][,tt.vec,drop=FALSE]
            lower.ests <- apply(mat,2,quantile,prob=c(alpha/2))
            upper.ests <- apply(mat,2,quantile,prob=c(1-(alpha/2)))
            for(ii in 1:length(tt.vec)){
                if(upper.ests[ii] < qq[1]){
                    point.qq <- find.quantile(subsamp,mat[,ii],twosided=TRUE)
                    anom.blocks[nrow(anom.blocks) + 1,]<- c("block",
                                                          ss,rr,tt.vec[ii],
                                                          point.qq)
                }else if(lower.ests[ii] > qq[2]){
                    point.qq <- find.quantile(subsamp,mat[,ii],twosided=TRUE)
                    anom.blocks[nrow(anom.blocks) + 1,]<- c("block",
                                                          ss,rr,tt.vec[ii],
                                                          point.qq)
                }
            }
        }
    }
    names(anom.blocks) <- c("Parameter","Sender","Receiver","Class",
                            "Intensity","Level")
    return(anom.blocks)
}


block.post.pred.test <- function(object.1, object.2, alpha=0.05,
                                 iter.draws = 10){


    if(missing(object.1)) stop("object must be provided")
    if(missing(object.2)) object.2 <- object.1
    tclasses <- 1:object.1$ee

    anom.block.mat <- NULL
    for(tt in tclasses){
        anom.block.mat <- rbind(anom.block.mat,
                                single.block.test(object.1=object.1,
                                                  object.2=object.2,
                                                  tclass=tt,alpha=alpha,
                                                  iter.draws=iter.draws))
    }

    return(anom.block.mat)
}




single.node.test <- function(object.1, object.2, tclass,
                             param=c("sender","receiver"),
                             alpha=0.05,iter.draws=10){
    ## browser()
    param <- match.arg(param)
    if(missing(object.1)) stop("object.1 must be provided")
    if(missing(object.2)) object.2 <- object.1
    if(tclass > object.1$ee) stop("tclass is greater than ee: ",
                                  object.1$ee)

    nn <- object.1$nn; TT <- object.1$TT

    samp <- param.post.predict(object=object.1,tclass=tclass,param=param,
                               iter.draws=iter.draws)
    if(param == "sender") node.mat <- object.2$chain$SS.mat
    else if(param == "receiver") node.mat <- object.2$chain$RR.mat

    anom.nodes <- data.frame(Parameter=character(),
                             Sender=numeric(),Receiver=numeric(),
                             Class=numeric(),
                             Intensity=numeric(),Level=numeric(),
                             stringsAsFactors=FALSE)

    for(node in 1:nn){

        qq <- quantile(samp[node,], c(alpha/2, 1 - (alpha/2)))
        tt.vec <- which(object.2$tmap == tclass)

        subsamp <- samp[node,]
        mat <- node.mat[node,,][,tt.vec,drop=FALSE]

        lower.ests <- apply(mat,2,quantile,prob=c(alpha/2))
        upper.ests <- apply(mat,2,quantile,prob=c(1-(alpha/2)))
        for(ii in 1:length(tt.vec)){
            if(upper.ests[ii] < qq[1]){
                point.qq <- find.quantile(subsamp,mat[,ii],twosided=TRUE)
                ## point.qq <- mean(samp[node,] < upper.ests[ii])
                anom.nodes[nrow(anom.nodes) + 1,]<- c(param,
                                                      node,node,tt.vec[ii],
                                                      point.qq)

            }else if(lower.ests[ii] > qq[2]){
                point.qq <- find.quantile(subsamp,mat[,ii],twosided=TRUE)
                ## point.qq <- mean(samp[node,] < lower.ests[ii])
                anom.nodes[nrow(anom.nodes) + 1,]<- c(param,
                                                      node,node,tt.vec[ii],
                                                      point.qq)
            }
        }
    }

    ## names(anom.nodes) <- c("Sender","Receiver","Class",
    ##                        "Intensity","Level")

    if(nrow(anom.nodes) > 0){
        if(param == "sender") anom.nodes$Receiver <- NA
        if(param == "receiver") anom.nodes$Sender <- NA
    }

    return(anom.nodes)
}


node.post.pred.test <- function(object.1, object.2, alpha=0.05,
                                iter.draws = 10,
                                param=c("sender","receiver")){

    param <- match.arg(param)

    if(missing(object.1)) stop("object must be provided")
    if(missing(object.2)) object.2 <- object.1
    tclasses <- 1:object.1$ee

    anom.mat <- NULL
    for(tt in tclasses){
        anom.mat <- rbind(anom.mat,
                          single.node.test(object.1=object.1,
                                           object.2=object.2,
                                           tclass=tt,alpha=alpha,
                                           iter.draws=iter.draws,
                                           param=param))
    }

    return(anom.mat)
}



#' @describeIn get.iter returns "dynsbm" model object for a specific iteration.
#' @export
get.iter.dynsbm.mcmc <- function(object,iter, ...){
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
}
