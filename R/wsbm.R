#####  sbm.R
#####  Functions Implementing SBM and MMSBM Fitting in C:
#####    sbm - Fits SBM Model
#####    mmsbm.C - Fits MMSBM Model
#####    mmsbm - Fits MMSBM choosing intelligent starting values
#####
##############################################################


######################################################################
###################  DATA GENERATION FUNCTION  #######################
######################################################################



#' Create Weighted Stochastic Block Model Object
#'
#' Generates a Model Object for a Weighted Stochastic Block Model
#'
#' @param nn number of nodes in the network
#' @param mmb block memberhsip vector
#' @param BB block propensity matrix
#' @param SS sender effect vector
#' @param RR receiver effect vector
#' @param mmb.prior prior distribution for multinomial distribution over
#' block memberships
#' @param self.ties if true self ties are allowed
#' @param hours scaling factor so that other factors can be interpreted
#' as hourly rate parameters.
#' @param gen if TRUE also samples a network from the model and includes it as a
#' component of the model object
#'
#' @return Returns a "wsbm" model object containing all the parameters needed
#' for sampling networks from the model.  If gen = TRUE, the object also
#' contains a network sampled from this parameterization.
#'
#' @examples
#' kk <- 3
#' nn <- kk^2 * 10
#' mmb <-  rep(1:kk,each=nn/kk)
#' BB <- array(1,c(kk,kk)); diag(BB) <- 5
#' SS <- RR <- rep(c(rep(5,nn/(kk^2)),rep(1,nn*(kk-1)/(kk^2))),kk)
#' set.seed(100)
#' dat <- wsbm(nn=nn,mmb=mmb,BB=BB,SS=SS,RR=RR,gen=TRUE)
#'
#' ### Plotting Intensity Matrix and Adjacency Matrix
#' par(mfrow=c(1,2))
#' plot(dat,main="Mean Intensity Matrix")
#' network.plot(dat,main="Adjacency Matrix")
#'
#' @seealso \code{\link{wsbm.fit}}
#' @export
wsbm <- function(nn,mmb,BB,SS=rep(1,nn),RR=rep(1,nn),
                 mmb.prior=c(.5,.5), self.ties=TRUE,hours=1,gen=FALSE){

    if(length(SS) != nn) stop("SS must have length n")
    if(length(RR) != nn) stop("RR must have length n")

    kk <- dim(BB)[1]
    if(missing(mmb)){
        if(length(mmb.prior) != kk) stop("mmb.prior must have length k")
        mmb <- sort(sample(kk,nn,replace=TRUE,prob=mmb.prior))
    }else{
        if(length(mmb) != nn) stop("mmb must have length n")
    }

    obj <- structure(list(mmb = mmb, BB=BB, SS=SS, RR=RR,
                          hours=hours, self.ties=self.ties,
                          nn=nn,kk=kk),
                     class="wsbm")

    if(gen){
        obj$net <- net.gen(obj)
    }
    return(obj)
}


#' @export
is.wsbm <- function(x) inherits(x,"wsbm")


#' Generate Weighted Stochastic Block Model Networks
#'
#' Generates networks from the weighted stochastic block model
#' with or without degree effects.
#'
#' @param x object of class "wsbm"
#'
#' @return Returns a matrix representing the adjacency matrix for a network
#' drawn from the given weighted stochastic block model
#'
#' @examples
#' kk <- 3
#' nn <- kk^2 * 10
#' mmb <-  rep(1:kk,each=nn/kk)
#' BB <- array(1,c(kk,kk)); diag(BB) <- 5
#' SS <- RR <- rep(c(rep(5,nn/(kk^2)),rep(1,nn*(kk-1)/(kk^2))),kk)
#' set.seed(100)
#' model <- wsbm(nn=nn,mmb=mmb,BB=BB,SS=SS,RR=RR)
#' net <- net.gen(model)
#'
#' @seealso \code{\link{wsbm}} \code{\link{wsbm.fit}}
#' @export
net.gen.wsbm <- function(x){

    PI <- mmb.to.PI(x$mmb)
    pmat <- PI %*% x$BB %*% t(PI)
    sr.mat <- as.matrix(x$SS) %*% t(as.matrix(x$RR))
    pmat <- pmat * sr.mat * x$hours
    net <- net.gen(pmat,self.ties=x$self.ties,d="poisson")

    ## net <- array(rpois(x$nn*x$nn,pmat),c(x$nn,x$nn))
    ## if(!(x$self.ties)){
    ##     diag(net) <- NA
    ## }

    return(net)
}


######################################################################
###################  MODEL FITTING FUNCTIONS  ########################
######################################################################

#' Weighted Stochastic Block Model Estimator
#'
#' Estimates parameters for the weighted stochastic block model using
#' an MCMC algorithm.
#'
#' @param net a matrix object representing the adjacency matrix
#' @param kk number of blocks
#' @param hours scaling parameter
#' @param self.ties if true assumes self ties are possible
#' @param priors priors for posterior inference.  See wsbm.priors for more
#' details
#' @param init.control control parameters for initialization routine.
#' If init.control contains a member named init.vals it is interpreted as a
#' list of parameters from which to start the MCMC chain.
#' @param mcmc.control control parameters for the MCMC algorithm
#' @param clean.out if true, removes MCMC chain from output and only
#' returns posterior means
#' @param verbose higher values correspond to more informative output as the
#' sampler runs
#'
#' @return Returns a wsbm object that has posterior means for each parameter
#' in the weighted directed degree corrected stochastic block model.  If
#' clean.out is FALSE, this object also contains an element 'chain' which
#' contains the thinned and burned in draws from the MCMC sampler.  The prior
#' parameters from the call to wsbm are also included.
#'
#' @seealso \code{\link{wsbm.priors}}
#' @export
wsbm.fit <- function(net, kk=3, hours=1, self.ties=TRUE,
                     priors=wsbm.priors(eta=rep(1/kk,kk)),
                     init.control=list(spectral.start=FALSE,
                                       multistart=0,
                                       multistart.total=20),
                     mcmc.control=list(total=1000,burn.in=1000,thin=10,
                                       extend.alpha=.001,extend.max=10,
                                       label.switch.mode="adhoc",
                                       label.switch.max=200,
                                       multi.impute=FALSE),
                     clean.out=FALSE, verbose=1){

    ##  Formatting Adjacency Matrix for C
    ## label.switch.mode <- match.arg(label.switch.mode)

    mcmc.control <- parse.mcmc.control(mcmc.control)
    init.control <- parse.init.control(init.control)

    ##  Loading Initial WSBM Object
    wsbm.init <- load.wsbm.init(init.control,priors,net,kk,hours,self.ties)
    wsbm.init$net <- net


    multi.int <- as.integer(ifelse(mcmc.control$multi.impute,1,0))
    net.clean <- net
    if(!self.ties){
        diag(net.clean) <- NA
    }
    net.clean[is.na(net)] <- -1
    nn <- nrow(net.clean)

    total <- mcmc.control$total
    thin <- mcmc.control$thin
    burn.in <- mcmc.control$burn.in

    flatMMB = integer(total * nn)
    flatBB = double(total * kk*kk)
    flatSS = double(total * nn)
    flatRR = double(total * nn)
    ll.vec <- double(total)
    HH.vec <- double(total*kk*nn)


    ##  Passing Initial Parameters
    start = 1
    flatMMB[1:nn] <- wsbm.init$mmb
    flatBB[1:(kk*kk)] <- wsbm.init$BB
    flatSS[1:nn] <- wsbm.init$SS
    flatRR[1:nn] <- wsbm.init$RR
    ll.vec[1] <- logLik(wsbm.init)
    HH.vec[1:(kk*nn)] <- mmb.to.PI(wsbm.init$mmb)


    extend.max <- as.integer(mcmc.control$extend.max)
    extend.qq <- as.double(qnorm(1 - mcmc.control$extend.alpha/2))


    ##  Calling C Implementation of MCMC
    #' @useDynLib lvnm wsbm_R
    out <- .C(wsbm_R,as.integer(total),as.integer(nn),as.integer(kk), ##3
              as.integer(net.clean), #4
              as.double(c(priors$sender.alpha,priors$sender.beta)), ## 5
              as.double(c(priors$receiver.alpha,priors$receiver.beta)), ## 6
              as.double(c(priors$block.alpha,priors$block.beta)), ## 7
              as.double(priors$eta), ## flatVec, 8
              flatBB, flatMMB, flatSS, flatRR, ## 12
              as.integer(burn.in), as.integer(thin), ##14
              as.integer(start),multi.int,ll.vec, ##17
              extend.max,extend.qq,HH.vec, ##20
              as.double(hours),as.integer(verbose))

    ##  Pulling the flat matrices from the C output
    BB.flat <- matrix(out[[9]],nrow=kk*kk)
    BB <- array(BB.flat,c(kk,kk,total))

    mmb <- array(out[[10]],c(nn,total))
    SS <- array(out[[11]],c(nn,total))
    RR <- array(out[[12]],c(nn,total))


    ll.vec <- as.vector(out[[17]]) - sum(lfactorial(net))

    HH.flat <- out[[20]]
    HH <- array(HH.flat,c(nn,kk,total))

    diag(net.clean) <- -1

    wsbm.chain.obj <- structure(list(mmb=mmb,BB=BB,
                                     SS=SS,RR=RR,HH=HH,
                                     logLik=ll.vec,
                                     hours=hours,self.ties=self.ties,
                                     nn=nn,kk=kk,
                                     mcmc.control=mcmc.control,
                                     init.control=init.control,priors=priors,
                                     clean=FALSE),
                                class="wsbm.chain")

    if(mcmc.control$label.switch.mode == "kl-loss"){
        wsbm.chain.obj <- switch.labels(wsbm.chain.obj,verbose=verbose,
                                        max.runs=mcmc.control$label.switch.max)
    }else{
        ##  Do nothing
    }

    ## Summarizing MCMC Chain
    wsbm.obj <- summary(wsbm.chain.obj)
    wsbm.obj$net <- net

    ## Calculating Log-Likelihood and Information Criteria
    wsbm.fit.obj <- wsbm.obj
    wsbm.fit.obj$chain <- wsbm.chain.obj
    class(wsbm.fit.obj) <- c("wsbm.fit","wsbm.mcmc","wsbm")

    if(clean.out){
        wsbm.fit.obj$chain$mmb <- NULL
        wsbm.fit.obj$chain$BB <- NULL
        wsbm.fit.obj$chain$SS <- NULL
        wsbm.fit.obj$chain$RR <- NULL
        wsbm.fit.obj$chain$HH <- NULL
        wsbm.fit.obj$chain$clean <- TRUE
    }

    return(wsbm.fit.obj)
}

#' @export
is.wsbm.fit <- function(x) inherits(x,"wsbm.fit")

#' @export
is.wsbm.mcmc <- function(x) inherits(x,"wsbm.mcmc")


#' Prior Parameters for wsbm Gamma Priors
#'
#' Generating prior paramerters for gamma distributions over Poisson rate
#' parameters
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
#' @seealso \code{\link{wsbm}}
#' @export
wsbm.priors <- function(eta,
                        block.alpha=1, block.beta=1,
                        sender.alpha=10, sender.beta=10,
                        receiver.alpha=10, receiver.beta=10){

    if(block.alpha < 0) stop("block.alpha is negative")
    if(block.beta < 0 ) stop("block.beta is negative")

    if(sender.alpha < 0) stop("sender.alpha is negative")
    if(sender.beta < 0 ) stop("sender.beta is negative")

    if(receiver.alpha < 0) stop("receiver.alpha is negative")
    if(receiver.beta < 0 ) stop("receiver.beta is negative")

    return(list(block.alpha=block.alpha,block.beta=block.beta,
                sender.alpha=sender.alpha,sender.beta=sender.beta,
                receiver.alpha=receiver.alpha,receiver.beta=receiver.beta,
                eta=eta))
}

######################################################################
#####################  PREDICTION FUNCTIONS  #########################
######################################################################


#' Predict Method for WSBM Fits
#'
#' Generates a matrix of predicted intensities for each node pair
#'
#' @param object a fitted object of the class "wsbm"
#'
#' @return Returns a matrix of predicted intensities.  These intensities are the
#' predicted values given the posterior mean parameter estimates
#'
#' @examples
#' data(commuter30)
#' fit.1 <- wsbm.fit(commuter30,kk=3,hours=1)
#' pred.mat <- predict(fit.1)
#'
#' @seealso \code{\link{wsbm.fit}}, \code{\link{wsbm}}
#' @export
predict.wsbm <- function(object,...){

    PI <- mmb.to.PI(object$mmb,kk=ncol(object$BB))
    pmat <- PI %*% object$BB %*% t(PI)
    if(!object$self.ties) diag(pmat) <- NA

    sr.mat <- as.matrix(object$SS) %*% t(as.matrix(object$RR))
    pmat <- pmat * sr.mat * object$hours
    return(pmat)
}


#' Posterior Predictive Distribution Method for WSBM Fits
#'
#' Generates a posterior predictive distribution for each tie in a network.
#'
#' @param x a fitted object of the class "wsbm".  Requires clean.out = FALSE
#' in the call to wsbm.
#'
#' @return Returns a 3 dimensional array of dimension n x n x total, where n
#' is the number of nodes in the network and total is the number of draws in
#' the MCMC chain generated using wsbm.
#'
#' @examples
#' data(commuter30)
#' fit.1 <- wsbm.fit(commuter30,kk=3,hours=1)
#'
#' pp.dist <- post.pred(fit.1)
#' density.ppd <- apply(post.pred,3,mean)
#' summary(density.ppd)
#'
#'
#' @seealso \code{\link{wsbm.fit}} \code{\link{predict.wsbm}}
#' @export
post.predict.wsbm.mcmc <- function(x, ...){
    if(x$chain$clean) stop("Use clean.out = TRUE to keep MCMC chains")
    nn <- x$nn
    total <- x$chain$mcmc.control$total

    pmat.post <- array(NA,c(nn,nn,total))
    for(ii in 1:total){
        pmat.post[,,ii] <- net.gen(get.iter(x,ii))
        ## tmp.obj <- with(x$chain,
        ##                 wsbm(nn=nn,BB=BB[,,ii],SS=SS[,ii],RR=RR[,ii],
        ##                      mmb=mmb[,ii],self.ties=self.ties,hours=hours))

        ## pmat.post[,,ii] <- net.gen(tmp.obj)

        ## with(x$chain,list(BB=BB[,,ii],SS=SS[,ii],RR=RR[,ii],
        ##                                  mmb=mmb[,ii],
        ##                                  self.ties=self.ties,hours=hours))
        ##     pmat.post[,,ii] <- data.gen.wsbm(tmp.obj)
    }
    return(pmat.post)
}


#' Predictive Matrix for WSBM Fits
#'
#' Generates predicted intensities for a WSBM model given a network.
#'
#' @param net a matrix object representing the adjacency matrix
#' @param kk number of blocks
#' @param hours scaling parameter
#' @param self.ties if true assumes self ties are possible
#' @param priors priors for posterior inference.  See wsbm.priors for more details
#' @param verbose higher values correspond to more informative output as the
#' @param ... other parameters to be passed to wsbm
#'
#' @return Returns a matrix object corresponding to the posterior mean intensities.
#'
#' @examples
#' data(commuter30)
#' met.3 <- wsbm.metric(commuter30,kk=3)
#'
#' @seealso \code{\link{wsbm}} \code{\link{predict.wsbm}}
#' @export
wsbm.metric <- function(net,kk=2,hours=1,self.ties=TRUE,
                        priors=wsbm.priors(eta=rep(1/kk,kk)),
                        verbose=0,...){

    wsbm.fit <- wsbm(net=net,kk=kk,hours=hours,self.ties=self.ties,
                     mcmc.control=mcmc.control,priors=priors,
                     verbose=verbose,clean.out=TRUE, ...)
    return(predict(wsbm.fit))
}


######################################################################
#####################  PLOTTING FUNCTIONS  #########################
######################################################################


#' Plotting Method for WSBM Fits
#'
#' Plots the expected mean tie probabilities for a given wsbm object.
#'
#' @param x object of class wsbm
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
#' data(commuter30)
#' fit.1 <- wsbm.fit(commuter30,kk=3,hours=1)
#' plot(fit.1)
#'
#' @seealso \code{\link{wsbm.fit}}
#' @export
plot.wsbm <- function(x,
                      node.order=c("default","membership"),
                      pal=grey((50:1)/50), ...){
    nn <- nrow(x$net)

    node.order <- match.arg(node.order)
    if(node.order == "default"){
        ord <- 1:nn
    }else if(node.order == "membership"){
        ord <- order(x$mmb)
    }

    pmat <- predict.wsbm(x)
    adj.image.plot(pmat,ord,pal,...)
}


#' @export
network.plot.wsbm <- function(x, pal=grey((50:1)/50), node.order=NULL,...){
    net <- x$net
    if(is.null(net))
        stop("x must have a net component.  Try calling wsbm with gen=TRUE")

    network.plot(net, pal=pal,node.order=node.order,...)

}





#' Diagnostic Plotting Method for WSBM fits
#'
#' Plots MCMC chains for various parameters in WSBM model.
#'
#' @param x object of class "wsbm"
#' @param type which parameters diagnostics to plot.  The default value "block"
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
#' fit.1 <- wsbm.fit(commuter30,kk=3,hours=1)
#'
#' diagnostic.plot(fit.1,param="block")
#' diagnostic.plot(fit.1,param="sender")
#' diagnostic.plot(fit.1,param="receiver")
#' diagnostic.plot(fit.1,param="ll")
#'
#' @seealso \code{\link{wsbm.fit}}
#' @export
diagnostic.plot.wsbm.mcmc <- function(x,
                                      param=c("block","sender","receiver","ll"),
                                      scale.ylim=FALSE, ...){
    if(x$chain$clean) stop("clean.out must be FALSE to produce diagnostics")

    param <- match.arg(param)
    if(param == "block"){
        block.diagnostic.wsbm.mcmc(x,scale.ylim=scale.ylim, ...)
    }else if(param == "sender"){
        sender.diagnostic.wsbm.mcmc(x, scale.ylim=scale.ylim, ...)
    }else if(param == "receiver"){
        receiver.diagnostic.wsbm.mcmc(x, scale.ylim=scale.ylim, ...)
    }else if(param == "ll"){
        ll.diagnostic.wsbm.mcmc(x, ...)
    }

}


block.diagnostic.wsbm.mcmc <- function(x,
                                       xlim=NULL,ylim=NULL,
                                       scale.ylim=FALSE,
                                       blocks=NULL,...){

    kk <- x$kk
    if(is.null(ylim) && scale.ylim) ylim <- range(x$chain$BB)
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
            plot(x$chain$BB[ii,jj,],ylim=ylim,xlim=xlim,type="l",
                 xlab="Iteration",ylab="",
                 main=paste0("B[",ii,",",jj,"]"),
                 ...)
        }
    }
    title(main="Block Effect Chains",outer=TRUE,cex.main=3)
    par(old.par)

}

sender.diagnostic.wsbm.mcmc <- function(x,
                                        xlim=NULL,ylim=NULL,
                                        scale.ylim=FALSE,
                                        nodes=NULL,...){

    node.diag.plot(x$chain$SS,
                   xlim=xlim,ylim=ylim,scale.ylim=scale.ylim,
                   main.prefix="S[",main.suffix="]",nodes=nodes)
    title(main="Sender Effect Chains",outer=TRUE,
          cex.main=2.5,line=-3)

}


receiver.diagnostic.wsbm.mcmc <- function(x,
                                          xlim=NULL,ylim=NULL,
                                          scale.ylim=FALSE,
                                          nodes=NULL,...){

    node.diag.plot(x$chain$RR,
                   xlim=xlim,ylim=ylim,scale.ylim=scale.ylim,
                   main.prefix="R[",main.suffix="]",nodes=nodes)
    title(main="Receiver Effect Chains",outer=TRUE,
          cex.main=2.5,line=-3)
}


ll.diagnostic.wsbm.mcmc <- function(x,...){

    ll.chain.plot(x,...)

}


node.boxplot <- function(mat, nodes=NULL, ...){

    nn <- dim(mat)[1]
    total <- dim(mat)[2]

    if(is.null(nodes)) nodes <- 1:nn
    par(cex.main=2,cex.lab=1.5,cex.axis=1.2)
    boxplot(t(mat[nodes,]),...)
}



#' Parameter Plotting Method for WSBM fits
#'
#' Plots boxplots of posterior samples for parameters in WSBM model.
#'
#' @param x object of class "wsbm"
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
#' fit.1 <- wsbm.fit(commuter30,kk=3,hours=1)
#'
#' param.plot(fit.1,param="block")
#' param.plot(fit.1,param="sender")
#' param.plot(fit.1,param="receiver")
#' param.plot(fit.1,param="ll")
#'
#' @seealso \code{\link{wsbm.fit}}
#' @export
param.plot.wsbm.mcmc <- function(x,
                                 param=c("block","sender","receiver"),
                                 ...){
    if(x$chain$clean) stop("clean.out must be FALSE to produce diagnostics")

    param <- match.arg(param)
    if(param == "block"){
        block.param.plot.wsbm.mcmc(x, ...)
    }else if(param == "sender"){
        sender.param.plot.wsbm.mcmc(x,  ...)
    }else if(param == "receiver"){
        receiver.param.plot.wsbm.mcmc(x,  ...)
    }

}

block.param.plot.wsbm.mcmc <- function(x, blocks=NULL,
                                       scale.ylim=TRUE,ylim=NULL,...){

    kk <- x$kk
    if(is.null(blocks)) blocks <- 1:kk

    if(dev.cur() == 1){
        len <- max(length(blocks)*2,7)
        dev.new(height=len,width=len)
    }
    if(is.null(ylim) && scale.ylim) ylim <- range(x$chain$BB[blocks,blocks,])


    old.par <- par(mfrow = c(length(blocks),length(blocks)),
                   cex.main=2,cex.lab=1.5,cex.axis=1.5,
                   mar=c(5,4,2,1),oma=c(0,1,3,1))

    for(ii in blocks){
        for(jj in blocks){
            boxplot(x$chain$BB[ii,jj,],ylim=ylim,
                    xlab=paste0("B[",ii,",",jj,"]"),ylab="Effect")
        }
    }
    title(main="Block Effect Posteriors",outer=TRUE,cex.main=3)
    par(old.par)
}

sender.param.plot.wsbm.mcmc <- function(x,
                                         nodes=NULL, ...){

    node.boxplot(x$chain$SS,nodes=nodes,
                 main="Sender Effect Posteriors",
                 xlab="Sending Node",ylab="Effect",...)
}

receiver.param.plot.wsbm.mcmc <- function(x,
                                         nodes=NULL, ...){

    node.boxplot(x$chain$RR,nodes=nodes,
                 main="Receiver Effect Posteriors",
                 xlab="Receiving Node",ylab="Effect",...)
}


######################################################################
##############  Likelihood and Informaiton Criteria  #################
######################################################################


#' Log-Likelihood method for WSBM Fits
#'
#' Returns log-likelihood for posterior mean estimates.
#'
#' @param x an object of class "wsbm"
#'
#' @return returns the log-likelihood of the data given the posterior mean
#' parameter estimates.
#'
#' @seealso \code{\link{wsbm}}
#' @export
logLik.wsbm <- function(x){
    if(is.null(x$net)) stop("no network present in object")
    pmat <- predict.wsbm(x)
    return(ll.pmat.pois(net=x$net,pmat=pmat,self.ties=x$self.ties))
}


#' AIC method for WSBM Fits
#'
#' Returns Akaike Information Criterion for WSBM Fit
#'
#' @param x an object of class "wsbm"
#'
#' @return returns the Akaike Information Criterion using the posterior mean
#' parameters as approximations to the MLE.  The number of parameters is taken
#' to be 2*(n-k) + k*k + n.
#'
#' @seealso \code{\link{wsbm}}
#' @export
AIC.wsbm <- function(x){
    nn <- length(x$SS)
    kk <- dim(x$BB)[2]
    edf <- 2*(nn - kk) + kk^2 + nn
    return(-2*logLik(x) + 2 * edf)
}


#' BIC method for WSBM Fits
#'
#' Returns Bayesian Information Criterion for WSBM Fit
#'
#' @param x an object of class "wsbm"
#'
#' @return returns the Bayesian Information Criterion using the posterior mean
#' parameters as approximations to the MLE.  The number of parameters is taken
#' to be 2*(n-k) + k*k + n and the sample size is taken to be n*n.
#'
#' @seealso \code{\link{wsbm}}
#' @export
BIC.wsbm <- function(x){
    nn <- length(x$SS)
    kk <- dim(x$BB)[2]
    edf <- 2*(nn - kk) + kk^2 + nn
    log.n <- 2 * log(nrow(x$net))
    return(-2*logLik(x) + log.n * edf)
}


#' DIC method for WSBM Fits
#'
#' Returns Deviance Information Criterion for WSBM Fit
#'
#' @param x an object of class "wsbm"
#'
#' @return returns the Deviance Information Criterion.
#'
#' @seealso \code{\link{wsbm}}
#' @export
DIC.wsbm <- function(x){
    if(!is.null(x$chain$logLik)){
        DIC.df <- (2*logLik(x) - 2*mean(x$chain$logLik))
        return(-2*logLik(x) + 2 * DIC.df)
    }else{
        return(NULL)
    }
}


######################################################################
#######################  DISPLAY FUNCTIONS  ##########################
######################################################################
#' @export
print.wsbm <- function(x,...){
    cat(format(x, ...), "\n")
}


#' @export
format.wsbm <- function(x, digits=4, max.width=78, ...){

    ## head.str <- "Posterior Mean Estimates for WSBM"

    mmb.str <- paste0("Block Membership:\n",
                      format.vector(x$mmb,digits=digits,max.width=max.width,...))

    blockmat.str <- apply(format(x$BB,digits=digits),1,paste,collapse=" ")
    blockmat.str <- paste(blockmat.str,collapse="\n")
    blockmat.str <- paste("Block Matrix:\n",blockmat.str,sep="")

    ss.str <- paste0("Sender Effects:\n",
                     format.vector(x$SS,digits=digits,max.width=max.width,...))

    rr.str <- paste0("Receiver Effects:\n",
                     format.vector(x$RR,digits=digits,max.width=max.width,...))

    ll.str <- paste("Log-Likelihood: ",format(logLik(x),digits=7))
    ic.str <- paste("BIC: ",format(BIC(x),digits=7),
                    "   AIC: ",format(AIC(x), digits=7),
                    "   DIC: ",format(DIC(x),digits=7))

    ## browser()
    return(paste(mmb.str,blockmat.str,
                 ss.str,rr.str,ll.str,ic.str,sep="\n\n"))

}



#' @export
print.wsbm.mcmc <- function(x,...){
    cat(format(x, ...), "\n")
}


#' @export
format.wsbm.mcmc <- function(x, digits=4, max.width=78, ...){

    head.str <- "Posterior Mean Estimates for WSBM"
    mcmc.str <- paste0("Number of Samples:  ",x$chain$mcmc.control$total,"  ",
                       "Thinning:  ",x$chain$mcmc.control$thin,"  ",
                       "Burn In:  ",x$chain$mcmc.control$burn.in,"\n")
    return(paste(head.str,mcmc.str,format.wsbm(x),sep="\n"))
}


#################################################################
######################  wsbm class functions  ###################
#################################################################
summary.wsbm.chain <- function(x,...){

    total <- x$mcmc.control$total
    nn <- x$nn; kk <- x$kk

    BB.hat <- apply(x$BB,c(1,2),mean)
    PI.mean <- t(apply(x$mmb,1,tabulate,nbins=kk) / total)
    mmb <- apply(PI.mean,1,which.max)

    SS.hat <- rowMeans(x$SS)
    RR.hat <- rowMeans(x$RR)

    wsbm.obj <- wsbm(nn=nn,BB=BB.hat,SS=SS.hat,RR=RR.hat,mmb=mmb,
                     self.ties=x$self.ties,hours=x$hours)

    ## structure(list(mmb=mmb,BB=BB.hat,
    ##                            SS=SS.hat,RR=RR.hat,PI.mean=PI.mean),
    ##                       class="wsbm")
    ## summ.obj$self.ties <- x$self.ties
    ## summ.obj$hours <- x$hours

    return(wsbm.obj)

}



load.wsbm.init <- function(init.control,priors,net,kk,hours,self.ties){

    nn <- nrow(net)
    init.vals <- init.control$init.vals

#####  Spectral Clustering Initialization
    if(is.null(init.vals) & init.control$spectral.start){
        warning("Spectral start is currently not allowed for WSBM")
        ## if(verbose > 0) message("Initializing with Spectral Clustering...")
        ## init.vals <- sbm.spectral(net=net,kk=kk,weighted=TRUE)
        ## if(verbose > 0) message("Initialization Complete.")
    }


    if(init.control$multistart > 0){
        ll.best <- -Inf
        multi.mcmc <- list(total=init.control$multistart.total,
                           burn.in=0,thin=1,extend.max=0,
                           label.switch.mode="adhoc")
        multi.init <- list(spectral.start=FALSE,multistart=0)
        for(ii in 1:init.control$multistart){


            ##  Update to new call specification
            fit.tmp <- wsbm.fit(net=net,kk=kk,hours=hours,self.ties=self.ties,
                                priors=priors,
                                mcmc.control=multi.mcmc,
                                init.control=multi.init,
                                clean.out=FALSE,verbose=0)

            obj.tmp <- get.iter.wsbm.mcmc(fit.tmp,
                                          init.control$multistart.total)

            ll.tmp <- logLik(obj.tmp)
            if(ll.tmp > ll.best){
                ll.best <- ll.tmp
                obj.best <- obj.tmp
            }
        }

        return(obj.best)
    }else{
        if(is.null(init.vals$BB)){
            init.vals$BB <- array(rbeta(kk^2,priors$block.alpha,
                                        priors$block.beta),
                                  c(kk,kk))
        }

        if(is.null(init.vals$SS) || is.null(init.vals$RR)){
            ss.init <- rowMeans(net,na.rm=TRUE)
            rr.init <- colMeans(net,na.rm=TRUE)

            if(any(ss.init == 0) || any(rr.init == 0)){
                ss.init <- ss.init + 0.001
                rr.init <- rr.init + 0.001
            }

            init.vals$SS <- ss.init * nn / sum(ss.init)
            init.vals$RR <- rr.init * nn / sum(rr.init)
        }

        if(is.null(init.vals$mmb)){
            init.vals$mmb <- sample(kk,nn,replace=TRUE,prob=priors$eta)
        }

        init.wsbm <- wsbm(nn=nn,BB=init.vals$BB,mmb=init.vals$mmb,
                          SS=init.vals$SS,RR=init.vals$RR,
                          hours=hours,self.ties=self.ties)

        return(init.wsbm)
    }

}



#' @export
get.iter.wsbm.mcmc <- function(object,iter){
    if(is.null(object$chain)){
        stop(paste("Object does not contain MCMC chain.",
                   "Call dynsbm with clean.out=FALSE"))
    }
    chain <- object$chain

    wsbm.obj <- wsbm(nn=object$nn,mmb=chain$mmb[,iter],
                     BB=chain$BB[,,iter],
                     SS=chain$SS[,iter], RR=chain$RR[,iter],
                     hours=object$hours,
                     self.ties=object$self.ties)
    wsbm.obj$net <- object$net

    return(wsbm.obj)

}




##########################################################
##################  ROTATION FUNCTIONS  ##################
##########################################################

## wsbm.load.init.vals <- function(init.vals,nn,kk){

##     if(is.null(init.vals$BB) | is.null(init.vals$PI)){
##         stop("init.vals must be a list containing BB and PI")
##     }
##     ## Checking BB
##     if(any(dim(init.vals$BB) != kk)) stop("BB must be a kk by kk matrix")
##     if(any(init.vals$BB < 0 | init.vals$BB > 1))
##         stop("BB must have values between 0 and 1")

##     ## Checking PI
##     if(any(dim(init.vals$PI)!=c(nn,kk))) stop("PI must be an nn by kk matrix")

##     BB.init <- double(kk^2)
##     PI.init <- double(kk*nn)
##     for(jj in 1:kk){
##         BB.init[((jj-1) * kk + 1):(jj*kk)] <- init.vals$BB[jj,]
##     }
##     for(jj in 1:nn){
##         PI.init[((jj-1)*kk + 1):(jj*kk)] <- init.vals$PI[jj,]
##     }
##     flatTable <- array(c(BB.init,PI.init),c(1,kk*(kk+nn)))
##     return(flatTable)
## }
