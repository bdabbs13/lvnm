#####  sbm.R
#####  Functions Implementing SBM and MMSBM Fitting in C:
#####    sbm - Fits SBM Model
#####    mmsbm.C - Fits MMSBM Model
#####    mmsbm - Fits MMSBM choosing intelligent starting values
#####
##############################################################
                                        #library(svd)
#####  Function to generate data from WSBM model
data.gen.wsbm <- function(nn=50,
                          BB=array(c(5,1,1,5),c(2,2)),
                          pp=c(.5,.5),
                          SS=rep(1,nn),RR=rep(1,nn),
                          mmb=NULL,PI=NULL,
                          self.ties=TRUE){

    kk <- length(pp)
    ll <- round(pp*nn)

    if(is.null(PI)){
        if(is.null(mmb)){
            mmb <- sort(sample(kk,nn,replace=TRUE,prob=pp))
        }else{
            if(length(mmb) != nn)
                stop("mmb must have length n")
        }
        PI <- mmb.to.PI(mmb)
    }else{
        if(ncol(PI) != ncol(BB) | nrow(PI) != nn)
            stop("PI must have dim ",nn," by ",ncol(BB))
    }

    pmat <- PI %*% BB %*% t(PI)
    sr.mat <- as.matrix(SS) %*% t(as.matrix(RR))
    pmat <- pmat * sr.mat
    net <- array(rpois(nn*nn,pmat),c(nn,nn))
    if(!self.ties){
        diag(net) <- NA
    }
    return(list(net=net,PI=PI,BB=BB,SS=SS,RR=RR,pmat=pmat))
}



wsbm.priors <- function(block.alpha=1, block.beta=1,
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
                receiver.alpha=receiver.alpha,receiver.beta=receiver.beta))
}

wsbm.convergence <- function(alpha=0.001, extend.max=20,shift.size=100){
    if(alpha < 0) stop("alpha is negative")
    if(extend.max < 0) stop("extend.max is negative")
    if(shift.size < 0) stop("shift.size is negative")

    return(list(alpha=alpha,extend.max = as.integer(extend.max),
                shift.size = as.integer(shift.size)))
}


###  Wrapper for C Implementation of WSBM
wsbm <- function(net, total=1000, burn.in=0, thin=1, kk=3, self.ties=TRUE,
                 priors=wsbm.priors(), eta=rep(1/kk,kk),
                 init.vals=NULL,spectral.start=FALSE,
                 clean.out=TRUE, verbose=0, max.runs=200,
                 label.switch.mode = c("adhoc","kl-loss"),
                 autoconverge=wsbm.convergence(),
                 multiImpute=FALSE,hours=1){

    ##  Formatting Adjacency Matrix for C
    label.switch.mode <- match.arg(label.switch.mode)
    multi.int <- as.integer(ifelse(multiImpute,1,0))
    net.clean <- net
    if(!self.ties){
        diag(net.clean) <- NA
    }
    net.clean[is.na(net)] <- -1
    nn <- nrow(net.clean)

    short.total <- total
    total <- short.total*thin + burn.in
    flatBB = double(short.total * kk*kk)
    flatMMB = integer(short.total * nn)
    flatSS = double(short.total * nn)
    flatRR = double(short.total * nn)
    ll.vec <- double(short.total)
    HH.vec <- double(short.total*kk*nn)

    if(is.null(init.vals) & spectral.start){
        if(verbose > 0) message("Initializing with Spectral Clustering...")
        init.vals <- sbm.spectral(net=net,kk=kk,weighted=TRUE)
        if(verbose > 0) message("Initialization Complete.")

    }

    if(is.null(init.vals$BB)){
        init.vals$BB <- array(rbeta(kk^2,priors$block.alpha,priors$block.beta),
                              c(kk,kk))
    }
    if(is.null(init.vals$mmb)){
        init.vals$mmb <- sample(kk,nn,replace=TRUE,prob=eta)
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
    if(is.null(init.vals$HH)){
        init.vals$HH <- mmb.to.PI(init.vals$mmb)
    }
    if(is.null(init.vals$logLik)){
        init.vals$logLik <- wsbm.log.like.net(net,BB=init.vals$BB,
                                              mmb=init.vals$mmb,
                                              SS=init.vals$SS, RR=init.vals$RR,
                                              self.ties=self.ties)
    }

    start = 1
    flatBB[1:(kk*kk)] <- init.vals$BB
    flatMMB[1:nn] <- init.vals$mmb
    flatSS[1:nn] <- init.vals$SS
    flatRR[1:nn] <- init.vals$RR
    ll.vec[1] <- init.vals$logLik
    HH.vec[1:(kk*nn)] <- init.vals$HH


    extend.max <- autoconverge$extend.max
    shift.size <- autoconverge$shift.size
    qq <- as.double(qnorm(1 - autoconverge$alpha/2))


    ##  Calling C Implementation of MCMC
    out <- .C("wsbm",as.integer(total),as.integer(nn),as.integer(kk), ##3
              as.integer(net.clean), #4
              as.double(c(priors$sender.alpha,priors$sender.beta)), ## 5
              as.double(c(priors$receiver.alpha,priors$receiver.beta)), ## 6
              as.double(c(priors$block.alpha,priors$block.beta)), ## 7
              as.double(eta), ## flatVec, 8
              flatBB, flatMMB, flatSS, flatRR, ## 12
              as.integer(burn.in), as.integer(thin), ##14
              as.integer(start),multi.int,ll.vec, ##17
              extend.max,shift.size,qq,HH.vec, ##21
              as.double(hours),as.integer(verbose))

                                        #  browser()
    ##  Pulling the flat matrices from the C output
    BB.flat <- matrix(out[[9]],nrow=kk*kk)
##    BB.flat <- apply(BB.flat,2,transpose.vector,nrow=kk)
    BB <- array(BB.flat,c(kk,kk,short.total))

    mmb <- array(out[[10]],c(nn,short.total))
    SS <- array(out[[11]],c(nn,short.total))
    RR <- array(out[[12]],c(nn,short.total))


    ll.vec <- as.vector(out[[17]])
    HH.flat <- out[[21]]
    HH <- array(HH.flat,c(nn,kk,short.total))

    pmat.mat <- NULL  ## I might want to compute this later...
    diag(net.clean) <- -1


    wsbm.out <- structure(list(BB=BB,mmb=mmb,
                               SS=SS,RR=RR,
                               net=net,logLik=ll.vec,hours=hours,
                               burn.in=burn.in,thin=thin,self.ties=self.ties,
                               pmat=pmat.mat,HH=HH),class="wsbm")

    if(label.switch.mode == "kl-loss"){
        wsbm.out <- switch.labels(wsbm.out,max.runs=max.runs,verbose=verbose)
    }else{
        ##  Do nothing
    }
                                        #  browser()
    ## Summarizing MCMC Chain
    wsbm.summ <- summary(wsbm.out)
    wsbm.summ$priors = priors

    if(clean.out){
        wsbm.summ$chain <- list(logLik=wsbm.out$logLik)
        wsbm.summ$clean <- TRUE
    }else{
        wsbm.summ$chain <- wsbm.out
        wsbm.summ$chain$net <- NULL
        wsbm.summ$clean <- FALSE
    }

    return(wsbm.summ)
}

##fast.mle <- function(mmb,net,nn=length(mmb),kk=max(mmb)){
#
#    net <- clean.net(net)
#    BB <- array(0.0,c(kk,kk))
#    out <- .C("getBlockMatMLE",as.integer(nn),as.integer(kk),
#              as.integer(as.vector(net)),as.double(as.vector(BB)),
#              as.integer(mmb))
#    BB <- array(out[[4]],c(kk,kk))
#    return(BB)
#}

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





#################################################################
######################  wsbm class functions  ####################
#################################################################


summary.wsbm <- function(object,...){
###  browser()
    total <- dim(object$BB)[3]; kk <- dim(object$BB)[1]

    BB.hat <- apply(object$BB,c(1,2),mean)
    ##  PI.mean <- apply(object$PI[,,my.sub],c(1,2),mean)
    PI.mean <- t(apply(object$mmb,1,tabulate,nbins=kk) / total)
    mmb <- apply(PI.mean,1,which.max)

    SS.hat <- rowMeans(object$SS)
    RR.hat <- rowMeans(object$RR)

    summ.obj <- structure(list(mmb=mmb,BB=BB.hat,
                                SS=SS.hat,RR=RR.hat,PI.mean=PI.mean),class="wsbm")
    summ.obj$self.ties <- object$self.ties
    summ.obj$hours <- object$hours

    summ.obj$pmat <- predict(summ.obj)
    summ.obj$net <- object$net
    ## Calculating DIC
    summ.obj$logLik <- with(summ.obj,wsbm.ll.pmat(net,pmat,
                                                  self.ties=self.ties))
    summ.obj$DIC <- 2*summ.obj$logLik - 4 * mean(object$logLik)
    summ.obj$burn.in <- object$burn.in; summ.obj$thin <- object$thin

    return(summ.obj)

}

predict.wsbm <- function(object,...){
    PI <- mmb.to.PI(object$mmb,kk=ncol(object$BB))
    pmat <- PI %*% object$BB %*% t(PI)
    if(!object$self.ties) diag(pmat) <- NA

    sr.mat <- as.matrix(object$SS) %*% t(as.matrix(object$RR))
    pmat <- pmat * sr.mat * object$hours
    return(pmat)
}


wsbm.metric <- function(graph,kk=2,total=1500,
                       thin=1,burn.in=500,verbose=0,...){

    ##    mode <- match.arg(mode)
    ##    if(mode == "mcmc"){
    wsbm.fit <- wsbm(total=total,net=graph,kk=kk,verbose=verbose,
                   thin=thin,burn.in=burn.in)
    ##    }
    return(wsbm.fit$pmat)
}



plot.net.wsbm <- function(wsbm.obj,
                         ord = order(wsbm.obj$mmb),...){
    nn <- nrow(wsbm.obj$net)
    image(1:nn,1:nn,wsbm.obj$net[ord[nn:1],ord],col=grey((50:1)/50),
          ylab="",xlab="",yaxt="n",xaxt="n",...)
}

plot.fit.wsbm <- function(wsbm.obj,
                          ord = order(wsbm.obj$mmb),...){
    nn <- nrow(wsbm.obj$net)
    image(1:nn,1:nn,wsbm.obj$pmat[ord[nn:1],ord],col=grey((50:1)/50),
          ylab="",xlab="",yaxt="n",xaxt="n",...)
}





##########################################################
##################  ROTATION FUNCTIONS  ##################
##########################################################


wsbm.log.like.net <- function(net,BB,mmb,SS,RR,PI,kk=max(mmb),
                              self.ties=TRUE,hours=1){
    ##browser()
    if(!self.ties){
        diag(net) <- NA
    }
    if(missing(PI)){
        if(missing(mmb)){
            stop("At least one of mmb or PI must be passed")
        }else{
            PI <- mmb.to.PI(mmb,kk=kk)
        }
    }
    summ.obj <- structure(list(mmb=mmb,BB=BB,SS=SS,RR=RR,hours=hours),
                          class="wsbm")
    summ.obj$self.ties <- self.ties

    pmat <- predict(summ.obj)
    return(sum(log(pmat)*net,na.rm=TRUE) - sum(pmat[!is.na(net)]))
}

wsbm.ll.pmat <- function(net,pmat,self.ties=TRUE){
    if(!self.ties){
        diag(net) <- NA
    }
    return(sum(log(pmat)*net,na.rm=TRUE) - sum(pmat[!is.na(net)]))
}

wsbm.marginal.ll.single <- function(net,BB,theta=rep(1/ncol(BB),ncol(BB))){
    nn <- ncol(net)
    kk <- ncol(BB)

    PI.t <- rmultinom(nn,1,theta)
    PP <- t(PI.t) %*% BB %*% PI.t
    diag(PP) <- NA
    ll <- sum(log(PP * net + (1-PP) * (1 - net)),na.rm=TRUE)
    return(ll)
}

wsbm.marginal.log.like.net <- function(net,BB,theta=rep(1/ncol(BB),ncol(BB)),
                                      iter.max=1e4){
    ll.vec <- replicate(iter.max,marginal.ll.single(net,BB,theta))
    ll.max <- max(ll.vec)
    ll.vec <- ll.vec - ll.max
    return(log(mean(exp(ll.vec))) + ll.max)
}



wsbm.load.init.vals <- function(init.vals,nn,kk){

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
