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




###  Wrapper for C Implementation of WSBM
wsbm <- function(total=1000,net,kk=3,verbose=0,init.vals=NULL,start=0,
                 priors=list(aa=1,bb=1,eta=rep(1/kk,kk)),clean.out=TRUE,
                 burn.in=0,thin=1,max.runs=200,spectral.start=FALSE,
                 label.switch.mode = c("adhoc","kl-loss"),
                 autoconverge=list(alpha=0.001,extend.max=20,shift.size=100),
                 multiImpute=FALSE,self.ties=TRUE){

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

    start = 0
    if(!is.null(init.vals)){
        if(is.null(init.vals$BB) || is.null(init.vals$mmb)){
            warning("Invalid init.vals, no BB or mmb attributed...")
        }else{
            if(is.null(init.vals$logLik)){
                init.vals$logLik <- wsbm.log.like.net(net,BB=init.vals$BB,
                                                      mmb=init.vals$mmb,
                                                      self.ties=self.ties)
            }
            if(is.null(init.vals$HH)){
                init.vals$HH <- mmb.to.PI(mmb)
            }
            start = 1
            flatBB[1:(kk*kk)] <- init.vals$BB
            flatMMB[1:nn] <- init.vals$mmb
            flatSS[1:nn] <- init.vals$SS
            flatRR[1:nn] <- init.vals$RR
            ll.vec[1] <- init.vals$logLik
            HH.vec[1:(kk*nn)] <- init.vals$HH
        }
    }


    if(is.null(autoconverge)){
        extend.max <- 0; shift.size <- 100
        qq <- 3
    }else{
        extend.max <- as.integer(autoconverge$extend.max)
        shift.size <- as.integer(autoconverge$shift.size)
        qq <- as.double(qnorm(1 - autoconverge$alpha/2))
    }

    ##  Calling C Implementation of MCMC
    out <- .C("wsbm",as.integer(total),as.integer(nn),as.integer(kk), ##3
              as.integer(t(net.clean)),as.double(c(priors$aa,priors$bb)), ## 5
              as.double(priors$eta), ## flatVec, 6
              flatBB, flatMMB, flatSS, flatRR, ## 10
              as.integer(burn.in), as.integer(thin), ##12
              as.integer(start),multi.int,ll.vec, ##15
              extend.max,shift.size,qq,HH.vec, ##19
              as.integer(verbose))

                                        #  browser()
    ##  Pulling the flat matrices from the C output
    BB.flat <- matrix(out[[7]],nrow=kk*kk)
    BB.flat <- apply(BB.flat,2,transpose.vector,nrow=kk)
    BB <- array(BB.flat,c(kk,kk,short.total))

    mmb <- array(out[[8]],c(nn,short.total))
    SS <- array(out[[9]],c(nn,short.total))
    RR <- array(out[[10]],c(nn,short.total))


    ll.vec <- as.vector(out[[15]])
    HH.flat <- out[[19]]
    HH <- array(HH.flat,c(nn,kk,short.total))

    pmat.mat <- NULL  ## I might want to compute this later...
    diag(net.clean) <- -1


    wsbm.out <- structure(list(BB=BB,mmb=mmb,
                               SS=SS,RR=RR,
                               net=net,logLik=ll.vec,
                               burn.in=burn.in,thin=thin,
                               pmat=pmat.mat,HH=HH),class="wsbm")

    if(label.switch.mode == "kl-loss"){
        wsbm.out <- switch.labels(wsbm.out,max.runs=max.runs,verbose=verbose)
    }else{
        ##  Do nothing
    }
                                        #  browser()
    wsbm.summ <- summary(wsbm.out,burn.in=0,thin=1)
    wsbm.summ$pmat <- predict(wsbm.summ)
    wsbm.summ$net <- wsbm.out$net
    ## Calculating DIC
    wsbm.summ$logLik <- with(wsbm.summ,wsbm.ll.pmat(net,pmat,
                                                    self.ties=self.ties))
    wsbm.summ$DIC <- 2*wsbm.summ$logLik - 4 * mean(ll.vec)
    wsbm.summ$burn.in <- burn.in; wsbm.summ$thin <- thin

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


summary.wsbm <- function(object,burn.in=0,thin=1,...){
                                        #  browser()
    total <- dim(object$BB)[3]; kk <- dim(object$BB)[1]
    my.sub <- seq(burn.in + thin, total, by = thin)

    BB.hat <- apply(object$BB[,,my.sub],c(1,2),mean)
    ##  PI.mean <- apply(object$PI[,,my.sub],c(1,2),mean)
    PI.mean <- t(apply(object$mmb,1,tabulate,nbins=kk) / total)
    mmb <- apply(PI.mean,1,which.max)

    return(structure(list(mmb=mmb,PI.mean=PI.mean,BB=BB.hat),class="wsbm"))

}

predict.wsbm <- function(object,...){
    PI <- mmb.to.PI(object$mmb,kk=ncol(object$BB))
    pmat <- PI %*% object$BB %*% t(PI)
    diag(pmat) <- NA
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


wsbm.log.like.net <- function(net,BB,mmb,PI,kk=max(mmb),self.ties=TRUE){
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
    pmat <- PI %*% BB %*% t(PI)
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
