#####  sbm.R
#####  Functions Implementing SBM and MMSBM Fitting in C:
#####    sbm - Fits SBM Model
#####    mmsbm.C - Fits MMSBM Model
#####    mmsbm - Fits MMSBM choosing intelligent starting values
#####
##############################################################
                                        #library(svd)

sbm.testing <- function(mat,kk=3){
    nn <- ncol(mat)
    out <- .C("sbmTesting",as.integer(nn),as.integer(kk),as.integer(mat))
}

rdirichlet <- function(n,alpha){
    l <- length(alpha)
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
    x/as.vector(sm)
}


#####  Function to generate data from SBM model
data.gen.sbm <- function(nn=50,
                         BB=array(c(.25,0.05,0.05,.25),c(2,2)),
                         pp=c(.5,.5),
                         mmb=NULL,PI=NULL){

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
    net <- array(rbinom(nn*nn,1,pmat),c(nn,nn))
    diag(net) <- NA
    return(list(net=net,PI=PI,BB=BB,pmat=pmat))
}




###  Wrapper for C Implementation of SBM
sbm <- function(total=1000,net,kk=3,verbose=0,init.vals=NULL,start=0,
                priors=list(aa=1,bb=1,eta=rep(1/kk,kk)),clean.out=TRUE,
                burn.in=0,thin=1,max.runs=200,spectral.start=FALSE,
                label.switch.mode = c("kl-loss","adhoc"),
                autoconverge=list(alpha=0.001,extend.max=20,shift.size=100),
                flatTable=NULL,ll.init=NULL,HH.init=NULL,multiImpute=FALSE){

    ##  Formatting Adjacency Matrix for C
    label.switch.mode <- match.arg(label.switch.mode)
    multi.int <- as.integer(ifelse(multiImpute,1,0))
    net.clean <- net
    diag(net.clean) <- -1
    net.clean[is.na(net)] <- -1
    nn <- nrow(net.clean)

    short.total <- total
    total <- short.total*thin + burn.in
    flatBB = double(short.total * kk*kk)
    flatMMB = integer(short.total * nn)
    ll.vec <- double(short.total)
    HH.vec <- double(short.total*kk*nn)

    if(is.null(init.vals) & spectral.start){
        if(verbose > 0) message("Initializing with Spectral Clustering...")
        init.vals <- sbm.spectral(net=net,kk=kk)
        if(verbose > 0) message("Initialization Complete.")

    }
    if(!is.null(init.vals)){
        ##  1234
        ##flatTable <- sbm.load.init.vals(init.vals,nn=nn,kk=kk)
        ll.init <- sbm.log.like.net(net.clean,init.vals$BB,PI=init.vals$PI)
        if(!is.null(init.vals$HH)){
            HH.init <- init.vals$HH
        }else{
            HH.init <- init.vals$PI
        }
    }

    start = 0
    ##  1234
    if(!is.null(flatTable)){
        if(is.null(ll.init)) stop("flatTables require ll.init vector")
        start = nrow(flatTable)
        total.start = start*(kk*(kk+nn))
        flatVec[1:total.start] = as.double(t(flatTable))
        ll.vec[1:start] <- ll.init
        HH.vec[1:(start *kk*nn)] <- as.double(t(HH.init))
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
    out <- .C("sbm",as.integer(total),as.integer(nn),as.integer(kk), ##3
              as.integer(t(net.clean)),as.double(c(priors$aa,priors$bb)), ## 5
              as.double(priors$eta), ## flatVec, 6
              flatBB, flatMMB, ## 8
              as.integer(burn.in), as.integer(thin), ##10
              as.integer(start),multi.int,ll.vec, ##13
              extend.max,shift.size,qq,HH.vec, ##17
              as.integer(verbose))

                                        #  browser()
    ##  Pulling the flat matrices from the C output
    ##full.mat <- matrix(out[[7]],ncol=kk*(kk+nn),byrow=TRUE)
    ##BB.flat <- t(as.matrix(full.mat[,1:kk^2]))
    BB.flat <- matrix(out[[7]],nrow=kk*kk)
    BB.flat <- apply(BB.flat,2,transpose.vector,nrow=kk)
    BB <- array(BB.flat,c(kk,kk,short.total))

    mmb <- array(out[[8]],c(nn,short.total))
    ##PI.flat <- t(as.matrix(full.mat[,-(1:kk^2)]))

    ll.vec <- as.vector(out[[13]])
##    HH.flat <- t(matrix(out[[17]],ncol=kk*nn,byrow=TRUE))
##    HH.flat <- apply(HH.flat,2,transpose.vector,nrow=nn)
##    HH <- array(HH.flat,c(nn,kk,short.total))
    HH.flat <- out[[17]]
    HH <- array(HH.flat,c(nn,kk,short.total))
    ##BB.flat.new <- apply(BB.flat.new,2,transpose.vector,nrow=kk)
    ##BB.new <- array(BB.flat.new,c(kk,kk,short.total))

    ##PI.flat <- apply(PI.flat,2,transpose.vector,nrow=nn)
    ##PI <- array(PI.flat,c(nn,kk,short.total))

    pmat.mat <- NULL  ## I might want to compute this later...
    diag(net.clean) <- -1


    sbm.out <- structure(list(BB=BB,mmb=mmb,##PI=PI,
                              net=net,logLik=ll.vec,
                              ##flat.mat=full.mat,
                              burn.in=burn.in,thin=thin,
                              pmat=pmat.mat,HH=HH),class="sbm")

    if(label.switch.mode == "kl-loss"){
        sbm.out <- switch.labels(sbm.out,max.runs=max.runs,verbose=verbose)
    }else{
        ##  Do nothing
    }
                                        #  browser()
    sbm.summ <- summary(sbm.out,burn.in=0,thin=1)
    sbm.summ$pmat <- predict(sbm.summ)
    sbm.summ$net <- sbm.out$net
    ## Calculating DIC
    sbm.summ$logLik <- with(sbm.summ,sbm.ll.pmat(net,pmat))
    sbm.summ$DIC <- 2*sbm.summ$logLik - 4 * mean(ll.vec)
    sbm.summ$burn.in <- burn.in; sbm.summ$thin <- thin

    if(clean.out){
        sbm.summ$chain <- list(logLik=sbm.out$logLik)
        sbm.summ$clean <- TRUE
    }else{
        sbm.summ$chain <- sbm.out
        sbm.summ$chain$net <- NULL
        sbm.summ$clean <- FALSE
    }

    return(sbm.summ)
}

PI.to.mmb <- function(PI){
    return(apply(PI,1,which.max))
}

mmb.to.PI <- function(mmb,kk=max(mmb)){
    nn <- length(mmb)
    PI <- array(0,c(nn,kk))
    PI[cbind(1:nn,mmb)] <- 1
    return(PI)
}

clean.net <- function(net){
    if(any(is.na(net))){
        net[is.na(net)] <- -1
    }
    return(net)
}

mean.net <- function(net){
    net[is.na(net)] <- mean(net,na.rm=TRUE)
    return(net)
}

fast.mle <- function(mmb,net,nn=length(mmb),kk=max(mmb)){

    net <- clean.net(net)
    BB <- array(0.0,c(kk,kk))
    out <- .C("getBlockMatMLE",as.integer(nn),as.integer(kk),
              as.integer(as.vector(net)),as.double(as.vector(BB)),
              as.integer(mmb))
    BB <- array(out[[4]],c(kk,kk))
    return(BB)
}

#####  OBSOLETE
get.BB.mle <- function(mmb,net,cols=1:ncol(net),kk=max(mmb)){
    nn <- length(mmb);
    PI <- array(0,c(nn,kk))
    PI[cbind(1:nn,mmb)] <- 1
    BB <- BB.tot <- array(0,c(kk,kk))

    for(ll in 1:kk){
        for(rr in 1:kk){
            mmb.2 <- PI[,ll,drop=FALSE] %*% t(PI[,rr,drop=FALSE])
            diag(mmb.2) <- NA
            mmb.2 <- mmb.2[,cols]
            BB.tot[ll,rr] <- max(sum(mmb.2,na.rm=TRUE),1)
            BB[ll,rr] <- sum(mmb.2*net,na.rm=TRUE)
        }
    }
    BB <- BB / BB.tot
    return(BB)

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


sbm.spectral <- function(net,kk=3,cols=1:ncol(net),fixed=FALSE,is.weighted=NULL,
                         svd.mode=c("default","trlan","propack",
                                    "lower","upper")){

    svd.mode <- match.arg(svd.mode)
    net <- net[,cols]
    nn <- ncol(net)
    if(is.null(is.weighted)){
        is.weighted <- any(net[!is.na(net)] < 0 | net[!is.na(net)] > 1)
        if(is.weighted){
            message("Network has values outside of 0 and 1.")
            message("Assuming network is weighted.")
        }
    }


    net.mean <- mean.net(net)
    ##  net.mean[is.na(net.mean)] <- mean(net.mean,na.rm=TRUE)

    if(svd.mode == "trlan"){
        eig <- svd::trlan.svd(net.mean,kk,opts=list(maxiter=1000))
        uu <- eig$u
    }else if(svd.mode == "propack"){
        eig <- svd::propack.svd(net.mean,kk)
        uu <- eig$u
    }else if(svd.mode == "default"){
        eig <- svd(net.mean,nu=kk,nv=kk)
        uu <- eig$u
    }else if(svd.mode == "lower"){
        net.sym <- symmetrize(net.mean,lower=TRUE)
        eig <- svd(net.sym,nu=kk,nv=kk)
        uu <- eig$u
    }else if(svd.mode == "upper"){
        net.sym <- symmetrize(net.mean,lower=FALSE)
        eig <- svd(net.sym,nu=kk,nv=kk)
        uu <- eig$u
    }

    mmb <- kmeans(uu,kk,iter.max=100,nstart=10)$cluster
    if(fixed) mmb <- sample(c(1:kk,sample(1:kk,size=nn-kk,replace=TRUE)))

    if(is.weighted){
        BB <- weighted.mle(mmb=mmb,net=net,kk=kk)
    }else{
        BB <- fast.mle(mmb=mmb,net=net,kk=kk)
    }
    spec.out <- list(BB=BB,mmb=mmb)

    spec.out$pmat <- predict.sbm(spec.out)
    spec.out$logLik <- sbm.ll.pmat(net,spec.out$pmat)
    spec.out$net <- net

    return(spec.out)
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


get.HH.init <- function(mmb,ee,kk=max(mmb)){
    HH <- mmb.to.PI(mmb)
    HH[HH == 1] <- 1 - ee
    HH[HH == 0] <- ee/(kk-1)
    return(HH)
}

get.random.params <- function(ee,kk){
    PI <- t(rmultinom(nn,1,rep(1,kk)/kk))
    array(runif(kk*kk,0,1),c(kk,kk))
    PI[PI==1] <- 1 - ee
    PI[PI==0] <- ee/(kk-1)
    return(list(BB=BB,PI=PI))
}


sbm.em <- function(net,kk=3,iter.max=1000,thresh=10e-4,verbose=0,
                   debug=FALSE,start=c("spectral","random","multi"),
                   b.min=0.001,ee=0.1,n.starts=100,calc.marginal.ll=FALSE,
                   mle.mode=c("fast","R"),
                   svd.mode=c("trlan","propack","default")){

    start <- match.arg(start)
    mle.mode <- match.arg(mle.mode)
    svd.mode <- match.arg(svd.mode)
    nn <- ncol(net)
    net.mean <- net
    net.mean[is.na(net.mean)] <- mean(net.mean,na.rm=TRUE)
    pi.prior <- rep(1/kk,kk)

    if(start == "multi"){
        em.out <- list(logLik.marginal = -Inf)
        for(ii in 1:n.starts){
            mmb <- sample(kk,nn,replace=TRUE)
            ##      obj <- get.random.params(ee,kk)
            ##      PI <- obj$PI
            ##BB <- obj$BB
            ##      mmb <- apply(PI,1,which.max)
            BB <- fast.mle(mmb,net)
            em.tmp <- sbmEM.C(net=net,kk=kk,BB=BB,mmb=mmb,pi.prior=pi.prior,
                              iter.max=iter.max,thresh=thresh,verbose=verbose,
                              debug=debug,calc.marginal.ll=FALSE,mle.mode=mle.mode)
            ##      if(em.out$logLik.marginal < em.tmp$logLik.marginal){
            if(em.out$logLik < em.tmp$logLik){
                em.out <- em.tmp
            }
        }
    }else{
        if(start == "spectral"){
            spec.fit <- sbm.spectral(net.mean,kk=kk,svd.mode=svd.mode)
            BB <- spec.fit$BB
            mmb <- spec.fit$mmb
            ##PI <- spec.fit$PI
            ##PI[PI==0] <- ee /(kk-1)
            ##PI[PI==1] <- 1 - ee
        }else if(start =="random"){
            mmb <- sample(kk,nn,replace=TRUE)
            BB <- fast.mle(mmb,net)
            ##obj <- get.random.params(ee,kk)
            ##BB <- obj$BB
            ##PI <- obj$PI
        }

        BB[BB==0] <- b.min

        em.out <- sbmEM.C(net=net,kk=kk,BB=BB,mmb=mmb,pi.prior=pi.prior,
                          iter.max=iter.max,thresh=thresh,verbose=verbose,
                          debug=debug,calc.marginal.ll=calc.marginal.ll)
    }
    em.out$pmat <- predict.sbm(em.out)
    ##em.out$mmb <- PI.to.mmb(em.out$PI)
    em.out$net <- net
    return(em.out)

}

sbmEM.C <- function(net,kk,BB,mmb,pi.prior,iter.max=1000,thresh=10e-4,
                    verbose=0,debug=FALSE,calc.marginal.ll=FALSE,ee=0.1,
                    mle.mode=c("fast","R")){
    ##  browser()
    mle.mode <- match.arg(mle.mode)
    nn <- ncol(net)
    kk <- ncol(BB)
    net.clean <- net;
    diag(net.clean) <- -1; net.clean[is.na(net.clean)] <- -1;
    ##mmb ##<- apply(PI,1,which.max)

    ##flatTable <- sbm.load.init.vals(list(BB=BB,PI=PI),nn=nn,kk=kk)
    flatBB <- as.double(t(BB))
    flatMMB <- as.integer(mmb)
    flatHH <- as.double(get.HH.init(mmb=mmb,ee=ee,kk=kk))
    ##flatVec <- as.double(t(flatTable))
    ll.init <- sbm.log.like.net(net.clean,BB,mmb=mmb,kk=kk)

    out <- .C("sbmEMout",as.integer(iter.max),as.integer(nn),as.integer(kk), #3
              as.integer(t(net.clean)),as.double(pi.prior), #5
              flatHH,flatBB,flatMMB,#8
              as.double(thresh), as.double(ll.init),as.integer(verbose))

    if(debug) browser()

    HH.flat <- out[[6]]
    BB.flat <- out[[7]]

    BB.post <- matrix(BB.flat,nrow=kk,byrow=FALSE) ##  This output is ignored...
    HH <- matrix(HH.flat,nrow=nn,byrow=FALSE)
    mmb <- PI.to.mmb(HH)

    if(mle.mode == "fast"){
        BB <- fast.mle(mmb=mmb,net=net,kk=kk)
    }else{
        BB <- get.BB.mle(mmb=mmb,net=net,kk=kk)
    }

    logLik.post <- out[[10]]
    logLik <- sbm.log.like.net(net=net,BB=BB,mmb=mmb,kk=kk)

    if(calc.marginal.ll){
        logLik.marginal <- sbm.marginal.log.like.net(net,BB,pi.prior)
    }else{
        logLik.marginal <- NULL
    }
    pi.prior <- out[[5]]

    return(list(BB=BB,mmb=mmb,PI=HH,pi.prior=pi.prior,BB.post=BB.post,
                logLik.post=logLik.post,logLik=logLik,
                logLik.marginal=logLik.marginal))

}




sbm.em.climb <- function(net,kk,BB,PI,pi.prior,iter.max=1000,thresh=10e-4,
                         verbose=FALSE,debug=FALSE,calc.marginal.ll){

    nn <- nrow(net)
    HH <- array(0,c(nn,kk))
    BB.tot <- array(0,dim(BB))
    if(debug) browser()
    BB.old <- BB
    for(ii in 1:iter.max){

        BB.old <- BB
###  E step
        HH <- 0 * HH
        for(rr in 1:kk){
            for(ss in 1:kk){
                ll.amat <- net * log(BB[rr,ss]) + (1 - net) * log(1 - BB[rr,ss])
                HH[,rr] = HH[,rr] + sapply(1:nn, function(x){sum(ll.amat[x,-x] * PI[-x,ss], na.rm=TRUE)})
                HH[,rr] = HH[,rr] + sapply(1:nn,function(x){sum(ll.amat[-x,x] * PI[-x,ss],na.rm=TRUE)})
            }
        }
        HH = HH - max(HH)
        HH.exp <- t(pi.prior * t(exp(HH)))
        PI = HH.exp / rowSums(HH.exp)

###  M step
        pi.prior = apply(PI, 2, sum) / nn

        for(ll in 1:kk){
            for(rr in 1:kk){
                mmb.2 <- PI[,ll,drop=FALSE] %*% t(PI[,rr,drop=FALSE])
                diag(mmb.2) <- NA
                mmb.2[is.na(net)] <- NA
                mmb.2 <- mmb.2
                BB.tot[ll,rr] <- max(sum(mmb.2,na.rm=TRUE),1)
                BB[ll,rr] <- sum(mmb.2*net,na.rm=TRUE)
            }
        }

        BB <- BB/BB.tot
        ##BB[BB==0] <- b.min

        delta <- sum((BB - BB.old)^2)
        if(delta < thresh){
            if(verbose) message("iter - ",ii," delta - ",delta)
            break;
        }
    }
    logLik <- sbm.log.like.net(net,BB,PI=PI)
    if(calc.marginal.ll){
        logLik.marginal <- sbm.marginal.log.like.net(net,BB,pi.prior)
    }else{
        logLik.marginal <- NULL
    }

    return(list(BB=BB,PI=PI,pi.prior=pi.prior,
                logLik=logLik,logLik.marginal=logLik.marginal))
}


#################################################################
###################  Modularity Maximization  ###################
#################################################################

sbm.modularity <- function(net){
    mod.obj <- mod.mat.init(net)
    ##mod.obj$mmb <- (mod.obj$ss + 1)/2 + 1
    if(all(mod.obj$ss == 1) | all(mod.obj$ss == -1)){
        mmb <- rep(1,ncol(net))
    }else{
        mmb.1 <- mod.mat.refine(mod.obj$DD,mod.obj$ss == 1,mod.obj$mm)
        mmb.2 <- mod.mat.refine(mod.obj$DD,mod.obj$ss == -1,mod.obj$mm)
        mmb <- mmb.1 + mmb.2
        mmb <- as.numeric(as.factor(mmb))
    }
    BB <- fast.mle(mmb,net=net)
    PI <- mmb.to.PI(mmb)

    sbm.obj <- list(BB=BB,mmb=mmb,PI=PI,net=net)
    sbm.obj$pmat <- predict.sbm(sbm.obj)
    sbm.obj$logLik <- sbm.log.like.net(net,BB,PI=PI)

    return(sbm.obj)
}

mod.mat.init <- function(net){
    mod.mat <- get.modularity.matrix(net)
    DD <- mod.mat + t(mod.mat)
    mm <- sum(net,na.rm=TRUE)

    mod.obj <- mod.mat.spec(DD,mm)
    mod.obj <- mod.fine.tune(mod.obj)
    ##  mod.obj$net <- net
    return(mod.obj)
}

mod.mat.refine <- function(DD,gg,mm){
    DD.sub <- DD[gg,gg]
    mod.sub <- mod.mat.split(DD.sub,mm)
    if(mod.sub$QQ > sqrt(.Machine$double.eps)){
        gg.1 <- gg.2 <- gg
        gg.1[gg.1] <- (mod.sub$ss == 1)
        mmb.1 <- mod.mat.refine(DD,gg.1,mm)
        gg.2[gg.2] <- (mod.sub$ss == -1)
        mmb.2 <- mod.mat.refine(DD,gg.2,mm)
        mmb <- mmb.1 + mmb.2
        return(mmb)
    }else{
        mmb <- rep(0,length(gg))
        mmb[gg] <- runif(1)
        return(mmb)
    }
}

mod.mat.spec <- function(DD,mm){
    decomp <- eigen(DD,symmetric=TRUE)
    ss <- 2 * as.numeric(decomp$vector[,1] > 0) - 1
    QQ <- (t(ss) %*% DD %*% ss) / (4 * mm)
    return(list(DD=DD,ss=ss,QQ=QQ,mm=mm))
}

mod.fine.tune <- function(mod.obj){
    nn <- nrow(mod.obj$DD)

    ss.cur <- ss.best <- mod.obj$ss
    QQ.best <- QQ.cur <- mod.obj$QQ
    change.vec <- get.change.vec(mod.obj$DD,ss.cur,mod.obj$mm)
    maximizer <- which.max(change.vec)
    changed <- c(maximizer)
    ss.cur[maximizer] <- ss.cur[maximizer] * -1
    QQ.cur <- QQ.cur + max(change.vec)
    if(QQ.cur > QQ.best){
        ss.best <- ss.cur
        QQ.best <- QQ.cur
    }

    for(ii in 2:nn){
        change.vec <- get.change.vec(mod.obj$DD,ss.cur,mod.obj$mm)
        order.vec <- order(change.vec,decreasing=TRUE)
        maximizer <- order.vec[sapply(order.vec,FUN=function(x) return(!any(changed == x)))][1]
        change.max <- change.vec[maximizer]
        QQ.cur <- QQ.cur + change.max

        changed <- c(changed,maximizer)
        ss.cur[maximizer] <- ss.cur[maximizer] * -1
        if(QQ.cur > QQ.best){
            ss.best <- ss.cur
            QQ.best <- QQ.cur
        }
    }
    if(QQ.best - mod.obj$QQ > sqrt(.Machine$double.eps)){
                                        #print("recurse")
        mod.obj$ss <- ss.best
        mod.obj$QQ <- QQ.best
        return(mod.fine.tune(mod.obj))
    }else{
        return(mod.obj)
    }
}

get.modularity.matrix <- function(net){
    diag(net) <- 0
    indegree <- colSums(net)
    outdegree <- rowSums(net)
    mm <- sum(net)

    mod.mat <- net - (indegree %*% t(outdegree) / mm)
    return(mod.mat)
}

get.change.vec <- function(DD,ss,mm){
    diag(DD) <- 0
    change.vec <- -1 * (DD %*% ss) * ss / mm
    return(as.vector(change.vec))
}

mod.mat.split <- function(DD,mm){
    DD <- DD - diag(rowSums(DD))
    mod.obj <- mod.mat.spec(DD,mm)
    mod.obj <- mod.fine.tune(mod.obj)
    delta.QQ <- t(mod.obj$ss) %*% DD %*% mod.obj$ss / (4*mm)
    return(mod.obj)
}

get.modularity <- function(net,mmb){
    mod.mat <- get.modularity.matrix(net)
    DD <- mod.mat + t(mod.mat)
    mm <- sum(net,na.rm=TRUE)
    kk <- max(mmb)
    mod <- 0
    for(ii in 1:kk){
        mod <- mod + sum(DD[mmb==ii,mmb==ii])
    }
    return(mod / (2*mm))
}



#################################################################
######################  sbm class functions  ####################
#################################################################


summary.sbm <- function(object,burn.in=0,thin=1,...){
                                        #  browser()
    total <- dim(object$BB)[3]; kk <- dim(object$BB)[1]
    my.sub <- seq(burn.in + thin, total, by = thin)

    BB.hat <- apply(object$BB[,,my.sub],c(1,2),mean)
    ##  PI.mean <- apply(object$PI[,,my.sub],c(1,2),mean)
    PI.mean <- t(apply(object$mmb,1,tabulate,nbins=kk) / total)
    mmb <- apply(PI.mean,1,which.max)

    return(structure(list(mmb=mmb,PI.mean=PI.mean,BB=BB.hat),class="sbm"))

}

predict.sbm <- function(object,...){
    PI <- mmb.to.PI(object$mmb,kk=ncol(object$BB))
    pmat <- PI %*% object$BB %*% t(PI)
    diag(pmat) <- NA
    return(pmat)
}


sbm.metric <- function(graph,kk=2,mode=c("mcmc","em","spectral"),total=1500,
                       thin=1,burn.in=500,verbose=0,...){

    mode <- match.arg(mode)
    if(mode == "mcmc"){
        sbm.fit <- sbm(total=total,net=graph,kk=kk,verbose=verbose,
                       thin=thin,burn.in=burn.in)
    }else if(mode == "em"){
        sbm.fit <- sbm.em(net=graph,kk=kk,...)
    }else if(mode == "spectral"){
        sbm.fit <- sbm.spectral(net=graph,kk=kk,...)
    }
    return(sbm.fit$pmat)
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


plot.net.sbm <- function(sbm.obj,
                         ord = order(sbm.obj$mmb),...){
    nn <- nrow(sbm.obj$net)
    image(1:nn,1:nn,sbm.obj$net[ord[nn:1],ord],col=grey((50:1)/50),
          ylab="",xlab="",yaxt="n",xaxt="n",...)
}

plot.fit.sbm <- function(sbm.obj,
                         ord = order(sbm.obj$mmb),...){
    nn <- nrow(sbm.obj$net)
    image(1:nn,1:nn,sbm.obj$pmat[ord[nn:1],ord],col=grey((50:1)/50),
          ylab="",xlab="",yaxt="n",xaxt="n",...)
}





##########################################################
##################  ROTATION FUNCTIONS  ##################
##########################################################


SBM.ID.rotation <- function(ID.labels, label.count=max(ID.labels)) {
                                        #ID.labels=c(3,4,4,1,1,2,5); label.count=max(ID.labels)

                                        #first pin down the keepers.
    assigned <- replaced <- rep(NA, label.count)
    replacements <- rep(NA, length(ID.labels))

    for (ii in 1:label.count) {
        if (is.na(replacements[ii])) {
            assigned[ID.labels[ii]] <- ii
            replaced[ii] <- ID.labels[ii]
            replacements[ID.labels==ID.labels[ii]] <- ii
        }
    }
    while (sum(is.na(assigned))>0) {
                                        #node <- min(which(is.na(replacements)))
        vacancy <- min(which(is.na(assigned)))
        newlabel <- min(which(is.na(replaced)))

        assigned[vacancy] <- newlabel
        replaced[newlabel] <- vacancy

        replacements[ID.labels==vacancy] <- newlabel
    }

    return (assigned)
}

SBM.rotate.bvector <- function(sbm.b.vector, rotation) {
    bvt <- cbind(make.edge.list.selfies(length(rotation)), sbm.b.vector)
    bvt[,1] <- rotation[bvt[,1]]; bvt[,2] <- rotation[bvt[,2]];
    bvt[,1:2] <- t(apply(rbind(bvt[,1:2]), 1, sort))
    bvt <- rbind(bvt[order(bvt[,1], bvt[,2]),])
    bvt[,3]
}

SBM.rotate.block <- function (sbm.block, rotation) {
    sbm.block[rotation, rotation]
}

rotate <- function (){
    rotation <- SBM.ID.rotation(membership, n.groups)
    membership <<- rotation[membership]
    block.matrix <<- SBM.rotate.block(block.matrix, rotation)
}



##########################################################
############  Useful Functions for MMSBMs ################
##########################################################
fake.data.check.sbm <- function(trials=100,draws=250,burn.in=1000,thin=10,kk=2,
                                nn=100,aa=1,bb=1,eta=rep(1/kk,kk),
                                verbose=0,FUN=function(BB){return(norm(BB,"F"))}){

    priors <- list(aa=aa,bb=bb,eta=eta)
    truth.mat <- array(NA,c(kk,kk,trials))#rep(NA,trials)
    est.mat <- array(NA,c(kk,kk,trials))
    qq.vec <- rep(NA,trials)
    for(trial in 1:trials){
        if(verbose > 0) message("Trial ",trial)
        ##  Generate parameters from prior
        BB <- array(rbeta(kk^2,priors$aa,priors$bb),c(kk,kk))
        truth.mat[,,trial] <- BB
        ##  Generate Data
        data.sample <- data.gen.sbm(nn=nn,BB=BB,pp=priors$eta)
        data.fit <- sbm(total=draws,net=data.sample$net,kk=kk,
                        burn.in=burn.in,thin=thin,clean.out=FALSE,
                        priors=priors)
        est.mat[,,trial] <- data.fit$BB
        fun.vec <- apply(data.fit$chain$BB,3,FUN)
        qq.vec[trial] <- mean(FUN(BB) < fun.vec)
    }

    ks.out <- ks.test(qq.vec,"punif")
    return(list(p.value=ks.out$p.value,truth.mat=truth.mat,est.mat=est.mat,
                qq.vec=qq.vec,ks.out=ks.out))

}

##input: ID.labels for a number of nodes, chosen to be the "dominant" ones.
##output: the permutation to switch labels to a simpler convention.
postprocess.SBM.IDs <- function(ID.labels, label.count=max(ID.labels)) {
    ##ID.labels=c(3,4,4,1,1,2,5); label.count=max(ID.labels)

    ##first pin down the keepers.
    assigned <- replaced <- rep(NA, label.count)
    replacements <- rep(NA, length(ID.labels))

    for (ii in 1:label.count) {
        if (is.na(replacements[ii])) {
            assigned[ID.labels[ii]] <- ii
            replaced[ii] <- ID.labels[ii]
            replacements[ID.labels==ID.labels[ii]] <- ii
        }
    }
    while (sum(is.na(assigned))>0) {
                                        #node <- min(which(is.na(replacements)))
        vacancy <- min(which(is.na(assigned)))
        newlabel <- min(which(is.na(replaced)))

        assigned[vacancy] <- newlabel
        replaced[newlabel] <- vacancy

        replacements[ID.labels==vacancy] <- newlabel
    }
    ##check replacements manually.

    return (assigned)

}


sbm.log.like.net <- function(net,BB,mmb,PI,kk=max(mmb)){
    diag(net) <- -1
    if(missing(PI)){
        if(missing(mmb)){
            stop("At least one of mmb or PI must be passed")
        }else{
            PI <- mmb.to.PI(mmb,kk=kk)
        }
    }
    pmat <- PI %*% BB %*% t(PI)
    return(sum(log(pmat[net==1])) + sum(log(1 - pmat[net==0])))
}

sbm.ll.pmat <- function(net,pmat){
    diag(net) <- -1
    return(sum(log(pmat[net==1])) + sum(log(1 - pmat[net==0])))
}

marginal.ll.single <- function(net,BB,theta=rep(1/ncol(BB),ncol(BB))){
    nn <- ncol(net)
    kk <- ncol(BB)

    PI.t <- rmultinom(nn,1,theta)
    PP <- t(PI.t) %*% BB %*% PI.t
    diag(PP) <- NA
    ll <- sum(log(PP * net + (1-PP) * (1 - net)),na.rm=TRUE)
    return(ll)
}

sbm.marginal.log.like.net <- function(net,BB,theta=rep(1/ncol(BB),ncol(BB)),
                                      iter.max=1e4){
    ll.vec <- replicate(iter.max,marginal.ll.single(net,BB,theta))
    ll.max <- max(ll.vec)
    ll.vec <- ll.vec - ll.max
    return(log(mean(exp(ll.vec))) + ll.max)
}



sbm.load.init.vals <- function(init.vals,nn,kk){

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

transpose.vector <- function(vec,nrow){
    return(as.vector(matrix(vec,nrow=nrow,byrow=TRUE)))
}


get.rotation.mat <- function(HH,max.runs=100,verbose=1){
                                        #  browser()
    total <- dim(HH)[3]
    nn <- dim(HH)[1]
    kk <- dim(HH)[2]

    ##  Step 1
    QQ <- apply(HH,c(1,2),mean)
    if(any(QQ == 1)){
        for(ii in 1:nrow(QQ)){
            if(any(QQ[ii,] == 1)){
                qq.1 <- which(QQ[ii,] == 1)
                QQ[ii,qq.1] <- 1 - .Machine$double.eps
                QQ[ii,-qq.1] <- .Machine$double.eps/(ncol(QQ) - 1)
            }
        }
    }

    QQ.old <- array(0,c(nn,kk))
    count <- 0

    rot.mat <- array(rep(diag(kk),total),c(kk,kk,total))
    CC <- array(NA,c(kk,kk))
    diff.vec <- rep(NA,max.runs)
    while(sum(abs(QQ.old - QQ)) > 1e-10){  #diff.vec[count] > 1e-10)!identical(QQ,QQ.old)){
        count <- count + 1
        if(count > max.runs) break

        ##  Step 2
        CC.mat <- get.cost.mat(HH,log(QQ))

        for(iter in 1:total){
            lp.out <- lpSolve::lp.assign(CC.mat[,,iter])
            rot.mat[,,iter] <- rot.mat[,,iter] %*% t(lp.out$solution)
            HH[,,iter] <- HH[,,iter] %*% t(lp.out$solution)
        }

        QQ.old <- QQ
        QQ <- apply(HH,c(1,2),mean)

        diff.vec[count] <- sum(abs(QQ.old - QQ))
    }
    if(verbose > 0){
        if(count-1 == max.runs){
            message("Warning:  Label-switching failed to converge after ",max.runs," iterations")
        }else{
            if(verbose > 1){
                message("Label-switching converged after ",count," iterations")
            }
        }
    }

    return(list(rot.mat=rot.mat,HH=HH))
}

mat.mult.left <- function(mat.1,mat.2) return(mat.2 %*% mat.1)
get.cost.mat <- function(HH,QQ.log){
    kk <- dim(HH)[2]
    total <- dim(HH)[3]
    bar <- array(apply(HH,3,mat.mult.left,mat.2=t(QQ.log)),c(kk,kk,total))
    foo <- array(rep(apply(HH[,,]*log(HH[,,]),c(2,3),sum),each=kk),c(kk,kk,total))
    return(foo - bar)
}

rotate.mmb <- function(mmb,rot.mat){
    change.vec <- apply(rot.mat,1,which.max)
    kk <- ncol(rot.mat)
    mmb.new <- rep(NA,length(mmb))
    for(ii in 1:kk){
        mmb.new[mmb == ii] <- change.vec[ii]
    }
    return(mmb.new)
}

switch.labels <- function(sbm.obj,max.runs=100,verbose=1){
                                        #  browser()
    rot.out <- get.rotation.mat(sbm.obj$HH,max.runs=max.runs,verbose=verbose)
    sbm.obj$HH <- rot.out$HH

    total = dim(sbm.obj$HH)[3]
    for(ii in 1:total){
        sbm.obj$BB[,,ii] <- t(rot.out$rot.mat[,,ii]) %*% sbm.obj$BB[,,ii] %*% rot.out$rot.mat[,,ii]
        ##sbm.obj$PI[,,ii] <- sbm.obj$PI[,,ii] %*% rot.out$rot.mat[,,ii]
        sbm.obj$mmb[,ii] <- rotate.mmb(sbm.obj$mmb[,ii],rot.out$rot.mat[,,ii])
    }

    return(sbm.obj)
}

diag.plot.sbm <- function(sbm.obj){

}



                                        #mmsbm.gen <- function(mmsbm.obj){
                                        #    p.mat <- predict(mmsbm.obj)
                                        #    nn <- nrow(p.mat)
                                        #    net <- array(NA,c(nn,nn))
                                        #    net[!is.na(p.mat)] <- rbinom(nn*(nn-1),1,p.mat[!is.na(p.mat)])
                                        #    return(net)
                                        #}



#####  Function to read mmsbm datafiles from C++ output
                                        #read.mmsbm <- function(file){
                                        #
#####  Reading in header information
                                        #    line1 <- as.vector(read.table(file,nrows=1,skip=0)[1,])
                                        #    total = line1[1,1]; nn = line1[1,2]; kk = line1[1,3];
                                        #
#####  Reading in flat matrix
                                        #    full.mat <- read.table(file,nrows=total,skip=1)
                                        #    BB.flat <- as.matrix(full.mat[,1:kk^2])
                                        #    PI.flat <- as.matrix(full.mat[,-(1:kk^2)])
                                        #
#####  Expanding Matrix
                                        #    BB <- array(NA,c(kk,kk,total))
                                        #    PI <- array(NA,c(nn,kk,total))
                                        #    for(ii in 1:total){
                                        #        for(jj in 1:kk){
                                        #            BB[jj,,ii] <- BB.flat[ii,((jj-1)*kk + 1):(jj*kk)]
                                        #        }
                                        #        for(jj in 1:nn){
                                        #            PI[jj,,ii] <- PI.flat[ii,((jj-1)*kk + 1):(jj*kk)]
                                        #        }
                                        #    }
                                        #    return(list(BB=BB,PI=PI,flat.mat=full.mat))
                                        #
                                        #}
