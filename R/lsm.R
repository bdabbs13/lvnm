#####  lsm-2-norm-3.R
#####
#####  Latent Space Network Models - Code for fitting these models
###################################################################

#library(truncnorm)
#library(ergm)

############################################
#####  EXTERNALLY AVAILABLE FUNCTIONS  #####
############################################

data.gen.lsm <- function(nn=50,pp = c(0.5,0.5),alpha=-0.5,
                         mean.list=list(c(-.5,-.5),c(1.25,1.25)),
                         sd.list=list(diag(.25,2),diag(.25,2)) ){

    nn.total <- nn
    nn <- nn*pp
    YY.gen <- array(0,c(nn.total,nn.total))

    ZZ.gen <- rmvnorm(nn[1],mean.list[[1]],sd.list[[1]])
    if(length(nn) > 1){
        for(ii in 2:length(nn)){
            ZZ.gen <- rbind(ZZ.gen,rmvnorm(nn[ii],mean.list[[ii]],sd.list[[ii]]))
        }
    }

    ZZ.gen <- t(ZZ.gen)

    dd <- rowSums(apply(ZZ.gen,1,FUN <- function(zz) diff.mat(zz,zz)))
    dd <- array(dd,dim(YY.gen))
    diag(dd) <- NA
    dd.vec <- dd[!is.na(dd)]

    YY.gen[!is.na(dd)] <- rbinom(nn.total*(nn.total-1),1,logit.inv(alpha - dd.vec))
    return(structure(list(YY=YY.gen,ZZ=t(ZZ.gen),alpha=alpha),
                     class="summary.lsm"))
}




lsm <- function(YY,total=100,kk=2,verbose=FALSE,
                tune.default = list(    MH=list(ZZ.sd=.5,alpha.sd=.1),
                    prior=list(ZZ.sd=10,alpha.sd=1.5),
                    init=list(ZZ.sd=5,alpha.sd=5)),
                burn.in=0,thin=1,clean.out=TRUE,
                flatTable=NULL,multiImpute=TRUE){

    multi.int <- as.integer(ifelse(multiImpute,1,0))
    YY.clean <- YY
    YY.clean[is.na(YY)] <- -1
    diag(YY.clean) <- 0
    nn <- nrow(YY)

    mhControl = with(tune.default,c(MH$alpha.sd,MH$ZZ.sd))
    alphaPrior = tune.default$prior$alpha.sd
    zzPrior = tune.default$prior$ZZ.sd
    initPrior = with(tune.default,c(init$alpha.sd,init$ZZ.sd))

    short.total <- total
    total <- short.total * thin + burn.in
    flatVec = double(short.total*(nn*kk+2))
    start = 0
    if(!is.null(flatTable)){
        start = nrow(flatTable)
        total.start = start*(nn*kk+1)
        flatVec[1:total.start] = as.double(t(flatTable))
    }

    out <- .C("lsm",as.integer(total),as.integer(nn),as.integer(kk),
              as.integer(YY.clean),
              as.double(mhControl),as.double(alphaPrior),as.double(zzPrior),
              as.double(initPrior),flatVec,
              as.integer(burn.in),as.integer(thin),as.integer(start),
              multi.int)

    full.mat <- matrix(out[[9]],ncol=nn*kk+2,byrow=TRUE)
    ll.vec <- as.vector(full.mat[,1])
    alpha.vec <- as.vector(full.mat[,2])
    ZZ.flat   <- as.matrix(full.mat[,-c(1,2)])

    ZZ <- array(NA,c(nn,kk,short.total))
    pmat.mat <- array(NA,c(nn,nn,short.total))

    diag(YY.clean) <- -1
    for(ii in 1:short.total){
        ZZ[,,ii] <- matrix(ZZ.flat[ii,],nrow=nn,byrow=TRUE)
        pmat.mat[,,ii] <- predict.lsm(list(ZZ=ZZ[,,ii],alpha=alpha.vec[ii]))
    }


    lsm.out <- structure(list(alpha=alpha.vec,ZZ=ZZ,YY=YY,
                              logLik=ll.vec,flat.mat=full.mat,
                              pmat=pmat.mat),
                         class="lsm")
    lsm.summ <- summary(lsm.out,burn.in=0,thin=1)
    lsm.summ$pmat <- predict(lsm.summ)
    lsm.summ$YY <- lsm.out$YY
    ##  Calculating DIC
    lsm.summ$logLik <- with(lsm.summ,lsm.logLik(YY,alpha,ZZ))
    lsm.summ$DIC <- 2*lsm.summ$logLik - 4*mean(ll.vec)
    lsm.summ$burn.in <- burn.in; lsm.summ$thin <- thin
    if(clean.out){
        lsm.summ$chain <- list(logLik=lsm.out$logLik)
        lsm.summ$clean <- TRUE
    }else{
        lsm.summ$chain <- lsm.out
        lsm.summ$chain$YY <- NULL
        lsm.summ$clean <- FALSE
    }

    return(lsm.summ)

}


summary.lsm <- function(object,burn.in=0,thin=1,...){
    nn <- nrow(object$YY)
    total <- length(object$alpha)

    my.seq <- seq(burn.in,total,by=thin)

    ZZ.means <- apply(object$ZZ,c(1,2),mean)
    alpha.mean <- mean(object$alpha[my.seq])

    return(structure(list(alpha=alpha.mean,ZZ=ZZ.means),
                     class="lsm"))
}


predict.lsm <- function(object,...){
    return(with(object,logit.inv(alpha - dist.mat(ZZ))))
}


lsm.metric <- function(graph,kk=2,total=100,
                       thin=10,burn.in=500,verbose=FALSE,
                       multiImpute=TRUE){

	lsm.fit <- lsm(total=total,graph,kk,verbose=verbose,
                       multiImpute=multiImpute,burn.in=burn.in,thin=thin)
#	lsm.summ <- summary(lsm.fit,thin=thin,burn.in=burn.in)
	return(lsm.fit$pmat)
}

lsm.logLik <- function(YY,alpha,ZZ){
    dd <- dist.mat(ZZ)
    YY.vec <- YY[!is.na(dd)]
    dd.vec <- dd[!is.na(dd)]
    ll <- 0
    for(ii in 1:length(YY.vec)){
        pp <- log.logit.inv(alpha - dd.vec[ii])
        ll <- ll + ifelse(YY.vec[ii] == 1,pp,1-pp)
    }
    return(ll)
    ##    return(sum(sapply(alpha - dd[!is.na(dd)],log.logit.inv)))
}




#####  UTILITY FUNCIONS  #####

logit.inv <- function(x) 1/(1+exp(-x))
log.logit.inv <- function(x) -1 * log(1+exp(-x))


diff.mat <- function(XX,YY){
    return(sapply(XX,FUN=	function(x,y){
        return((x-y)^2)
    },y=YY))
}

dist.mat <- function(ZZ){
    ZZ <- t(ZZ)
    dd <- rowSums(apply(ZZ,1,FUN <- function(zz) diff.mat(zz,zz)))
    dd <- array(dd,c(ncol(ZZ),ncol(ZZ)))
    diag(dd) <- NA
    return(dd)
}


rmvnorm <- function(n, mean = rep(0,nrow(sigma)),sigma = diag(length(mean))){
    if(!isSymmetric(sigma)){
        stop("sigma must be a symmetric matrix")
    }

    ev <- eigen(sigma,symmetric=TRUE)
    if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
        warning("sigma is numerically not positive definite")
    }
    retval <- ev$vectors %*% diag(sqrt(ev$values),length(ev$values)) %*% t(ev$vectors)
    retval <- matrix(rnorm(n*ncol(sigma)),nrow=n, byrow=FALSE) %*% retval
    retval <- sweep(retval,2,mean,"+")
    return(retval)
}


