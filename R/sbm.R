#####  sbm.R
#####  Functions Implementing SBM and MMSBM Fitting in C:
#####    sbm - Fits SBM Model
#####    mmsbm.C - Fits MMSBM Model
#####    mmsbm - Fits MMSBM choosing intelligent starting values
#####
##############################################################


rdirichlet <- function(n,alpha){
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  x/as.vector(sm)
}

#####  Function to generate data from SBM model
data.gen.sbm <- function(nn=50,
                         BB=array(c(.25,0.05,0.05,.25),c(2,2)),
                         pp=c(.5,.5)){

  kk <- length(pp)
  ll <- round(pp*nn)

  block.vec <- sort(sample(kk,nn,replace=TRUE,prob=pp))
  PI <- array(0,c(nn,kk))
  for(ii in 1:nn){
    PI[ii,block.vec[ii]] <- 1
  }

  p.mat <- PI %*% BB %*% t(PI)
  YY <- array(rbinom(nn*nn,1,p.mat),c(nn,nn))
  diag(YY) <- NA
  return(list(YY=YY,PI=PI,BB=BB,p.mat=p.mat))
}




###  Wrapper for C Implementation of SBM
sbm <- function(total=1000,YY,kk=3,verbose=0,init.vals=NULL,start=0,
                priors=list(aa=1,bb=1,eta=rep(1/kk,kk)),clean.out=TRUE,
                burn.in=0,thin=1,
                autoconverge=list(alpha=0.001,extend.max=20,shift.size=100),
                flatTable=NULL,ll.init=NULL,multiImpute=TRUE){

  ##  Formatting Adjacency Matrix for C
  multi.int <- as.integer(ifelse(multiImpute,1,0))
  YY.clean <- YY
  diag(YY.clean) <- -1
  YY.clean[is.na(YY)] <- -1
  nn <- nrow(YY.clean)

  short.total <- total
  total <- short.total*thin + burn.in
  flatVec = double(short.total*(kk*(kk+nn)))
  ll.vec <- double(short.total)

  if(!is.null(init.vals)){
    flatTable <- sbm.load.init.vals(init.vals)
    ll.init <- sbm.log.like.YY(YY.clean,init.vals$BB,init.vals$PI)
  }

  start = 0
  if(!is.null(flatTable)){
    if(is.null(ll.init)) stop("flatTables require ll.init vector")
    start = nrow(flatTable)
    total.start = start*(kk*(kk+nn))
    flatVec[1:total.start] = as.double(t(flatTable))
    ll.vec[1:start] <- ll.init
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
  out <- .C("sbm",as.integer(total),as.integer(nn),as.integer(kk),
            as.integer(t(YY.clean)),as.double(c(priors$aa,priors$bb)),
            as.double(priors$eta), flatVec,
            as.integer(burn.in), as.integer(thin),
            as.integer(start),multi.int,ll.vec,
            extend.max,shift.size,qq,as.integer(verbose))

  ##  Pulling the flat matrices from the C output
  full.mat <- matrix(out[[7]],ncol=kk*(kk+nn),byrow=TRUE)

  BB.flat <- as.matrix(full.mat[,1:kk^2])
  PI.flat <- as.matrix(full.mat[,-(1:kk^2)])
  ll.vec <- as.vector(out[[12]])

  ##  Formatting data for output
  BB <- array(NA,c(kk,kk,short.total))
  PI <- array(NA,c(nn,kk,short.total))
  pmat.mat <- array(NA,c(nn,nn,short.total))
  ll.vec.old <- rep(NA,short.total)
  diag(YY.clean) <- -1
  for(ii in 1:short.total){

    ##  3D PI Matrix
    for(jj in 1:nn){
      PI[jj,,ii] <- PI.flat[ii,((jj-1)*kk + 1):(jj*kk)]
    }
    ##  3D BB Matrix
    for(jj in 1:kk){
      BB[jj,,ii] <- BB.flat[ii,((jj-1)*kk + 1):(jj*kk)]
    }

    ##  Log-Likelihood Vector

    ll.vec.old[ii] <- sbm.log.like.YY(YY.clean,BB[,,ii],PI[,,ii])
    pmat.mat[,,ii] <- predict.mmsbm(list(PI=PI[,,ii],BB=BB[,,ii]))
  }

  sbm.out <- structure(list(BB=BB,PI=PI,YY=YY,logLik=ll.vec,
                            logLik.old=ll.vec.old,flat.mat=full.mat,
                            burn.in=burn.in,thin=thin,
                            pmat=pmat.mat),class="sbm")
  sbm.summ <- summary(sbm.out,burn.in=0,thin=1)
  sbm.summ$pmat <- predict(sbm.summ)
  sbm.summ$YY <- sbm.out$YY
  ## Calculating DIC
  sbm.summ$logLik <- with(sbm.summ,sbm.log.like.YY(YY,BB,PI))
  sbm.summ$DIC <- 2*sbm.summ$logLik - 4 * mean(ll.vec.old)
  sbm.summ$burn.in <- burn.in; sbm.summ$thin <- thin
  if(clean.out){
    sbm.summ$chain <- list(logLik=sbm.out$logLik)
    sbm.summ$clean <- TRUE
  }else{
    sbm.summ$chain <- sbm.out
    sbm.summ$chain$YY <- NULL
    sbm.summ$clean <- FALSE
  }

  return(sbm.summ)
}


sbm.spectral <- function(YY,kk=3,cols=1:ncol(YY),mode="short"){
  if(any(is.na(diag(YY)))){
    diag(YY) <- 0
  }
  nn <- ncol(YY)
  YY <- YY[,cols]
  eig <- svd(YY,kk,kk)
  mmb <- kmeans(eig$u,kk,iter.max=100,nstart=10)$cluster

  BB <- BB.tot <- array(0,c(kk,kk))

  if(mode == "short"){
    PI <- array(0,c(nn,kk))
    PI[cbind(1:nn,mmb)] <- 1
    for(ll in 1:kk){
      for(rr in 1:kk){
        mmb.2 <- PI[,ll,drop=FALSE] %*% t(PI[,rr,drop=FALSE])
        diag(mmb.2) <- NA
        mmb.2 <- mmb.2[,cols]
        BB.tot[ll,rr] <- sum(mmb.2,na.rm=TRUE)
        BB[ll,rr] <- sum(mmb.2*YY,na.rm=TRUE)
      }
    }
  }

  if(mode == "long"){
    for(ii in 1:nrow(YY)){
      for(jj in 1:ncol(YY)){
        if(ii != jj){
          BB.tot[mmb[ii],mmb[cols[jj]]] <- BB.tot[mmb[ii],mmb[cols[jj]]] + 1
          BB[mmb[ii],mmb[cols[jj]]] <- BB[mmb[ii],mmb[cols[jj]]] + YY[ii,jj]
        }
      }
    }
  }

  BB <- BB / BB.tot
  PI <- array(0,c(nn,kk))
  PI[cbind(1:nn,mmb)] <- 1
  logLik <- sbm.log.like.YY(YY,BB,PI)

  return(list(PI=PI,BB=BB,logLik=logLik))
}

get.random.params <- function(ee){
  BB <- array(runif(kk*kk,0,1),c(kk,kk))
  PI <- t(rmultinom(nn,1,rep(1,kk)/kk))
  PI[PI==1] <- 1 - ee
  PI[PI==0] <- ee/(kk-1)
  return(list(BB=BB,PI=PI))
}


sbm.em <- function(YY,kk=3,iter.max=1000,thresh=10e-4,verbose=FALSE,
                   debug=FALSE,start=c("spectral","random","multi"),
                   b.min=0.001,ee=0.1,n.starts=100,calc.marginal.ll=FALSE){
  start <- match.arg(start)
  nn <- ncol(YY)
  YY.na <- YY
  YY.na[is.na(YY.na)] <- mean(YY.na,na.rm=TRUE)
  pi.prior <- rep(1/kk,kk)

  if(start == "multi"){
    em.out <- list(logLik.marginal = -Inf)
    for(ii in 1:n.starts){
      obj <- get.random.params(ee)
      BB <- obj$BB
      PI <- obj$PI
      em.tmp <- sbm.em.climb(YY=YY,kk=kk,BB=BB,PI=PI,pi.prior=pi.prior,
                             iter.max=iter.max,thresh=thresh,verbose=verbose,
                             debug=debug,calc.marginal.ll=TRUE)
      if(em.out$logLik.marginal < em.tmp$logLik.marginal){
        em.out <- em.tmp
      }
    }

    return(em.out)
  }else{
    if(start == "spectral"){
      spec.fit <- sbm.spectral(YY.na,kk=kk)
      BB <- spec.fit$BB
      PI <- array(ee/(kk-1),c(nn,kk))
      PI[cbind(1:nn,spec.fit$mmb)] <- 1 - ee
    }else if(start =="random"){
      obj <- get.random.params(ee)
      BB <- obj$BB
      PI <- obj$PI
    }

    BB[BB==0] <- b.min
    em.out <- sbm.em.climb(YY=YY,kk=kk,BB=BB,PI=PI,pi.prior=pi.prior,
                           iter.max=iter.max,thresh=thresh,verbose=verbose,
                           debug=debug,calc.marginal.ll=calc.marginal.ll)
    return(em.out)
  }
}
sbm.em.climb <- function(YY,kk,BB,PI,pi.prior,iter.max=1000,thresh=10e-4,
                         verbose=FALSE,debug=FALSE,calc.marginal.ll){
  HH <- array(0,c(nn,kk))
  BB.tot <- BB*0
  if(debug) browser()
  for(ii in 1:iter.max){
    BB.old <- BB
###  E step
    HH <- 0 * HH
    for(rr in 1:kk){
      for(ss in 1:kk){
        ll.amat <- YY * log(BB[rr,ss]) + (1 - YY) * log(1 - BB[rr,ss])
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
        mmb.2[is.na(YY)] <- NA
        mmb.2 <- mmb.2
        BB.tot[ll,rr] <- sum(mmb.2,na.rm=TRUE)
        BB[ll,rr] <- sum(mmb.2*YY,na.rm=TRUE)
      }
    }
    BB <- BB/BB.tot
                                        #BB[BB==0] <- b.min

    delta <- sum((BB - BB.old)^2)
    if(delta < thresh){
      if(verbose) message("iter - ",ii," delta - ",delta)
      break;
    }
  }
  logLik <- sbm.log.like.YY(YY,BB,PI)
  if(calc.marginal.ll){
    logLik.marginal <- sbm.marginal.log.like.YY(YY,BB,pi.prior)
  }else{
    logLik.marginal <- NULL
  }

  return(list(BB=BB,PI=PI,pi.prior=pi.prior,
              logLik=logLik,logLik.marginal=logLik.marginal))
}


summary.sbm <- function(object,burn.in=0,thin=1,...){
  total <- dim(object$PI)[3]; kk <- dim(object$BB)[1]
  my.sub <- seq(burn.in + thin, total, by = thin)

  BB.hat <- apply(object$BB[,,my.sub],c(1,2),mean)
  PI.mean <- apply(object$PI[,,my.sub],c(1,2),mean)
  PI.hat <- diag(kk)[apply(PI.mean,1,which.max),]

  return(structure(list(PI=PI.hat,BB=BB.hat),class="sbm"))

}

predict.sbm <- function(object,...){
  p.mat <- object$PI %*% object$BB %*% t(object$PI)
  diag(p.mat) <- NA
  return(p.mat)
}


sbm.metric <- function(graph,kk=2,total=1500,
                       thin=1,burn.in=500,verbose=0){
  sbm.fit <- sbm(total=total,YY=graph,kk=kk,verbose=verbose,
                 thin=thin,burn.in=burn.in)
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
                         ord = order(apply(sbm.obj$PI,1,which.max)),...){
  nn <- nrow(sbm.obj$YY)
  image(1:nn,1:nn,sbm.obj$YY[ord[nn:1],ord],col=grey((50:1)/50),
        ylab="",xlab="",yaxt="n",xaxt="n",...)
}

plot.fit.sbm <- function(sbm.obj,
                         ord = order(apply(sbm.obj$PI,1,which.max)),...){
  nn <- nrow(sbm.obj$YY)
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
fake.data.check <- function(trials=10,draws=100,burn.in=1000,thin=50,
                            FUN=function(bb) abs(bb[1,1] - bb[2,2])){

  truth.vec <- rep(NA,trials)
  est.mat <- array(NA,c(draws,trials))
  BB <- array(NA,c(2,2))
  for(trial in 1:trials){
    ##  Generate parameters from prior
    for(ii in 1:2){
      for(jj in 1:2){
        BB[ii,jj] <- rbeta(1,1,4)
      }
    }


    ##  Generate Data
    data.sample <- data.gen.mmsbm(nn=100,BB=BB)
    total = burn.in + thin*draws
    data.fit <- mmsbm(total=total,YY=data.sample$YY,kk=2)

    truth.vec[trial] <- FUN(BB)
    my.seq <- seq(burn.in + thin , total, by = thin)
    est.mat[,trial] <- apply(data.fit$BB[,,my.seq],3,FUN)

  }

  return(list(truth=truth.vec,est=est.mat))

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


sbm.log.like.YY <- function(YY,BB,PI){
  diag(YY) <- -1
  p.mat <- PI %*% BB %*% t(PI)
  return(sum(log(p.mat[YY==1])) + sum(log(1 - p.mat[YY==0])))
}

marginal.ll.single <- function(YY,BB,theta=rep(1/ncol(BB),ncol(BB))){
  nn <- ncol(YY)
  kk <- ncol(BB)

  PI.t <- rmultinom(nn,1,theta)
  PP <- t(PI.t) %*% BB %*% PI.t
  diag(PP) <- NA
  ll <- sum(log(PP * YY + (1-PP) * (1 - YY)),na.rm=TRUE)
  return(ll)
}

sbm.marginal.log.like.YY <- function(YY,BB,theta=rep(1/ncol(BB),ncol(BB)),
                                     iter.max=1e4){
  ll.vec <- replicate(iter.max,marginal.ll.single(YY,BB,theta))
  ll.max <- max(ll.vec)
  ll.vec <- ll.vec - ll.max
  return(log(mean(exp(ll.vec))) + ll.max)
}



sbm.load.init.vals <- function(init.vals){

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







                                        #mmsbm.gen <- function(mmsbm.obj){
                                        #    p.mat <- predict(mmsbm.obj)
                                        #    nn <- nrow(p.mat)
                                        #    YY <- array(NA,c(nn,nn))
                                        #    YY[!is.na(p.mat)] <- rbinom(nn*(nn-1),1,p.mat[!is.na(p.mat)])
                                        #    return(YY)
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
