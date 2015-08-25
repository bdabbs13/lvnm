#####  mmsbm.R
#####  Functions Implementing SBM and MMSBM Fitting in C:
#####    sbm - Fits SBM Model
#####    mmsbm.C - Fits MMSBM Model
#####    mmsbm - Fits MMSBM choosing intelligent starting values
#####
##############################################################


#####  Function to generate data from MMSBM model
data.gen.mmsbm <- function(nn=50,
                           BB=array(c(.25,0.05,0.05,.25),c(2,2)),
                           PI.alpha=c(.25,.25)){

  PI <- rdirichlet(nn,alpha=PI.alpha)
  p.mat <- PI %*% BB %*% t(PI)
  YY <- array(rbinom(nn*nn,1,p.mat),c(nn,nn))
  diag(YY) <- NA
  return(list(YY=YY,PI=PI,BB=BB,p.mat=p.mat))
}




#####  Wrapper function implementing sbm starts near the mode
mmsbm <- function(YY,total=200,kk=3,mode=c("multistart","random"),
                  priors=list(aa=1,bb=1,alpha=rep(.25,kk)),init.vals=NULL,
                  max.starts=40,start.length=50,steps=2,clean.out=TRUE,
                  burn.in=0,thin=1,flatTable=NULL,ll.init=NULL,
                  autoconverge=list(alpha=0.001,extend.max=10,shift.size=100),
                  verbose=0,multiImpute=TRUE){

  mode <- match.arg(mode)
  flatTable <- NULL; ll.init <- NULL

  ##  Finding Initial Chain Near Mode
  if(mode == "multistart"){
    initial.obj <- sbm.mode.search(max.starts=max.starts,
                                   start.length=start.length,
                                   steps=steps,
                                   YY=YY,kk=kk,verbose=verbose,
                                   multiImpute=multiImpute)
    flatTable <- initial.obj$flat.mat
    init.vals <- NULL  #NEED TO ADD ACTUAL FUNCTIONALITY HERE!!
  }

  ##  Running Full Chain from Initial Chain
  if(verbose>0){
    print("Running Chain...")
  }

  mmsbm.out <- mmsbm.C(total=total,YY=YY,kk=kk,verbose=verbose,
                       burn.in=burn.in,thin=thin,init.vals=init.vals,
                       flatTable=flatTable,ll.init=ll.init,
                       autoconverge=autoconverge,
                       priors=priors, multiImpute=multiImpute)
  if(verbose>1) browser()
  if(verbose>0) print("Finished Chain")
                                        #    mmsbm.out <- postprocess.mmsbm(mmsbm.out)
  mmsbm.summ <- summary(mmsbm.out,burn.in=0,thin=1)
  mmsbm.summ$pmat <- predict(mmsbm.summ)
  mmsbm.summ$YY <- mmsbm.out$YY
  ##  Calculating DIC
  mmsbm.summ$logLik <- with(mmsbm.summ,mmsbm.log.like.YY(YY,BB,PI))
  mmsbm.summ$DIC <- 2*mmsbm.summ$logLik - 4 * mean(mmsbm.out$logLik)
  mmsbm.summ$burn.in <- burn.in; mmsbm.summ$thin <- thin

  if(clean.out){
    mmsbm.summ$chain <- list(logLik=mmsbm.out$logLik)
    mmsbm.summ$clean <- TRUE
  }else{
    mmsbm.summ$chain <- mmsbm.out
    mmsbm.summ$chain$YY <- NULL
    mmsbm.summ$clean <- FALSE
  }

  return(mmsbm.summ)
}


summary.mmsbm <- function(object,burn.in=0,thin=1,...){

  total <- dim(object$PI)[3]
  my.sub <- seq(burn.in+thin,total,by=thin)

  PI.hat <- t(apply(object$PI[,,my.sub],1,rowMeans))
  BB.hat <- t(apply(object$BB[,,my.sub],1,rowMeans))

  return(structure(list(PI=PI.hat,BB=BB.hat),class="mmsbm"))
}

predict.mmsbm <- function(object,...){
  p.mat <- object$PI %*% object$BB %*% t(object$PI)
  diag(p.mat) <- NA
  return(p.mat)
}

mmsbm.metric <- function(graph,kk=3,total=2000,
                         thin=1,burn.in=1000,verbose=0){
  mmsbm.fit <- mmsbm(total=total,YY=graph,kk=kk,verbose=verbose,
                     thin=thin,burn.in=burn.in)
  return(mmsbm.fit$pmat)
}

plot.net.mmsbm <- function(mmsbm.obj,
                         ord = order(apply(mmsbm.obj$PI,1,which.max)),...){
  nn <- nrow(mmsbm.obj$YY)
  image(1:nn,1:nn,mmsbm.obj$YY[ord[nn:1],ord],col=grey((50:1)/50),
        ylab="",xlab="",yaxt="n",xaxt="n",...)
}

plot.fit.mmsbm <- function(mmsbm.obj,
                         ord = order(apply(mmsbm.obj$PI,1,which.max)),...){
  nn <- nrow(mmsbm.obj$YY)
  image(1:nn,1:nn,mmsbm.obj$pmat[ord[nn:1],ord],col=grey((50:1)/50),
        ylab="",xlab="",yaxt="n",xaxt="n",...)
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

postprocess.mmsbm <- function(mmsbm.obj,burn.in=0,thin=1){

  nn <- dim(mmsbm.obj$PI)[1]
  kk <- dim(mmsbm.obj$PI)[2]
  total <- dim(mmsbm.obj$PI)[3]

  my.seq <- seq(burn.in + thin,total,by=thin)

  PP.summ <- apply(mmsbm.obj$PI[,,my.seq],c(1,2),mean)
  IDs <- apply(PP.summ,1,which.max)

  rot.vec <- postprocess.SBM.IDs(IDs)
  rot.mat <- array(0,c(kk,kk))
  for(ii in 1:kk){
    rot.mat[ii,rot.vec[ii]] <- 1
  }
  rotate <- function(mat,rot.mat){
    solve(rot.mat) %*% mat %*% rot.mat
  }
  mult <- function(mat,rot.mat){
    mat %*% rot.mat
  }

  mmsbm.obj$BB <- array(apply(mmsbm.obj$BB,3,rotate,rot.mat=rot.mat),
                        c(kk,kk,total))

  mmsbm.obj$PI <- array(apply(mmsbm.obj$PI,3,mult,rot.mat=rot.mat),
                        c(nn,kk,total))

  return(mmsbm.obj)

}


#####  Uses a large number of relatively short chains to find an
#####  initial value in the mode.
#####  Note:  This mode should also be close to the mmsbm mode.
sbm.mode.search <- function(max.starts=40,start.length=50,YY,kk=2,verbose=0,
                            priors=list(aa=1,bb=1,eta=rep(1/kk,kk)),steps=2,
                            multiImpute=TRUE){

  if(verbose>0) message("Trying Initial Values...")
  start.out <- NULL

  ##  Matrix holding likelihood values
  ll.mat <- array(-Inf,c(start.length,max.starts))

  max.found <- FALSE
  iter <- 0
  ##  Continue generating chains until a significant maximum is found
  while(!max.found & iter < max.starts){
    iter <- iter + 1
    if(verbose>1){
      message("Start ", iter)
    }

    ##  Generating SBM chain
    start.out[[iter]] <- sbm(total=start.length,YY=YY,kk=kk,
                             verbose=verbose-1,priors=priors,clean.out=FALSE,
                             multiImpute=multiImpute)$chain
    ##  Updating Likelihood Matrix
    ll.mat[,iter] <- start.out[[iter]]$logLik

    ##  Checking for Maximum value
    if(iter > 9){
      max.val = max(ll.mat[start.length,1:iter])
      mean.val = mean(ll.mat[start.length,1:iter])
      sd.val = sd(ll.mat[start.length,1:iter])
      ##  Check if max is steps stand. dev. greater than mean
      max.found <- max.val  >  (mean.val + steps*sd.val)
    }

  }

  ##  Notify if Maximum was Found
  if(verbose>0){
    if(iter < max.starts){
      message("Mode found after ",iter," iterations.")
    }else{
      message("Mode not found after ",max.starts,
              " iterations.  Continuing with best result")
    }
  }
  ##  Returning best start
  best.start <- start.out[[which.max(ll.mat[start.length,1:iter])]]

  return(best.start)

}

mmsbm.C <- function(total=10,YY,kk=3,verbose=0,
                    burn.in=0,thin=1,init.vals=NULL,
                    priors=list(aa=1,bb=1,alpha=rep(.25,kk)),
                    autoconverge=list(alpha=0.001,extend.max=10,shift.size=100),
                    flatTable=NULL,ll.init=NULL,multiImpute=TRUE){

  ##  Formatting Adjacency Matrix for C
  multi.int <- as.integer(ifelse(multiImpute,1,0))
  YY.clean <- YY
  diag(YY.clean) <- 0
  YY.clean[is.na(YY)] <- -1
  nn <- nrow(YY.clean)

  short.total <- total
  total <- short.total * thin + burn.in
  flatVec = double(short.total*(kk*(kk+nn)))
  ll.vec <- double(short.total)

  if(!is.null(init.vals)){
    flatTable <- sbm.load.init.vals(init.vals)
    ll.init <- mmsbm.log.like.YY(YY.clean,init.vals$BB,init.vals$PI)
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
  out <- .C("mmsbm",as.integer(total),as.integer(nn),as.integer(kk),
            as.integer(YY.clean),as.double(c(priors$aa,priors$bb)),
            as.double(priors$alpha), flatVec,
            as.integer(burn.in),as.integer(thin),
            as.integer(start),multi.int,ll.vec,
            extend.max,shift.size,qq,as.integer(verbose))

  ##  Pulling the flat matrices from the C output
  full.mat <- matrix(out[[7]],ncol=kk*(kk+nn),byrow=TRUE)

  BB.flat <- as.matrix(full.mat[,1:kk^2])
  PI.flat <- as.matrix(full.mat[,-(1:kk^2)])

  ##  Formatting data for output
  BB <- array(NA,c(kk,kk,short.total))
  PI <- array(NA,c(nn,kk,short.total))
  pmat.mat <- array(NA,c(nn,nn,short.total))
  ll.vec.old <- rep(NA,short.total)
  ll.vec <- as.vector(out[[12]])
  diag(YY.clean) <- -1
  for(ii in 1:short.total){

    ##  3D BB Matrix
    for(jj in 1:kk){
      BB[jj,,ii] <- BB.flat[ii,((jj-1)*kk + 1):(jj*kk)]
    }
    ##  3D PI Matrix
    for(jj in 1:nn){
      PI[jj,,ii] <- PI.flat[ii,((jj-1)*kk + 1):(jj*kk)]
    }
    ##  Log-Likelihood Vector
    ll.vec.old[ii] <- mmsbm.log.like.YY(YY.clean,BB[,,ii],PI[,,ii])
    pmat.mat[,,ii] <- predict.mmsbm(list(PI=PI[,,ii],BB=BB[,,ii]))
  }

  return(structure(list(BB=BB,PI=PI,YY=YY,logLik=ll.vec,logLik.old=ll.vec.old,flat.mat=full.mat,pmat=pmat.mat),class="mmsbm"))
}


mmsbm.log.like.YY <- function(YY,BB,PI){
  diag(YY) <- -1
  p.mat <- PI %*% BB %*% t(PI)
  return(sum(log(p.mat[YY==1])) + sum(log(1 - p.mat[YY==0])))
}

