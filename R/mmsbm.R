#####  cid-C.R
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

#####  Function to generate data from MMSBM model
data.gen.sbm <- function(nn=50,
                         BB=array(c(.25,0.05,0.05,.25),c(2,2)),
                         pp=c(.5,.5)){

  kk <- length(pp)
  ll <- round(pp*nn)
#  total <- 1
#  for(ii in 1:(kk-1)){
#    block.vec[total:(total+ll[ii])] <- ii
#    total = total + ll[ii]
#  }

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




                                        #  Wrapper for C Implementation of SBM
sbm <- function(total=10,YY,kk=3,verbose=FALSE,init.vals=NULL,start=0,
                priors=list(aa=1,bb=1,eta=rep(1/kk,kk)),clean.out=TRUE,
                burn.in=0,thin=1,
                flatTable=NULL,multiImpute=TRUE){

  ##  Formatting Adjacency Matrix for C
  multi.int <- as.integer(ifelse(multiImpute,1,0))
  YY.clean <- YY
  diag(YY.clean) <- 0
  YY.clean[is.na(YY)] <- -1
  nn <- nrow(YY.clean)

  short.total <- total
  total <- short.total*thin + burn.in
  flatVec = double(short.total*(kk*(kk+nn)))
  start = 0
  if(!is.null(flatTable)){
    start = nrow(flatTable)
    total.start = start*(kk*(kk+nn))
    flatVec[1:total.start] = as.double(t(flatTable))
  }
  ##  Calling C Implementation of MCMC
  out <- .C("sbm",as.integer(total),as.integer(nn),as.integer(kk),
            as.integer(YY.clean),as.double(c(priors$aa,priors$bb)),
            as.double(priors$eta), flatVec,
            as.integer(burn.in), as.integer(thin),
            as.integer(start),multi.int)

  ##  Pulling the flat matrices from the C output
  full.mat <- matrix(out[[7]],ncol=kk*(kk+nn),byrow=TRUE)

  BB.flat <- as.matrix(full.mat[,1:kk^2])
  PI.flat <- as.matrix(full.mat[,-(1:kk^2)])

  ##  Formatting data for output
  BB <- array(NA,c(kk,kk,short.total))
  PI <- array(NA,c(nn,kk,short.total))
  pmat.mat <- array(NA,c(nn,nn,short.total))
  ll.vec <- rep(NA,short.total)
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

##### Rotating Blocks
#    old.memb <- apply(PI[,,ii],1,which.max)
#    rotation <- SBM.ID.rotation(old.memb,kk)
#    new.memb <- rotation[old.memb]
#    for(jj in 1:nn){
#      PI[jj,,ii] <- 0
#      PI[jj,new.memb[jj],ii] <- 1
#    }
#    BB[,,ii] <- BB[rotation,rotation,ii]

    ##  Log-Likelihood Vector

    ll.vec[ii] <- mmsbm.log.like.YY(YY.clean,BB[,,ii],PI[,,ii])
    pmat.mat[,,ii] <- predict.mmsbm(list(PI=PI[,,ii],BB=BB[,,ii]))
  }

  sbm.out <- structure(list(BB=BB,PI=PI,YY=YY,
                            logLik=ll.vec,flat.mat=full.mat,
                            burn.in=burn.in,thin=thin,
                            pmat=pmat.mat),class="sbm")
  sbm.summ <- summary(sbm.out,burn.in=0,thin=1)
  sbm.summ$pmat <- predict(sbm.summ)
  sbm.summ$YY <- sbm.out$YY
  ## Calculating DIC
  sbm.summ$logLik <- with(sbm.summ,mmsbm.log.like.YY(YY,BB,PI))
  sbm.summ$DIC <- 2*sbm.summ$logLik - 4 * mean(ll.vec)
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


#####  Wrapper function implementing sbm starts near the mode
mmsbm <- function(YY,total=200,kk=3,mode=c("multistart","random"),
                  priors=list(aa=1,bb=1,alpha=rep(.25,kk)),
                  max.starts=40,start.length=50,steps=2,clean.out=TRUE,
                  burn.in=0,thin=1,
                  verbose=FALSE,multiImpute=TRUE){

  mode <- match.arg(mode)
  flat.mat <- NULL
  ##  Finding Initial Chain Near Mode
  if(mode == "multistart"){
    initial.obj <- sbm.mode.search(max.starts=max.starts,
                                   start.length=start.length,
                                   steps=steps,
                                   YY=YY,kk=kk,verbose=verbose,
                                   multiImpute=multiImpute)
    flat.mat <- initial.obj$flat.mat
  }

  ##  Running Full Chain from Initial Chain
  if(verbose){
    print("Running Chain...")
  }

  mmsbm.out <- mmsbm.C(total=total,YY=YY,kk=kk,verbose=verbose,
                       burn.in=burn.in,thin=thin,
                       priors=priors,flatTable=flat.mat,
                       multiImpute=multiImpute)
  if(verbose) browser()
  if(verbose) print("Finished Chain")
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

summary.sbm <- function(object,burn.in=0,thin=1,...){
  total <- dim(object$PI)[3]; kk <- dim(object$BB)[1]
  my.sub <- seq(burn.in + thin, total, by = thin)

  BB.hat <- apply(object$BB[,,my.sub],c(1,2),mean)
  PI.mean <- apply(object$PI[,,my.sub],c(1,2),mean)
  PI.hat <- diag(kk)[apply(PI.mean,1,which.max),]

  return(structure(list(PI=PI.hat,BB=BB.hat),class="sbm"))

}

predict.mmsbm <- function(object,...){
  p.mat <- object$PI %*% object$BB %*% t(object$PI)
  diag(p.mat) <- NA
  return(p.mat)
}

predict.sbm <- predict.mmsbm


mmsbm.metric <- function(graph,kk=3,total=2000,
                         thin=1,burn.in=1000,verbose=FALSE){
  mmsbm.fit <- mmsbm(total=total,YY=graph,kk=kk,verbose=verbose,
                     thin=thin,burn.in=burn.in)
                                        #	mmsbm.summ <- summary(mmsbm.fit,thin=thin,burn.in=burn.in)
  return(mmsbm.fit$pmat)
}

sbm.metric <- function(graph,kk=2,total=1500,
                       thin=1,burn.in=500,verbose=FALSE){
  sbm.fit <- sbm(total=total,YY=graph,kk=kk,verbose=verbose,
                 thin=thin,burn.in=burn.in)
                                        #    sbm.summ <- summary(sbm.fit,thin=thin,burn.in=burn.in)
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


plot.net.sbm <- function(sbm.obj,ord = order(apply(sbm.obj$PI,1,which.max)),...){
###    ord <- order(apply(sbm.obj$PI,1,which.max))
  nn <- nrow(sbm.obj$YY)
###    image(1:nn,1:nn,sbm.obj$YY[ord,ord],col=grey((50:1)/50))
  image(1:nn,1:nn,sbm.obj$YY[ord[nn:1],ord],col=grey((50:1)/50),
        ylab="",xlab="",yaxt="n",xaxt="n",...)
}

plot.fit.sbm <- function(sbm.obj,ord = order(apply(sbm.obj$PI,1,which.max)),...){
###    ord <- order(apply(sbm.obj$PI,1,which.max))
  nn <- nrow(sbm.obj$YY)
###    image(1:nn,1:nn,sbm.obj$pmat[ord,ord],col=grey((50:1)/50))
  image(1:nn,1:nn,sbm.obj$pmat[ord[nn:1],ord],col=grey((50:1)/50),
        ylab="",xlab="",yaxt="n",xaxt="n",...)
}

plot.ll.sbm <- function(sbm.obj,add=FALSE){
  if(add){
    lines(sbm.obj$chain$logLik,type="l")
  }else{
    plot(sbm.obj$chain$logLik,type="l")
  }
}





##########################################################
##################  INTERNAL FUNCTIONS  ##################
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

                                        #input: ID.labels for a number of nodes, chosen to be the "dominant" ones.
                                        #output: the permutation to switch labels to a simpler convention.
postprocess.SBM.IDs <- function(ID.labels, label.count=max(ID.labels)) {
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
                                        #check replacements manually.

  return (assigned)

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
sbm.mode.search <- function(max.starts=40,start.length=50,YY,kk=2,verbose=FALSE,
                            priors=list(aa=1,bb=1,eta=rep(1/kk,kk)),steps=2,
                            multiImpute=TRUE){

  if(verbose) message("Trying Initial Values...")
  start.out <- NULL

  ##  Matrix holding likelihood values
  ll.mat <- array(-Inf,c(start.length,max.starts))

  max.found <- FALSE
  iter <- 0
  ##  Continue generating chains until a significant maximum is found
  while(!max.found & iter < max.starts){
    iter <- iter + 1
    if(verbose){
      message("Start ", iter)
    }

    ##  Generating SBM chain
    start.out[[iter]] <- sbm(total=start.length,YY=YY,kk=kk,
                             verbose=verbose,priors=priors,clean.out=FALSE,
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
  if(verbose){
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

mmsbm.C <- function(total=10,YY,kk=3,verbose=FALSE,
                    burn.in=0,thin=1,init.vals=NULL,
                    priors=list(aa=1,bb=1,alpha=rep(.25,kk)),
                    flatTable=NULL,multiImpute=TRUE){

  ##  Formatting Adjacency Matrix for C
  multi.int <- as.integer(ifelse(multiImpute,1,0))
  YY.clean <- YY
  diag(YY.clean) <- 0
  YY.clean[is.na(YY)] <- -1
  nn <- nrow(YY.clean)

  short.total <- total
  total <- short.total * thin + burn.in
  flatVec = double(short.total*(kk*(kk+nn)))
  start = 0
  if(!is.null(flatTable)){
    start = nrow(flatTable)
    total.start = start*(kk*(kk+nn))
    flatVec[1:total.start] = as.double(t(flatTable))
  }

  ##  Calling C Implementation of MCMC
  out <- .C("mmsbm",as.integer(total),as.integer(nn),as.integer(kk),
            as.integer(YY.clean),as.double(c(priors$aa,priors$bb)),
            as.double(priors$alpha), flatVec,
            as.integer(burn.in),as.integer(thin),
            as.integer(start),multi.int)

  ##  Pulling the flat matrices from the C output
  full.mat <- matrix(out[[7]],ncol=kk*(kk+nn),byrow=TRUE)

  BB.flat <- as.matrix(full.mat[,1:kk^2])
  PI.flat <- as.matrix(full.mat[,-(1:kk^2)])

  ##  Formatting data for output
  BB <- array(NA,c(kk,kk,short.total))
  PI <- array(NA,c(nn,kk,short.total))
  pmat.mat <- array(NA,c(nn,nn,short.total))
  ll.vec <- rep(NA,short.total)
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
    ll.vec[ii] <- mmsbm.log.like.YY(YY.clean,BB[,,ii],PI[,,ii])
    pmat.mat[,,ii] <- predict.mmsbm(list(PI=PI[,,ii],BB=BB[,,ii]))
  }

  return(structure(list(BB=BB,PI=PI,YY=YY,logLik=ll.vec,flat.mat=full.mat,pmat=pmat.mat),class="mmsbm"))
}


mmsbm.log.like.YY <- function(YY,BB,PI){
  diag(YY) <- 0
  p.mat <- PI %*% BB %*% t(PI)
  return(sum(log(p.mat[YY==1])) + sum(log(1 - p.mat[YY==0])))
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
