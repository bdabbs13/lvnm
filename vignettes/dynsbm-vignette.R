library(lvnm)
data(dc60)

plot(dc60,marginal=TRUE)
plot(dc60,marginal=FALSE)
network.plot(dc60)

##  Setting up priors
avg <- 1; size <- 1e-3; eps <- 0.1
pp <- (avg - eps)^size; qq <- (avg + eps)*size; rr <- ss <- size

avg <- 1; size <- 1e-2; eps <- 0.1
pp.med <- (avg - eps)^size; qq.med <- (avg + eps)*size; rr.med <- ss.med <- size

kk <- 3
flat.priors <- dynsbm.priors(eta=rep(1/kk,kk),
                             c(pp,qq,rr,ss),
                             c(pp.med,qq.med,rr.med,ss.med),
                             c(pp.med,qq.med,rr.med,ss.med))

dyn.test <- dynsbm.fit(dc60$net.mat,kk=3,
                       tmap=dc60$tmap,hours.vec=dc60$hours.vec,
                       mcmc.control=list(total=1000,burn.in=100,thin=10),
                       priors=flat.priors,
                       init.control=list(multistart=10),
                       verbose=2)
dyn.test

plot(dyn.test,marginal=TRUE)
plot(dyn.test,marginal=FALSE)

diagnostic.plot.dynsbm.mcmc(dyn.test,t=3,marginal=FALSE,type="b")
diagnostic.plot.dynsbm.mcmc(dyn.test,t=3,marginal=FALSE,type="s")
diagnostic.plot.dynsbm.mcmc(dyn.test,t=3,marginal=FALSE,type="r")
diagnostic.plot.dynsbm.mcmc(dyn.test,t=3,marginal=FALSE,type="l")

diagnostic.plot.dynsbm.mcmc(dyn.test,t=2,marginal=TRUE,type="b")
diagnostic.plot.dynsbm.mcmc(dyn.test,t=2,marginal=TRUE,type="s")
diagnostic.plot.dynsbm.mcmc(dyn.test,t=2,marginal=TRUE,type="r")
diagnostic.plot.dynsbm.mcmc(dyn.test,t=2,marginal=TRUE,type="l")


param.plot(dyn.test,tclass=2,type="b")
param.plot(dyn.test,tclass=2,type="s")
param.plot(dyn.test,tclass=2,type="r")

block.prior.plot.dynsbm.mcmc(dyn.test)




#####  Generating Example Network
## nn <- 60; kk <- 3; ee <- 2; TT <- 6; tmap <- c(1,2,1,2,1,2)
## mmb <- sort(rep(1:3,length.out=nn))
## noise <- 1; signal <- 3
## BB.prior <- array(noise,c(kk,kk,2,ee))
## for(tt in 1:ee) diag(BB.prior[,,1,tt]) <- signal
## SS.prior <- RR.prior <- array(10*noise,c(nn,2,ee))



## seed <- 123456
## set.seed(seed)
## dc60 <- dynsbm(nn,TT,tmap,hours.vec=c(1,1),
##                   BB.prior=BB.prior,mmb=mmb,
##                   SS.prior=SS.prior,RR.prior=RR.prior,
##                   normalize=TRUE,gen=TRUE)
## save(dc60,file="dc60.RData")
