library(lvnm)
data(commuter30)

out.wsbm <- wsbm.fit(net=commuter30$net,kk=3)
out.wsbm

old.par <- par(mfrow=c(1,2))
network.plot(commuter30$net,main="True Network")
plot(out.wsbm,main="WSBM Estimate")
par(old.par)

fit.pred <- predict(out.wsbm)
sqrt(mean((fit.pred - commuter30$net)^2))
sqrt(mean((fit.pred - predict(commuter30))^2))

diagnostic.plot(out.wsbm,param="b") # Block Effects
diagnostic.plot(out.wsbm,param="s") # Sender Effects
diagnostic.plot(out.wsbm,param="r") # Receiver Effects
diagnostic.plot(out.wsbm,param="l") # Log-Likelihood

param.plot(out.wsbm,param="b")
param.plot(out.wsbm,param="s")
param.plot(out.wsbm,param="r")


out.wsbm.2 <- wsbm.fit(net=commuter30$net,kk=2)
out.wsbm.4 <- wsbm.fit(net=commuter30$net,kk=4)
AIC(out.wsbm.2); AIC(out.wsbm); AIC(out.wsbm.4)
BIC(out.wsbm.2); BIC(out.wsbm); BIC(out.wsbm.4)



#####  Experimenting with Sender/Receiver Effects
set.seed(123456)
nn <- 60; kk <- 3; noise <- 1; signal <- 2.5
BB <- array(noise,c(kk,kk)); diag(BB) <- signal
mmb <- sort(rep(1:kk,length.out=nn))

## Generating a simulated network
net.2 <- wsbm(nn=nn,BB=BB,mmb=mmb,
              SS=c(3,5,10,rgamma(nn-3,10,10)),
              RR=rgamma(nn,100,100),gen=TRUE)
old.par <- par(mfrow=c(1,2))
plot(1:nn,net.2$SS,ylim=c(0,10)); plot(1:nn,net.2$RR,ylim=c(0,10))
par(old.par)

#####  Fitting model with weak priors on sender/receiver effects
fit.base <- wsbm.fit(net.2$net,kk=3,verbose=0,
                     prior=wsbm.priors(eta=rep(1/kk,kk),
                                       sender.alpha=1e1,sender.beta=1e1,
                                       receiver.alpha=1e1,receiver.beta=1e1),
                     mcmc.control=list(total=1000,burn.in=1000,thin=10),
                     init.control=list(multistart=100))
## Memberships are Correct
fit.base$mmb


#####  Fitting model with strong priors on sender/receiver effects
fit.strong <- wsbm.fit(net=net.2$net,kk=3,hours=1,self.ties=T,
                       prior=wsbm.priors(eta=rep(1/kk,kk),
                                         sender.alpha=1e6,sender.beta=1e6,
                                         receiver.alpha=1e6,receiver.beta=1e6),
                       mcmc.control=list(total=1000,burn.in=1000,thin=10),
                       init.control=list(multistart=100))
## Nodes 2 and 3 are assigned to a new block
## Blocks 2 and 3 are combined
fit.strong$mmb

#####  Adding a fourth block with strong priors on sender/receiver effects
fit.strong.4 <- wsbm.fit(net=net.2$net,kk=4,
                         prior=wsbm.priors(eta=rep(1/4,4),
                                           sender.alpha=1e6,sender.beta=1e6,
                                           receiver.alpha=1e6,receiver.beta=1e6),
                         mcmc.control=list(total=1000,burn.in=1000,thin=10),
                         init.control=list(multistart=100))
##  Nodes 2 and 3 are still assigned to a new block
fit.strong.4$mmb

old.par <- par(mfrow=c(3,2))
plot(fit.strong,zlim=c(0,35))
network.plot(net.2$net)

plot(fit.base,zlim=c(0,35))
network.plot(net.2$net)

plot(fit.strong.4,zlim=c(0,35))
network.plot(net.2$net)


#####  Comparing Model Fits
c(DIC(fit.strong),DIC(fit.base),DIC(fit.strong.4))
c(AIC(fit.strong),AIC(fit.base),AIC(fit.strong.4))
c(BIC(fit.strong),BIC(fit.base),BIC(fit.strong.4))




#####  Experimenting with No Sender/Receiver Effects
set.seed(123456)
nn <- 60; kk <- 3; noise <- 1; signal <- 2.5
BB <- array(noise,c(kk,kk)); diag(BB) <- signal
mmb <- sort(rep(1:kk,length.out=nn))

net.3 <- wsbm(nn=nn,BB=BB,mmb=mmb,
              SS=rep(1,nn),RR=rep(1,nn),gen=TRUE)

old.par <- par(mfrow=c(1,2))
plot(net.3$SS,1:nn,xlim=c(0,10)); plot(net.3$RR,1:nn,xlim=c(0,10))
par(old.par)

#####  Improper control for degree effects
fit3.strong <- wsbm.fit(net=net.3$net,kk=3,hours=1,self.ties=T,
                        prior=wsbm.priors(eta=rep(1/kk,kk),
                                          sender.alpha=1e6,sender.beta=1e6,
                                          receiver.alpha=1e6,receiver.beta=1e6),
                        mcmc.control=list(total=1000,burn.in=1000,thin=10),
                        init.control=list(multistart=100))
fit3.strong$mmb


fit3.base <- wsbm.fit(net.3$net,kk=3,verbose=0,
                      prior=wsbm.priors(eta=rep(1/kk,kk),
                                        sender.alpha=1e1,sender.beta=1e1,
                                        receiver.alpha=1e1,receiver.beta=1e1),
                      mcmc.control=list(total=1000,burn.in=1000,thin=10),
                      init.control=list(multistart=100))
fit3.base$mmb


fit3.strong.4 <- wsbm.fit(net=net.3$net,kk=4,
                          prior=wsbm.priors(eta=rep(1/4,4),
                                            sender.alpha=1e6,sender.beta=1e6,
                                            receiver.alpha=1e6,receiver.beta=1e6),
                          mcmc.control=list(total=1000,burn.in=1000,thin=10),
                          init.control=list(multistart=100))
fit3.strong.4$mmb

old.par <- par(mfrow=c(3,2))
plot(fit3.strong,zlim=c(0,35))
network.plot(net.3$net)

plot(fit3.base,zlim=c(0,35))
network.plot(net.3$net)

plot(fit3.strong.4,zlim=c(0,35))
network.plot(net.3$net)


#####  Comparing Model Fits
c(DIC(fit3.strong),DIC(fit3.base),DIC(fit3.strong.4))
c(AIC(fit3.strong),AIC(fit3.base),AIC(fit3.strong.4))
c(BIC(fit3.strong),BIC(fit3.base),BIC(fit3.strong.4))




#####  Generation of commuter30
## library(lvnm)
## nn <- 30; kk <- 3; noise <- 1; signal <- 3
## BB <- array(noise,c(kk,kk))
## BB[cbind(1:3,c(2,3,1))] <- signal

## diag(BB) <- signal
## mmb <- sort(rep(1:kk,length.out=nn))

## commuter30 <- data.gen.wsbm(nn=nn,BB=BB,mmb=mmb,
##                             SS=rgamma(nn,10,10),
##                             RR=rgamma(nn,10,10))
## save(commuter30,file="commuter30.RData")
