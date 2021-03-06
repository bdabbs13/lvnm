\name{mmsbm}
\alias{mmsbm}
\title{Mixed Membership Stochastic Block Model}
\description{Fit a Mixed Membership Stochastic Block Model.}
\usage{
mmsbm(YY, total=200, kk=3, mode=c("multistart","random"),verbose=FALSE,priors=list(aa=1,bb=1,alpha=rep(.25,kk)), max.starts=40,start.length=50,steps=2)
%\method{mmsbm.wrapper}{default}(x, y, ...)
}
\arguments{
  \item{YY}{ Adjacency matrix for the network to model }
  \item{total}{ number of steps for MCMC chian. }
  \item{kk}{ number of blocks for the MMSBM }
  \item{mode}{multistart uses multiple sbm starts to find an initial
    mode, random uses a truly random starting value.}
  \item{verbose}{TRUE causes additional output to be printed while runnin}
  \item{priors}{controls the priors for the MCMC algorithm}
  \item{max.starts}{ number of starts to try. }
  \item{start.length}{ length of initial starts. }
  \item{steps}{controls mode detection for initial start}
}
\value{
  An object of class \code{mmsbm}, basically a list including elements
  \item{BB}{ 3 dimensional array of tie matrix values }
  \item{PI}{ 3 dimensional array of block membership values }
}
\author{Beau Dabbs}
\examples{
data(latent50)

##  Running MCMC Chain
mmsbm.fit <- mmsbm(total=1500,max.starts=20,YY=latent50)
plot(mmsbm.fit$logLik)

##  Summarizing Results
mmsbm.summ <- summary(mmsbm.fit)
mmsbm.summ$BB

##  Creating Matrix of Edge Probabilities
pred.mat <- predict(mmsbm.summ)


}
\keyword{networks}