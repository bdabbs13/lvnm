\name{sbm}
\alias{sbm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Stochastic Block Model
%%  ~~function to do ... ~~
}
\description{This function fits a stochastic block model (SBM) to
  network data using an MCMC algorithm.  The function is able to handle missing data.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sbm(total = 10, YY, kk = 3, verbose = FALSE, init.vals = NULL, start = 0, priors = list(aa = 1, bb = 1, eta = rep(1/kk, kk)), flatTable = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{total}{number of iterations for MCMC algorithm
%%     ~~Describe \code{total} here~~
}
  \item{YY}{adjacency matrix for network
%%     ~~Describe \code{YY} here~~
}
  \item{kk}{number of blocks
%%     ~~Describe \code{kk} here~~
}
  \item{verbose}{setting this value to TRUE provides more feedback while running
%%     ~~Describe \code{verbose} here~~
}
  \item{init.vals}{obsolete
%%     ~~Describe \code{init.vals} here~~
}
  \item{start}{number of initial values in flatTable
%%     ~~Describe \code{start} here~~
}
  \item{priors}{a list of priors which control the MCMC algorithm
%%     ~~Describe \code{priors} here~~
}
  \item{flatTable}{an initial chain for the MCMC algorithm to continue
%%     ~~Describe \code{flatTable} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{The output is an \code{sbm} object which has a summary method.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{Beau Dabbs
%%  ~~who you are~~
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{mmsbm}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(latent50)

##  Running MCMC Algorithm
sbm.fit <- sbm(total=1500,YY=latent50,kk=2)
plot(sbm.fit$logLik)

##  Getting summary of chain.
sbm.summ <- summary(sbm.fit,burn.in=500,thin=10)
sbm.summ$BB

##  Predicting probability of edges
pred.mat <- predict(sbm.summ)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{networks}

