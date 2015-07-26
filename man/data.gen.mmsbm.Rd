\name{data.gen.mmsbm}
\alias{data.gen.mmsbm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MMSBM Data Generation
%%  ~~function to do ... ~~
}
\description{Data generation function for a Mixed Membership Stochastic
  Block Model Network.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
data.gen.mmsbm(nn = 50, BB = array(c(0.25, 0.05, 0.05, 0.25), c(2, 2)), PI.alpha = c(0.25, 0.25))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nn}{Number of nodes in network to be generated
%%     ~~Describe \code{nn} here~~
}
  \item{BB}{Matrix of tie probabilities between blocks.  Must be a
    square matrix
%%     ~~Describe \code{BB} here~~
}
  \item{PI.alpha}{Parameter for dirichlet prior which will generated
    block membership vectors
%%     ~~Describe \code{PI.alpha} here~~
k
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{Returns a list containing four objects:
  \item{YY}{Adjacency matrix for simulated network}
  \item{PI}{Membership Vectors}
  \item{BB}{Matrix of tie probabilities between blocks}
  \item{p.mat}{Matrix of tie probabilities for each pair of nodes}
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
##  Block Tie Probability Matrix
BB <- array(c(.3,.1,.1,.3),c(2,2))
alpha <- c(.5,.5)  ## Dirichlet Prior Parameter

##  Generating Data
dat.mmsbm <- data.gen.mmsbm(nn=100,BB=BB,PI.alpha=alpha)

##  Adjacency Matrix
dat.mmsbm$YY

##  Membership Vectors
dat.mmsbm$PI

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{networks}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
