\name{data.gen.sbm}
\alias{data.gen.sbm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{SBM Data Generations
%%  ~~function to do ... ~~
}
\description{Data generation function for a single membership Stochastic
  Block Model network.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
data.gen.sbm(nn = 50, BB = array(c(0.25, 0.05, 0.05, 0.25), c(2, 2)), pp = c(0.5, 0.5))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nn}{Number of nodes in network to be generated
%%     ~~Describe \code{nn} here~~
}
  \item{BB}{Matrix of tie probabilities between blocks.  Must be a
    square matrix.
%%     ~~Describe \code{BB} here~~
}
  \item{pp}{Proportion of memberships in each block.
%%     ~~Describe \code{pp} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{Returns a list containing four objects:
  \item{YY}{Adjacency matrix for simulated network}
  \item{PI}{Membership Vectors}
  \item{BB}{matrix of tie probabilities between each block}
  \item{p.mat}{Matrix of tie probabilities for each pair of nodes}
  %    return(list(YY = YY, PI = PI, BB = BB, p.mat = p.mat))

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

\seealso{\code{\link{sbm}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##  Matrix of Block Tie Probabilities
BB <- array(c(.3,.1,.1,.3),c(2,2))
pp <- c(0.5,0.5)  ##  Equal Block memberships

##  Generating Data:
dat.sbm <- data.gen.sbm(nn=100,BB=BB,pp=pp)

##  Adjacency Matrix
dat.sbm$YY

##  Membership Vectors
dat.sbm$PI
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{networks}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
