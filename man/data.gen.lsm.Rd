\name{data.gen.lsm}
\alias{data.gen.lsm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{LSM Data Generation
%%  ~~function to do ... ~~
}
\description{Data generation function for a Latent Space Network Model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
data.gen.lsm(nn = 50, pp = c(0.5,0.5), alpha = -0.5, mean.list = list(c(-0.5, -0.5), c(1.25, 1.25)), sd.list = list(diag(0.25, 2), diag(0.25, 2)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nn}{Number of nodes in network to be generated
%%     ~~Describe \code{nn} here~~
}
\item{pp}{Percentage of nodes in each group}
  
  \item{alpha}{Intercept term for latent space model
%%     ~~Describe \code{alpha} here~~
}
  \item{mean.list}{List of means for multivariate normals
%%     ~~Describe \code{mean.list} here~~
}
  \item{sd.list}{List of variance matrices for multivariate normals
%%     ~~Describe \code{sd.list} here~~
}
}
\details{Generates data from a 2-Dimensional Latent Space Model.  this
  is done by creating draws from multiple 2-Dimensional normal
  distributions for the latent space positions.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{YY}{Adjaceny matrix for simulated network}
  \item{ZZ}{Matrix of Latent Space Positions}
  %    return(list(YY = YY.gen, ZZ = ZZ.gen))
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

\seealso{\code{\link{lsm}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Mean and Variance for first Normal Distribution
mu.1 <- c(0,0)
sigma.1 <- array(c(1,.25,.25,1),c(2,2))
## Mean and Variance for second Normal Distribution
mu.2 <- c(2,1)
sigma.2 <- diag(.5,2)

##  Combining into lists
mu.list <- list(mu.1,mu.2)
sigma.list <- list(sigma.1,sigma.2)

dat.lsm <- data.gen.lsm(nn=100,mean.list=mu.list,sd.list=sigma.list)

## Adjacency Matrix
dat.lsm$YY

## Latent Space Vectors
dat.lsm$ZZ
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{networks}

