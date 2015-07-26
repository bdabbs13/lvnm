\name{summary.lsm}
\alias{summary.lsm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarizing Latent Space Network Model
%%  ~~function to do ... ~~
}
\description{'summary' method for clas 'lsm'
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{summary}{lsm}(object, burn.in = 0, thin = 1, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class 'lsm', usually a result of a call to lsm
%%     ~~Describe \code{lsm.obj} here~~
}
  \item{burn.in}{number of iterations to burn in the chain
%%     ~~Describe \code{burn.in} here~~
}
  \item{thin}{amount to thin the chain
%%     ~~Describe \code{thin} here~~
}
\item{...}{further arguments passed to or from different methods.}
}
\details{This method summarizes the MCMC chain from the output of lsm.
  The chain is first burned in and then thinned by the amounts given.
  mean values are then computed on the reduced chain.
%%  ~~ If necessary, more details than the description above ~~
}
\value{The result is a list with two items:
  \item{BB}{Block Matrix}
  \item{PI}{Membership Matrix}
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

\seealso{code{\link{lsm}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{networks}

