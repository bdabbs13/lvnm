\name{predict.sbm}
\alias{predict.sbm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Predict Method for Stochastic Block Models
  %%  ~~function to do ... ~~
}
\description{Returns edges probabilities given an sbm summary object.
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{predict}{sbm}(object, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object of class sbm
    %%     ~~Describe \code{sbm.summ} here~~
  }
  \item{...}{ further arguments passed to or from other methods.
  }
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{A matrix of edge probabilities.  The diagonal values are
  defaulted to NA.
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

\seealso{\code{\link{sbm}}, \code{\link{summary.sbm}}
  %% ~~objects to See Also as \code{\link{help}}, ~~~
}
%\examples{

%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{networks}

