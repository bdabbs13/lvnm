#####  gamma.prior.R
#####
#####  Functions for sampling from the conjugate prior to the gamma
#####  distribution.
######################################################################



#' The Gamma Conjugate Prior Distribution
#'
#' Density and random generation for the conjugate prior to the Gamma
#'   distribution.
#'
#' @param x vector of length 2 at which to evaluate density
#' @param p parameter representing product of samples from a Gamma distribution
#' @param q parameter representing sum of samples from a Gamma distribution
#' @param r parameter representing samples used to calculate product
#' @param s parameter representing samples used to calculate sum
#' @param log if TRUE, probabilities P are returned as log(P)
#'
#' @details Returns the density for points drawn from the bivariate conjugate
#'   gamma prior.
#'
#' @export
dgamma.prior <- function(x,p=.99,q=1.01,r=1,s=r,log=FALSE){
    alpha <- x[1]; beta <- x[2]
    l.num <- (alpha-1)*log(p) - beta*q + alpha*s*log(beta)
    l.den <- r*lgamma(alpha)

    if(log){
        return(l.num - l.den)
    }else{
        return(exp(l.num - l.den))
    }
}

#' @rdname dgamma.prior
#' @param n number of samples to draw
#' @param ... additional parameters for conjugate Gamma prior sampler
#' @export
rgamma.prior <- function(n,p=0.99,q=1.01,r=1,s=r, ...){
    return(rgamma.prior.C(n,p,q,r,s,...))
}


rgamma.prior.C <- function(n, p=0.99, q=1.01, r=1, s=1,
                           mcmc.control=list(thin=100,burn.in=10000,
                                             alpha.init=1,beta.init=1,
                                             prop.sd=1.0)){

    thin <- mcmc.control$thin
    burn.in <- mcmc.control$burn.in
    alpha.init <- mcmc.control$alpha.init
    beta.init <- mcmc.control$beta.init
    prop.sd <- mcmc.control$prop.sd

    ##  Allocating Memory for C++
    alpha.vec <- beta.vec <- rep(0,n)
    ##    print(prop.sd)
    #' @useDynLib lvnm RGammaPriorNorm
    out <- .C(RGammaPriorNorm,
              as.integer(n),as.integer(thin),as.integer(burn.in),
              as.double(alpha.init),as.double(beta.init),
              as.double(p),as.double(q),as.double(r),as.double(s),
              as.double(prop.sd),as.double(alpha.vec),as.double(beta.vec))
    return(cbind(out[[11]],out[[12]]))
}

