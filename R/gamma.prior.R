#####  gamma.prior.R
#####
#####  Functions for sampling from the conjugate prior to the gamma
#####  distribution.
######################################################################

#' @export
rgamma.prior <- function(n=100,p=0.99,q=1.01,r=1,s=1,
                         thin=100,burn.in=10000,
                         alpha.init = 1, beta.init = 1, prop.sd=1.0,
                         method=c("norm","mv","alpha-beta")){

    method <- match.arg(method)
    ##  Allocating Memory for C++
    alpha.vec <- beta.vec <- rep(0,n)
    ##    print(prop.sd)
    if(method == "alpha-beta"){
        out <- .C("RGammaPrior",
                  as.integer(n),as.integer(thin),as.integer(burn.in),
                  as.double(alpha.init),as.double(beta.init),
                  as.double(p),as.double(q),as.double(r),as.double(s),
                  as.double(prop.sd),as.double(alpha.vec),as.double(beta.vec))
    }else if(method == "mv"){
        out <- .C("RGammaPriorMV",
                  as.integer(n),as.integer(thin),as.integer(burn.in),
                  as.double(alpha.init),as.double(beta.init),
                  as.double(p),as.double(q),as.double(r),as.double(s),
                  as.double(prop.sd),as.double(alpha.vec),as.double(beta.vec))
    }else if(method == "norm"){
        out <- .C("RGammaPriorNorm",
                  as.integer(n),as.integer(thin),as.integer(burn.in),
                  as.double(alpha.init),as.double(beta.init),
                  as.double(p),as.double(q),as.double(r),as.double(s),
                  as.double(prop.sd),as.double(alpha.vec),as.double(beta.vec))
    }
    return(cbind(out[[11]],out[[12]]))
}



#' @export
dgamma.prior <- function(x,p=1,q=2,r=2,s=2,log=FALSE){
    alpha <- x[1]; beta <- x[2]
    l.num <- (alpha-1)*log(p) - beta*q + alpha*s*log(beta)
    l.den <- r*lgamma(alpha)

    if(log){
        return(l.num - l.den)
    }else{
        return(exp(l.num - l.den))
    }
}
