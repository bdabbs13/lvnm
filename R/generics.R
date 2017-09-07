#####  generics.R
#####  Declarations for new generic functions
#####
##############################################################



#' Network Plot
#'
#' Plots a Heatmap Representation of an Adjacency Matrix
#'
#' @param x network object
#' @param ...
#'
#' @details This is a generic function.
#'
#' @seealso \code{\link{network.plot.wsbm}},
#' \code{\link{network.plot.dynsbm}}
#'
#' @export
network.plot <- function(x, ...){
    UseMethod("network.plot",x)
}

#' Network Generation Function
#'
#' Generates Adjacency Matrix Corresponding to a Model Object
#'
#' @param x model object
#' @param ... additional arguments for specific methods
#'
#' @details This is a generic funciton.  See specific method definitions for
#' additional information.
#'
#' @seealso  \code{\link{net.gen.wsbm}},
#' \code{\link{net.gen.dynsbm}}
#' @export
net.gen <- function(x, ...) {
    UseMethod("net.gen",x)
}


#' Diagnostic Plotting Function
#'
#' Plots Diagnostics for Various Model Fits
#'
#' @param x model object
#' @param ... additional arguments for specific methods
#'
#' @details This is a generic funciton.  See specific method definitions for
#' additional information.
#'
#' @seealso \code{\link{diagnostic.plot.wsbm.mcmc}},
#' \code{\link{diagnostic.plot.dynsbm.mcmc}}
#' @export
diagnostic.plot <- function(x, ...) {
    UseMethod("diagnostic.plot",x)
}


#' Parameter Plotting Function
#'
#' Plots Parameter Summaries for Various Model Fits
#'
#' @param x model object
#' @param ... additional arguments for specific methods
#'
#' @details This is a generic funciton.  See specific method definitions for
#' additional information.
#'
#' @seealso \code{\link{param.plot.wsbm.mcmc}}
#' \code{\link{param.plot.dynsbm.mcmc}}
#' @export
param.plot <- function(x, ...) {
    UseMethod("param.plot",x)
}



#' Posterior Predictive Distribution Estimator
#'
#' Estimates Posterior Predictive Distribution using MCMC Chain Output
#'
#' @param x model object which contains output from an MCMC chain
#' @param ... additional arguments for specific methods
#'
#' @details This is a generic function.  See specific method definitions for
#' additional information.
#'
#' @seealso \code{\link{post.predict.wsbm.mcmc}}
#' \code{\link{post.predict.dynsbm.mcmc}}
#' @export
post.predict <- function(x, ...) {
    UseMethod("post.predict",x)
}


#' Posterior Predictive Distribution Estimator
#'
#' Estimates Posterior Predictive Distribution using MCMC Chain Output
#'
#' @param x model object which contains output from an MCMC chain
#' @param ... additional arguments for specific methods
#'
#' @details This is a generic function.  See specific method definitions for
#' additional information.
#'
#' @seealso \code{\link{param.post.predict.dynsbm.mcmc}}
#'
#' @export
param.post.predict <- function(x, ...) {
    UseMethod("param.post.predict",x)
}



#' Get Model Object for an MCMC Iteration
#'
#' Returns Model Object from MCMC Iteration
#'
#' @param x mdl.mcmc object
#' @param iter iteration to return
#' @param ...
#'
#' @details This is a generic function.  Given an object of type
#' @export
get.iter <- function(x, iter, ...){
    UseMethod("get.iter",x)
}


#' Calculate Deviance Information Criterion
#'
#' \code{DIC} returns the Deviance Information Criterion for a Model Object
#'
#' @param x object for which to calculate DIC
#' @param ... additional parameters
#'
#' @details When comparing models using DIC, the smaller the value, the better
#' the fit.
#'
#' This is a generic function: methods can be defined for it directly
#' or via the \code{\link{Summary}} group generic. For this to work properly,
#' the arguments \code{...} should be unnamed, and dispatch is on the first
#' argument.
#'
#' @return Returns the Deviance Information Criterion.
#'
#' @seealso \code{\link{AIC}} \code{\link{BIC}}
#' @export
DIC <- function(x, ...) {
    UseMethod("DIC",x)
}

