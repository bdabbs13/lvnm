#####  generics.R
#####  Declarations for new generic functions
#####
##############################################################



#' Network Plot
#'
#' Plots a Heatmap Representation of an Adjacency Matrix
#'
#' @param x network object
#' @param ... additional graphical parameters
#'
#' @details Plots a heatmap representation of the adjacency matrix for a
#' network object.
#'
#' @seealso \code{\link{wsbm}},
#' \code{\link{dynsbm}}
#'
#' @export
network.plot <- function(x, ...){
    UseMethod("network.plot",x)
}

#' Network Generation Function
#'
#' Generates Adjacency Matrix Corresponding to a Model Object
#'
#' @param object model object
#' @param ... additional arguments for specific methods
#'
#' @details This is a generic funciton.  See specific method definitions for
#' additional information.
#'
#' @seealso  \code{\link{net.gen.wsbm}},
#' \code{\link{net.gen.dynsbm}}
#' @export
net.gen <- function(object, ...) {
    UseMethod("net.gen",object)
}


#' Diagnostic Plotting Function
#'
#' Plots Diagnostics for Various Model Fits
#'
#' @param object model object
#' @param ... additional arguments for plotting methods
#'
#' @details Plots diagnostics for various parameter estimations for a given
#'   model object.
#'
#' @seealso \code{\link{wsbm.fit}}, \code{\link{dynsbm.fit}}
#' @export
diagnostic.plot <- function(object, ...) {
    UseMethod("diagnostic.plot",object)
}


#' Parameter Plotting Function
#'
#' Plots Parameter Summaries for Various Model Fits
#'
#' @param object model object
#' plots the Log-Likelihood for each step of the chain.
#' @param ... other graphical parameters to pass to the plotting functions
#'
#' @details Displays boxplots for various sets of parameters for the model
#'   object.
#'
#' @seealso \code{\link{wsbm.fit}}, \code{\link{dynsbm.fit}}
#' @export
param.plot <- function(object, ...) {
    UseMethod("param.plot",object)
}



#' Posterior Predictive Distribution Estimator
#'
#' Estimates Posterior Predictive Distribution using MCMC Chain Output
#'
#' @param object model object which contains output from an MCMC chain
#' @param ... additional arguments for specific methods
#'
#' @details This is a generic function.  See specific method definitions for
#' additional information.
#'
#' @seealso \code{\link{post.predict.wsbm.mcmc}}
#' \code{\link{post.predict.dynsbm.mcmc}}
#' @export
post.predict <- function(object, ...) {
    UseMethod("post.predict",object)
}


#' Posterior Predictive Distribution Estimator
#'
#' Estimates Posterior Predictive Distribution using MCMC Chain Output
#'
#' @param object model object which contains output from an MCMC chain
#' @param ... additional arguments for specific methods
#'
#' @details This is a generic function.  See specific method definitions for
#' additional information.
#'
#' @seealso \code{\link{param.post.predict.dynsbm.mcmc}}
#'
#' @export
param.post.predict <- function(object, ...) {
    UseMethod("param.post.predict",object)
}



#' Get Model Object for an MCMC Iteration
#'
#' Returns Model Object from MCMC Iteration
#'
#' @param object mdl.mcmc object
#' @param iter iteration to return
#' @param ... additional parameters
#'
#' @return Returns an object of the class for which the given MCMC object
#'   is designed to sample.
#' @details Assuming the MCMC sampler is iteratively drawing parameters from
#'   the joint distribution for a given model class, this function will
#'   return the sample from a given iteration of the correct class.  This
#'   model object can then be used as a starting point for extending a chain,
#'   or can be analyzed by any of the functions that allow you to interact
#'   with that model object.
#'
#' @seealso \code{\link{wsbm.fit}}, \code{\link{dynsbm.fit}}
#' @export
get.iter <- function(object, iter, ...){
    UseMethod("get.iter",object)
}


#' Calculate Deviance Information Criterion
#'
#' \code{DIC} returns the Deviance Information Criterion for a Model Object
#'
#' @param object object for which to calculate DIC
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
DIC <- function(object, ...) {
    UseMethod("DIC",object)
}




#' Calculate Bayesian Information Criterion
#'
#' \code{BIC} returns the Bayesian Information Criterion for a Model Object
#'
#' @param object object for which to calculate BIC
#' @param ... additional parameters
#'
#' @details When comparing models using BIC, the smaller the value, the better
#' the fit.
#'
#' This is a generic function: methods can be defined for it directly
#' or via the \code{\link{Summary}} group generic. For this to work properly,
#' the arguments \code{...} should be unnamed, and dispatch is on the first
#' argument.
#'
#' @return Returns the Bayesian Information Criterion.
#'
#' @seealso \code{\link{AIC}} \code{\link{DIC}}
## BIC <- function(object, ...) {
##     UseMethod("BIC",x)
## }

