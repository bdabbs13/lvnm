### data.R
### Documentation for datasets in the lvnm package
###
######################################################################


#' Simulated Commuter Network
#'
#' A dataset containing the model parameters and simulated network
#' from a weighted stochastic block model.
#'
#' @format An object of class "wsbm"
#'
"commuter30"


#' Simulated Dynamic Commuter Network
#'
#' A dataset containing the model parameters and simulated networks
#' from a dynamic weighted stochastic block model.
#'
#' @format An object of class "dynsbm":
#' \describe{
#'   \item{TT=10}{Number of networks in the dataset}
#'   \item{ee=2}{Number of equivalence classes}
#'   \item{nn=60}{Number of nodes}
#'   \item{mmb}{The network contains 3 communities of equal size}
#' }
"dc60"
