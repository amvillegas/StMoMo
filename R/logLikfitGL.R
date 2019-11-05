
#'Log-Likelihood of a fitGL object
#'
#'Returns the log-likelihood of a  Stochastic Mortality Model fitted with 
#'group regularised penalty.
#'
#'@param object an object of class \code{fitGL} representing a 
#'Stochastic Mortality Model fitted to some data using the group lasso function
#'\code{\link{glStMoMo}}.
#'@param ... other arguments.
#'
#'@return The log-likelihood of the fitted model.
#'@export
logLik.grpfitStMoMo <- function (object, ...) {
  structure(object$loglik, df = object$npar, nobs=object$nobs, class = "logLik")
}