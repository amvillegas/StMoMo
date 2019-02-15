#' Remove Intercept from Coefficients
#' 
#' @param object  an object of class \code{"StMoMo"} defining the stochastic
#' mortality model, with only parametric age-period and age-cohort terms.
#' @param tempfitttedCoef estimated coefficients from \code{grpreg} using
#' \code{extractCoefficientsFromGnm}.
#' @param intercept estimated intercept from \code{grpreg}.
#' @param ages vector of ages corresponding to rows of \code{Dxt} and \code{Ext}.
#' @keywords internal
glRemoveIntercept <- function(object, tempfittedCoef, intercept, ages) {
  
  # Include intercept in ax if possible
  if (object$staticAgeFun == "TRUE") {
    tempfittedCoef$ax <- tempfittedCoef$ax + intercept
  } 
  
  # Include intercept in gc if possible
  else if ("1" %in% object$cohortAgeFun) {
    tempfittedCoef$gc <- tempfittedCoef$gc + intercept
  } 
  
  # Include intercept in kt if possible
  else if ("1" %in% object$periodAgeFun) {
    
    for (i in 1:object$N) {
      if (length(unique(tempfittedCoef$bx[,i])) == 1) {
        tempfittedCoef$kt[i,] <- tempfittedCoef$kt[i,] + intercept
      }
    }
  } 
  
  # Force model to include ax to include intercept
  else {
    ax <- rep(intercept, length(ages))
    names(ax) <- ages
    tempfittedCoef$ax <- ax
    warning("Forced inclusion of static age term, ax, to account for intercept")
  }
  
  return(tempfittedCoef)
  
}