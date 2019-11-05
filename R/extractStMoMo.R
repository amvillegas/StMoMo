#' Extract a Stochastic Mortality Model fit from group regularised fitting
#' 
#' Extract a Stochastic Mortality Model fit from group regularised fitting.
#' This is useful to get the fit for a single lambda in the regularisation path
#' of a group regularised fit obtained with function \code{\link{grpfit}}. 
#' The output is a \code{"fitStMoMo"} object which can then be used with other 
#' StMoMo functions such as \code{\link{plot.fitStMoMo}}, 
#' \code{\link{forecast.fitStMoMo}} and \code{\link{simulate.fitStMoMo}}
#' 
#' @param object an object of class \code{"grpfitStMoMo"} with the regularisation
#' path of the fitted Stochastic Mortality Model
#' @param  k integer index indicating which position (lambda) in the regularisation 
#' path is being extracted
#' @param simplify optional logical value indicating if the model should be simplified
#' to only include model terms which are different from 0. 
#' 
#' @return An list of class \code{"fitStMoMo"} with a fitted Stochastic Mortality Model.
#' See \code{\link{fit.StMoMo}} for the details of the elements in the list.
#' 
#' @export
extractStMoMo <- function(object, k, simplify = TRUE){
  
  if (class(object) != "grpfitStMoMo") {
    stop("Argument object needs to be of class grpfitStMoMo.")
  }
  
  if (k < 1 || k > length(object$lambda)) {
    stop(paste0("The index k=", k, "is out of bounds. Maximum number is", 
                max(object$lambda)))
  }
  
  out <- object$beta[[k]]
  model <- object$model
  textFormula <- model$textFormula
  if (simplify == TRUE){
    model <- object$model
    
    #Simplify alpha
    staticAgeFun <- model$staticAgeFun
    ax <- NULL
    if(staticAgeFun){
      if (sum(abs(out$ax)) > 0) {
        ax <- out$ax
        names(ax) <- object$ages
      } else {
        textFormula <- sub(" [+] a[[]x[]]", "", textFormula)
        textFormula <- sub("a[[]x[]]", "", textFormula)
        staticAgeFun <- FALSE
      }
    }
    out$ax <- ax
    
    #simplify age-period terms 
    bx <- NULL
    kt <- NULL
    periodAgeFun <- NULL
    if (model$N>0) {
      nonzero <- rowSums(abs(out$kt)) > 0
      N <- sum(nonzero)
      if(N > 0){
        kt <- matrix(out$kt[nonzero, ], N, ncol(out$kt))
        bx <- matrix(out$bx[, nonzero], nrow(out$bx), N)
        row.names(kt) <- colnames(bx) <- (1:model$N)[nonzero]
        colnames(kt) <- object$years
        rownames(bx) <- object$ages
        periodAgeFun <- model$periodAgeFun[nonzero]
      }
      for (i in 1:model$N){
        if (!nonzero[i]) {
          textFormula <- sub(paste0(" [+] f", i, "[[]x[]] k", i, "[[]t[]]"), "", textFormula)
          textFormula <- sub(paste0("f", i, "[[]x[]] k", i, "[[]t[]]"), "", textFormula)
          textFormula <- sub(paste0(" [+] k", i, "[[]t[]]"), "", textFormula)
          textFormula <- sub(paste0("k", i, "[[]t[]]"), "", textFormula)
        }
      }
    }
    out$kt <- kt
    out$bx <- bx
    
    #simplify age-cohort terms 
    b0x <- NULL
    gc <- NULL
    cohortAgeFun <- model$cohortAgeFun
    if (!is.null(cohortAgeFun)){
      if (sum(abs(out$gc)) > 0) {
        b0x <- out$b0x
        names(b0x) <- object$ages
        gc <- out$gc
        names(gc) <- object$cohorts
      } else{
        textFormula <- sub(" [+] f0[[]x[]] g[[]t[-]x[]]", "", textFormula)
        textFormula <- sub("f0[[]x[]] g[[]t[-]x[]]", "", textFormula)
        textFormula <- sub("g[[]t[-]x[]]", "", textFormula)
        textFormula <- sub(" [+] g[[]t[-]x[]]", "", textFormula)
        cohortAgeFun <- NULL
      }
    }
    out$b0x <- b0x
    out$gc <- gc
    
    #Create new model
    model <- StMoMo(link = model$link, staticAgeFun = staticAgeFun, periodAgeFun = periodAgeFun, 
                    cohortAgeFun = cohortAgeFun)
    model$textFormula <- textFormula
  }
  out$model <- model
  
  #Add the other components of a fitStMoMo object
  out$data <- object$data
  out$Ext <- object$Ext
  out$Dxt <- object$Dxt
  out$oxt <- object$Dxt*0
  out$wxt <- object$wxt
  out$ages <- object$ages
  out$years <- object$years
  out$cohorts <- object$cohorts
  out$loglik <- object$loglik[k]
  out$deviance <- object$deviance[k]
  out$npar <- object$npar[k]
  out$nobs <- object$nobs
  out$call <- match.call()
  class(out) <- "fitStMoMo"
  out
}


