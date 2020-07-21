

#' Plot the regularisation path of a Stochasic Mortality Model
#' 
#' Plots the regularisation path of a Stochasic Mortality Model 
#' fitted with function \code{grpfit}.
#' 
#' @usage 
#' \method{plot}{grpfitStMoMo}(x, label = FALSE, include.ax = FALSE,  ...)
#' 
#' @param x an object of class \code{"grpfitStMoMo"}.
#' @param label If \code{TRUE}, annotates the plot with text labels in the right margin 
#' describing which age/period or cohort terms the lines belong o
#' @param include.ax If \code{FALSE} the the norm of each group of parameters is calculated
#' relative to the predictor witout \eqn{\alpha_x} and the line related to \eqn{\alpha_x} is
#' not plotted. 
#' @param ... additional arguments to control graphical appearance.
#' See \code{\link[graphics]{plot}}. 
#' 
#' @export 
#' @method plot grpfitStMoMo
plot.grpfitStMoMo <- function(x, label = FALSE, include.ax = FALSE, ...) {
  
  xlab <- expression(Log(lambda))
  ylab <- expression("||X"*beta[j]*"||/||X"*beta*"||")
  ages <- x$ages
  years <- x$years
  nAges <- length(ages)
  nYears <- length(years)
  
  
  nlambda <- length(x$lambda)
  N <- x$model$N
  
  anorm <- rep(0,nlambda)
  bnorm <- matrix(0, nrow = nlambda, ncol = N)
  gnorm <- rep(0,nlambda)
  
  for (i in 1:nlambda){
    
    
    mxtNorm <- sqrt(sum(predictLink(ax = x$beta[[i]]$ax, bx = x$beta[[i]]$bx, 
                           kt = x$beta[[i]]$kt, b0x = x$beta[[i]]$b0x, 
                           gc = x$beta[[i]]$gc, oxt = NULL, 
                           ages = ages, years = years)^2))
    
    if (x$model$staticAgeFun && include.ax){
      
      anorm[i] <- sqrt(sum(predictLink(ax = x$beta[[i]]$ax, bx = NULL, kt = NULL, 
                                      b0x = NULL, gc = NULL, oxt = NULL, 
                                      ages = ages, years = years)^2))/mxtNorm
      if (is.nan(anorm[i])) anorm[i] <- 0
      
    }
    if (N > 0){
      for (j in 1:N) {
        bx <- matrix(x$beta[[i]]$bx[,j], nrow =nAges, ncol = 1)
        kt <- matrix(x$beta[[i]]$kt[j, ], nrow = 1, ncol = nYears)
        bnorm[i, j] <- sqrt(sum(predictLink(ax = NULL, bx = bx, kt = kt, 
                                         b0x = NULL, gc = NULL, oxt = NULL, 
                                         ages = ages, years = years)^2))/mxtNorm
        if (is.nan(bnorm[i, j])) bnorm[i, j] <- 0
        
      }  
    }
    
    if (!is.null(x$model$cohortAgeFun)){
      gnorm[i] <- sqrt(sum(predictLink(ax = NULL, bx = NULL, kt = NULL, 
                                       b0x = x$beta[[i]]$b0x, gc = x$beta[[i]]$gc, oxt = NULL, 
                                       ages = ages, years = years)^2))/mxtNorm
      if (is.nan(gnorm[i])) gnorm[i] <- 0
    }
    
  }
  
  
  allnorm <- matrix(NA, nrow = nlambda, ncol = 1 + N + 1)
  
  if(include.ax) {
    allnorm[, 1] <- anorm
  } 
  
  if (x$model$N>0) {
    for (i in 1:x$model$N){
      if (any(bnorm[, i] > 0, na.rm = TRUE))
        allnorm[, 1 + i] <- bnorm[, i]
    }
  }
  
  
  
  allnorm[, N + 2]  <- gnorm
  indexNoNA <- !is.na(colSums(allnorm))
  
  xtext <- min(log(x$lambda)) - (max(log(x$lambda)) - min(log(x$lambda))) * 0.03
  
  matplot(x = log(x$lambda), y = allnorm[, indexNoNA], 
          xlab = xlab, ylab = ylab, type  = "l", 
          xlim = c(xtext,
                   max(log(x$lambda))), ...)
  
  
  #Add annotation
  
  if(include.ax && indexNoNA[1]) {
    text(xtext, allnorm[nlambda, 1], expression(alpha))
  }
  if (x$model$N>0) {
    for (i in 1:x$model$N){
      if (indexNoNA[i + 1])
        text(xtext, allnorm[nlambda, i + 1], substitute(f^{"("*k*")"}, list(k = i)))
    }
  }
  if(indexNoNA[N+2]) {
    text(xtext, allnorm[nlambda, N + 2], expression(gamma))
  }
  
}
