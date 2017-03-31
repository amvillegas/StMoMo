#' Predict method for Stochastic Mortality Models fits
#' 
#' Obtain predictions from a Stochastic Mortality Model
#' fit. 
#' 
#' This function evaluates 
#' \deqn{\hat{\eta}_{xt} = o_{xt} + \alpha_x + 
#' \sum_{i=1}^N \beta_x^{(i)}\kappa_t^{(i)} + \beta_x^{(0)}\gamma_{t-x}}
#' for a fitted Stochastic Mortality model.
#' In producing a prediction the static age function, \eqn{\alpha_x}, and the
#' age-modulating parameters, \eqn{\beta_x^{(i)}, i=0, ..., N}, are taken from
#' the fitted model in \code{object} while the period indexes, 
#' \eqn{\kappa_t^{(i)}, i=1,..., N}, and cohort index, \eqn{\gamma_{t-x}}, 
#' are taken from the function arguments.  
#' 
#' This function can be useful, for instance, in producing forecasts of 
#' mortality rates using time series models different to those available 
#' in \code{\link{forecast.fitStMoMo}} (See examples below). 
#' 
#' @param object an object of class \code{"fitStMoMo"} with the fitted 
#' parameters of a stochastic mortality model.
#' 
#' @param years vector of years for which a prediction is required.
#' 
#' @param kt matrix of values of the period indexes to use for the prediction. 
#' If the model has any age-period term this argument needs to be provided and 
#' the number of rows in \code{kt} must be equal to the number of age-period 
#' terms in the model and the number of columns in \code{kt} must correspond 
#' to the length of \code{years}. If the Stochastic Mortality Model doesn't 
#' have any age-period terms this argument is ignored and needs not be 
#' provided.
#'  
#' @param gc vector of values of the cohort indexes to use for the prediction. 
#' If the model has a cohort effect this argument needs to be provided. 
#' In this case the length of \code{gc} must be equal to the number of cohorts 
#' for which a prediction is being produced, namely, 
#' \code{length(object$ages) + length(years) - 1}. If the Stochastic Mortality 
#' Model doesn't have a cohort effect  this argument is ignored and needs not 
#' be provided. 
#' 
#' @param oxt optional matrix/vector or scalar of known offset to be used in 
#' the prediction.
#' @param type the type of the predicted values that should be returned. The 
#' alternatives are \code{"link"}(default) and \code{"rates"}.
#' @param ... other arguments.
#'   
#' @return A matrix with the predicted values.
#'   
#' @seealso \code{\link{forecast.fitStMoMo}}   
#'   
#' @examples
#' 
#' \dontrun{
#' ##M6 Forecast using VARIMA(1,1) model
#' library(MTS)
#' 
#' # fit m6
#' years <- EWMaleData$years
#' ages.fit <- 55:89
#' M6 <- m6(link = "log")
#' M6fit <- fit(M6, data = EWMaleData, ages.fit = ages.fit)
#' 
#' # Forecast kt using VARIMA(1,1) model from MTS
#' h <- 50
#' kt.M6 <- t(M6fit$kt) 
#' kt.M6.diff <- apply(kt.M6, 2, diff)
#' fit.kt.M6.11 <- VARMA(kt.M6.diff, p = 1, q = 1)
#' pred.ktdiff.M6.11 <- VARMApred(fit.kt.M6.11, h = h)
#' pred.kt.M6.11 <- apply(rbind(tail(kt.M6, n = 1),
#'                              pred.ktdiff.M6.11$pred), 
#'                        2, cumsum)[-1, ]
#' 
#' # set row names
#' years.forecast <- seq(tail(years, 1) + 1, length.out = h)
#' rownames(pred.kt.M6.11) <- years.forecast
#' 
#' # plot kt1
#' plot(x = c(years, years.forecast),
#'      y = c(kt.M6[, 1], pred.kt.M6.11[, 1]),
#'      col = rep(c("black", "red"), times = c(length(years), h)),
#'      xlab = "time",
#'      ylab = "k1")
#' 
# plot kt2
#' plot(x = c(years, years.forecast),
#'      y = c(kt.M6[, 2], pred.kt.M6.11[, 2]),
#'      col = rep(c("black", "red"), times = c(length(years), h)),
#'      xlab = "time",
#'      ylab = "k2")
#' 
#' # forecast cohort effect
#' # the following cohorts are required:
#' # from 2012 - 89 = 1923
#' # to 2061 - 55 = 2006
#' pred.gc.M6 <- forecast(auto.arima(M6fit$gc, max.d = 1), h = h)
#' 
#' # use predict to get rates
#' pred.qxt.M6.11 <- predict(object = M6fit,
#'                           years = years.forecast,
#'                           kt = t(pred.kt.M6.11),
#'                           gc = c(tail(M6fit$gc, 34), pred.gc.M6$mean),
#'                           type = "rates")
#' 
#' qxthatM6 <- fitted(M6fit, type = "rates")
#' 
#' # plot mortality profile at age 60, 70 and 80
#' matplot(1961 : 2061,
#'         t(cbind(qxthatM6, pred.qxt.M6.11)[c("60", "70", "80"), ]),
#'         type = "l", col = "black", xlab = "years", ylab = "rates",
#'         lwd = 1.5)
#' }
#' @export 
predict.fitStMoMo <- function(object, years, kt = NULL, gc = NULL, oxt = NULL, 
                              type = c("link", "rates"), ...) {
  type <- match.arg(type)
  ages <- object$ages
  nAges <- length(ages)
  nYears <- length(years)
  cohorts <- (years[1] - ages[nAges]):(years[nYears] - ages[1])
  nCohorts <- length(cohorts)
  
  #Check inputs
  if (object$model$N > 0) {
    kt <- as.matrix(kt)
    if (ncol(kt) == 1) 
      kt <- t(as.matrix(kt))
    if (ncol(kt) != nYears) {
      stop( "Mismatch between the dimension of kt and the number of years.")
    }
    if (nrow(kt) != object$model$N) {
      stop("Mismatch between the dimension of kt and the number of age-period 
           terms in the model.")
    }
  } else {
    if (!is.null(kt)) {
      warning("kt argument ignored as the model doesn't have any age-period
              terms.")
    }
    kt <- NULL
  }
  if (!is.null(object$model$cohortAgeFun)) {
    gc <- as.vector(gc)
    if (length(gc) != nCohorts) {
      stop(paste("Mismatch between the length of gc and the number of cohorts 
                  required. Required number of cohorts:", nCohorts))
    }
  } else {
    if (!is.null(gc)) {
      warning("gc argument ignored as the model doesn't have an age-cohort
              term.")
     }
    gc <- NULL
  }
  if (is.null(oxt)) {
    oxt <- matrix(0, nrow = nAges, ncol = nYears)
    rownames(oxt) <- ages
    colnames(oxt) <- years        
  } else {
    oxt <- matrix(oxt, nrow = nAges, ncol = nYears)
    rownames(oxt) <- ages
    colnames(oxt) <- years             
  }   
  
  #compute prediction
  link <- predictLink(ax = object$ax, bx = object$bx, kt = kt, 
                      b0x = object$b0x, gc = gc, oxt = oxt, 
                      ages = ages, years = years)
  rates <- switch(object$model$link, log = exp(link), logit = invlogit(link))
  switch(type, rates = rates, link = link)  
}