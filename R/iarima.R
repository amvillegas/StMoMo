#' Fit independent arima series to a multivariate time series
#' 
#' Fits independent arima series to \code{x}, a multivariate 
#' time series.
#' 
#' 
#' @param x numeric matrix with a multivariate time series.
#' Series are arranged in rows with columns representing time.
#' @param order an optional matrix with one row per time series
#' specifying the ARIMA models: for the ith row the three 
#' components \eqn{(p, d, q)} are the AR order, the degree of 
#' differencing, and the MA order. If absent the arima models 
#' are fitted using \code{\link[forecast]{auto.arima}}.
#' @param include.constant an optional vector of logical values 
#' indicating if the ARIMA model for the ith series should include a 
#' constant value. The default is \code{TRUE}. This parameter is ignored
#' if \code{order} is \code{NULL}. 
#' @param ... additional parameters for \code{\link[forecast]{auto.arima}}
#' 
#' @return an object of class \code{"iarima"} with components:
#' \item{models}{a list with the arima models fitted to each time 
#' series.}
#' \item{x}{the original time series.}
#' 
#' @details 
#' The fitting of the ARIMA models for each time series is done with 
#' function \code{\link[forecast]{Arima}} from package 
#' \pkg{forecast}. See the latter function for further details on 
#' input arguments \code{kt.order} and \code{kt.include.constant}.
#' 
#' @export
iarima <- function(x, order = NULL, include.constant = TRUE, ...) {  
  x <- as.matrix(x)
  if (ncol(x) == 1L) 
    x <- t(x)
  N <- nrow(x)
  
  if(!is.null(order)){
    order <- as.matrix(order)
    if (ncol(order) == 1L) 
      order <- t(order)
    if (ncol(order) != 3)
      stop("order must have three columns representing (p,d,q)")
    if (nrow(order) == 1) 
      order <- matrix(order, nrow = N, ncol = 3, byrow = TRUE)
    if (nrow(order) != N)
      stop("number of rows in order does not match number of series.")
    if (length(include.constant) == 1) 
      include.constant <- rep(include.constant, N)
    if (length(include.constant) != N) 
      stop("Length of include.constant does not match number of series")
  }
  
  models <- list()
  if(is.null(order)){
    for (i in 1:N)  {
      models[[i]] <- forecast::auto.arima(x[i, ], ...)    
    }
  } else {
    for (i in 1:N)  {
      models[[i]] <- forecast::Arima(x[i, ], order = order[i, ], 
                           include.constant = include.constant[i])    
    }
  }
  structure(list(models = models, x = x), class = "iarima")
}


#' Forecast independent arima series
#' 
#' Returns forecasts and other information for a group of
#' independent arima series.
#' 
#' @param object an object of class \code{"iarima"}.
#' @param h Number of periods for forecasting.
#' @param level confidence level for prediction intervals.
#' @param fan if \code{TRUE}, level is set to \code{seq(50, 99, by = 1)}. 
#' This is suitable for fan plots.
#' @param ... other arguments.
#' 
#' @return An object of class \code{"iarimaForecast"} with components:
#' \item{model}{a list containing information about the fitted arima models.}
#' \item{mean}{ array with the central forecast.}
#' \item{lower}{ three dimensional array with lower limits for prediction 
#'  intervals.}
#' \item{upper}{ three dimensional array with upper limits for prediction 
#'  intervals.}
#'  \item{level}{ the confidence values associated with the prediction 
#'  intervals.}
#'  @export  
forecast.iarima <- function(object, h = 10, level = c(80,95), fan = FALSE, ...) {
  
  x <- object$x
  nYear <- ncol(x)
  N <- nrow(x)  
  yearsFor <- (as.numeric(colnames(x)[nYear]) + 1):(as.numeric(colnames(x)[nYear]) + h)
  
  mean <- array(NA, dim = c(N, h), dimnames = list(rownames(x), yearsFor))
  
  if (fan) 
    level <- seq(51, 99, by = 3)
  else {
    if (min(level) > 0 & max(level) < 1) 
      level <- 100 * level
    else if (min(level) < 0 | max(level) > 99.99) 
      stop("Confidence limit out of range")
  }
  nconf <- length(level)
  lower <- upper <- array(NA, c(N, h, nconf), 
                          dimnames = list(rownames(x), yearsFor, 
                                          paste(level, "%", sep = "")))
  
  for(i in 1:N){
    kt.for <- forecast::forecast(object$models[[i]], h = h, level = level)
    mean[i, ] <- kt.for$mean
    lower[i, ,] <- kt.for$lower 
    upper[i, ,] <- kt.for$upper
  }
  
  structure(list(model = object, level = level, mean = mean, lower = lower, 
                 upper = upper), class = "iarimaForecast")    
}


#' Simulate independent arima series
#' 
#' Returns one simulated path of the group of independent 
#' arima series in \code{object}.
#' 
#' @param object An object of class \code{"iarima"}.
#' @param nsim number of periods for the simulated series.
#' @param seed either \code{NULL} or an integer that will be used in a 
#' call to \code{\link{set.seed}} before simulating the time series. 
#' The default, \code{NULL} will not change the random generator state.
#' @param ... other arguments.
#' 
#' @export
simulate.iarima <- function(object, nsim = 10, seed = NULL, ...) {
  #Hack to remove notes in CRAN check
  x <- NULL
  if (!exists(".Random.seed", envir = .GlobalEnv)) 
    runif(1)
  if (is.null(seed)) 
    RNGstate <- .Random.seed
  else {
    R.seed <- .Random.seed
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  
  sim <- t(simplify2array(lapply(object$models, FUN = simulate, nsim = nsim)))
  rownames(sim) <- rownames(object$x)
  x <- object$x
  nYear <- ncol(x)
  yearsSim <- (as.numeric(colnames(x)[nYear]) + 1):(as.numeric(colnames(x)[nYear]) + nsim)
  colnames(sim) <- yearsSim
  sim
}
