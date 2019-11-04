
#' Compute the link for a given mortality models
#' @keywords internal
predictLink <- function(ax, bx, kt, b0x, gc, oxt, ages, years) {
  
  nAges <- length(ages)
  nYears <- length(years)
  cohorts <- (years[1] - ages[nAges]):(years[nYears] - ages[1])
  nCohorts <- length(cohorts)
    
  #offset
  if (is.null(oxt)) oxt <- array(0, dim = c(nAges, nYears)) 
  
  link <- oxt  
  #static table
  if (!is.null(ax)) {    
    link <- link + array(ax, dim = c(nAges, nYears))  
  }
  #age-period cohort terms
  if (!is.null(bx)) {    
    N <- dim(bx)[2]  
    for (i in 1:N) {
      link <- link + array(bx[, i], dim = c(nAges, nYears)) * 
              t(array(kt[i, ], dim = c(nYears, nAges)))
    }
  }
  
  #age-cohort term
  if (!is.null(gc)) {    
    ft <- gl(nYears, nAges) 
    ftx <- factor(as.numeric(ft) + seq(nAges - 1, 0))
    link <- link + array(b0x, dim = c(nAges, nYears)) * 
                  array(gc[ftx], dim = c(nAges, nYears))
  }  
  dimnames(link)<- list(ages, years)  
  link  
}


#' Compute Poisson loglikelihod
#' 
#' @param obs observed number of deaths
#' @param fit fitted number of deaths
#' @param weight weights given to each observation#' 
#' @keywords internal
computeLogLikPoisson <- function(obs, fit, weight) {
  ind <- (weight > 0)
  res <- array(NA, dim = dim(weight))
  res[ind] <- weight[ind] * (obs[ind] * log(fit[ind]) - 
                          fit[ind] - lfactorial(obs[ind]))
  sum(res, na.rm=TRUE)
}

#' Compute Binomial loglikelihod
#' 
#' @param obs observed rates
#' @param fit fitted rates
#' @param exposure observed exposure
#' @param weight weights given to each observation 
#' @keywords internal
computeLogLikBinomial <- function(obs, fit, exposure, weight) {
  ind <- (weight > 0)
  res <- array(NA, dim = dim(weight))
  res[ind] <- weight[ind] * (exposure[ind] * (obs[ind] * log(fit[ind]) 
              + (1 - obs[ind]) * log(1 - fit[ind])) +
                   lchoose(round(exposure[ind]), round(obs[ind] * exposure[ind])))
  sum(res, na.rm = TRUE)
}


#' Compute Gaussian loglikelihod
#' 
#' @param obs observed log rates
#' @param fit fitted log rates
#' @param sigma2 estimated vairance
#' @param weight weights given to each observation 
#' @keywords internal
computeLogLikGaussian <- function(obs, fit, weight) {
  ind <- (weight > 0)
  res <- array(NA, dim = dim(weight))
  res[ind] <- weight[ind] * (obs[ind] - fit[ind])^2
  sigma2 <- sum(res[ind]) / sum(weight[ind])
  res[ind] <- weight[ind] * (-res[ind]/(2*sigma2) - log(sqrt(2*pi*sigma2)))
  sum(res, na.rm = TRUE)
}


#' Compute Poisson deviance
#' 
#' @param obs observed number of deaths
#' @param fit fitted number of deaths
#' @param weight weights given to each observation#' 
#' @keywords internal
computeDeviancePoisson <- function(obs, fit, weight) {
  ind <- (weight > 0)
  dev <- array(NA, dim = dim(weight))
  dev[ind] <- 2 * weight[ind] * (obs[ind] * log(obs[ind] / fit[ind])
                                 - (obs[ind] - fit[ind]))
  sum(dev, na.rm = TRUE)
}

#' Binomial deviance
#' @param obs observed rates
#' @param fit fitted rates
#' @param exposure observed exposure
#' @param weight weights given to each observation#' 
#' @keywords internal
computeDevianceBinomial <- function(obs, fit, exposure, weight) {
  ind <- (weight > 0)
  dev <- array(NA, dim = dim(weight))
  
  dev[ind] <- 2 * weight[ind] * exposure[ind] * 
             (obs[ind] * log(obs[ind] / fit[ind]) + 
        (1 - obs[ind]) * log((1 - obs[ind]) / (1 - fit[ind])))
  sum(dev, na.rm = TRUE)
}


#' Gaussian deviance
#' @param obs observed log rates
#' @param fit fitted log rates
#' @param sigma2 estimated vairance
#' @param weight weights given to each observation 
#' @keywords internal
computeDevianceGaussian <- function(obs, fit, weight) {
  ind <- (weight > 0)
  dev <- array(NA, dim = dim(weight))
  
  dev[ind] <- weight[ind] * (obs[ind] - fit[ind])^2 
  sum(dev, na.rm = TRUE)
}


#'Logit function
#' \code{logit} computes the logit function
#' @keywords internal
logit <- function(x) log(x / (1 - x))

#' Inverse Logit function
#' \code{invlogit} computes the inverse logit function
#' @keywords internal
invlogit <- function(x) {
  y <- exp(x)
  y / (y + 1)
}

#' Extract a lighter version of a fitted Stochastic Mortality Model
#' 
#' Obtain a lighter version of a fitted Stochastic Mortality Model
#' with the essential information for plotting, forecasting, 
#' and simulation.
#' 
#' @param object an object of class \code{"fitStMoMo"} with the fitted 
#' parameters of a stochastic mortality model.
#' 
#' @return A list with class \code{"fitStMoMo"} with components
#'   
#'   \item{model}{ The \code{StMoMo} defining the fitted stochastic 
#'   mortality model.}
#'   
#'   \item{ax}{ Vector with the fitted values of the static age function
#'   \eqn{\alpha_x}. If the model does not have a static age function or 
#'   failed to fit this is set to \code{NULL}.}
#'     
#'   \item{bx}{ Matrix with the values of the period age-modulating functions
#'   \eqn{\beta_x^{(i)}, i=1, ..., N}. If the \eqn{i}-th age-modulating 
#'   function is non-parametric (e.g. as in the Lee-Carter model) 
#'   \code{bx[, i]} contains the estimated values. If the model does not have
#'   any age-period terms (i.e. \eqn{N=0}) or failed to fit this is set to
#'   \code{NULL}.}
#'   
#'   \item{kt}{ Matrix with the values of the fitted period indexes
#'   \eqn{\kappa_t^{(i)}, i=1, ..., N}. \code{kt[i, ]} contains the estimated
#'   values of the \eqn{i}-th period index. If the model does not have any 
#'   age-period terms (i.e. \eqn{N=0}) or failed to fit this is set to 
#'   \code{NULL}.}
#'   
#'   \item{b0x}{ Vector with the values of the cohort age-modulating function
#'   \eqn{\beta_x^{(0)}}. If the age-modulating function is non-parametric 
#'   \code{b0x} contains the estimated values. If the model does not have a 
#'   cohort effect or failed to fit this is set to \code{NULL}.}
#'     
#'   \item{gc}{ Vector with the fitted cohort index \eqn{\gamma_{c}}. If the
#'   model does not have a cohort effect or failed to fit this is set to 
#'   \code{NULL}.}
#'   
#'   \item{Dxt}{ Matrix of deaths used in the fitting.}
#'   
#'   \item{Ext}{ Matrix of exposures used in the fitting.}
#'   
#'   \item{oxt}{ Matrix of known offset values used in the fitting.}
#'   
#'   \item{wxt}{ Matrix of 0-1 weights used in the fitting.}
#'   
#'   \item{ages}{ Vector of ages in the data.}
#'   
#'   \item{years}{ Vector of years in the data.}
#'   
#'   \item{cohorts}{ Vector of cohorts in the data.}
#'   
#' @keywords internal
getMinimalFitStMoMo <- function(object) {  
  structure(list(model = object$model, ax = object$ax, bx = object$bx, 
                 kt = object$kt, b0x = object$b0x, gc = object$gc, 
                 oxt = object$oxt, ages = object$ages, years = object$years, 
                 cohorts = object$cohorts), class = "fitStMoMo")  
}


#' Copy of unexported function arima.string from forecast
#' @keywords internal
arima.string <- function (object, padding = FALSE) 
{
  order <- object$arma[c(1, 6, 2, 3, 7, 4, 5)]
  result <- paste("ARIMA(", order[1], ",", order[2], ",", order[3], 
                  ")", sep = "")
  if (order[7] > 1 & sum(order[4:6]) > 0) 
    result <- paste(result, "(", order[4], ",", order[5], 
                    ",", order[6], ")[", order[7], "]", sep = "")
  if (is.element("constant", names(object$coef)) | is.element("intercept", 
                                                              names(object$coef))) 
    result <- paste(result, "with non-zero mean")
  else if (is.element("drift", names(object$coef))) 
    result <- paste(result, "with drift        ")
  else if (order[2] == 0 & order[5] == 0) 
    result <- paste(result, "with zero mean    ")
  else result <- paste(result, "                  ")
  if (!padding) 
    result <- gsub("[ ]*$", "", result)
  return(result)
}




#' Remove levels from kt 
#' @keywords internal
constRemoveLevelKt <- function(ax, bx, kt){
  K <- dim(kt)[1]
  for (i in 1:K){
    phi <- kt[i, 1] 
    kt[i, ] <- kt[i, ] - phi[1]
    ax <- ax + phi[1] * bx[, i]
  }
  list(ax = ax, bx = bx, kt = kt)
}

