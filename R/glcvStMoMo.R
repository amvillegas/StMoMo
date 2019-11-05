#' Cross validation for Stochastic Mortality Model fitted 
#' with Group Regularised Penalties
#' 
#' Perform cross validation by period for a Stochastic Mortality Model
#' being fitted using the group regularised penalties.
#' 
#' @param object an object of class \code{"StMoMo"} defining the 
#' stochastic mortality model.
#' 
#' @param h number of years for forecasting horizon.
#' 
#' @param lambda optional vector of lambda values for group lasso 
#' fitting. If not specified, a sequence of values of length nlambda 
#' is computed automatically, equally spaced on the log scale.
#' 
#' @param nlambda the number of lambda values (if \code{lambda}) is
#' left unspecificed. Default is 50.
#' 
#' @param data an optional object of type StMoMoData containing 
#' information on deaths and exposures to be used for training the model. 
#' This is typically created with function \code{\link{StMoMoData}}. 
#' If this is not provided then the training data is taken from 
#' arguments, \code{Dxt}, \code{Ext}, \code{ages}, \code{years}.  
#' 
#' @param Dxt optional matrix of deaths data.
#' 
#' @param Ext optional matrix of observed exposures of the same 
#' dimension of \code{Dxt}.
#' 
#' @param ages optional vector of ages corresponding to rows of 
#' \code{Dxt} and \code{Ext}. 
#' 
#' @param years optional vector of years corresponding to rows of 
#' \code{Dxt} and \code{Ext}. 
#' 
#' @param ages.train optional vector of ages to include in the 
#' training. Must be a subset of \code{ages}. 
#' 
#' @param years.train optional vector of years to include in the 
#' training. Must be a subset of \code{years}.
#' 
#' @param index optional vector describing the grouping of the coefficients. If
#' there are coefficients to be included in the model without being penalised, assign
#' them to group 0 (or "0"). If this is not provided, the grouping is automatically
#' determined and all groups are penalised.
#'  
#' @param returnY a logical value. If \code{TRUE}, \code{cv.grpStMoMo}
#' returns a list of length \code{nlambda}, containing matrices of the fitted mortality rates 
#' for the cross validation folds for each value of \code{lambda}. 
#' 
#' @param type the type of the predicted values that the cross 
#' validation is performed with respect to. The alternatives are 
#' \code{"rates"} and \code{"logrates"}. If \code{"rates"}, \code{cvStMoMo}
#' returns the mean squared error using the deviation between
#' realised and predicted mortality rates. If \code{"logrates"}, 
#' \code{cvStMoMo} returns the mean squared error using the 
#' deviation between realised and predicted log mortality rates. \code{lambda.min}, 
#' \code{min}, \code{sd}, and \code{se} also depend on \code{type}.
#'     
#' @param verbose a logical value. If \code{TRUE} progress 
#' indicators are printed as the model is fitted. Set 
#' \code{verbose = FALSE} to silent the fitting and avoid 
#' progress messages.
#' 
#' @return A list with class \code{"cv.grpStMoMo"} with components:
#'   
#'   \item{model}{ the object of class \code{"StMoMo"} defining 
#'   the fitted stochastic mortality model.}
#'   
#'   \item{data}{ StMoMoData object provided for training the model.}
#'   
#'   \item{Dxt}{ matrix of deaths used in the training.}
#'   
#'   \item{Ext}{ matrix of exposures used in the training.}
#'   
#'   \item{qxt}{ matrix of mortality rates (\eqn{D_{x,t}/E_{x,t}}) 
#'   used in the training.} 
#'    
#'   \item{Lqxt}{ matrix of log mortality rates (\eqn{log(q_{x,t})}) 
#'   used in the training.} 
#'   
#'   \item{lambda}{ vector of lambda values used in the fitting.}
#'   
#'   \item{cv.rates}{ if \code{returnY=TRUE}, matrix of predicted 
#'   mortality rates.} 
#'   
#'   \item{cv.mse}{ if \code{type="rates"}, vector mean squared errors using predicted mortality rates and realised mortality rates. If \code{type="logrates"}, vector of mean squared errors using predicted log mortality rates and realised log mortality rates.}  
#'   
#'   \item{sd}{ standard deviation for cross validation.}
#'   
#'   \item{se}{ standard error for cross validation.}
#'   
#'   \item{min}{ the index of lambda corresponding to lambda.min.}
#'   
#'   \item{lambda.min}{ the value of lambda with the minimum cross-validation error
#'   (\code{cv.mse}).}
#'  
#' @examples
#' 
#' # One-year forecasting horizon based on log-rates using 25 lambda values
#' APCcvgrp <- cv.grpStMoMo(apc(link = "log-Gaussian"), h = 1, data = EWMaleData, ages.train = 55:89, 
#'                      type = "logrates", nlambda = 25)
#' 
#' # Ten-year forecasting using arguments Dxt, Ext, ages, years to pass fitting data, using pre-defined lambda grid
#' lambda <- c(0.5, 0.4, 0.3, 0.2, 0.1)
#' APCcvgrp <- cv.grpStMoMo(apc(link = "log-Gaussian"), h = 10, Dxt = EWMaleData$Dxt, 
#'                      Ext = EWMaleData$Ext, ages = EWMaleData$ages, years = EWMaleData$years, 
#'                      ages.train = 55:89, type = "logrates", lambda = lambda)
#'             
#' @export 
cv.grpStMoMo <- function(object,  h = 1, lambda = NULL, nlambda = 50, data = NULL, Dxt = NULL, Ext = NULL, 
                       ages.train = NULL, years.train = NULL, ages = NULL, years = NULL, index = NULL, 
                       returnY = FALSE, type = c("rates", "logrates"), verbose = TRUE) {
  
  
  type <- match.arg(type)  
  # Determine fitting ages and years are specified
  if(!is.null(data)) {
    if (class(data) != "StMoMoData") {
      stop("Argument data needs to be of class StMoMoData.")
    }
    ages <- data$ages
    years <- data$years
  } else {
    if (is.null(Dxt) || is.null(Ext)) {
      stop("Either argument data or arguments Dxt and Ext need to be provided.")
    }
    if (is.null(ages)) ages <- 1:nrow(Dxt)
    if (is.null(years)) years <- 1:ncol(Dxt)
  }
  if (is.null(ages.train)) ages.train <- ages
  if (is.null(years.train)) years.train <- years
  
  # Reasonability check between length of data set and forecasting horizon 
  if (h > length(years.train)/2) {
    warning("Forecasting horizon is too large for data set supplied.")
  } 
  
  # Determine number of test sets
  folds <- length(years.train) - h
  
  # Determine lambda grid
  if (is.null(lambda)) {
    fitTemp <- grpfit(object, nlambda = nlambda, data = data, Dxt = Dxt, Ext = Ext, 
                        ages = ages, years = years, ages.fit = ages.train, years.fit = years.train, verbose = FALSE) 
    lambda <- fitTemp$lambda
  } else {
    fitTemp <- grpfit(object, lambda = lambda, data = data, Dxt = Dxt, Ext = Ext, 
                        ages = ages, years = years, ages.fit = ages.train, years.fit = years.train, verbose = FALSE) 
    nlambda <- length(lambda)
  }
  
  # Initialise out-of-sample predicted values
  cv.rates <- replicate(nlambda, genWeightMat(ages = ages.train, years = years.train)*0, simplify=FALSE)
  cv.mse <- rep(0, nlambda)
  cv.Lmse <- rep(0, nlambda)
  
  # Perform cross validation iterations
  for (i in 2:(length(years.train)-h+1)) {
    
    if (verbose) {
     cat("StMoMo: Fitting test set", i-1, "of", folds, "\n")
    }
    
    # Define training set weight matrix
    wxtT <- genWeightMat(ages = ages.train, years = years.train)
    wxtT[,(i:(i+h-1))] <- 0
    
    # Define test set weight matrix
    wxtF <- genWeightMat(ages = ages.train, years = years.train)
    wxtF[,((i+h-1):(i+h-1))] <- 0
    
    # Fit using grpfit
    fit <- grpfit(object, lambda = lambda, data = data, Dxt = Dxt, Ext = Ext, 
                    ages = ages, years = years, ages.fit = ages.train, years.fit = years.train, 
                    wxt = wxtT, verbose = FALSE) 
    
    # Compute test set errors

    for (n in 1:nlambda) {
      
      # Extract and impute kappa parameters
      N <- object$N
      
      if (N > 0) {
        kt <- vector()
        for (j in 1:N) {
          k <- fit$beta[[n]]$kt[j,]
          if (min(k, na.rm = TRUE) == max(k, na.rm = TRUE)) {
            drift <- 0
          } else {
            arima.k <- tryCatch( Arima(k, order = c(0,1,0), include.drift = TRUE), error = function( err ) FALSE, warning = function( err ) FALSE )
            if(!is.logical(arima.k)) {
              drift <- arima.k$coef
            } else {
              drift <- NA
              warning(paste("Time series model did not converge for test set",  i-1, "of", folds, "for lambda[", n, "]\n", sep = ""))
            }
          }
          for (p in 2:length(k)) {
            if (is.na(k[p])){
              k[p] <- k[p-1] + drift
            }
          }
          kt <- rbind(kt, k)
        }
      }
    
      # Compute fitted values (training and test sets)
      kt <- as.matrix(kt)
      gc <- as.vector(fit$beta[[n]]$gc)
      Lqhat <- predictLink(ax = fit$beta[[n]]$ax, bx = fit$beta[[n]]$bx, kt = kt, b0x = fit$beta[[n]]$b0x, 
                           gc = gc, oxt = NULL, ages = ages.train, years = years.train)  
      qhat <- exp(Lqhat)
      
      # Update rates matrix
      qhat[is.na(qhat)] <- 0
      cv.rates[[n]] <- cv.rates[[n]] + (qhat * (1 - wxtF))
      
    }
    
  }

  # Prepare data output
  data <- fitTemp$data
  Dxt <- fitTemp$Dxt
  Ext <- fitTemp$Ext
  qxt <- Dxt/Ext
  Lqxt <- log(qxt)
    
  # Prepare cv prediction output
  for (j in 1:nlambda) {
    cv.rates[[j]][cv.rates[[j]] == 0] <- NA
    cv.Lrates <- log(cv.rates[[j]])
    
    # Prepare cv error output
    cv.error <- (cv.rates[[j]] - qxt)^2
    cv.mse[j] <- mean(cv.error[!is.na(cv.error)])
    
    cv.Lerror <- (cv.Lrates - Lqxt)^2
    cv.Lmse[j] <- mean(cv.Lerror[!is.na(cv.Lerror)])
  }

  # Determine optimal lambda
  min <- which.min(cv.mse)
  lambda.min <- lambda[min]
  
  Lmin <- which.min(cv.Lmse)
  lambda.Lmin <- lambda[Lmin]
  
  # Determine standard errors at optimal lambda
  cv.error.lambda.min <- (cv.rates[[min]] - qxt)^2
  cv.lambda.min <- colMeans(cv.error.lambda.min, na.rm = TRUE)
  sd <- sd(cv.lambda.min, na.rm = TRUE)
  se <- sd/sqrt(folds)
  
  cv.Lerror.lambda.Lmin <- (log(cv.rates[[Lmin]]) - Lqxt)^2
  cv.lambda.Lmin <- colMeans(cv.Lerror.lambda.Lmin, na.rm = TRUE)
  Lsd <- sd(cv.lambda.Lmin, na.rm = TRUE)
  Lse <- Lsd/sqrt(folds)
  
  # Return output
  out <- list(model = object, data = data, Dxt = Dxt, Ext = Ext, rates = qxt, logrates = Lqxt, 
              lambda = lambda, call = match.call())
  
  if (returnY == TRUE) {
    out$cv.rates = cv.rates
  } 
  
  if (type == "rates") {
    out$cv.mse = cv.mse 
    out$min = min
    out$lambda.min = lambda.min
    out$sd = sd
    out$se = se
  } else if (type == "logrates") {
    out$cv.mse = cv.Lmse
    out$min = Lmin
    out$lambda.min = lambda.Lmin
    out$sd = Lsd
    out$se = Lse
  } 
 
  class(out) <- "cv.grpStMoMo"
  return(out)
}