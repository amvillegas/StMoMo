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
#' @param penalise.ax optional logical value indicating whether the \eqn{\alpha_x}
#' should be penalsied in the fitting. Default is \code{FALSE} so \eqn{\alpha_x} is
#' not penalised.
#' 
#' @param penalise.kt optional logical vector of length \eqn{N} indicating whether 
#' the \eqn{i}-th period term, \eqn{\kappa^{(i)}_t}, should be penalised. By default
#' all period term are penalised.
#' 
#' @param penalise.gc optional logical value indicating whether the cohort effect
#' \eqn{\gamma_{c}} should be penalised in the fitting. By default it is not penalised.
#'  
#' @param returnY a logical value. If \code{TRUE}, \code{cv.grpStMoMo}
#' returns a list of length \code{nlambda}, containing matrices of the fitted mortality rates 
#' for the cross validation folds for each value of \code{lambda}. 
#' 
#' @param type the type of the predicted values that the cross 
#' validation is performed with respect to. The alternatives are 
#' \code{"rates"} and \code{"logrates"}. If \code{"rates"}, \code{cv.grpStMoMo}
#' returns the mean squared error using the deviation between
#' realised and predicted mortality rates. If \code{"logrates"}, 
#' \code{cv.grpStMoMo} returns the mean squared error using the 
#' deviation between realised and predicted log mortality rates. \code{lambda.min}, 
#' \code{min} and \code{se} also depend on \code{type}.
#'     
#' @param all.horizon a logical value. If \code{TRUE}, 
#'  the evaluation is done pooling all horizons
#'  ranging from 1 to \code{h}. If \code{FALSE}, it is done only for horizon
#' \code{h}.
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
#'   \item{cv.rates}{ if \code{returnY=TRUE}, list with matrices of predicted 
#'   mortality rates at horizon \code{h} if \code{all.horizon = FALSE}. 
#'   If \code{all.horizon = TRUE} then it is list of three dimensional arrays
#'   where the first dimension indicates each horizon from 1 to \code{h}} 
#'   
#'   \item{cv.mse}{ if \code{type="rates"}, mean squared error using 
#'   predicted mortality rates and realised mortality rates. If 
#'   \code{type="logrates"}, mean squared error using predicted log 
#'   mortality rates and realised log mortality rates. 
#'   If \code{all.horizon = TRUE} this MSE is pooled across all horizon
#'   from 1 to \code{h}.}  
#'   
#'   \item{cv.mse.h}{ Only present if \code{all.horizon=TRUE}.
#'   if \code{type="rates"}, vector with mean squared error 
#'   using     predicted mortality rates and realised mortality rates for each 
#'   horizon. If \code{type="logrates"}, vector with mean squared error 
#'   using predicted log mortality rates and realised log mortality rates 
#'   for each horizon.}
#'   
#'   \item{se}{ the estimated standard error for \code{cv.mse}.}
#'   
#'   \item{min}{ the index of lambda corresponding to \code{lambda.min}.}
#'   
#'   \item{lambda.min}{ the value of lambda with the minimum cross-validation error
#'   (\code{cv.mse}).}
#'   
#'   \item{lambda.1.se}{ the value of lambda such that error is within 1 
#'   standard error of the minimum.}
#'   
#'   \item{min.1se}{ the index of lambda corresponding to \code{lambda.1se}.}
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
cv.grpStMoMo <- function(object,  h = 1, lambda = NULL, nlambda = 50, data = NULL, Dxt = NULL, 
                         Ext = NULL, ages.train = NULL, years.train = NULL, ages = NULL, years = NULL, 
                         penalise.ax = FALSE, penalise.kt = rep(TRUE, object$N), penalise.gc = TRUE,
                         returnY = FALSE, type = c("logrates", "rates"), 
                         all.horizon = FALSE,
                         verbose = TRUE) {
  
  
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
                      ages = ages, years = years, ages.fit = ages.train, years.fit = years.train, 
                      penalise.ax = penalise.ax, penalise.kt = penalise.kt, penalise.gc = penalise.gc,
                      verbose = FALSE) 
    lambda <- fitTemp$lambda
  } else {
    fitTemp <- grpfit(object, lambda = lambda, data = data, Dxt = Dxt, Ext = Ext, 
                      ages = ages, years = years, ages.fit = ages.train, years.fit = years.train, 
                      penalise.ax = penalise.ax, penalise.kt = penalise.kt, penalise.gc = penalise.gc,
                      verbose = FALSE) 
    nlambda <- length(lambda)
  }
  
  # Initialise out-of-sample predicted values
  if (all.horizon == FALSE){
    cv.rates <- replicate(nlambda, genWeightMat(ages = ages.train, years = years.train)*0, simplify=FALSE)
    
  } else {
    cv.rates.j <- genWeightMat(ages = ages.train, years = years.train)*0
    cv.rates <- array(NA, dim = c(h, nrow(cv.rates.j), ncol(cv.rates.j)), 
                      dimnames = list(1:h, 
                                      rownames(cv.rates.j), 
                                      colnames(cv.rates.j)))
    for (j in 1:h) {
      cv.rates[j, , ] <- cv.rates.j
    }
    
    cv.rates_list <- vector("list", nlambda) 
    for (n in 1:nlambda) {
      cv.rates_list[[n]] <- cv.rates
    }  
    cv.rates <- cv.rates_list
    
  }
  
  cv.mse <- rep(0, nlambda)
  cv.se <- rep(0, nlambda)
  cv.Lmse <- rep(0, nlambda)
  cv.Lse <- rep(0, nlambda)
  
  if(all.horizon == TRUE) {
    cv.mse.h <- replicate(nlambda, rep(NA, h), simplify=FALSE)
    cv.Lmse.h <- replicate(nlambda, rep(NA, h), simplify=FALSE)
  }  
  
  # Perform cross validation iterations
  for (i in 2:(length(years.train)-h+1)) {
    
    if (verbose) {
     cat("StMoMo: Fitting test set", i-1, "of", folds, "\n")
    }
    
    # Define training set weight matrix
    wxtT <- genWeightMat(ages = ages.train, years = years.train)
    wxtT[,(i:(i+h-1))] <- 0
    
    # Define test set weight matrix
    if(all.horizon == FALSE){
      wxtF <- genWeightMat(ages = ages.train, years = years.train)
      wxtF[,((i+h-1):(i+h-1))] <- 0
    } else{
      wxtF.j <- genWeightMat(ages = ages.train, years = years.train)
      wxtF <- array(NA, dim = c(h, nrow(wxtF.j), ncol(wxtF.j)), 
                    dimnames =  list(1:h, 
                                     rownames(wxtF.j), 
                                     colnames(wxtF.j)))
      for (j in 1:h) {
        wxtF[j, , ] <- wxtF.j
        wxtF[j, ,((i+j-1):(i+j-1))] <- 0
        
      }
    }
    
    # Fit using grpfit
    fit <- grpfit(object, lambda = lambda, data = data, Dxt = Dxt, Ext = Ext,  ages = ages, 
                  years = years, ages.fit = ages.train, years.fit = years.train,  wxt = wxtT, 
                  penalise.ax = penalise.ax, penalise.kt = penalise.kt, penalise.gc = penalise.gc,
                  verbose = FALSE) 
    
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
            arima.k <- tryCatch( mrwd(k), error = function( err ) FALSE, warning = function( err ) FALSE )
            if(!is.logical(arima.k)) {
              drift <- arima.k$drift
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
      
      if(all.horizon == FALSE){
        cv.rates[[n]] <- cv.rates[[n]] + (qhat * (1 - wxtF))
      } else{
        for (j in 1:h) {
          cv.rates[[n]][j, ,] <- cv.rates[[n]][j, ,] + (qhat * (1 - wxtF[j, ,]))
        }
      }
      
      
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
    if (all.horizon == FALSE) {
      cv.error <- (cv.rates[[j]] - qxt)^2
      cv.mse[j] <- mean(cv.error[!is.na(cv.error)])
      cv <- colMeans(cv.error, na.rm = TRUE)   
      sd <- sd(cv, na.rm = TRUE)
      cv.se[j] <- sd/sqrt(folds)
      
      cv.Lerror <- (cv.Lrates - Lqxt)^2
      cv.Lmse[j] <- mean(cv.Lerror[!is.na(cv.Lerror)])
      cv.L <- colMeans(cv.Lerror, na.rm = TRUE)
      Lsd <- sd(cv.L, na.rm = TRUE)
      cv.Lse[j] <- Lsd/sqrt(folds)
    } else {
      cv.error.h <- cv.Lerror.h <- NA*cv.rates[[j]]
      for (k in 1:h) {
        cv.error.h[k, , ] <- (cv.rates[[j]][k, , ] - qxt)^2
        cv.error <- cv.error.h[k, , ]
        cv.mse.h[[j]][k] <- mean(cv.error, na.rm = TRUE)
        cv.Lerror.h[k, , ] <- (cv.Lrates[k, , ] - Lqxt)^2
        cv.Lerror <- cv.Lerror.h[k, , ]
        cv.Lmse.h[[j]][k] <- mean(cv.Lerror, na.rm = TRUE)
      }
      cv.mse[j] <- mean(cv.mse.h[[j]])
      cv.Lmse[j] <- mean(cv.Lmse.h[[j]])
      
      cv <- cv.L <- rep(0, folds)
      for (i in 2:(length(years.train)-h+1)) {
        for (k in 1:h) {
          cv[i-1] <- cv[i-1] + mean(cv.error.h[k, , k + i - 1], na.rm = TRUE)/h
          cv.L[i-1] <- cv.L[i-1] + mean(cv.Lerror.h[k, , k + i - 1], na.rm = TRUE)/h 
        }
      }
      sd <- sd(cv, na.rm = TRUE)
      cv.se[j] <- sd/sqrt(folds)
      Lsd <- sd(cv.L, na.rm = TRUE)
      cv.Lse[j] <- Lsd/sqrt(folds)
      
    }
    
    
  }

  # Determine optimal lambda
  min <- which.min(cv.mse)
  lambda.min <- lambda[min]
  min.1se <- min(which(cv.mse <= cv.mse[min]+cv.se[min]))
  lambda.1se <- lambda[min.1se]
  
  
  Lmin <- which.min(cv.Lmse)
  lambda.Lmin <- lambda[Lmin]
  min.L1se <- min(which(cv.Lmse <= cv.Lmse[Lmin]+cv.Lse[Lmin]))
  lambda.L1se <- lambda[min.L1se]
  
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
    out$se = cv.se
    out$min.1se = min.1se
    out$lambda.1se = lambda.1se
    if (all.horizon) {
      out$cv.mse.h <- cv.mse.h
    }
  } else if (type == "logrates") {
    out$cv.mse = cv.Lmse
    out$min = Lmin
    out$lambda.min = lambda.Lmin
    out$se = cv.Lse
    out$min.1se = min.L1se
    out$lambda.1se = lambda.L1se
    if (all.horizon) {
      out$cv.mse.h <- cv.Lmse.h
    }
  } 
 
  class(out) <- "cv.grpStMoMo"
  return(out)
}


#' Plot the cross-validation curve produced by cv.grpreg
#' 
#' Plots the cross-validation currve, upper and lower stamdard errpr curves,
#' as a function of the \code{lambda} values used.
#' 
#' @usage 
#' \method{plot}{cv.grpStMoMo}(x, sign.lambda = 1, ...)
#' 
#' @param x an object of class \code{"cv.grpStMoMo"} with the cross-validation results
#' of the group regularised fitting  of a stochastic mortality model.
#' @param sign.lambda Either plot against \code{log(lambda)} (default) or its
#' negative if \code{sign.lambda=-1}.
#' @param ... additional arguments to control graphical appearance.
#' See \code{\link[graphics]{plot}}.
#' 
#' @seealso \code{cv.grpStMoMo}.
#' 
#' @examples 
#' 
#' 
#' APCcvgrp <- cv.grpStMoMo(apc(link = "log-Gaussian"), h = 1, data = EWMaleData, 
#'                          ages.train = 55:89, nlambda = 25)
#' 
#' plot(APCcvgrp)
#' 
#' 
#' #Long computing times
#' \dontrun{
#' #Create big model
#' strikes <- seq(25,85,5)
#' bModel <- StMoMo(link = "log-Gaussian", staticAgeFun = TRUE, 
#'                  periodAgeFun = c("1", genPoly(2:10), genCall(strikes),
#'                                                    genPut(strikes)),
#'                  cohortAgeFun = "1")
#'                  
#' #Cross-validation for 20 year ahead forecast                  
#' EWcv.20 <- cv.grpStMoMo(bModel, h = 20, data = EWMaleData, 
#'                         ages.train = 20:89)
#'                         
#' plot(EWcv.20)
#' }                       
#' 
#' 
#' @export 
#' @method plot cv.grpStMoMo
plot.cv.grpStMoMo <- function(x, sign.lambda = 1, ...) {
  
  xlab <- expression(Log(lambda))
  if(sign.lambda<0)xlab=paste("-",xlab,sep="")
  ylab <- "Mean-Squared Error"
  upper <- x$cv.mse+x$se
  lower <- x$cv.mse-x$se
  
  plot(x = sign.lambda*log(x$lambda), y = x$cv.mse, ylim = range(upper, lower),
       xlab = xlab, ylab = ylab, type = "n", ...)
  
  error.bars(sign.lambda*log(x$lambda),upper,lower,width=0.01,
             col="darkgrey")
  
  points(sign.lambda*log(x$lambda),x$cv.mse,pch=20,
         col="red")
  abline(v=sign.lambda*log(x$lambda.min),lty=3)
  abline(v=sign.lambda*log(x$lambda.1se),lty=3)
  
  
  
}


#' Plot error bars 
#' @keywords internal
error.bars <-
  function(x, upper, lower, width = 0.02, ...)
  {
    xlim <- range(x)
    barw <- diff(xlim) * width
    segments(x, upper, x, lower, ...)
    segments(x - barw, upper, x + barw, upper, ...)
    segments(x - barw, lower, x + barw, lower, ...)
    range(upper, lower)
  }
