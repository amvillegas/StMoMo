#' Cross validation for Stochastic Mortality Model
#' 
#' Perform cross validation by period for a Stochastic Mortality Model.
#' 
#' @param object an object of class \code{"StMoMo"} defining the 
#' stochastic mortality model.
#' 
#' @param h number of years for forecasting horizon.
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
#' @param returnY a logical value. If \code{TRUE}, \code{cv.StMoMo}
#' returns the fitted mortality rates for the cross validation
#' folds. 
#' 
#' @param type the type of the predicted values that the cross 
#' validation is performed with respect to. The alternatives are 
#' \code{"rates"} and \code{"logrates"}. If \code{"rates"}, \code{cv.StMoMo}
#' returns the mean squared error using the deviation between
#' realised and predicted mortality rates. If \code{"logrates"}, 
#' \code{cv.StMoMo} returns the mean squared error using the 
#' deviation between realised and predicted log mortality rates.    
#' 
#' @param verbose a logical value. If \code{TRUE} progress 
#' indicators are printed as the model is fitted. Set 
#' \code{verbose = FALSE} to silent the fitting and avoid 
#' progress messages.
#' 
#' @return A list with class \code{"cv.StMoMo"} with components:
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
#'   \item{rates}{ matrix of mortality rates (\eqn{D_{x,t}/E_{x,t}}) 
#'   used in the training.} 
#'    
#'   \item{logrates}{ matrix of log mortality rates (\eqn{log(D_{x,t}/E_{x,t})}) 
#'   used in the training.} 
#'   
#'   \item{cv.rates}{ if \code{returnY=TRUE}, matrix of predicted 
#'   mortality rates.} 
#'   
#'   \item{cv.mse}{ if \code{type="rates"}, mean squared error using 
#'   predicted mortality rates and realised mortality rates. If 
#'   \code{type="logrates"}, mean squared error using predicted log 
#'   mortality rates and realised log mortality rates.}  
#'  
#' @examples
#' 
#' # Lee-Carter model only for one year forecasting horizon based on log-rates
#' LCcv <- cv.StMoMo(lc(), h = 1, data = EWMaleData, ages.train = 55:89, type = "logrates")
#' 
#' # APC model using arguments Dxt, Ext, ages, years to pass fitting data, based on rates
#' APCcv <- cv.StMoMo(apc(), h = 10, Dxt = EWMaleData$Dxt, Ext = EWMaleData$Ext, 
#' ages = EWMaleData$ages, years = EWMaleData$years, ages.train = 55:89, type = "rates")
#'             
#' @export 
cv.StMoMo <- function(object,  h = NULL, data = NULL, Dxt = NULL, Ext = NULL, 
                     ages.train = NULL, years.train = NULL, ages = NULL, years = NULL, 
                     returnY = FALSE, type = c("logrates", "rates"), verbose = TRUE) {
  
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
  
  # Initialise out-of-sample predicted values
  cv.rates <- genWeightMat(ages = ages.train, years = years.train)*0
  
  # Determine number of test sets
  folds <- length(years.train) - h
  
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

     # Fit using fit.StMoMo
     fit <- fit(object, data = data, Dxt = Dxt, Ext = Ext, ages = ages, years = years, 
                ages.fit = ages.train, years.fit = years.train, wxt = wxtT, verbose = FALSE)

     # Extract and impute kappa parameters
     N <- object$N

     if (N > 0) {
       kt <- vector()
       for (j in 1:N) {
         k <- fit[["kt"]][j,]
         if (min(k, na.rm = TRUE) == max(k, na.rm = TRUE)){
           drift <- 0
         } else {
           arima.k <- tryCatch( mrwd(k), error = function( err ) FALSE, 
                                warning = function( err ) FALSE )
           if(!is.logical(arima.k)) {
             drift <- arima.k$drift
           } else {
             drift <- NA
             warning(paste("Time series model did not converge for test set",  i-1, "of", folds, "\n", sep = ""))
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
     if (is.null(object$cohortAgeFun)) {
       qhat <- predict(fit, years = years.train, kt = kt, type = "rates")
     } else {
       qhat <- predict(fit, years = years.train, kt = kt, gc = fit[["gc"]], type = "rates")
     }
    
    qhat[is.na(qhat)] <- 0
    cv.rates <- cv.rates + (qhat * (1 - wxtF))
    
    }

   # Prepare data output
   data <- fit$data
   Dxt <- fit$Dxt
   Ext <- fit$Ext
   qxt <- Dxt/Ext
   Lqxt <- log(qxt)
   
   # Prepare cv prediction output
   cv.rates[cv.rates == 0] <- NA
   cv.Lrates <- log(cv.rates)

   # Prepare cv error output
   cv.error <- (cv.rates - qxt)^2
   cv.mse <- mean(cv.error[!is.na(cv.error)])
   
   cv.Lerror <- (cv.Lrates - Lqxt)^2
   cv.Lmse <- mean(cv.Lerror[!is.na(cv.Lerror)])
   
   # Return output
   out <- list(model = object, data = data, Dxt = Dxt, Ext = Ext, rates = qxt, 
               logrates = Lqxt, call = match.call())
   
   if (returnY == TRUE) {
     out$cv.rates = cv.rates
   } 
   
   if (type == "rates") {
     out$cv.mse = cv.mse 
   } else if (type == "logrates") {
     out$cv.mse = cv.Lmse
   } 
   
   class(out) <- "cv.StMoMo"
   return(out)
      
 }