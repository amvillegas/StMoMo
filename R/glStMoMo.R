#' Fit a Group Lasso Stochastic Mortality Model
#' 
#' Fit a Group Lasso Stochastic Mortality Model to a given data set. The fitting is done using package \code{grpreg}.
#' 
#' Fitting is done using function \code{\link[grpreg]{grpreg}} within package 
#' \code{grpreg}. This is achieved by minimising the residual sum of squares subject 
#' to the minimax concave penalty (MCP). This function is not yet supported for  
#' logit-Binomial or log-Poisson models. Ages and years in the data should be of 
#' type numeric. Data points with zero exposure are assigned a zero weight 
#' and are ignored in the fitting process. Similarly, \code{NA} are assigned a
#' zero weight and ignored in the fitting process.
#' 
#' @param object an object of class \code{"StMoMo"} defining the stochastic
#' mortality model, with only parametric age-period and age-cohort terms.
#' 
#' @param lambda a user supplied sequence of lambda values. If left unspecified,
#' the function will compute a grid of lambda values of length \code{nlambda} 
#' that range uniformly on the log scale. 
#' 
#' @param nlambda number of lambda values. Default is 50.
#' 
#' @param index optional vector describing the grouping of the coefficients. If
#' there are coefficients to be included in the model without being penalised, assign
#' them to group 0 (or "0"). If this is not provided, the grouping is automatically
#' determined and all groups are penalised.
#' 
#' @param data an optional object of type StMoMoData containing information on
#' deaths and exposures to be used for fitting the model. This is typically created 
#' with  function \code{\link{StMoMoData}}. If this is not provided then the fitting 
#' data is taken from arguments, \code{Dxt}, \code{Ext}, \code{ages}, \code{years}.
#' 
#' @param Dxt optional matrix of deaths data.
#' 
#' @param Ext optional matrix of observed exposures of the same dimension of 
#' \code{Dxt}. 
#' 
#' @param ages optional vector of ages corresponding to rows of \code{Dxt} and 
#' \code{Ext}. 
#' 
#' @param years optional vector of years corresponding to rows of \code{Dxt} and 
#' \code{Ext}. 
#' 
#' @param ages.fit optional vector of ages to include in the fit. Must be a 
#' subset of \code{ages}. 
#' 
#' @param years.fit optional vector of years to include in the fit. Must be a 
#' subset of \code{years}. 
#' 
#' @param wxt optional matrix of 0-1 weights to be used in the fitting process. 
#' This can be used, for instance, to zero weight some cohorts in the data.
#' See \code{\link{genWeightMat}} which is a helper function for defining 
#' weighting matrices.
#' 
#' @param verbose a logical value. If \code{TRUE} progress indicators are 
#' printed as the model is fitted. Set \code{verbose = FALSE} to silent the 
#' fitting and avoid progress messages.
#' 
#' @return A list with class \code{"glStMoMo"} with components:
#'   
#'   \item{model}{ the object of class \code{"StMoMo"} defining the fitted 
#'   stochastic mortality model.}
#'   
#'   \item{beta}{ list of fitted coefficients for each value of lambda. Each
#'   entry contains \code{ax}, \code{kt}, and \code{gc}}. 
#'   
#'   \item{ax}{ vector with the fitted values of the static age function 
#'   \eqn{\alpha_x}. If the model does not have a static age function or 
#'   failed to fit this is set to \code{NULL}.}
#'   
#'   \item{kt}{ matrix with the values of the fitted period indexes 
#'   \eqn{\kappa_t^{(i)}, i=1, ..., N}. \code{kt[i, ]} contains the estimated 
#'   values of the \eqn{i}-th period index. If the model does not have any 
#'   age-period terms (i.e. \eqn{N=0}) or failed to fit this is set to 
#'   \code{NULL}.}
#'     
#'   \item{gc}{ vector with the fitted cohort index \eqn{\gamma_{c}}.
#'   If the model does not have a cohort effect or failed to fit this is set
#'   to \code{NULL}.}
#'   
#'   \item{data}{ StMoMoData object provided for fitting the model.}
#'   
#'   \item{Dxt}{ matrix of deaths used in the fitting.}
#'   
#'   \item{Ext}{ matrix of exposures used in the fitting.}
#'   
#'   \item{wxt}{ matrix of 0-1 weights used in the fitting.}
#'   
#'   \item{ages}{ vector of ages used in the fitting.}
#'   
#'   \item{years}{ vector of years used in the fitting.}
#'   
#'   \item{cohorts}{ vector of cohorts used in the fitting.}
#'   
#'   \item{loglik}{ vector of log-likelihoods of fitted model.}
#'   
#'   \item{lambda}{ vector of lambda values used in the fitting.}
#'     
#' @examples
#' glCBD <- glStMoMo(CBD, data = EWMaleData, ages.fit = 20:89, years.fit = 1965:2010)    
#' 
#' @export 
glStMoMo <- function(object, lambda = NULL, nlambda = 50, data = NULL, Dxt = NULL, Ext = NULL, ages = NULL, years = NULL, ages.fit = NULL, years.fit = NULL, wxt = NULL, index = NULL, verbose = TRUE) {
  
  # Group lasso not yet supported for Binomial/Poisson models, and parametric models

  if (object$link == "logit" || object$link == "log")
    stop("Group lasso not yet supported for logit-Binomial or log-Poisson models.\n")
  if ("NP" %in% object$periodAgeFun || "NP" %in% object$cohortAgeFun)
    stop("Group lasso not yet supported for parametric models.\n")
  
  # Select data from data or from Dxt, Ext, ages, years
  
  if(!is.null(data)) {
    if (class(data) != "StMoMoData")
      stop("Argument data needs to be of class StMoMoData.")
    if (data$type == "initial")
      warning("Normal model fitted to initial exposure data\n")
    Dxt <- data$Dxt
    Ext <- data$Ext
    ages <- data$ages
    years <- data$years
  } else {
    if (is.null(Dxt) || is.null(Ext))
      stop("Either argument data or arguments Dxt and Ext need to be provided.")
    if (is.null(ages)) ages <- 1:nrow(Dxt)
    if (is.null(years)) years <- 1:ncol(Dxt)
    data <- structure(list(Dxt = Dxt, Ext = Ext, ages = ages, years = years, type = "central", series = "unknown", label = "unknown"), class = "StMoMoData")
  }
  if (is.null(ages.fit)) ages.fit <- ages
  if (is.null(years.fit)) years.fit <- years
  
  # Construct fitting data
  
  nAges <- length(ages)
  nYears <- length(years)    
  Dxt <- as.matrix(Dxt)
  if (nrow(Dxt) != nAges ||  ncol(Dxt) != nYears) {
    stop("Mismatch between the dimension of Dxt and the number of years or ages")
  }
  rownames(Dxt) <- ages
  colnames(Dxt) <- years
  
  Ext <- as.matrix(Ext)
  if (nrow(Ext) != nAges ||  ncol(Ext) != nYears) {
    stop("Mismatch between the dimension of Ext and the number of years or ages")
  }
  rownames(Ext) <- ages
  colnames(Ext) <- years  

  # Extract the specific ages and years for fitting 
  
  if (length(ages.fit) != length(which(ages.fit %in% ages))) {
    stop("ages.fit is not a subset of ages")
  }
  if (length(years.fit) != length(which(years.fit %in% years))) {
    stop("years.fit is not a subset of years")
  }  
  Dxt <- Dxt[which(ages %in% ages.fit), which(years %in% years.fit)]
  Ext <- Ext[which(ages %in% ages.fit), which(years %in% years.fit)]  
  ages <- ages.fit
  years <- years.fit
  
  # Construct fitting data
  
  nAges <- length(ages)
  nYears <- length(years)  
  cohorts <- (years[1] - ages[nAges]):(years[nYears] - ages[1])
  nCohorts <- length(cohorts)  
  
  fitDataD <- (reshape2::melt(Dxt, value.name = "D", varnames = c("x", "t")))
  fitDataE <- (reshape2::melt(Ext, value.name = "E", varnames = c("x", "t")))
  fitData <- merge(fitDataD, fitDataE)      
  fitData <- transform(fitData, c = t - x, response = log(D/E))      
  
  if (is.null(wxt)) {
    wxt <- matrix(1, nrow = nAges, ncol = nYears)
  } else {
    wxt <- as.matrix(wxt)
    if (nrow(wxt) != nAges ||  ncol(wxt) != nYears) {
      stop("Mismatch between the dimension of the weight matrix wxt and the number of fitting years or ages")
    }    
  }
  rownames(wxt) <- ages
  colnames(wxt) <- years
  if (any(Ext <= 0, na.rm = TRUE)) { #Non-positive exposures
    indExt <- (Ext <= 0)
    wxt[indExt] <- 0
    warning(paste("StMoMo: ", sum(indExt), " data points have non-positive exposures and have been zero weighted\n", sep = ""))
  }
  if (any(is.na(Dxt / Ext))) { #Missing values
    indqxt <- is.na(Dxt / Ext)
    wxt[indqxt] <- 0
    warning(paste("StMoMo: ", sum(indqxt), " missing values which have been zero weighted\n", sep = ""))
  }
  
  fitDataW <- (reshape2::melt(wxt, value.name = "w", varnames = c("x", "t")))
  fitData <- merge(fitData, fitDataW)
  
  # Data for age-parametric terms
  N <- object$N
  if (N > 0) {
    for (i in 1:N) {      
      if (is.function(object$periodAgeFun[[i]])) {
        f <- function(x) object$periodAgeFun[[i]](x, ages)
        fitData$B__ <- apply(as.array(fitData$x), MARGIN = 1, FUN = f)
        names(fitData)[names(fitData) == "B__"] <- paste("B", i, sep = "")
      }      
    }
  }
  
  # Data for age-cohort term  
  if (is.function(object$cohortAgeFun)) {
    f <- function(x) object$cohortAgeFun(x, ages)
    fitData$B0 <- apply(as.array(fitData$x), MARGIN = 1, FUN = f)    
  } 
  
  # Remove from the data ages, years or cohorts with 0 weight 

  wxTemp <- aggregate(data = fitData, w ~ x, FUN = sum)  
  zeroWeigthAges <- as.character(wxTemp$x[which((wxTemp$w <= 0))])
  wtTemp <- aggregate(data = fitData, w ~ t, FUN = sum)  
  zeroWeigthYears <- as.character(wtTemp$t[which((wtTemp$w <= 0))])  
  wcTemp <- aggregate(data = fitData, w ~ c, FUN = sum)  
  zeroWeigthCohorts <- as.character(wcTemp$c[which((wcTemp$w <= 0))])      
  fitData <- subset(fitData, w > 0)
  
  if (verbose) {
    if (length(zeroWeigthAges) > 0) {
      cat("StMoMo: The following ages have been zero weigthed:", 
          zeroWeigthAges, "\n")
    }
    if (length(zeroWeigthYears) > 0) {
      cat("StMoMo: The following years have been zero weigthed:", 
          zeroWeigthYears, "\n")
    }
    if (length(zeroWeigthCohorts) > 0) {
      cat("StMoMo: The following cohorts have been zero weigthed:", 
          zeroWeigthCohorts, "\n")
    }    
  }  

  # Determine design matrix and index 

  designmatrix <- gnm(formula = as.formula(object$gnmFormula), data = fitData, family = gaussian, method = "model.matrix")
  if (is.null(index)) {
    index <- attributes(designmatrix)$assign
  }
  
  # Fit using grpreg
  
  if (verbose) cat("StMoMo: Start fitting with grpreg\n")
  
  if (is.null(lambda)) {
    fit <- grpreg(designmatrix, fitData$response, group = index, penalty = "grMCP", family = "gaussian", nlambda = nlambda)
  } else {
    fit <- grpreg(designmatrix, fitData$response, group = index, lambda = lambda, penalty = "grMCP", family = "gaussian")
    if (identical(lambda, fit$lambda)==FALSE) {
      warning(paste("Model fitted with lambda =", fit$lambda, "instead of", lambda, "\n"))
    }
  }
  
  if (verbose) cat("StMoMo: Finish fitting with grpreg\n")  
  
  ### Prepare output
  
  # Estimated coefficients
  
  lambda = fit$lambda
  beta <- vector("list", length(lambda)) 
  
  for (i in 1:length(lambda)) {
    
    # Extract un-adjusted coefficients
    coefGnmModel <- fit$beta[,i]
    tempfittedCoef <- extractCoefficientsFromGnm(object, fit$beta[,i], ages, years, cohorts, zeroWeigthAges, zeroWeigthYears, zeroWeigthCohorts)
    
    # Extract intercept
    intercept <- fit$beta[1,i]
    
    # Adjust coefficients for intercept
    fittedCoef <- glRemoveIntercept(object, tempfittedCoef, intercept, ages)
    
    # Update beta
    beta[[i]] <- fittedCoef
    
  }
  
  # Log likelihood
  
  loglik <- logLik(fit)

  out <- list(model = object, data = data, Dxt = Dxt, Ext = Ext, wxt = wxt, 
              ages = ages, years = years, cohorts = cohorts,
              lambda = lambda, index = index, beta = beta, loglik = loglik) 
  
  class(out) <- "fitGL"
  return(out)  

}