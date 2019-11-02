#' Create a Renshaw and Haberman (Lee-Carter with cohorts) mortality model
#' 
#' Utility function to initialise a \code{StMoMo} object representing a 
#' Renshaw and Haberman (Lee-Carter with cohorts) mortality model introduced
#' in Renshaw and Haberman (2006).
#' 
#' The created model is either a log-Poisson, a  logit-Binomial, or a log-Gaussian
#'  version of the Renshaw and Haberman model which has predictor structure   
#' \deqn{\eta_{xt} = \alpha_x + \beta^{(1)}_x\kappa_t + \beta^{(0)} \gamma_{t-x}.}
#' or
#' \deqn{\eta_{xt} = \alpha_x + \beta^{(1)}_x\kappa_t + \gamma_{t-x}.}
#' depending on the value of argument \code{cohortAgeFun}.
#'   
#' To ensure identifiability the following constraints are imposed
#' \deqn{\sum_t\kappa_t = 0, \sum_x\beta^{(1)}_x = 1, \sum_c\gamma_c = 0}
#' plus
#' \deqn{\sum_x\beta^{(0)}_x = 1}
#' if \code{cohortAgeFun = "NP"}
#' 
#' In addition, if \code{approxConst=TRUE} then the approximate 
#' identifiability constraint
#' \deqn{\sum_c (c-\bar{c})\gamma_c = 0}
#' is applied to improve the stability and robustness of the model 
#' (see Hunt and Villegas (2015)).
#'
#' By default \eqn{\beta^{(0)}_x = 1} as this model has shown to be more
#' stable (see Haberman and Renshaw (2011) and Hunt and Villegas (2015)).
#' 
#' @inheritParams StMoMo
#' @param cohortAgeFun defines the cohort age modulating parameter 
#'   \eqn{\beta_x^{(0)}}. It can take values: \code{"NP"} for a non-parametric age
#'   term or \code{"1"} for \eqn{\beta_x^{(0)}=1} (the default). 
#'   
#' @param approxConst defines if the approximate identifiability constraint of 
#' Hunt and Villegas (2015) is applied or not. If \code{TRUE}, the output object 
#' is of class \code{rh} and subsequent model fitting is performed with 
#' \code{\link{fit.rh}}. If \code{FALSE}, the output object is of class 
#' \code{StMoMo} and subsequent model fitting is performed with 
#' \code{\link{fit.StMoMo}}. This functionality is not yet implemented for the
#' \code{"log-Gaussian"} link.
#' 
#' @return An object of class \code{"StMoMo"} or \code{"rh"}.
#' 
#' @seealso \code{\link{fit.rh}}, \code{\link{StMoMo}}, \link{lc}, \link{apc}
#'  
#' @references
#'
#' Haberman, S., & Renshaw, A. (2011). A comparative study of parametric 
#' mortality projection models. Insurance: Mathematics and Economics, 
#' 48(1), 35-55. 
#' 
#' Hunt, A., & Villegas, A. M. (2015). Robustness and convergence in the 
#' Lee-Carter model with cohorts. Insurance: Mathematics and Economics, 
#' 64, 186-202. 
#' 
#' Renshaw, A. E., & Haberman, S. (2006). A cohort-based extension to the 
#' Lee-Carter model for mortality reduction factors. 
#' Insurance: Mathematics and Economics, 38(3), 556-570.
#' 
#' @examples 
#' 
#' LCfit <-  fit(lc(), data = EWMaleData, ages.fit = 55:89)
#' wxt <- genWeightMat(55:89,  EWMaleData$years, clip = 3)
#' RHfit <- fit(rh(), data = EWMaleData, ages.fit = 55:89, wxt = wxt, 
#'              start.ax = LCfit$ax, start.bx = LCfit$bx, start.kt = LCfit$kt)
#' plot(RHfit)
#' 
#' #Impose approximate constraint as in Hunt and Villegas (2015)    
#' \dontrun{
#' RHapprox <- rh(approxConst = TRUE)
#' RHapproxfit <- fit(RHapprox, data = EWMaleData, ages.fit = 55:89, 
#'                     wxt = wxt)
#' plot(RHapproxfit) 
#' }
#' 
#' 
#' @export
rh <- function(link = c("log", "logit", "log-Gaussian"), cohortAgeFun = c("1", "NP"), 
               approxConst = FALSE) {
  link <- match.arg(link)
  cohortAgeFun <- match.arg(cohortAgeFun)
  constRHgeneral <- function(ax, bx, kt, b0x, gc, wxt, ages, cohortAgeFun, approxConst) {
    #\sum k[t] = 0
    c1 <- mean(kt[1, ], na.rm = TRUE) 
    ax <- ax + c1 * bx[, 1]
    kt[1, ] <- kt[1, ] - c1
    #\sum b[x, 1] = 0
    c2 <- sum(bx[, 1], na.rm = TRUE)
    bx[, 1] <- bx[, 1] / c2
    kt[1, ] <- kt[1, ] * c2
    #\sum g[c] = 0
    c3 <- mean(gc, na.rm = TRUE)
    ax <- ax + c3 * b0x
    gc <- gc - c3
    #\sum b0[x] = 1
    if (cohortAgeFun == "NP") {
      c4 <- sum(b0x, na.rm = TRUE)
      b0x <- b0x / c4
      gc <- gc * c4
    }
    
    list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)
  }
  constRH <- function(ax, bx, kt, b0x, gc, wxt, ages) {
    constRHgeneral(ax, bx, kt, b0x, gc, wxt, ages, cohortAgeFun, approxConst) 
  }
  out <- StMoMo(link = link, staticAgeFun = TRUE, periodAgeFun = "NP",
         cohortAgeFun = cohortAgeFun, constFun = constRH)
  if (!is.logical(approxConst)) {
    stop( "approxConst should be a logical variable")
  }
  if (approxConst && link == "log-Gaussian") {
    stop( "approxConst fitting not yet available for the log-Gaussian link")
  }
  out$approxConst <- approxConst
  if(approxConst)  class(out) <- c("rh", "StMoMo")
  out
}



#' Fit a Renshaw and Haberman (Lee-Carter with cohorts) mortality model
#' 
#' Fit a Renshaw and Haberman (Lee-Carter with cohorts) mortality model
#' using the iterative Newton-Raphson procedure presented in Algorithm 1
#' of Hunt and Villegas (2015). This approach helps solve the 
#' well-known robustness and converges issues of the Lee-Carter model 
#' with cohort-effects.
#' 
#' @inheritParams fit.StMoMo
#' @param object an object of class \code{"rh"} created with function 
#' \code{\link{rh}}.
#' @param tolerance a positive numeric value specifying the tolerance 
#' level for convergence.
#' @param iterMax  a positive integer specifying the maximum number of 
#' iterations to perform.
#' @param ... arguments to be passed to or from other methods.
#' 
#' @return 
#' 
#'   \item{model}{ the object of class \code{"rh"} defining the fitted 
#'   stochastic mortality model.}
#'   
#'   \item{ax}{ vector with the fitted values of the static age function 
#'   \eqn{\alpha_x}. If the model does not have a static age function or 
#'   failed to fit this is set to \code{NULL}.}
#'     
#'   \item{bx}{ matrix with the values of the period age-modulating functions 
#'   \eqn{\beta_x^{(i)}, i=1, ..., N}. If the \eqn{i}-th age-modulating 
#'   function is non-parametric (e.g. as in the Lee-Carter model) 
#'   \code{bx[, i]} contains the estimated values. If the model does not have 
#'   any age-period terms (i.e. \eqn{N=0}) or failed to fit this is set to 
#'   \code{NULL}.}
#'   
#'   \item{kt}{ matrix with the values of the fitted period indexes 
#'   \eqn{\kappa_t^{(i)}, i=1, ..., N}. \code{kt[i, ]} contains the estimated 
#'   values of the \eqn{i}-th period index. If the model does not have any 
#'   age-period terms (i.e. \eqn{N=0}) or failed to fit this is set to 
#'   \code{NULL}.}
#'   \item{b0x}{ vector with the values of the cohort age-modulating function 
#'   \eqn{\beta_x^{(0)}}. If the age-modulating function is non-parametric 
#'   \code{b0x} contains the estimated values. If the model does not have a 
#'   cohort effect or failed to fit this is set to \code{NULL}.}
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
#'   \item{oxt}{ matrix of known offset values used in the fitting.}
#'   
#'   \item{wxt}{ matrix of 0-1 weights used in the fitting.}
#'   
#'   \item{ages}{ vector of ages used in the fitting.}
#'   
#'   \item{years}{ vector of years used in the fitting.}
#'   
#'   \item{cohorts}{ vector of cohorts used in the fitting.}
#'   
#'   \item{fittingModel}{ output from the iterative fitting algorithm.}
#'   
#'   \item{loglik}{ log-likelihood of the model. If the fitting failed to 
#'   converge this is set to \code{NULL}.}
#'   
#'   \item{deviance}{ deviance of the model. If the fitting failed to 
#'   converge this is set to \code{NULL}.}
#'  
#'   \item{npar}{ effective number of parameters in the model. If the fitting
#'   failed to converge this is set to \code{NULL}.}
#'    
#'    \item{nobs}{ number of observations in the model fit. If the fitting
#'    failed to converge this is set to \code{NULL}.}
#'
#'    \item{fail}{ \code{TRUE} if a model could not be fitted and 
#'    \code{FALSE} otherwise.}    
#'            
#'    \item{conv}{ \code{TRUE} if the model fitting converged and 
#'    \code{FALSE} if it didn't.}
#' 
#' @references 
#' 
#' Hunt, A., & Villegas, A. M. (2015). Robustness and convergence in the 
#' Lee-Carter model with cohorts. Insurance: Mathematics and Economics, 
#' 64, 186-202. 
#' 
#' @examples 
#'
#' LCfit <-  fit(lc(), data = EWMaleData, ages.fit = 55:89)
#' wxt <- genWeightMat(55:89,  EWMaleData$years, clip = 3)
#' RHfit <- fit(rh(), data = EWMaleData, ages.fit = 55:89, 
#'              wxt = wxt, start.ax = LCfit$ax,
#'              start.bx = LCfit$bx, start.kt = LCfit$kt)
#' plot(RHfit)
#'  
#' #Impose approximate constraint as in Hunt and Villegas (2015)    
#' \dontrun{
#' RHapprox <- rh(approxConst = TRUE)
#' RHapproxfit <- fit(RHapprox, data = EWMaleData, ages.fit = 55:89, 
#'                     wxt = wxt)
#' plot(RHapproxfit) 
#' }
#' 
#' @export 
fit.rh <- function(object, data = NULL, Dxt = NULL, Ext = NULL,
                       ages = NULL, years = NULL, ages.fit = NULL, 
                       years.fit = NULL, oxt = NULL, wxt = NULL, 
                       start.ax = NULL, start.bx = NULL, start.kt = NULL,
                       start.b0x = NULL, start.gc = NULL, verbose = TRUE,
                       tolerance = 1e-04, iterMax = 10000,
                       ...) {
  
  #Hack to remove notes in CRAN check
  x <- NULL
  
  # Select data from data or from Dxt, Ext, ages, years
  
  if(!is.null(data)) {
    if (class(data) != "StMoMoData")
      stop("Argument data needs to be of class StMoMoData.")
    Dxt <- data$Dxt
    Ext <- data$Ext
    ages <- data$ages
    years <- data$years
  } else {
    if (is.null(Dxt) || is.null(Ext))
      stop("Either argument data or arguments Dxt and Ext need to be provided.")
    if (is.null(ages)) ages <- 1:nrow(Dxt)
    if (is.null(years)) years <- 1:ncol(Dxt)
    data <- structure(list(Dxt = Dxt, Ext = Ext, ages = ages, years = years, 
                           type = ifelse(object$link == "log", "central", "initial"), 
                           series = "unknown", label = "unknown"), class = "StMoMoData")
  }
  if (is.null(ages.fit)) ages.fit <- ages
  if (is.null(years.fit)) years.fit <- years
  
  
  # Construct fitting data
  
  nAges <- length(ages)
  nYears <- length(years)    
  Dxt <- as.matrix(Dxt)
  if (nrow(Dxt) != nAges ||  ncol(Dxt) != nYears) {
    stop( "Mismatch between the dimension of Dxt and the 
          number of years or ages")
  }
  rownames(Dxt) <- ages
  colnames(Dxt) <- years
  
  Ext <- as.matrix(Ext)
  if (nrow(Ext) != nAges ||  ncol(Ext) != nYears) {
    stop( "Mismatch between the dimension of Ext and the 
          number of years or ages")
  }
  rownames(Ext) <- ages
  colnames(Ext) <- years  
  if (data$type == "central" && object$link == "logit")
    warning( "logit-Binomial model fitted to central exposure data\n")
  if (data$type == "initial" && object$link == "log")
    warning( "log-Poisson model fitted to initial exposure data\n")
  
  #Extract the specific ages and years for fitting 
  if (length(ages.fit) != length(which(ages.fit %in% ages))) {
    stop( "ages.fit is not a subset of ages")
  }
  if (length(years.fit) != length(which(years.fit %in% years))) {
    stop( "years.fit is not a subset of years")
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
  
  if (is.null(oxt)) {
    oxt <- matrix(0, nrow = nAges, ncol = nYears)
    rownames(oxt) <- ages
    colnames(oxt) <- years
  } else {
    oxt <- matrix(oxt, nrow = nAges, ncol = nYears)
    rownames(oxt) <- ages
    colnames(oxt) <- years
  }
  
  if (is.null(wxt)) {
    wxt <- matrix(1, nrow = nAges, ncol = nYears)
  } else {
    wxt <- as.matrix(wxt)
    if ( nrow(wxt) != nAges ||  ncol(wxt) != nYears) {
      stop( "Mismatch between the dimension of the weigth matrix wxt and the 
            number of fitting years or ages")
    }    
  } 
  rownames(wxt) <- ages
  colnames(wxt) <- years
  if (any(Ext <= 0)) { #Non-positive exposures
    indExt <- (Ext <= 0)
    wxt[indExt] <- 0
    warning(paste("StMoMo: ", sum(indExt), " data points have 
                  non-positive exposures and have been zero weighted\n", 
                  sep = ""))
  }
  if (any(is.na(Dxt / Ext))) { #Missing values
    indqxt <- is.na(Dxt / Ext)
    wxt[indqxt] <- 0
    warning(paste("StMoMo: ", sum(indqxt), 
                  " missing values which have been zero weighted\n", sep = ""))
  }
  
  # Identify the data ages, years or cohorts with 0 weight 
  wxtDF <- (reshape2::melt(wxt, value.name = "w", varnames = c("x", "t")))
  wxtDF <- transform(wxtDF, c = t - x)
  wxTemp <- aggregate(data = wxtDF, w ~ x, FUN = sum)  
  zeroWeigthAges <- as.character(wxTemp$x[which((wxTemp$w <= 0))])
  indX <- which((!(as.character(ages) %in% zeroWeigthAges)))
  wtTemp <- aggregate(data = wxtDF, w ~ t, FUN = sum)  
  zeroWeigthYears <- as.character(wtTemp$t[which((wtTemp$w <= 0))])  
  indT <- which((!(as.character(years) %in% zeroWeigthYears)))
  wcTemp <- aggregate(data = wxtDF, w ~ c, FUN = sum)  
  zeroWeigthCohorts <- as.character(wcTemp$c[which((wcTemp$w <= 0))])      
  indC <- which((!(as.character(cohorts) %in% zeroWeigthCohorts)))
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
  
  
  ### Set starting values
  if (is.null(start.ax)) {
    if (object$link == "log") { 
      logmxt <- log(Dxt / Ext)
      start.ax <- rowSums((logmxt-oxt) * wxt) / rowSums(wxt)
    }else{
      logoddsxt <- log(Dxt / (Ext - Dxt))
      start.ax <-  rowSums((logoddsxt - oxt) * wxt) / rowSums(wxt)
    }
  } else {
    if (length(start.ax) != nAges) {
      stop( "Mismatch between the number of ages and start.ax")
    }
  }
  names(start.ax) <- ages
  start.ax[zeroWeigthAges] <- NA
  
  if (is.null(start.bx)) {
    start.bx <- array(1 / nAges, c(nAges, 1), dimnames = list(ages, 1))
  } else {
    if (nrow(start.bx) != nAges) {
      stop( "Mismatch between the number of ages and start.bx.")
    }
    if (ncol(start.bx) != 1) {
      stop( "Mismatch between the number of age/period terms and start.bx.")
    }
  }
  dimnames(start.bx) <- list(ages, 1)
  start.bx[zeroWeigthAges, ] <- NA
  
  if (is.null(start.kt)) {
    start.kt <- array(0, c(1, nYears), dimnames = list(1, years))
  } else {
    if (ncol(start.kt) != nYears) {
      stop( "Mismatch between the number of years and start.kt.")
    }
    if (nrow(start.kt) != 1) {
      stop( "Mismatch between the number of age/period terms and start.kt.")
    }
  }
  dimnames(start.kt) <- list(1, years)
  start.kt[, zeroWeigthYears] <- NA
    
  if (is.null(start.b0x) || object$cohortAgeFun == "1" ) {
    start.b0x <- rep(1, nAges)
  } else {
    if (length(start.b0x) != nAges && object$cohortAgeFun == "NP") {
      stop( "Mismatch between the number of ages and start.b0x.")
    }
  }
  names(start.b0x) <- ages
  if (object$cohortAgeFun != "1" )  start.b0x[zeroWeigthAges] <- NA
  
  if (is.null(start.gc)) {
    start.gc <- rep(0, nCohorts)
  } else {
    if (length(start.gc) != nCohorts) {
      stop( "Mismatch between the number of cohorts and start.gc.")
    }
  }
  names(start.gc) <- cohorts
  start.gc[zeroWeigthCohorts] <- NA
  
  
  out <- list(model = object, ax = start.ax, bx = start.bx, 
              kt = start.kt, b0x = start.b0x, gc = start.gc, 
              data = data, Dxt = Dxt, Ext = Ext, oxt = oxt , 
              wxt = wxt, ages = ages, years = years, 
              cohorts = cohorts, fittingModel = NULL, 
              loglik = NA, deviance = NA,   
              npar = (2 + (object$cohortAgeFun == "NP")) * 
                (nAges - length(zeroWeigthAges)) + 
                (nYears - length(zeroWeigthYears)) + 
                (nCohorts - length(zeroWeigthCohorts)) - 
                (3 + object$approxConst + (object$cohortAgeFun == "NP")),
              nobs = sum(wxt),  conv = NA, fail = FALSE, 
              call  = match.call())   
  
  class(out) <- c("fitrh", "fitStMoMo")
  
  
  ### Do maximisation iterations
  
  #auxiliary variables for itertions
  matbx <- array(NA, dim = c(nAges, nYears))
  matbx2 <- array(NA, dim = c(nAges, nYears))
  matkt <- array(NA, dim = c(nYears, nAges))
  matkt2 <- array(NA, dim = c(nYears, nAges))
  matb0x <- array(NA, dim = c(nAges, nYears))
  matb0x2 <- array(NA, dim = c(nAges, nYears))
  diagonal <- col(Ext) - row(Ext)
  sumDiags <- function(A) sapply(split(A, diagonal), FUN = sum, na.rm = TRUE)
  matgc <- array(NA, dim = c(nAges, nYears))
  matgc2 <- array(NA, dim = c(nAges, nYears))
  gcIndex <- as.numeric(gl(nYears, nAges)) + seq(nAges - 1, 0)
  xbar <- 1:nAges
  names(xbar) <- ages
  xbar[zeroWeigthAges] <- NA
  xbar <- xbar - mean(xbar, na.rm = TRUE)
  tbar <- 1:nYears
  names(tbar) <- years
  tbar[zeroWeigthYears] <- NA
  tbar <- tbar - mean(tbar, na.rm = TRUE)
  sumtbar2 <- sum(tbar ^ 2, na.rm = TRUE);
  cbar <- 1:nCohorts
  names(cbar) <- cohorts
  cbar[zeroWeigthCohorts] <- NA
  cbar <- cbar - mean(cbar, na.rm = TRUE)
  sumcbar2 <- sum(cbar ^ 2, na.rm = TRUE);
  
  #initialisation
  mxt <- Dxt / Ext
  mhatxt <- fitted.fitStMoMo(out, type = "rates")
  Dhatxt <- Ext * mhatxt
  out$loglik <- ifelse(object$link == "log",
                       computeLogLikPoisson(obs = Dxt, fit = Dhatxt,
                                            weight = wxt),
                       computeLogLikBinomial(obs = mxt, fit = mhatxt,
                                             exposure = Ext, weight = wxt))
  
  absErrL  <- 10
  iter <- 0
  if (verbose) {
    cat("StMoMo: Start fitting with iterative Newton-Raphson\n")
    cat("Running iterations")
  }
  
  
  while(absErrL > tolerance & iter < iterMax){
    #update ax
    sumNum <- rowSums(wxt * (Dxt - Dhatxt), na.rm = TRUE)
    sumDen <- switch(object$link, 
                     log = rowSums(wxt * Dhatxt, na.rm = TRUE),
                     logit = rowSums(wxt * Dhatxt * (1 - mhatxt), 
                                     na.rm = TRUE))
    out$ax[indX] <- out$ax[indX] + sumNum[indX] / sumDen[indX]
    
    #update kt
    mhatxt <- fitted.fitStMoMo(out, type = "rates")
    Dhatxt <- Ext * mhatxt
    matbx[] <- out$bx[, 1]
    matbx2[] <- out$bx[, 1] ^ 2
    sumNum <- colSums(wxt * (Dxt - Dhatxt) * matbx, na.rm = TRUE)
    sumDen <- switch(object$link, 
                     log = colSums(wxt * Dhatxt * matbx2, na.rm = TRUE),
                     logit = colSums(wxt * Dhatxt * matbx2 * (1 - mhatxt), 
                                     na.rm = TRUE))
    out$kt[1, indT] <- out$kt[1, indT] + sumNum[indT] / sumDen[indT]
    
    #update bx
    mhatxt <- fitted.fitStMoMo(out, type = "rates")
    Dhatxt <- Ext * mhatxt
    matkt[] <- out$kt[1, ]
    matkt2[] <- out$kt[1, ] ^ 2
    sumNum <- rowSums(wxt * (Dxt - Dhatxt) * t(matkt), na.rm = TRUE)
    sumDen <- switch(object$link,
                     log = rowSums(wxt * Dhatxt * t(matkt2), na.rm = TRUE),
                     logit = rowSums(wxt * Dhatxt * t(matkt2) * (1 - mhatxt), 
                                     na.rm = TRUE))
    out$bx[indX, 1] <- out$bx[indX, 1] + sumNum[indX] / sumDen[indX]
    
    #update gc
    mhatxt <- fitted.fitStMoMo(out, type = "rates")
    Dhatxt <- Ext * mhatxt
    if (object$cohortAgeFun == "NP"){
      matb0x[] <- out$b0x
      matb0x2[] <- out$b0x ^ 2
      sumNum <- sumDiags(wxt * (Dxt - Dhatxt) * matb0x)
      sumDen <- switch(object$link,
                       log = sumDiags(wxt * Dhatxt * matb0x2),
                       logit = sumDiags(wxt * Dhatxt * matb0x2 * (1 - mhatxt)))
    } else {
      sumNum <- sumDiags(wxt * (Dxt - Dhatxt))
      sumDen <- switch(object$link,
                       log = sumDiags(wxt * Dhatxt),
                       logit = sumDiags(wxt * Dhatxt * (1 - mhatxt)))
    }
    out$gc[indC] <- out$gc[indC] + sumNum[indC] / sumDen[indC]
    
    #update b0x
    if (object$cohortAgeFun == "NP"){
      mhatxt <- fitted.fitStMoMo(out, type = "rates")
      Dhatxt <- Ext * mhatxt
      matgc[] <- out$gc[gcIndex]
      matgc2[] <- (out$gc ^ 2)[gcIndex]
      sumNum <- rowSums(wxt * (Dxt - Dhatxt) * matgc, na.rm = TRUE)
      sumDen <- switch(object$link,
                       log = rowSums(wxt * Dhatxt * matgc2, na.rm = TRUE),
                       logit = colSums(wxt * Dhatxt * matgc2 * (1 - mhatxt), na.rm = TRUE))
      out$b0x[indX] <- out$b0x[indX] + sumNum[indX] / sumDen[indX]
      mhatxt <- fitted.fitStMoMo(out, type = "rates")
      Dhatxt <- Ext * mhatxt
    }
    
    #apply constraints
    constPar<- object$constFun(out$ax, out$bx, out$kt, 
                               out$b0x, out$gc, wxt, ages)
    out$ax <- constPar$ax
    out$bx <- constPar$bx
    out$kt <- constPar$kt
    out$b0x <- constPar$b0x
    out$gc <- constPar$gc
    
    #Apply approx constraint
    if (object$approxConst) {
      K <- sum(tbar[indT] * out$kt[1, indT]) / sumtbar2
      e <- -sum(cbar[indC] * out$gc[indC]) / sumcbar2
      out$ax[indX] <- out$ax[indX] + e * out$b0x[indX] * xbar[indX]
      out$bx[indX, 1] <- K / (K - e) * out$bx[indX, 1] - e / (K - e) * out$b0x[indX]
      out$kt[, indT] <- (K - e) / K * out$kt[, indT]
      out$gc[indC] <- out$gc[indC] + e * cbar[indC]
    }
    
    #update iteration
    oldLoglik <- out$loglik
    mhatxt <- fitted.fitStMoMo(out, type = "rates")
    Dhatxt <- Ext * mhatxt
    out$loglik <- ifelse(object$link == "log",
                         computeLogLikPoisson(obs = Dxt, fit = Dhatxt,
                                              weight = wxt),
                         computeLogLikBinomial(obs = mxt, fit = mhatxt,
                                               exposure = Ext, weight = wxt))
    absErrL <- abs(out$loglik - oldLoglik)    
    iter <- iter + 1    
    if (verbose ) cat(".")
    
    
  }
  
  if (verbose ){
    cat('\nDone\n')
    cat('StMoMo: Finish fitting with iterative Newton-Raphson\n')
  }
  
  
  ### Prepare output
  out$deviance <- ifelse(object$link == "log",
                       computeDeviancePoisson(obs = Dxt, fit = Dhatxt,
                                            weight = wxt),
                       computeDevianceBinomial(obs = mxt, fit = mhatxt,
                                             exposure = Ext, weight = wxt))
  out$conv <- (iter < iterMax)
  if (!out$conv) warning( "Iterative Newton-Raphson procedure did not converge.
  See fittingModel for information about the last iteration.\n")
  if (is.na(out$loglik)) {
    warning( "Iterative Newton-Raphson procedure failed.
  See fittingModel for information about the last iteration.\n")
    out$fail <- TRUE
    out$conv <- FALSE
    
  }
  
  out$fittingModel <- list(tolerance = tolerance, iterMax = iterMax, iter = iter, absErr = absErrL) 
  
  out
  
}

