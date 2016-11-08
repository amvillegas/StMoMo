
#'Generate weight matrix
#'
#'Generates a weight matrix given a group of ages and years
#'and a set of cohorts which are to be given zero weight. This
#'is useful for excluding some data points when fitting a 
#'Stochastic Mortality Model (see \code{\link{fit.StMoMo}}).
#'
#'@param ages vector of ages.
#'@param years vector of years.
#'@param clip number of cohorts in the boundary to assign a zero 
#'weight. This can be be used to zero weigh some of the first and 
#'last cohorts in the data.
#'@param zeroCohorts other cohort for which a zero weight is to be assigned.
#'   
#'@return A 0-1 matrix with 0 for the zero-weighed cohorts.
#'
#'@seealso \code{\link{fit.StMoMo}}
#'
#'@examples
#' #Zero-weight the first three and last three cohorts
#' wxt1 <- genWeightMat(55:89,  EWMaleData$years, clip = 3)
#' APCfit1 <- fit(apc(), data = EWMaleData, ages.fit = 55:89, 
#'                wxt = wxt1)
#' plot(APCfit1, parametricbx = FALSE, nCol = 3)
#' 
#' #Also Zero-weight the 1886 cohort
#' wxt2 <- genWeightMat(55:89,  EWMaleData$years, clip = 3, 
#'                      zeroCohorts = 1886)
#' APCfit2 <- fit(apc(), data = EWMaleData, ages.fit = 55:89, 
#'                wxt = wxt2)
#' plot(APCfit2, parametricbx = FALSE, nCol = 3)
#'
#'@export
genWeightMat <- function(ages, years, clip = 0, zeroCohorts = NULL) {
  nAges <- length(ages)
  nYears <- length(years)  
  cohorts <- (years[1] - ages[nAges]):(years[nYears] - ages[1])
  nCohorts <- length(cohorts)  
  wxt <- matrix(1, nrow = nAges, ncol = nYears)
  rownames(wxt) <- ages
  colnames(wxt) <- years
  if (clip > 0) {
    zeroCohorts <- c(zeroCohorts, head(cohorts, clip), tail(cohorts, clip))
  }
  for (c in zeroCohorts){
    h <- c - cohorts[1] + 1 - nAges
    if (h <= 0){
      col <- 1
      row <- -h + 1
    } else {
      row <- 1
      col <- h + 1
    }
    while (col <= nYears && row <= nAges) {
      wxt[row, col] <- 0
      row <- row + 1
      col <- col + 1          
    }  
  }
  wxt
}


#'Extract cohort from an age-period array
#'
#'Extract cohorts from an age-period array. This is useful to
#'construct a life table or to perform actuarial/demographic 
#'calculations on a cohort basis using the output of several functions
#'in \code{StMoMo}.
#'
#'@param A an age-period array with a demographic quantity. This array 
#'can have two or more dimensions, with the first dimension being the age and 
#'the second dimension being the period (calendar year). Note that the names 
#'of these two dimension are taken to represent the possible ages and periods
#'in the array.
#'   
#'@param age optional age for defining the cohort to be extracted. If 
#'argument \code{age} is provided (and argument \code{cohort} is not) then the 
#'extracted cohort corresponds to those born in the year \code{period-age}.
#'@param period optional period (calendar year) for defining the cohort to be 
#'extracted. If argument \code{period} is provided (and argument \code{cohort} 
#'is not) then the extracted cohort corresponds to those born in the year 
#'\code{period-age}.
#'@param cohort optional cohort to be extracted. If this argument is provided
#'then arguments \code{age} and \code{period} are ignored.
#'
#'@return If the the input array is two dimensional the the output is a 
#'a vector with the quantity along the cohort. Otherwise if \code{A} is an 
#'N-dimensional  array the output is an (N-1)-dimensional array with the first 
#'dimension representing the cohort.
#'
#'
#'@seealso \code{\link{fitted.fitStMoMo}}, \code{\link{forecast.fitStMoMo}},
#'\code{\link{simulate.fitStMoMo}}, \code{\link{simulate.bootStMoMo}}
#'
#'@examples
#' LCfit <- fit(lc(), data = EWMaleData, ages.fit = 55:89)
#' #Plot forecast mortality rates for the 1950 cohort
#' LCfor <- forecast(LCfit)
#' plot(55:61, extractCohort(fitted(LCfit, type = "rates"), cohort = 1950), 
#'      type = "l", log = "y", xlab = "age", ylab = "Mortality rate", 
#'      main = "Mortality rates for the 1950 cohort", 
#'      xlim = c(55,89), ylim = c(0.005, 0.12))
#' lines(62:89, extractCohort(LCfor$rates, cohort = 1950), lty = 2, col = "blue")
#' 
#' 
#' #Plot 10 simulated sets of mortality rates for the cohort 
#' # aged 60 in year 2010 (i.e., the 1950 cohort)
#' LCsim <- simulate(LCfit, nsim = 10)
#' mSim <- extractCohort(LCsim$rates, age = 60, period = 2010)
#' plot(55:61, extractCohort(fitted(LCfit, type = "rates"), cohort = 1950), 
#'      type = "l", log = "y", xlab = "age", ylab = "Mortality rate", 
#'      main = "Mortality rates for the 1950 cohort", 
#'      xlim = c(55,89), ylim = c(0.005, 0.12))
#' matlines(62:89, mSim, lty = 2)
#' 
#'@export
extractCohort <- function(A, age = as.numeric(dimnames(A)[[1]][1]), 
                          period = as.numeric(dimnames(A)[[2]][1]), 
                          cohort = period - age){
  if (!is.array(A)) stop( "extractCohort only defined for arrays.")
  
  ages <- as.numeric(dimnames(A)[[1]])   
  years <- as.numeric(dimnames(A)[[2]])
  nAges <- length(ages)
  nYears <- length(years)  
  cohorts <- (years[1] - ages[nAges]):(years[nYears] - ages[1])
  
  if (!(cohort %in% cohorts)) 
    stop("The required cohort is not part of the data.")
  
  #Find roe and column indexes corresponding to the cohort
  diag <- cohort - (years[1]-ages[1])
  if (diag >= 0){
    i0 <- 1
    j0 <- diag + 1
  } else {
    i0 <- -diag + 1
    j0 <- 1
  }
  nc <- min(nAges-i0, nYears - j0)
  
  #Extract cohort
  nDims <- length(dim(A))
  if (nDims == 2){
    out <- A[cbind(i0:(i0 + nc),  j0:(j0 + nc))]
    names(out) <- i0:(i0 + nc) + ages[1]-1
  } else {
    out <- apply(A, 3:nDims, function(x) x[cbind(i0:(i0 + nc),  
                                                 j0:(j0 + nc))])
    dimnames(out)[[1]] <- i0:(i0 + nc) + ages[1]-1
  }
  out
}


#' Create StMoMoData object from demogdata object 
#' 
#' Create StMoMoData object suitable for fitting a Stochastic Mortality 
#' Model using function \code{\link{fit.StMoMo}}. 
#' 
#' @param  data demogdata object of type "mortality". It is either the
#' output from functions \code{\link[demography]{read.demogdata}} or 
#' \code{\link[demography]{hmd.mx}} of package \code{demography}.
#' @param series name of series within \code{data} to use.
#' @param type the type of exposure that should be included in the 
#' output. The alternatives are \code{"central"} (default) and 
#' \code{"initial"}.  \code{"central"} exposures are suitable for fitting
#' models under a log-Poisson framework while \code{"initial"} exposures
#' are suitable under a logit-Binomial framework. 
#' 
#' @return A list with class \code{"StMoMoData"} with components:
#'   
#'   \item{Dxt}{ matrix of deaths data.}   
#'   
#'   \item{Ext}{ matrix of observed exposures.} 
#'   
#'   \item{ages}{ vector of ages corresponding to rows of \code{Dxt} and 
#'   \code{Ext}.} 
#'   
#'   \item{years}{ vector of years corresponding to rows of \code{Dxt} and 
#'   \code{Ext}.} 
#'   
#'   \item{type}{ the type of exposure in the data.}
#'   
#'   \item{series}{ name of the extracted series.}
#'   
#'   \item{label}{ label of the data.}
#' 
#' @examples 
#' \dontrun{
#' library(demography)
#' NZdata <- hmd.mx(country = "NZL_NP", username = username, password = password, 
#' label = "New Zealand")
#' NZStMoMo <- StMoMoData(NZdata, series = "male")
#' summary(NZStMoMo)
#' }
#' 
#' @export
StMoMoData <- function(data, series = names(data$rate)[1], 
                       type = c("central", "initial")){
  
  if (class(data) != "demogdata")
    stop("Argument data needs to be of class demogdata.")
  if (data$type != "mortality")
    stop("Argument data needs to be of class demogdata and type mortality.")
  type <- match.arg(type)
  
  Ext <- data$pop[[series]]
  Dxt <- data$pop[[series]]  * data$rate[[series]]
  if (type == "initial") Ext <- Ext + 0.5 * Dxt
  
  rownames(Ext) <- rownames(Dxt) <- data$age
  
  structure(list(Dxt = Dxt, Ext = Ext, ages = data$age, 
                 years = data$year, type = type, 
                 series = series, label = data$label), 
            class = "StMoMoData")
    
}

#' @export
print.StMoMoData <- function(x, ...)
{
  cat(paste("Mortality data for", x$label))
  cat("\n    Series: ", x$series)
  cat(paste("\n    Years:", min(x$years), "-", max(x$years)))
  cat(paste("\n    Ages: ", min(x$ages), "-", max(x$ages)))
  cat("\n    Exposure: ", x$type, "\n")
}

#' @export
summary.StMoMoData <- function(object, ...)
{
  print(object)
}

#' Transform StMoMoData from central to initial exposures
#' 
#' Transform StMoMoData from central to initial exposures. Initial exposures
#' are computed by adding one half of the deaths to the central exposures.
#' 
#' @param data StMoMoData object of type "central" created with function 
#' \code{\link{StMoMoData}}.
#' 
#' @return A StMoMoData object of type "initial". 
#' 
#' @seealso \code{\link{initial2central}}
#' 
#' @export
central2initial <- function(data){
  if (class(data) != "StMoMoData")
    stop("Argument data needs to be of class StMoMoData.")
  if (data$type != "central")
    stop("Argument data needs to be of class StMoMoData and type central.")
  data$Ext <- data$Ext + 0.5 * data$Dxt
  data$type <- "initial"
  data
}

#' Transform StMoMoData from initial to central exposures
#' 
#' Transform StMoMoData from initial to central exposures. Central exposures
#' are computed by substracting one half of the deaths from the initial exposures.
#' 
#' @param data StMoMoData object of type "initial" created with function 
#' \code{\link{StMoMoData}}.
#' 
#' @return A StMoMoData object of type "central". 
#' 
#' @seealso \code{\link{central2initial}}
#' 
#' @export
initial2central <- function(data){
  if (class(data) != "StMoMoData")
    stop("Argument data needs to be of class StMoMoData.")
  if (data$type != "initial")
    stop("Argument data needs to be of class StMoMoData and type initial.")
  data$Ext <- data$Ext + 0.5 * data$Dxt
  data$type <- "central"
  data
}