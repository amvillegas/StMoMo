#' Create a Lee-Carter model
#' 
#' Utility function to initialise a \code{StMoMo} object representing a Lee-Carter model.
#' 
#' The created model is either a log-Poisson (see Brouhns et al (2002)) or a 
#' logit-Binomial version of the Lee-Carter model which has predictor structure   
#' \deqn{\eta_{xt} = \alpha_x + \beta_x\kappa_t.}
#' To ensure identifiability one of the  following constraints is imposed
#' \deqn{\sum\kappa_t = 0,\,\kappa_1 = 0,\, \kappa_n = 0}
#' depending on the value of \code{const}, and
#' \deqn{\sum\beta_x = 1.}
#' 
#' @inheritParams StMoMo
#' @param const defines the constraint to impose to the period index of the model
#' to ensure identifiablity. The alternatives are \code{"sum"}(default),  \code{"last"}
#' and \code{"first"} which apply constraints \eqn{\sum\kappa_t = 0}, \eqn{\kappa_n = 0}
#' and \eqn{\kappa_1 = 0}, respectively.
#' 
#' @return An object of class "StMoMo".
#' 
#' @seealso \link{StMoMo}
#'  
#' @references
#' Brouhns, N., Denuit, M., & Vermunt, J. K. (2002). A Poisson log-bilinear regression 
#' approach to the construction of projected lifetables. Insurance: Mathematics and 
#' Economics, 31(3), 373-393.
#' 
#' Lee, R. D., & Carter, L. R. (1992). Modeling and forecasting U.S. mortality. 
#' Journal of the American Statistical Association, 87(419), 659-671. 
#' 
#' @examples
#' 
#' #sum(kt) = 0 and log link
#' LC1 <- lc()
#' LCfit1<-fit(LC1, Dxt = EWMaleData$Dxt,Ext = EWMaleData$Ext, 
#'             ages = EWMaleData$ages, years = EWMaleData$years)
#' plot(LCfit1)
#' 
#' #kt[1]= 0 and log link
#' LC2 <- lc(const = "first")
#' LCfit2<-fit(LC2, Dxt = EWMaleData$Dxt,Ext = EWMaleData$Ext, 
#'             ages = EWMaleData$ages, years = EWMaleData$years)
#' plot(LCfit2)
#' 
#' #kt[n]= 0 and logit link
#' LC3 <- lc("logit", "last")
#' LCfit3<-fit(LC3, Dxt = EWMaleData$Dxt,Ext = EWMaleData$Ext, 
#'             ages = EWMaleData$ages, years = EWMaleData$years)
#' plot(LCfit3)
#' 
#' @export
lc <- function(link = c("log", "logit"), const = c("sum", "last", "first")){
  link <- match.arg(link)
  const <- match.arg(const)  
  constLC <- function(ax, bx, kt, b0x, gc, wxt, ages){    
    c1 <- switch(const, sum = mean(kt[1, ], na.rm = TRUE),
                 first = kt[1, 1], last = tail(kt[1, ], 1))
    ax <- ax + c1 * bx[, 1]
    kt[1, ] <- kt[1, ] - c1
    c2 <- sum(bx[, 1], na.rm = TRUE)
    bx[, 1] <- bx[, 1] / c2
    kt[1, ] <- kt[1, ] * c2
    list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)
  }
  StMoMo(link = link, staticAgeFun = TRUE, periodAgeFun = "NP", constFun = constLC)    
}

#' Create a Cairns-Blake-Dowds mortality model
#' 
#' Utility function to initialise a \code{StMoMo} object representing a 
#' Cairns-Blake-Dowds mortality model.
#' 
#' The created model is either a logit-Binomial or a log-Poisson version of the 
#' Cairns-Blake-Dowd mortality moddel which has predictor structure 
#' \deqn{\eta_{xt} = \kappa_t^{(1)} + (x-\bar{x})\kappa_t^{(2)}}
#' 
#' @param link defines the link function and error distibution associated with 
#'   the mortality model. \code{"log"} would assume that deaths follow a Poisson
#'   distribution and use a log link while \code{"logit"} would assume that 
#'   deaths follow a Binomial distribution and a logit link. Note that the default 
#'   is the logit link.
#' @return An object of class "StMoMo".
#' 
#' @seealso \link{StMoMo}
#'  
#' @references
#' Cairns, A. J. G., Blake, D., & Dowd, K. (2006). A Two-Factor Model for Stochastic 
#' Mortality with Parameter Uncertainty: Theory and Calibration. 
#' Journal of Risk and Insurance, 73(4), 687-718.
#' 
#' @examples
#' 
#' CBD <- cbd()
#' CBDfit <- fit(CBD, Dxt = EWMaleData$Dxt, Ext = EWMaleData$Ext, 
#'             ages = EWMaleData$ages, years = EWMaleData$years,
#'             ages.fit = 55:89)
#' plot(CBDfit, parametricbx = FALSE)
#' 
#' @export
cbd <- function(link = c("logit", "log")){
  link <- match.arg(link)
  f1 <- function(x,ages) x - mean(ages)
  StMoMo(link = link, staticAgeFun = FALSE, periodAgeFun=c("1", f1))
}


#' Create an Age-Period-Cohort mortality model
#' 
#' Utility function to initialise a \code{StMoMo} object representing a 
#' Age-Period-Cohort mortality model.
#' 
#' The created model is eithera log-Poisson or a logit-Binomial version of the 
#' classical age-period-cohort mortality model which has predictor structure 
#' \deqn{\eta_{xt} = \alpha_x + \kappa_t + \gamma_{t-x}.}
#' 
#' To ensure identifiability we follow Cairns et al (2009) and impose constraints 
#' \deqn{\sum \gamma_c = 0}  and  \deqn{\sum c\gamma_c = 0}.
#' 
#' @inheritParams StMoMo
#' @return An object of class "StMoMo".
#' 
#' @seealso \link{StMoMo}
#' 
#' @references
#' 
#' Cairns, A. J. G., Blake, D., Dowd, K., Coughlan, G. D., Epstein, D., Ong, A., 
#' & Balevich, I. (2009). A quantitative comparison of stochastic mortality models using 
#' data from England and Wales and the United States. 
#' North American Actuarial Journal, 13(1), 1-35.
#' 
#' @examples
#' 
#' APC <- apc()
#' wxt <- genWeightMat(EWMaleData$ages,  EWMaleData$years, clip = 3)
#' APCfit <- fit(APC, Dxt = EWMaleData$Dxt, Ext = EWMaleData$Ext, 
#'               ages = EWMaleData$ages, years = EWMaleData$years, 
#'               wxt = wxt)
#' plot(APCfit, parametricbx = FALSE, nCol = 3)
#' 
#' @export
apc <- function(link = c("log", "logit")){
  link <- match.arg(link)
  constAPC <- function(ax, bx, kt, b0x, gc, wxt, ages){    
    nYears <- dim(wxt)[2]  
    x <- ages  
    t <- 1:nYears
    c <- (1 - tail(ages, 1)):(nYears - ages[1])    
    #\sum g(c)=0  and  \sum cg(c)=0   
    indC <- !is.na(gc)
    phiReg <- lm(gc ~ 1 + c, data.frame(gc = gc[indC], c = c[indC]))  
    phi <- coef(phiReg)   
    gc[indC] <- residuals(phiReg)  
    ax <- ax + phi[1] - phi[2] *x 
    kt <- kt + phi[2] * t   
    #\sum k(t)=0 
    c1 <- mean(kt, na.rm = TRUE)
    kt <- kt - c1
    ax <- ax + c1    
    list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc) 
  }
  StMoMo(link = link, staticAgeFun = TRUE, periodAgeFun = "1",
                cohortAgeFun = "1", constFun = constAPC)

}