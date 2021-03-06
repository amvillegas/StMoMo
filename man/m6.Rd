% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standardModels.R
\name{m6}
\alias{m6}
\title{Create an M6 type extension of the Cairns-Blake-Dowd mortality model}
\usage{
m6(link = c("logit", "log"))
}
\arguments{
\item{link}{defines the link function and random component associated with 
the mortality model. \code{"log"} would assume that deaths follow a 
Poisson distribution and use a log link while \code{"logit"} would 
assume that deaths follow a Binomial distribution and a logit link. 
Note that the default is the logit link.}
}
\value{
An object of class \code{"StMoMo"}.
}
\description{
Utility function to initialise a \code{StMoMo} object representing the 
M6 (CBD with cohorts) extension of the Cairns-Blake-Dowd mortality model 
introduced in Cairns et al (2009).
}
\details{
The created model is either a logit-Binomial or a log-Poisson version of the 
M6 model which has predictor structure 
\deqn{\eta_{xt} = \kappa_t^{(1)} + (x-\bar{x})\kappa_t^{(2)} + \gamma_{t-x},} 
where \eqn{\bar{x}} is the average age in the data.

Identifiability of the model is accomplished by applying parameters 
constraints
\deqn{\sum_c\gamma_c = 0, \sum_c c\gamma_c = 0}
which ensure that the cohort effect fluctuates around zero and has no 
linear trend. These constraints are applied using the strategy discussed 
in Appendix A of Haberman and Renshaw (2011).
}
\examples{

M6 <- m6()
wxt <- genWeightMat(55:89,  EWMaleData$years, clip = 3)
M6fit <- fit(M6, data = central2initial(EWMaleData), ages.fit = 55:89)
plot(M6fit, parametricbx = FALSE)

}
\references{
Cairns, A. J. G., Blake, D., Dowd, K., Coughlan, G. D., Epstein, D., 
Ong, A., & Balevich, I. (2009). A quantitative comparison of stochastic 
mortality models using data from England and Wales and the United States. 
North American Actuarial Journal, 13(1), 1-35.

Haberman, S., & Renshaw, A. (2011). A comparative study of parametric 
mortality projection models. Insurance: Mathematics and Economics, 
48(1), 35-55.
}
\seealso{
\code{\link{StMoMo}}, \code{\link{central2initial}}, 
\code{\link{cbd}}, \code{\link{m7}}, \code{\link{m8}}
}
