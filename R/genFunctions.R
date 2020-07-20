#'Generate call functions
#'
#'Generates a list of call functions of order \eqn{j} given a vector of strikes.
#'For each strike \eqn{k}, the functions are of the form
#' \deqn{f(x) = \max(x-k,0)^j} 
#'This is useful for constructing a large abstract Stochastic
#'Mortality Model to be fitted using group regularisation.
#'
#'@param strikes vector of strikes.
#'@param order optional scalar indicating the order \eqn{j} of the call function.
#'   
#'@return A list of call functions for the given vector of strikes.
#'
#'@export
genCall <- function(strikes, order = 1) {
  order <- order[1]
  vector_of_functions = NULL
  call <- function(x, ages, k) {
    ifelse((x-k)<0, 0, (x-k)^order)
  }
  for (i in strikes) {
    vector_of_functions = c(vector_of_functions, eval(substitute(function(x, ages) call(x=x, ages=ages, k=i), 
                                                                 list(i=i))))
  }
  return(vector_of_functions)
}

#'Generate put functions
#'
#'Generates a list of put functions of order \eqn{j} given a vector of strikes.
#'For each strike \eqn{k}, the functions are of the form
#' \deqn{f(x) = \max(k-x,0)^j}  
#'This is useful for constructing a large abstract Stochastic
#'Mortality Model to be fitted using group regularisation.
#'
#'@param strikes vector of strikes.
#'@param order optional scalar indicating the order \eqn{j} of the call function.
#'   
#'@return A list of put functions for the given vector of strikes.
#'
#'@export
genPut <- function(strikes, order = 1) {
  order <- order[1]
  vector_of_functions = NULL
  put <- function(x, ages, k) {
    ifelse((k-x)<0, 0, (k-x))
  }
  for (i in strikes) {
    vector_of_functions = c(vector_of_functions, eval(substitute(function(x, ages) put(x=x, ages=ages, k=i), list(i=i))))
  }
  return(vector_of_functions)
}

#'Generate above indicator functions
#'
#'Generates a list of above indicator functions given a vector of strikes.
#'For each strike \eqn{k}, the functions are of the form
#' \deqn{f(x) = I_{x>k}} 
#'This is useful for constructing a large abstract Stochastic
#'Mortality Model to be fitted using group regularisation.
#'
#'@param strikes vector of strikes.
#'   
#'@return A list of above indicator functions for the given vector of strikes.
#'
#'@export
genAbove <- function(strikes) {
  vector_of_functions = NULL
  above <- function(x, ages, k) {
    ifelse((x-k)>0, 1, 0)
  }
  for (i in strikes) {
    vector_of_functions = c(vector_of_functions, eval(substitute(function(x, ages) above(x=x, ages=ages, k=i), list(i=i))))
  }
  return(vector_of_functions)
}

#'Generate below indicator functions
#'
#'Generates a list of above indicator functions given a vector of strikes.
#'For each strike \eqn{k}, the functions are of the form
#' \deqn{f(x) = I_{x<k}} 
#'This is useful for constructing a large abstract Stochastic
#'Mortality Model to be fitted using group regularisation.
#'
#'@param strikes vector of strikes.
#'   
#'@return A list of below indicator functions for the given vector of strikes.
#'
#'@export
genBelow <- function(strikes) {
  vector_of_functions = NULL
  below <- function(x, ages, k) {
    ifelse((k-x)>0, 1, 0)
  }
  for (i in strikes) {
    vector_of_functions = c(vector_of_functions, eval(substitute(function(x, ages) below(x=x, ages=ages, k=i), list(i=i))))
  }
  return(vector_of_functions)
}

#'Generate polynomial functions
#'
#'Generates a list of polynomial functions given a vector of exponents.
#'For each exponent \eqn{j}, the functions are of the form
#' \deqn{f(x) = (x-\bar{x})^j} 
#'This is useful for constructing a large abstract Stochastic
#'Mortality Model to be fitted using group regularisation.
#'
#'@param exp vector of exponents.
#'   
#'@return A list of polynomial functions for the given vector of exponents.
#'
#'@export
genPoly <- function(exp) {
  vector_of_functions = NULL
  poly <- function(x, ages, j) {
    (x - mean(ages))^j
  }
  for (i in exp) {
    vector_of_functions = c(vector_of_functions, eval(substitute(function(x, ages) poly(x=x, ages=ages, j=i), list(i=i))))
  }
  return(vector_of_functions)
}