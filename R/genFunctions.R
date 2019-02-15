#'Generate call functions
#'
#'Generates a list of call functions given a vector of strikes.
#'This is useful for constructing a large abstract Stochastic
#'Mortality Model to be fitted using the group lasso.
#'
#'@param strikes vector of strikes.
#'   
#'@return A list of call functions for the given vector of strikes.
#'
#'@export
genCall <- function(strikes) {
  vector_of_functions = NULL
  call <- function(x, ages, k) {
    ifelse((x-k)<0, 0, (x-k))
  }
  for (i in strikes) {
    vector_of_functions = c(vector_of_functions, eval(substitute(function(x, ages) call(x=x, ages=ages, k=i), list(i=i))))
  }
  return(vector_of_functions)
}

#'Generate put functions
#'
#'Generates a list of put functions given a vector of strikes.
#'This is useful for constructing a large abstract Stochastic
#'Mortality Model to be fitted using the group lasso.
#'
#'@param strikes vector of strikes.
#'   
#'@return A list of put functions for the given vector of strikes.
#'
#'@export
genPut <- function(strikes) {
  vector_of_functions = NULL
  put <- function(x, ages, k) {
    ifelse((k-x)<0, 0, (k-x))
  }
  for (i in strikes) {
    vector_of_functions = c(vector_of_functions, eval(substitute(function(x, ages) put(x=x, ages=ages, k=i), list(i=i))))
  }
  return(vector_of_functions)
}

#'Generate call squared functions
#'
#'Generates a list of call squared functions given a vector of strikes.
#'This is useful for constructing a large abstract Stochastic
#'Mortality Model to be fitted using the group lasso.
#'
#'@param strikes vector of strikes.
#'   
#'@return A list of call squared functions for the given vector of strikes.
#'
#'@export
genCallSq <- function(strikes) {
  vector_of_functions = NULL
  callSq <- function(x, ages, k) {
    ifelse((x-k)<0, 0, (x-k)^2)
  }
  for (i in strikes) {
    vector_of_functions = c(vector_of_functions, eval(substitute(function(x, ages) callSq(x=x, ages=ages, k=i), list(i=i))))
  }
  return(vector_of_functions)
}

#'Generate put squared functions
#'
#'Generates a list of put squared functions given a vector of strikes.
#'This is useful for constructing a large abstract Stochastic
#'Mortality Model to be fitted using the group lasso.
#'
#'@param strikes vector of strikes.
#'   
#'@return A list of put squared functions for the given vector of strikes.
#'
#'@export
genPutSq <- function(strikes) {
  vector_of_functions = NULL
  putSq <- function(x, ages, k) {
    ifelse((k-x)<0, 0, (k-x)^2)
  }
  for (i in strikes) {
    vector_of_functions = c(vector_of_functions, eval(substitute(function(x, ages) putSq(x=x, ages=ages, k=i), list(i=i))))
  }
  return(vector_of_functions)
}

#'Generate above indicator functions
#'
#'Generates a list of above indicator functions given a vector of strikes.
#'This is useful for constructing a large abstract Stochastic
#'Mortality Model to be fitted using the group lasso.
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
#'This is useful for constructing a large abstract Stochastic
#'Mortality Model to be fitted using the group lasso.
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
#'This is useful for constructing a large abstract Stochastic
#'Mortality Model to be fitted using the group lasso.
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