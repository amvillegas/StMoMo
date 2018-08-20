# Mon Aug 20 16:02:25 2018 --------- Marius D. Pascariu ---

#' Hack gnm::Mult weird dependency
#' The way this function is called in StMoMo functions makes it necesary 
#' to have the gnm dependecy. A workaround is redefine it in the package and 
#' export it from here, and import the gnm package. 
#' @inheritParams gnm::Mult
#' @keywords internal
#' @source \code{\link[gnm]{Mult}} or from here: 
#' \url{https://github.com/hturner/gnm/blob/master/R/Mult.R}
#' @keywords internal
#' @export
Mult <- function(..., inst = NULL){
  if ("multiplicity" %in% names(match.call()[-1]))
    stop("multiplicity argument of Mult has been replaced by",
         "\"inst\" argument.")
  dots <- match.call(expand.dots = FALSE)[["..."]]
  list(predictors = dots,
       term = function(predLabels, ...) {
         paste("(", paste(predLabels, collapse = ")*("), ")", sep = "")
       },
       call = as.expression(match.call()),
       match = seq(dots))
}
class(Mult) <- "nonlin"
