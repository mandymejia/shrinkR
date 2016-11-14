#' @title Row Variances
#'
#' @description Computes variance of each row of a matrix.
#' @param x a matrix
#' @param na.rm If TRUE (default), \code{NA} values are stripped before
#' computation.  If FALSE, \code{NA} values are included in computation.
#' @export
#' @return A vector of numerics of length \code{nrow(x)}
#'
rowVars <- function (x, na.rm = TRUE) 
{
  sqr = function(x) x * x
  n = rowSums(!is.na(x))
  n[n <= 1] = NA
  return(rowSums(sqr(x - rowMeans(x,na.rm = na.rm)), na.rm = na.rm)/(n - 1))
}