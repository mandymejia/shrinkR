#' @title Row Variances
#'
#' @description Computes variance of each row of a matrix.
#' @param x a matrix
#' @export
#' @return A vector of numerics, same length as the number of rows of the 
#' matrix
#' @examples \dontrun{
#'
#'}
rowVars <- function (x, na.rm = TRUE) 
{
  sqr = function(x) x * x
  n = rowSums(!is.na(x))
  n[n <= 1] = NA
  return(rowSums(sqr(x - rowMeans(x,na.rm = na.rm)), na.rm = na.rm)/(n - 1))
}