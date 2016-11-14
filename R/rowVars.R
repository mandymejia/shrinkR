#' @title Computes variance of each row of a matrix.
#'
#' @description 
#' @param x a matrix
#' @details 
#' @export
#' @return 
#' @importFrom 
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