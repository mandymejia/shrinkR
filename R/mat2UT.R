#' @title Matrix to Upper Triangular Vector
#'
#' @description Returns the vectorized upper triangle of a square matrix
#' @param x A square matrix
#' 
#' @export
#' @return The vectorized upper triangle of x.  
#'
mat2UT <- function(x){

  UT <- x[upper.tri(x)]  
  return(UT)
}