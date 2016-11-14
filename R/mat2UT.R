#' @title Matrix to Upper Triangular
#'
#' @description Returns the vectorized upper triangle of a square matrix
#' @param x A square matrix
#' 
#' @export
#' @return The vectorized upper triangle of x.  
#' @examples \dontrun{
#'
#'}
mat2UT <- function(x){

  UT <- x[upper.tri(x)]  
  return(UT)
}