#' @title Returns the vectorized upper triangle of a square matrix
#'
#' @description 
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