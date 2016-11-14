#' @title Returns the symmatric square matrix from a vector containing the upper 
#' triangular elements
#'
#' @description 
#' @param x A vector containing the upper triangular elements of a square, symmetric 
#' matrix.
#' @param diag A scalar equal to the diagonal values of the matrix to be returned.
#' 
#' @export
#' @return A symmetric matrix with the values of x in the upper and lower triangles 
#' and the value \code{diag} on the diagonals.
#' @examples \dontrun{
#'
#'}
UT2mat <- function(x, diag=1){

  if(!is.vector(x) | !is.numeric(x)) stop('x must be a numeric vector')
  if(!is.numeric(diag) | length(diag)>1) stop('diag must be numeric and of length 1')
  
  M <- length(x)
  V <- (1+sqrt(8*M+1))/2
  if(round(V) != V) stop('Length of x not equal to V(V-1)/2, for some integer V.')
  
  mat <- matrix(0, nrow=V, ncol=V)
  mat[upper.tri(mat)] <- x
  mat <- mat + t(mat)
  diag(mat) <- diag
  
  return(mat)
}