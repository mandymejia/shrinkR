#' @title Computes estimates of quantity(ies) of interest from subsamples of time series for use in shrinkIt
#'
#' @description This function computes 4 estimates of the quantity(ies) of interest from the 
#' time series of a single subject.  The resulting estimates should be vectorized and concatenated into 
#' 4 matrices (each subject is a column) for input to shrinkIt.
#' @param Y A \eqn{T\times K} matrix, where T is the number of time points in the time series and K is the number of 
#' variables (e.g. voxels, regions) observed at each time point.
#' @param FUN A function that can be applied to a TxK matrix to obtain an estimate of the quantity of 
#' interest (e.g. \code{FUN=cor} returns the KxK pairwise correlation matrix)
#' @param block.size The number of time points forming each block for estimation of sampling variance 
#' (see Details).  Default is 5 time points.
#' @param block.space The number of time points separating consecutive blocks (see Details).  Default is
#' one time point.
#' @details The shrinkIt function takes in four estimates of the quantity of interest, \eqn{x_i}, for 
#' each subject \eqn{i=1,\dots,n}, in the form of 4 \eqn{M\times n} matrices, where M is the length of the vectorized 
#' quantity of interest (e.g. number of elements in the upper triangle of the \eqn{K\times K} correlation or pairwise 
#' correlation matrix).  This function subsamples the time series of the ith subject to obtain the \eqn{i}th 
#' column of these four matrices.
#' 
#' The four inputs to shrinkIt are W_part1, W_part2, W_odd and W_even.  W_odd and W_even are used to 
#' estimate the sampling variance associated with the estimation of the quantity of interest.  W_part1 
#' and W_part2 are used to estimate the intrasession variance of the signal of interest, which serves as 
#' a proxy for intersession variance.  For details on the computation of both variance components, see 
#' shrinkIt.
#' 
#' Let \eqn{Y_i} be the \eqn{T\times K} timeseries of \eqn{K} variables for subject \eqn{i}.
#'
#' The \eqn{i}th column of W_part1 (the first argument to shrinkIt) is the estimate of \eqn{x_i} produced 
#' from the first \eqn{T'} time points of \eqn{Y_i}, where \eqn{T'\leq T/2}.  The \eqn{i}th column of 
#' W_part2 (the second argument to shrinkIt) is the estimate of \eqn{x_i} produced from the last \eqn{T'} 
#' time points of \eqn{Y_i}.
#' 
#' The \eqn{i}th column of W_odd (the third argument to shrinkIt) is the estimate of \eqn{x_i} produced 
#' from the odd blocks of length block.size of \eqn{Y_i}. The \eqn{i}th column of W_even (the fourth 
#' argument to shrinkIt) is the estimate of \eqn{x_i} produced from the even blocks of length block.size of 
#' \eqn{Y_i}.  Blocks of \code{block.size} consecutive time points, separated by \code{block.space} time 
#' points, are used instead of individual time points to reduce dependence between estimates produced from 
#' the timeseries formed by odd and even blocks, respectively.
#' 
#' 
#' @export
#' @return A list containing four elements, each a vector, to be input into the \eqn{i}th columns of W_part1, 
#' W_part2, W_odd and W_even, respectively.  
#' @examples \dontrun{
#'
#'}
getSubEstimates <- function(Y, FUN, block.size=5, block.space=1){
  
  ntime <- nrow(Y)
  inds <- 1:ntime
  
  # Identify blocks of length block.size, separated by spaces of length block.space
  end <- max(which(inds %% (block.size+block.space) == block.size)) #identify end of last full block
  blocks <- floor(((1:end)-1)/(block.size+block.space)) + 1 #identify each block
  blanks <- which((1:end) %% (block.size+block.space) == 0) #identify blanks at end of each block
  blocks[blanks] <- NA #label spaces as NA

  # Identify odd and even blocks
  is.even <- function(x) { x %% 2 == 0 }
  is.odd <- function(x) { x %% 2 != 0 }
  inds_odd <- which(is.odd(blocks))
  inds_even <- which(is.even(blocks))
  inds_odd <- inds_odd[1:length(inds_even)] #may have more odd blocks than even.  if so, shorten.
  
  # Identify first and second "halves" of timeseries (will have a gap in the middle if block.space > 0)
  ntime2 <- length(inds_even)
  inds1 <- 1:ntime2
  inds2 <- (ntime-ntime2+1):ntime
  
  # Apply FUN to each sub-timeseries
  FUN <- match.fun(FUN)
  w_part1 <- FUN(Y[inds1,])
  w_part2 <- FUN(Y[inds2,])
  w_odd <- FUN(Y[inds_odd,])
  w_even <- FUN(Y[inds_even,])
  
  result <- list(w_part1, w_part2, w_odd, w_even)
  names(result) <- c('w_part1','w_part2','w_odd','w_even')
  return(result)
}