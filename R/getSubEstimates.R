#' @title Computes Estimates from Sub Time Series
#'
#' @description This function computes 5 estimates of the quantity(ies) of 
#' interest from the time series of a single subject.  The resulting estimates 
#' should be vectorized and concatenated into 5 matrices across subjects for 
#' input to shrinkIt (see \code{link{shrinkIt}}).
#' @param Y A \eqn{T\times K} matrix, where T is the number of time points in 
#' the time series and K is the number of variables (e.g. voxels, regions) 
#' observed at each time point.
#' @param FUN A function that can be applied to a TxK matrix to obtain an 
#' estimate of the quantity of interest (e.g. \code{FUN=cor} returns the KxK 
#' pairwise correlation matrix)
#' @param block.size The number of time points forming each block for estimation 
#' of sampling variance (see Details).  Default is 5.
#' @param block.space The number of time points separating consecutive blocks 
#' (see Details).  Default is 1.
#' @details 
#' 
#' The main arguments to \code{\link{shrinkIt}} are W_all, W_part1, W_part2, 
#' W_odd and W_even, each a matrix size \eqn{M\times n}, where M is the length 
#' of the vectorized quantity of interest (e.g. number of elements in the upper 
#' triangle of a \eqn{K\times K} correlation matrix).  W_odd and W_even are 
#' used to estimate sampling variance, while W_part1 and W_part2 are used to 
#' estimate intrasession variance. For more details, see \code{\link{shrinkIt}}.  
#' This function subsamples the time series of a single subject \eqn{i} to 
#' obtain the \eqn{i}th columns of W_all, W_part1, W_part2, W_odd and W_even,
#' respectively.  
#' 
#' Let \eqn{Y_i} be the \eqn{T\times K} timeseries of \eqn{K} variables for 
#' subject \eqn{i}.  The \eqn{i}th column of \code{W_all} is the estimate of 
#' \eqn{x_i} produced from the full time series \eqn{Y_i}.  The \eqn{i}th column 
#' of \code{W_part1} is the estimate of \eqn{x_i} produced from the first 
#' \eqn{T'} time points of \eqn{Y_i}, where \eqn{T'\leq T/2} is the length of 
#' the even block time series (see below).  The \eqn{i}th column of 
#' \code{W_part2} is the estimate of \eqn{x_i} produced from the last \eqn{T'} 
#' time points of \eqn{Y_i}. 
#' 
#' The \eqn{i}th column of \code{W_odd} is the estimate of \eqn{x_i} produced 
#' from the odd blocks (of length \code{block.size}) of \eqn{Y_i}. The \eqn{i}th 
#' column of \code{W_even} is the estimate of \eqn{x_i} produced from the even 
#' blocks (of length \code{block.size}) of \eqn{Y_i}.  Blocks are defined as 
#' \code{block.size} consecutive time points, separated by \code{block.space} 
#' time points.  These are used instead of individual time points to reduce 
#' dependence between estimates produced from the odd and even timeseries.
#' 
#' @export
#' @return A list containing five elements, each a vector, to be input into the 
#' \eqn{i}th columns of \code{W_all}, \code{W_part1}, \code{W_part2}, 
#' \code{W_odd} and \code{W_even}, respectively.  
#'
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
  w_all <- FUN(Y)
  w_part1 <- FUN(Y[inds1,])
  w_part2 <- FUN(Y[inds2,])
  w_odd <- FUN(Y[inds_odd,])
  w_even <- FUN(Y[inds_even,])
  
  result <- list(w_all, w_part1, w_part2, w_odd, w_even)
  names(result) <- c('w_all','w_part1','w_part2','w_odd','w_even')
  return(result)
}