#' @title Performs empirical Bayes shrinkage of quantity(ies) of interest from stationary or nonstationary 
#' time series data
#'
#' @description This function performs empirical Bayes shrinkage towards the group mean of subject-level 
#' estimates of some summary statistic of time series data (e.g. functional connectivity matrices produced
#' from resting-state fMRI data).
#' @param W_part1 An M-by-n matrix, where the ith column contains estimates of the M quantities of interest
#' produced from the first part of the time series of subject i (see getSubEstimates for details)
#' @param W_part2 An M-by-n matrix, where the ith column contains estimates of the M quantities of interest
#' produced from the second part of the time series of subject i (see getSubEstimates for details)
#' @param W_odd An M-by-n matrix, where the ith column contains estimates of the M quantities of interest
#' produced from the odd blocks of the time series of subject i (see getSubEstimates for details)
#' @param W_even An M-by-n matrix, where the ith column contains estimates of the M quantities of interest
#' produced from the even blocks of the time series of subject i (see getSubEstimates for details)
#' @details This function performs empirical Bayes shrinkage towards the group 
#' mean of subject-level 
#' estimates of some summary statistic of stationary or nonstationary 
#' time series data (e.g. functional connectivity 
#' matrices produced from resting-state fMRI data).  If the quantity of 
#' interest is a voxel pair-wise correlation, 
#' similarity or distance matrix for use in clustering, shrinkage estimates 
#' can be used in place of raw subject-level 
#' estimates to produce more reliable subject-level brain parcellations.  
#' The shrinkage parameter \eqn{\lambda}
#' ranges from 0 (no shrinkage of subject-level estimates) to 1 (complete 
#' shrinkage, so subject-level 
#' estimates are replaced with the group average), and is determined by the 
#' relationship between within-subject 
#' variance and between-subject variance.  The shrinkage parameters 
#' \eqn{\lambda_m}, \eqn{m = 1,\dots,M}, are 
#' computed separately for every quantity of interest (e.g. every element 
#' in the upper triangle of a correlation matrix),
#' but are shared across subjects. The shrinkage estimate of the value 
#' \eqn{X_i(m)} for subject \eqn{i=1,\dots,n} 
#' and quantity \eqn{m = 1,\dots,M} is
#' 
#'          \eqn{\tilde{X}_i(m) = \lambda_m(\bar{W}(m)) + [1-\lambda_m] W_i(m)}.
#'          
#' The "raw" observation \eqn{W_i} is the average of the \eqn{i}th column of 
#' W_part1 and W_part2, and \eqn{\bar{W}(m)}
#' is the average of \eqn{W_i} across all subjects \eqn{i=1,\dots,n}.
#'          
#' The shrinkage parameter \eqn{\lambda_m} represents the optimal trade-off 
#' between signal and noise.  More reliable 
#' subject-level estimates will receive less shrinkage, while less reliable 
#' subject-level estimates will receive more 
#' shrinkage.  \eqn{\lambda_m} is computed from the data as
#' 
#'        \eqn{\lambda_m = \sigma^2_{w/in}/(\sigma^2_{w/in}+\sigma^2_{b/wn})}.
#' 
#' The within-subject variance \eqn{\sigma^2_{w/in}} is composed of two components: sampling variance, or error due to 
#' random variation in the data, and intrasession variance in the signal over time.  The second component is used
#' as a proxy for intersession variance, which cannot be computed in the absence of multiple manifestations of the
#' time series (e.g. multiple fMRI sessions).  The between-subject variance \eqn{\sigma^2_{b/wn}} represents
#' the variance between subjects in the \emph{signal} of interest and is equal to the total variance in the observed
#' estimates (those based on the full time series) minus the total within-subject variance.
#' 
#' 
#' 
#' @export
#' @return A list containing (1) the shrunken estimates (an m-by-n matrix), (2) the shrinkage 
#' parameter lambda (a vector of length m), (3) the noise variance estimate (a scalar), and (4) 
#' the signal variance estimate (a vector of length m).  
#' @examples \dontrun{
#'
#'}
shrinkIt <- function(W_part1, W_part2, W_odd, W_even, estimate.only = FALSE){
  
  ## Perform Checks
  if(!all.equal(dim(W_part1),dim(W_part2)) | 
     !all.equal(dim(W_part1),dim(W_even)) | 
     !all.equal(dim(W_part1),dim(W_odd))) stop("Dimensions of inputs do not match")
  if(!is.numeric(W_part1) | !is.numeric(W_part2) | !is.numeric(W_even) | !is.numeric(W_odd)) stop("Inputs must be numeric")
  
  dims <- dim(W_part1)
  n <- dims[2] #number of subjects
  m <- dims[1] #number of observations per subject
  
  ## Compute SAMPLING Variance (varU)
  
  D <- W_even - W_odd
  varU <- (1/4)*rowVars(D)
  
  ## Compute SIGNAL Variance (varZ) -- close to zero for stationary time series
  
  D <- W_part1 - W_part2
  varD <- rowVars(D)
  varZ <- (1/2)*(varD - 4*varU)
  
  ## Compute WITHIN-SUBJECT Variance (varZ+varU)
  
  varWITHIN <- varU + varZ
  varWITHIN[varWITHIN < 0] <- 0
  
  ## Compute TOTAL Variance
  
  W <- (W_part1+W_part2)/2 #proxy for estimate from full timeseries
  varTOT <- rowVars(W)
  
  ## Compute SHRINKAGE PARAMETER (lambda) - the higher the within-subject variance relative to the between-subject variance, the greater shrinkage towards the group mean
  
  lambda <- varWITHIN/varTOT
  lambda[lambda > 1] <- 1
  lambda[lambda < 0] <- 0
  lambda[varTOT==0] <- 0 #for estimates with zero variance (e.g. diagonals of correlation matrix, which are always zero)

  ## Perform Shrinkage
  
  #Wbar and lambda are vectors and will get recycled by column 
  Wbar <- rowMeans(W)
  Wshrink <- lambda*Wbar + (1 - lambda)*W
  
  if(estimate.only){
    result <- Wshrink
  }
  else {
    result <- list(Wshrink, lambda, varTOT, varWITHIN, varU, varZ)
    names(result) <- c('shrinkage.estimates', 'shrinkage.parameter', 'total.var', 'within.subject.var', 'sampling.var', 'intrasession.var')
  }
  return(result)
}