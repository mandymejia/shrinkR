#' @title Performs shrinkage for subject-level parcellation using resting-state fMRI
#'
#' @description This function performs shrinkage towards the group mean of subject-level 
#' observations of any voxel-by-voxel similarity or distance matrix.  This can be used as a 
#' pre-processing step before performing clustering in order to improve reliability of subject-
#' level parcellations.  If there are M variables (e.g. the elements in the upper triangle of the 
#' voxel-by-voxel similarity/distance matrix) observed for each subject, then a different 
#' shrinkage parameter will be 
#' computed for each variables The shrinkage parameter ranges from 0 (no shrinkage) to 1 (complete 
#' shrinkage), and is determined by the relationship between within-subject (noise) variance and 
#' between-subject (signal) variance.  Both variance components and the shrinkage parameter are
#' estimated from the data.
#' @param X1 An m-by-n matrix, where the ith column contains the m observed values of subject i 
#' @param X2 An m-by-n matrix, where the ith column contains a second set of the m observed values 
#' of subject i (see Description for details)
#' @param len A scalar equal to the total amount of scan time (in minutes) collected for each
#' subject.  This value will be used to adjust the noise variance estimate, since it will tend
#' to be over-estimated by splitting the data to compute X1 and X2.
#' @details The shrinkage estimate of value m for subject i is given by
#' 
#'          X*_i(m) = λ_m(X_bar(m))+[1-λ_m]X_i(m).
#'          
#' The "raw" observation X_i is the average of X1_i (the first observation vector for subject i) 
#' and X2_i (the second observation vector for subject i).  
#'          
#' The shrinkage parameter λ(m) for variable m is computed as
#' 
#'        λ(m) = σ^2_u/(σ^2_u+σ^2_x(m)).
#' 
#' Two observations (X1 and X2) are required to compute the noise variance σ^2_u.  If a total of 
#' T minutes of scan time is available for each subject, then X1 should be computed using the 
#' first T/2 minutes, and X2 should be computed using the last T/2 minutes.  The difference 
#' between X1 and X2 is used to estimate the noise variance.  However, since only T/2 minutes of
#' scan time are used to compute X1 and X2, the resulting estimate will be appropriate for 
#' scans of length T/2 minutes.  Therefore, we apply a correction, which is determined by the 
#' total scan time T.  The details of this correction are described in Mejia et al. (in press).  
#' The noise variance is estimated globally (across all observed values m), but the signal 
#' variance is estimated separately for each variable m.  The shrinkage parameter λ_m is therefore
#' different for each observed value m.
#' 
#' @export
#' @return A list containing (1) the shrunken estimates (an m-by-n matrix), (2) the shrinkage 
#' parameter lambda (a vector of length m), (3) the noise variance estimate (a scalar), and (4) 
#' the signal variance estimate (a vector of length m).  
#' @importFrom matrixStats rowVars
#' @examples \dontrun{
#'
#'}
shrinkIt <- function(X1, X2, len){
  
  dims <- dim(X1)
  n <- dims[2] #number of subjects
  m <- dims[1] #number of observations per subject
  
  ## Perform Checks
  if(!all.equal(dims,dim(X2))) stop("Dimensions of X1 and X2 do not match")
  if(length(len) != 1) stop("len should be a scalar")
  if(!is.numeric(len)) stop("len must be numeric")
  if(!is.numeric(X1) | !is.numeric(X2)) stop("X1 and X2 must be numeric")
  
  ## Compute Noise Variance
  
  Var.u <- mean(1/2*rowVars(X2-X1))
  theta <- 0.590 + 0.129*log(len)
  Var.u <- Var.u*theta
  
  ## Compute Total Variance
  
  X <- (X1+X2)/2
  Var.w <- rowVars(X)
  
  ## Compute Signal Variance

  Var.x <- Var.w - Var.u
  Var.x[Var.x<0] <- 0
  
  ## Compute Shrinkage Parameter 
  
  #the higher the noise variance, the more shrinkage towards the group mean
  #if noiseVar="individual" or "scaled", Var.u is a matrix and Var.x is a vector
  #in this case, the vector Var.x will get recycled by column and the result will be a matrix
  lambda <- Var.u/(Var.x+Var.u) 

  ## Perform Shrinkage
  
  #group.mean is a vector and will get recycled by column 
  #if noiseVar="common" or "global", lambda is a vector and will get recycled by column
  X.bar <- rowMeans(X)
  X.shrink <- lambda*X.bar + (1 - lambda)*X
  
  result <- list(X.shrink, lambda, Var.u, Var.x)
  names(result) <- c("X.shrink", "lambda", "Var.u", "Var.x")
  return(result)
}