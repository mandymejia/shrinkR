#' @title Perform shrinkage of subject-level values towards group mean
#'
#' @description This function performs shrinkage of subject values towards the group mean.
#' A common (group) noise variance, individual noise variance, or individually-scaled noise
#' variance may be computed.  The noise variance is computed for each variable 1 to m. 
#' Alternatively, a pre-specified global noise variance may be used.
#' @param observation The observation matrix (see description) on which shrinkage will be performed
#' @param obs1 Observation matrix (see description) of the first observation for each subject
#' @param obs1 Observation matrix (see description) of the second observation for each subject
#' @param noiseVar Either a character string equal to "common", "individual" or "scaled",
#' in which case the noise variance will be computed, or a numeric, non-negative scalar 
#' to be used as the global noise variance.
#' @details Each observation matrix is m-by-n, where m is the number of variables observed
#' for each subject, and n is the number of subjects.  The variables should match across
#' subjects.  For example, in an fMRI context, m might be the number of elements in 
#' the upper triangle of each subject's V-by-V correlation matrix.  If noise variance is to 
#' be computed, obs1 and obs2 are used to compute the noise variance.  If a global 
#' noise variance value is provided, obs1 and obs2 will not be used.
#' @export
#' @return A list containing the shrunken observation matrix, shrinage parameter (lambda) values, 
#' noise variance values, and signal variance values.  The shrinkage parameter values range
#' between 0 and 1, with 1 representing complete shrinkage towards the group mean and
#' 0 representing no shrinkage towards the group mean.
#' @importFrom genefilter rowVars
#' @examples \dontrun{
#'
#'}
shrinkIt <- function(observation, obs1=NULL, obs2=NULL, noiseVar="common"){
  
  dims <- dim(observation)
    
  ## Perform Checks
  
  #check that obs1 and obs2 provided and dimensions of all observation matrices match
  if(noiseVar %in% c("common","individual","scaled")){
    if(is.null(obs1) | is.null(obs2)) stop("Both obs1 and obs2 should be provided in order to compute noise variance")
    if(dims != dim(obs1) | dims != dim(obs2)) stop("Dimensions do not match")
  }
  
  #check that noise Var is "common", "individual", "scaled", or a positive scalar
  if(is.numeric(noiseVar)){
    if(length(noiseVar) == 1){
      if(noiseVar < 0) stop("noiseVar cannot be negative")
    } else {
      stop("noiseVar is numeric but is not of length 1")
    }
  } else if(!(noiseVar %in% c("common","individual","scaled"))){
    stop("noiseVar parameter should be \"common\", \"individual\", \"scaled\", or a numeric scalar")
  }  
  
  n <- dims[2] #number of subjects
  m <- dims[1] #number of observations per subject
  
  ## Compute Noise Variance
  
  if(is.numeric(noiseVar)){
    Var.u <- Var.u.avg <- rep(noiseVar, M)
  } else if(noiseVar=="common") {
    Var.u <- Var.u.avg <- 1/2*rowVars(obs2-obs1)
  } else if(noiseVar=="individual"){
    Var.u <- 1/2*(obs2-obs1)^2
    Var.u.avg <- rowMeans(Var.u)
  } else if(noiseVar=="scaled"){
    MSE.1scan <- colMeans((obs2-obs1)^2)
    scale <- MSE.1scan/mean(MSE.1scan)
    #compute variance components and shrinkage parameter
    Var.u.avg <- 1/2*rowVars(obs2-obs1)
    Var.u <- matrix(Var.u.avg, nrow=m, ncol=n) * matrix(scale, nrow=m, ncol=n, byrow=TRUE)
  }
  
  ## Compute Total Variance
  
  if(is.numeric(noiseVar)){
    Var.w <- rowVars(observation)
  } else {
    Var.w1 <- rowVars(obs1)
    Var.w2 <- rowVars(obs2)
    Var.w <- (Var.w1 + Var.w2)/2
  }
  
  ## Compute Signal Variance

  if(noiseVar=="individual"){
    Var.x <- (n-1)/n * (Var.w - Var.u.avg)
  } else {
    Var.x <- Var.w - Var.u.avg
  }
  Var.x[Var.x<0] <- 0
  
  ## Compute Shrinkage Parameter 
  
  #the higher the noise variance, the more shrinkage towards the group mean
  #if noiseVar="individual" or "scaled", Var.u is a matrix and Var.x is a vector
  #in this case, the vector Var.x will get recycled by column and the result will be a matrix
  lambda <- Var.u/(Var.x+Var.u) 

  ## Perform Shrinkage
  
  #group.mean is a vector and will get recycled by column 
  #if noiseVar="common" or "global", lambda is a vector and will get recycled by column
  group.mean <- rowMeans(observation)
  obs.shrink <- lambda*group.mean + (1 - lambda)*observation
  
  result <- list(obs.shrink, lambda, Var.u, Var.x)
  names(result) <- c("obs.shrink", "lambda", "Var.u", "Var.x")
  return(result)
}