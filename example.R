library(corrplot)

##############################################################
# Load timeseries data 
##############################################################

# dat is a list, each element is the TxV time series for a single subject (10 subjects)
# T = 4800 (4 runs of 1200 volumes each), V = 25 (regions from an ICA)
# Will use the first run to perform shrinkage and the third run to evaluate predictive accuracy
# Third run is used since it is from a future visit

data("dat", envir = environment()) #loads dat object
n <- length(dat) #number of subjects

##############################################################
# Compute verious estimates of RSFC
##############################################################

sub_estimates <- list() #store estimates produced from sub-timeseries
raw_estimates1 <- list() #store raw estimates from visit 1 timeseries
raw_estimates2 <- list() #store raw estimates from visit 2 timeseries (for performance evaluation)
for(i in 1:n){

  #visit 1
  dat_i <- dat[[i]][1:1200,]
  raw_estimates1[[i]] <- cor(dat_i)
  list_i <- getSubEstimates(dat_i, cor) #each element is one of the sub-timeseries estimates
  vecs_i <- sapply(list_i, mat2UT) #each column is the UT of one of the sub-timeseries estimates
    
  #initialize Mxn matrices to store sub-timeseries estimates for all subjects
  if(i==1){
      M <- nrow(vecs_i)
      sub_estimates[[1]] <- matrix(nrow=M, ncol=n)
      sub_estimates[[2]] <- matrix(nrow=M, ncol=n)
      sub_estimates[[3]] <- matrix(nrow=M, ncol=n)
      sub_estimates[[4]] <- matrix(nrow=M, ncol=n)
      names(sub_estimates) <- names(list_i)
  }
  
  #input subject i's estimates into sub_estimates[[1:4]]
  sub_estimates[[1]][,i] <- vecs_i[,1]
  sub_estimates[[2]][,i] <- vecs_i[,2]
  sub_estimates[[3]][,i] <- vecs_i[,3]
  sub_estimates[[4]][,i] <- vecs_i[,4]
  
  #visit 2
  dat_i <- dat[[i]][2401:3600,]
  raw_estimates2[[i]] <- cor(dat_i)
  
}

 
#####################################################
# Perform shrinkage
#####################################################

result <- shrinkIt(W_part1 = sub_estimates$w_part1, 
                   W_part2 = sub_estimates$w_part2,
                   W_odd = sub_estimates$w_odd,
                   W_even = sub_estimates$w_even)


#####################################################
# Transform estimates back to matrices
#####################################################

shrinkage_estimates <- list()
MSE_raw <- rep(0, n)
MSE_shrink <- rep(0, n)

for(i in 1:n){
  
  MSE_raw[i] <- mean((mat2UT(raw_estimates1[[i]]) - mat2UT(raw_estimates2[[i]]))^2)
  MSE_shrink[i] <- mean((result$shrinkage.estimates[,i] - mat2UT(raw_estimates2[[i]]))^2)
  shrinkage_estimates[[i]] <- UT2mat(result$shrinkage.estimates[,i])
  
}

#####################################################
# Visualize matrix of shrinkage parameters
#####################################################

corrplot(UT2mat(result$shrinkage.parameter, diag=0), method='color', tl.pos="n", cl.lim=c(0,1))


#####################################################
# Visualize estimates
#####################################################

#each row is a subject
#column 1: raw estimates from visit 1
#column 2: shrinkage estimates from visit 1
#column 3: raw estimates from visit 2
par(mfrow = c(5,3), mar=c(0,0,0,0))
for(i in 1:5){
  diag(raw_estimates1[[i]]) <- NA
  diag(shrinkage_estimates[[i]]) <- NA
  diag(raw_estimates2[[i]]) <- NA
  corrplot(raw_estimates1[[i]], method='color', tl.pos="n", cl.pos="n", na.label=" ")
  corrplot(shrinkage_estimates[[i]], method='color', tl.pos="n", cl.pos="n", na.label=" ")
  corrplot(raw_estimates2[[i]], method='color', tl.pos="n", cl.pos="n", na.label=" ")
}

