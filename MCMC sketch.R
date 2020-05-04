dim <- 75 #dimension of square matrix
n <- length(seismic$V1) #length of l-vector
M <- 1000 #number of iteration
MCMC_output <- matrix(0,nrow=n,ncol=M)
sigma <- #standard deviation
mean <- #mean, different for p0 and p1
beta <- #estimated value of beta
current_l <- # some initial vector
  
for (j in (1:M)) { #MCMC-sampling (M steps)
  current_l <- one_step(n,current_l,mean,sigma,beta) #updates the grid
  MCMC_output[,j] <- current_l #stores current value in column j
}

one_step <- function(n,current_l,mean,sigma,beta) {
  for (node in (1:n)) { # Will be running slowly, iterates through entire grid
    p_di# Calculate p(d_i) as function of l_i, return vector-valued (p0, p1)
    # Count number of neighbours? Write a function?
    p_l0 <- #effect from  neighbours
    p_l1 <- #effect from neighbours
    p0_unscaled <- p_l0*p_di[1]
    p1_unscaled <- p_l1*p_di[2]
    p0 <- p0_unscaled(p0_unscaled+p1_unscaled) # scale the probability
    draw <- runif(1)
    if (draw < p0) {
      current_l[i] = 0
    }
    else {
      current_l[i] = 1
    }
  }
  return (current_l)
}