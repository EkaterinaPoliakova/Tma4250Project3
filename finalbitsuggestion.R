obs_matrix <- matrix (0,nrow=5000,ncol=n)
obs_matrix[i,] <- #realization[45001:50000,] stores values at converged stage
mean_posterior <- colMeans(obs_matrix) #n-dimensional mean vector
#Then ggplot mean_posterior

variance_posterior <- var(obs_matrix) #variance of each node
#Then ggplot variance_posterior





