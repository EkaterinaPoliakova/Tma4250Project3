complit_vec = as.vector(t(complit_mat))
dim <- 75 #dimension of square matrix
n <- length(seismic$V1) #length of l-vector
M <- 1000 #number of iteration
MCMC_output <- matrix(0,nrow=n,ncol=M)
convergence_vector <- vector(0,length=M)
sigma <- 0.06 #standard deviation
mean_vector <-  c(0.02,0.08) #mean, different for p0 and p1
beta_0 <- 1 #intital guess
possible_neigh <- matrix(nrow=n,ncol=8)
for (i in (1:n)) {
  possible_neigh[i,] <- c(i-1,i+1,i-dim,i+dim,i+1+dim,i+1-dim,i+dim-1,i-dim-1)
}
current_l <- complit_vec #initial value in the grid

####### HELP FUNCTIONS######

accept_change <- function(p_new,p_old) {
  accepted <- FALSE
  exponent <- min(0, p_new-p_old)
  prob <- exp(exponent)
  draw <- runif(1)
  if (draw < prob) {
    accepted <- TRUE
  }
  return (accepted)
} #accepts change

find_neighbours <- function(i,dim,n,possible_neigh) {
  valid_neigh <- vector()
  for (j in (1:8)) {
    if (possible_neigh[j] > 0 && possible_neigh[j] < n) { #checks if in domain
      valid_neigh <- append(valid_neigh,possible_neigh[j])
    }
    return (valid_neigh)
  }
}
#finds only valid neighbours (checks if in domain)

evaluate_neigbours <- function(i,dim,n,possible_neigh,current_l) {
  valid_neigh <- find_neighbours(i,dim,n,possible_neigh)
  neighbour_sum <- 0
  for (j in (1:length(valid_neigh))) {
    neighbour_sum <- neighbour_sum + current_l[j]
  }
  return (list(neighbour_sum,length(valid_neigh)-neigbour_sum))
} 
#returns number of 1 and 0-neighbours

p_d_l <- function(current_l,n,mean_vector,sigma,beta) {
  p_d_l <- matrix(nrow=n,ncol=2)
  for (i in 1:n) {
    p_di = dnorm(current_l[i], mean = mean_vector, sd = sigma)# return vector-valued (p0, p1)
    neig_eval <- evaluate_neigbours(i,dim,n,possible_neigh[i,],current_l)
    # evaluate_neighbours returns vector with number with 1 and 0 neighbours
    n0 <- neigh_eval[2]
    n1 <- neig_eval[1]
    p_l0 <- beta^n0 #effect from neighbours
    p_l1 <- beta^n1 #effect from neighbours
    p0_unscaled <- p_l0*p_di[1] #multiply contributions
    p1_unscaled <- p_l1*p_di[2]
    p0 <- p0_unscaled/(p0_unscaled+p1_unscaled) # scale the probability
    p1 <- 1 - p0
    p_d_l[i,] <- c(p0,p1)
  }
  return (p_d_l)
} 
#returns matrix with probablities for cases 0 and 1

logbeta <- function(beta) {
  totalsum <- 0
  for (i in (1:n)) {  #adds contributions node-wise
    number_neigh <- evaluate_neigbours(i,dim,n,possible_neigh,current_l)
    totalsum <- totalsum + log(p_d_l[i,2])+log(p_d_l[i,1]) +
      (number_neigh[1] + number_neigh[2])*log(beta) #number of neighbours
  }
  return (totalsum)
}
#returns loglikelihood, to be optimized

##### MCMC-STEP FUNCTION####

one_step <- function(node_change,current_l,mean,sigma,beta) {
  number_of_changes <- 0 #keeps track of convergence
  for (i in (1:n)) {
    p_di = dnorm(current_l[i], mean = mean_vector, sd = sigma)# return vector-valued (p0, p1)
    neig_eval <- evaluate_neigbours(i,dim,n,possible_neigh[i,],current_l)
    # evaluate_neighbours returns vector with number with 1 and 0 neighbours
    n0 <- neigh_eval[2]
    n1 <- neig_eval[1]
    p_l0 <- beta^n0 #effect from neighbours
    p_l1 <- beta^n1 #effect from neighbours
    p0_unscaled <- p_l0*p_di[1] #multiply contributions
    p1_unscaled <- p_l1*p_di[2]
    p0 <- p0_unscaled/(p0_unscaled+p1_unscaled) # scale the probability
    p1 <- 1-p0
    #Is p_old the p for current state, and p_new the p for the other state?
    change <- accept_change(p_new,p_old) #boolean return, is change accepted?
    if (change) {
      # CODE BELOW THIS SHOULD BE CHANGED
      draw <- runif(1) #draw random number between 0 and 1
      new_value <- (draw > p0) #TRUE=1, FALSE = 0
      if (new_value) {
        if (current_l[node]!= 1) {
          current_l[node] <- 1
          number_of_changes <- number_of_changes + 1
        }
      }
      else {
        if (current_l[node]!= 0) {
          current_l[node] <- 0
          number_of_changes <- number_of_changes + 1
        }
      }
    }
  }
  return (list(current_l,number_of_changes))
}
######PROCEDURE#######

beta <- optim(beta_0, logbeta, method = 'Brent', upper = 10, lower = 0,)#estimated value of beta

for (j in (1:M)) { #MCMC-sampling (M steps)
  node_change <-round(n*runif(1)) #draws random node
  step <- one_step(node_change,current_l,mean,sigma,beta) #perform one step
  current_l <- step[1] #updates the grid
  convergence_vector[j] <- step[2] #stores the number of changes
  MCMC_output[,j] <- current_l #stores current value in column
}

plot(rep(1,M),convergence_vector) #plots convergence
# last output grid is stored in MCMC_output[,M]

## TODO: Add plotting of grid after MCMC-sampling


