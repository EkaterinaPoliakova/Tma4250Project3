library(ggplot2)
library(reshape2)
seismic=read.delim("seismic.txt",sep=" ",header = FALSE)
complit=read.delim("complit.txt",sep=" ",header = FALSE)
complit_matrix <- as.matrix(complit)
complit_vec = as.vector((complit_matrix))
dim <- 66 #dimension of square matrix
n <- length(complit_vec) #length of l-vector
M <- 50 #number of iteration
MCMC_output <- matrix(0,nrow=n,ncol=M) 
convergence_vector <- replicate(M,0) 
sigma <- 0.06 #standard deviation
mean_vector <-  c(0.02,0.08) #mean, different for p0 and p1
beta_0 <- 1 #intital guess
possible_neigh <- matrix(nrow=n,ncol=8)
for (i in (1:n)) {
  possible_neigh[i,] <- c(i-dim-1,i-dim,i-dim+1,i-1,i+1,i+dim-1,i+dim,i+dim+1)
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

find_neighbours <- function(i,dim,n) {
  valid_neigh <- vector()
  valid <- rep(TRUE,8)
  possible <- possible_neigh[i,]
  if (i < dim+1) { #bottom
    valid[1:3] <-FALSE
  }
  if (i > dim^dim-dim) { #top
    valid[6:8] <- FALSE
  }
  if (i %% dim == 1) { #left
    valid[1] <- FALSE
    valid[4] <- FALSE
    valid[6] <- FALSE
  }
  if (i %% dim == 0) { #right
    valid[3] <- FALSE
    valid[5] <- FALSE
    valid[7] <- FALSE
  }
  for (j in (1:8)) {
    if (valid[j]) {
      valid_neigh <- append(valid_neigh,possible[j])
    }
  }
  return (valid_neigh)
}
#finds only valid neighbours (currently only checks if in domain, can be expanded)

evaluate_neigbours <- function(i,dim,n,current_l) {
  valid_neigh <- find_neighbours(i,dim,n)
  neighbour_sum <- 0
  for (j in (1:length(valid_neigh))) {
    neighbour_sum <- neighbour_sum + current_l[j]
  }
  return (c(neighbour_sum,length(valid_neigh)-neighbour_sum))
} 
#returns number of 1 and 0-neighbours

pdl_func <- function(current_l,n,mean_vector,sigma,beta) {
  p_d_l <- matrix(nrow=n,ncol=2)
  for (i in (1:n)) {
    p_di = dnorm(current_l[i], mean = mean_vector, sd = sigma)# return vector-valued (p0, p1)
    neigh_eval <- evaluate_neigbours(i,dim,n,current_l)
    # evaluate_neighbours returns vector with number with 1 and 0 neighbours
    p_l0 <- beta^neigh_eval[2] #effect from neighbours
    p_l1 <- beta^neigh_eval[1] #effect from neighbours
    p0_unscaled <- p_l0*p_di[1] #multiply contributions
    p1_unscaled <- p_l1*p_di[2]
    p0 <- p0_unscaled/(p0_unscaled+p1_unscaled) # scale the probabilityb
    p1 <- 1 - p0
    p_d_l[i,] <- c(p0,p1)
  }
  p_d_l[1,] <- p_d_l [2,] #for now, something is wrong...
  return (p_d_l)
} 
#returns matrix with probablities for cases 0 and 1

pdl_onenode <- function(current_l,n,mean_vector,sigma,beta,node,p_d_l) {
  p_di = dnorm(current_l[node], mean = mean_vector, sd = sigma)# return vector-valued (p0, p1)
  neigh_eval <- evaluate_neigbours(node,dim,n,current_l)
  # evaluate_neighbours returns vector with number with 1 and 0 neighbours
  p_l0 <- beta^neigh_eval[2] #effect from neighbours
  p_l1 <- beta^neigh_eval[1] #effect from neighbours
  p0_unscaled <- p_l0*p_di[1] #multiply contributions
  p1_unscaled <- p_l1*p_di[2]
  p0 <- p0_unscaled/(p0_unscaled+p1_unscaled) # scale the probabilityb
  p1 <- 1 - p0
  p_d_l[node,] <- c(p0,p1)
  return (p_d_l)
}
#only calulates and updates p for changed node

logbeta <- function(beta) {
  p_d_l <- pdl_func(current_l,n,mean_vector,sigma,beta)
  totalsum <- 0
  for (i in (1:n)) {  #adds contributions node-wise
    number_neigh <- evaluate_neigbours(i,dim,n,current_l) #number of 1 and 0 neigbours
    if (current_l[i]==1) {
      totalsum <- totalsum + log(p_d_l[i,2]) + number_neigh[1]*log(beta)
    }
    else {
      totalsum <- totalsum + log(p_d_l[i,1]) + number_neigh[2]*log(beta)
    } 
  }
  return (totalsum)
}
#returns loglikelihood, to be optimized

find_p_old <- function(current_l,p_d_l,n) {
  output <- vector(length=n)
  for (i in 1:n) {
    if (current_l[i]==0) {
      output[i] <- log(p_d_l[i,1]) #log prob of state = 0
    }
    else {
      output[i] <- log(p_d_l[i,2]) #log prob of state = 1
    }
  }
  return (output)
}

##### MCMC-STEP FUNCTION####

one_step <- function(node_change,current_l,mean,sigma,beta) {
  number_of_changes <- 0 #keeps track of convergence
  p_d_l <- pdl_func(current_l,n,mean_vector,sigma,beta) #updates current probabilities
  p_old_vec <- find_p_old(current_l,p_d_l,n)
  p_old <- sum(p_old_vec)
  for (node in (1:n)) { # all n nodes
    p0 <- p_d_l[node,1] #current probability for 0
    p1 <- p_d_l[node,2] #current probability for 1
    sand <- (current_l[node] == 1) #boolean
    if (sand) {
      p_new <- p_old - p_old_vec[node] + log(p0) #prob for change
    }
    else {
      p_new <- p_old - p_old_vec[node] + log(p1) #prob for change
    }
    change <- accept_change(p_new,p_old) #boolean return, is change accepted?
    if (change) {
      number_of_changes <- number_of_changes + 1 #update number of changes
      if (sand) {
          current_l[node] <- 0 #change value to 0
         }
      else {
          current_l[node] <- 1 #change value to 1
      }
      p_d_l <- pdl_onenode(current_l,n,mean_vector,sigma,beta,node,p_d_l) #updates current probabilities
      p_old_vec <- find_p_old(current_l,p_d_l,n)
      p_old <- sum(p_old_vec)
    }
  }
  my_list <- list("grid" = current_l, "changes" = number_of_changes)
  return (my_list)
}
######PROCEDURE#######
current_l <- complit_vec
beta <- optim(par=10^3, logbeta, method = 'Brent', upper = 10^15, 
              lower = 10^(-5))
beta <- beta$par#estimated value of beta
beta
# THIS ONE IS WRONG

current_l <- complit_vec #initial value in the grid
convergence_vector <- replicate(M,0) 
for (j in (1:M)) { #MCMC-sampling (M steps)
  print("Step")
  print(j)
  node_change <-round(n*runif(1)) #draws random node
  step <- one_step(node_change,current_l,mean,sigma,0.8) #perform one step
  current_l <- step$grid #updates the grid
  convergence_vector[j] <- step$changes
  #MCMC_output[,j] <- current_l #stores current value in column
}
plot(seq(1, M),convergence_vector,type='l',main="Changes in each iteration"
     , xlab = "Iteration number", ylab="Changes") #plots convergence
final_l <- current_l # saves converges output


## TODO: Add plotting of grid after MCMC-sampling
final_df = melt(vapply(seq(1, dim), rep, rep(1,dim), times=dim)) #add new column
final_df$value = final_l
colnames(final_df) <- c("x","y","Value")
myplot =
  ggplot(final_df, aes(x=x, y=y, fill=Value)) +
  theme_minimal() +
  geom_raster() +
  scale_fill_gradient2()
print(myplot)

