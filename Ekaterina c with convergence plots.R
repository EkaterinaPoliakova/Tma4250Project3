library(ggplot2)
library(reshape2)
library(purrr) #for binomial sample
seismic=read.delim("seismic.dat",sep=" ",header = FALSE)
complit=read.delim("complit.dat",sep=" ",header = FALSE)
complit_matrix <- as.matrix(complit)
complit_vec = as.vector((complit_matrix))
dim <- sqrt(length(seismic$V1)) #dimension of square matrix
n <- length(seismic$V1) #length of l-vector
M <- 50 #number of iteration
MCMC_output <- matrix(0,nrow=n,ncol=M) 
convergence_vector <- replicate(M,0) 
sigma <- 0.06 #standard deviation
mean_vector <-  c(0.02,0.08) #mean, different for p0 and p1
beta_0 <- 1 #intital guess
Compute_possible_neighours<-function(n)
{
dim=round(sqrt(n))  
possible_neigh <- matrix(nrow=n,ncol=8)
  for (i in (1:n)) {
    if ( i== 1) { #lowerleft corner
      possible_neigh[i,] <- c(2,1+dim,2+dim,dim,dim*2,dim*dim,dim*dim-dim+1,dim*dim-dim+2)
    } else
    if ( i== dim) { #lowerright corner
      possible_neigh[i,] <- c(1,1+dim,dim-1,dim*2-1,dim*2,dim*dim,dim*dim-1,dim*dim-dim+1)
    } else
    if ( i== dim*dim-dim+1) { #upperleft corner
      possible_neigh[i,] <- c(1,2,dim,dim*dim,dim*dim-dim,dim*dim-dim+2,dim*dim-dim*2+2,dim*dim-dim*2+1)
    } else
    if ( i== dim*dim) { #upperright corner
      possible_neigh[i,] <- c(1,dim-1,dim,dim*dim-dim,dim*dim-dim-1,dim*dim-dim*2+1,dim*dim-dim+1,dim*dim-1)
    } else
    if (i < dim+1) { #bottom
      possible_neigh[i,] <- c(dim*dim+i-dim-1,dim*dim+i-dim,dim*dim+i-dim+1,i-1,i+1,i+dim-1,i+dim,i+dim+1)
    }  else
    if (i > dim*dim-dim) { #top
      possible_neigh[i,] <- c(i-dim-1,i-dim,i-dim+1,i-1,i+1,i-dim*dim+dim-1,i-dim*dim+dim,i-dim*dim+dim+1)
    } else
    if (i %% dim == 1) { #left
      possible_neigh[i,] <- c(i-1,i-dim,i-dim+1,i-1+dim,i+1,i+dim*2-1,i+dim,i+dim+1)
    } else
    if (i %% dim == 0) { #right
      possible_neigh[i,] <- c(i-dim-1,i-dim,i-dim*2+1,i-1,i-dim+1,i+dim-1,i+dim,i+1)
    }
    else
    possible_neigh[i,] <- c(i-dim-1,i-dim,i-dim+1,i-1,i+1,i+dim-1,i+dim,i+dim+1)
  }
return(possible_neigh)  
}

#####we use one of simulations from b as initial value for lithology
set.seed(1) #reproducable results
N <- length(seismic$V1)
p0_gaussian = dnorm(seismic$V1, mean=0.02, sd=0.06) #Gaussian sample
p1_gaussian = dnorm(seismic$V1, mean=0.08, sd=0.06)
pd = 0.5 * (p0_gaussian + p1_gaussian) #Equal probability for 0 and 1
prob = 0.5 * p1_gaussian / pd
realizations = replicate(6,(as.integer(rbernoulli(N,prob)))) #N binomial samples with p = prob
#current_l <- complit_vec #initial value in the grid
current_l <-realizations[,1]
current_s <-as.vector(as.matrix(seismic$V1))
####### HELP FUNCTIONS######

accept_change <- function(p_new,p_old) {
  accepted <- FALSE
  exponent <- min(0, p_new-p_old)
  prob <- exp(exponent)
  draw <- runif(1)
  if (is.numeric(prob))
  if (draw < prob) {
    accepted <- TRUE
  }
  return (accepted)
} #accepts change

find_neighbours <- function(i,dim,n) {
  return (possible_neigh[i,])
}
#finds only valid neighbours (currently only checks if in domain, can be expanded)

evaluate_neigbours <- function(i,dim,n,current_l) {
  valid_neigh <- find_neighbours(i,dim,n)
  neighbour_sum <- 0
  for (j in (1:length(valid_neigh))) {
    neighbour_sum <- neighbour_sum + current_l[valid_neigh[j]]
  }
  return (c(neighbour_sum,length(valid_neigh)-neighbour_sum)) #first value gives number of 1s
} 
#returns number of 1 and 0-neighbours

pdl_func <- function(current_l,n,mean_vector,sigma,beta) {
  p_d_l <- matrix(nrow=n,ncol=2)
  for (i in (1:n)) {
    p_di = dnorm(current_s[i], mean = mean_vector, sd = sigma)# return vector-valued (p0, p1)
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
 # p_d_l[1,] <- p_d_l [2,] #for now, something is wrong...
  return (p_d_l)
} 
#returns matrix with probablities for cases 0 and 1

pdl_onenode <- function(current_l,n,mean_vector,sigma,beta,node,p_d_l) {
  nodes_to_change=c(node,find_neighbours(node,dim,n))
  p_di1 = dnorm(current_s[nodes_to_change], mean = mean_vector[1], sd = sigma)# return vector-valued (p0, p1)
  p_di2 = dnorm(current_s[nodes_to_change], mean = mean_vector[2], sd = sigma)# return vector-valued (p0, p1)
    for (j in length(nodes_to_change))
  {  
  neigh_eval <- evaluate_neigbours(j,dim,n,current_l)
  # evaluate_neighbours returns vector with number with 1 and 0 neighbours
  p_l0 <- beta^neigh_eval[2] #effect from neighbours
  p_l1 <- beta^neigh_eval[1] #effect from neighbours
  p0_unscaled <- p_l0*p_di1[j] #multiply contributions
  p1_unscaled <- p_l1*p_di2[j]
  p0 <- p0_unscaled/(p0_unscaled+p1_unscaled) # scale the probabilityb
  p1 <- 1 - p0
  p_d_l[nodes_to_change[j],] <- c(p0,p1)
  }
  return (p_d_l)
}
#only calulates and updates p for changed node


pdl_no_seismic_data <- function(current_l,beta) {
  n=length(current_l)
  p_d_l <- matrix(nrow=n,ncol=2)
  dim=sqrt(n)
  for (i in (1:n)) {
    neigh_eval <- evaluate_neigbours(i,dim,n,current_l)
    # evaluate_neighbours returns vector with number with 1 and 0 neighbours
    p0_unscaled <-  beta^neigh_eval[2] #effect from neighbours
    p1_unscaled <-  beta^neigh_eval[1] #effect from neighbours
    p0 <- p0_unscaled/(p0_unscaled+p1_unscaled) # scale the probabilityb
    p1 <- 1 - p0
    p_d_l[i,] <- c(p0,p1)
  }
  # p_d_l[1,] <- p_d_l [2,] #for now, something is wrong...
  return (p_d_l)
} 

possible_neigh=Compute_possible_neighours(length(complit_vec))

logbeta <- function(beta) {
  p_d_l <- pdl_no_seismic_data(complit_vec,beta)
  totalsum<-sum(log(p_d_l[,1])*(1-complit_vec))+sum(log(p_d_l[,2])*(complit_vec))
  return (totalsum)
}

#log_pseudolik <- function(beta,current_l) {
#  p_d_l <- pdl_func(current_l,length(complit_vec),mean_vector,sigma,beta)
#  totalsum<-sum(log(p_d_l[,1])*(1-complit_vec))+sum(log(p_d_l[,2])*(complit_vec))
#  return (c(totalsum,p_d_l[,2]))
#}
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
  p_d_l <- pdl_func(current_l,n,mean_vector,sigma,beta)
  p_old <- sum(log(p_d_l[,1])*(1-current_l))+sum(log(p_d_l[,2])*(current_l))
  current_l_candidate=current_l
  neighb=evaluate_neigbours(node_change,dim,n,current_l)
  Seism_eff= dnorm(current_s[node_change], mean = mean_vector, sd = sigma)
  p_new1_unscaled=Seism_eff[2]*beta^neighb[1]
  p_new0_unscaled=Seism_eff[1]*beta^neighb[2]
  node_new=rbinom(1,1,p_new1_unscaled/(p_new1_unscaled+p_new0_unscaled))
  current_l_candidate[node_change]=node_new
  p_d_l_cand=pdl_onenode(current_l_candidate,n,mean_vector,sigma,beta,node_change,p_d_l)
  p_new<- sum(log(p_d_l_cand[,1])*(1- current_l_candidate))+sum(log(p_d_l_cand[,2])*( current_l_candidate))
    change <- accept_change(p_new,p_old) #boolean return, is change accepted?
    if (change) {
      current_l[node_change] <-node_new
      }
  return (current_l)
}
######PROCEDURE#######
#current_l <- complit_vec
#beta <- optim(par=10^3, logbeta, method = 'Brent', upper = 10, 
     #         lower = 10^(-5))
#beta <- beta$par#estimated value of beta
#beta
possible_neigh=Compute_possible_neighours(length(complit_vec))
betas=seq(from=0.1,to=15,by=0.1)#rough maximizing
ybetas=rep(0,length(betas))
for (i in 1:length(betas))
  ybetas[i]=logbeta(betas[i])
plot(betas,ybetas,type='l') #check that it looks reasonably

#finer maximizing
betas=seq(from=betas[which.max(ybetas)]-0.1,to=betas[which.max(ybetas)]+0.1,by=0.001)
ybetas=rep(0,length(betas))
for (i in 1:length(betas))
  ybetas[i]=logbeta(betas[i])
beta=betas[which.max(ybetas)]#2.325
obs_matrix <- matrix (0,nrow=6,ncol=n)

for (ii in 1:6)
{
current_l <-realizations[,ii] #initial value in the grid
possible_neigh=Compute_possible_neighours(length(current_l))
n=length(current_l)

#convergence_vector <- replicate(M,0) 
M=70000
Sand_proportion=matrix(nrow=M,ncol=6,data=0)
for (j in (1:M)) { #MCMC-sampling (M steps)
  if(j %% 1000 ==0) print(paste("Step",j))
  node_change <-round(n*runif(1)+0.5) #draws random node, avoid "node Nr 0"
  current_l <- one_step(node_change,current_l,mean,sigma,beta)  #updates the grid
  Sand_proportion[j,ii] <- 1-sum(current_l)/n
}
xlit=rep((1:75),rep(75,75))
ylit=rep((1:75),75)
complit.df =data.frame(lythology=current_l,x=xlit,y=ylit)
lab=c("shale","sand")
newplot = #should change the visual look of this
  ggplot(complit.df, main="Text", aes(x=x, y=y)) + 
  coord_fixed(ratio = 1)+
  theme_minimal() +
  geom_tile(aes(fill=as.factor(lythology)))+
  scale_fill_manual(name = "Lythology",breaks = waiver(),labels=rev(lab),  values=c(0.1,1.1),  
                    na.value="black")
pdf(paste("realisation",ii,".pdf"), width = 4, height = 4)
print(newplot)
dev.off()

pdf(paste("Convergence_plot",ii,".pdf"), width = 4, height = 4)
plot(1:M,Sand_proportion[1:M,ii],type='l',xlab="Markov step",ylab="Proportion of sand",main=paste("Realization",i))
dev.off()
obs_matrix[ii,]  <- current_l

}
for (i in 1:6)
pdf(paste("Convergence_plot",i,".pdf"), width = 4, height = 4)
plot(1:M,Sand_proportion[1:M,i],type='l',xlab="Markov step",ylab="Proportion of sand",main=paste("Realization",i))
dev.off()

mean_posterior <- rowMeans(obs_matrix) #n-dimensional mean vector
#Then ggplot mean_posterior

variance_posterior <- mean_posterior*(1-mean_posterior) #variance of each node
#Then ggplot variance_posterior
