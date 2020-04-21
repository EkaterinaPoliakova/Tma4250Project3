seismic=read.delim("seismic.txt",sep=" ",header = FALSE)
complit=read.delim("complit.txt",sep=" ",header = FALSE)
library(ggplot2)
library(reshape2)
library(purrr) #for binomial sample

################## OPPGÅVE A #################
seismic.df = melt(vapply(seq(1, 75), rep, rep(1,75), times=75)) #add new column
seismic.df$value = seismic$V1
colnames(seismic.df) <- c("x","y","Value")
myplot =
  ggplot(seismic.df, aes(x=x, y=y, fill=Value)) +
  theme_minimal() +
  geom_raster() +
  scale_fill_gradient2()
print(myplot)



############### OPPGÅVE B #############
set.seed(1) #reproducable results
N <- length(seismic$V1)
p0_gaussian = dnorm(seismic$V1, mean=0.02, sd=0.06) #Gaussian sample
p1_gaussian = dnorm(seismic$V1, mean=0.08, sd=0.06)
pd = 0.5 * (p0_gaussian + p1_gaussian) #Equal probability for 0 and 1
prob = 0.5 * p1_gaussian / pd
realizations = replicate(6,round(rbernoulli(N,prob),5)) #N binomial samples with p = prob

for (i in (1:6)) {
  seismic.df$real= realizations[,i]
  newplot = #should change the visual look of this
    ggplot(seismic.df, main="Text", aes(x=x, y=y, fill=real)) +
    theme_minimal() +
    geom_raster() +
    scale_fill_gradient2()
  print(newplot)
}
seismic.df$mean = prob #p
seismic.df$var = (c(rep(1,N))-prob)*prob #p(1-p)

meanplot =
  ggplot(seismic.df, aes(x=x, y=y, fill=mean)) +
  theme_minimal() +
  geom_raster() +
  scale_fill_gradient2()
print(meanplot)

varplot =
  ggplot(seismic.df, aes(x=x, y=y, fill=var)) +
  theme_minimal() +
  geom_raster() +
  scale_fill_gradient2()
print(varplot)