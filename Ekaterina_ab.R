seismic=read.delim("seismic.dat",sep=" ",header = FALSE)
complit=read.delim("complit.dat",sep=" ",header = FALSE)
library(ggplot2)
library(reshape2)
library(purrr) #for binomial sample

################## OPPG?VE A #################
seismic.df = melt(vapply(seq(1, 75), rep, rep(1,75), times=75)) #add new column
seismic.df$value = seismic$V1
colnames(seismic.df) <- c("x","y","Value")
myplot =
  ggplot(seismic.df, aes(x=x, y=y, fill=Value)) +
  coord_fixed(ratio = 1)+
  theme_minimal() +
  geom_raster() +
  scale_fill_gradient2()
print(myplot)


xlit=rep((1:66),rep(66,66))
ylit=rep((1:66),66)
complit.df =data.frame(lythology=as.vector(as.matrix(complit)),x=xlit,y=ylit)
View(complit.df)


lab=c("shale","sand")
#plotting lythology such way 
#works under Ubuntu and RStudio version 1.2.1355
#but doesn't work under windows 10
  newplot = #should change the visual look of this
    ggplot(complit.df, main="Text", aes(x=x, y=y)) + 
    coord_fixed(ratio = 1)+
    theme_minimal() +
    geom_tile(aes(fill=as.factor(lythology)))+
    scale_fill_manual(name = "Lythology",breaks = waiver(),labels=rev(lab),  values=c(0.1,1.1),  
                      na.value="black")
  print(newplot)




############### OPPG?VE B #############
set.seed(1) #reproducable results
N <- length(seismic$V1)
p0_gaussian = dnorm(seismic$V1, mean=0.02, sd=0.06) #Gaussian sample
p1_gaussian = dnorm(seismic$V1, mean=0.08, sd=0.06)
pd = 0.5 * (p0_gaussian + p1_gaussian) #Equal probability for 0 and 1
prob = 0.5 * p1_gaussian / pd
realizations = replicate(6,(as.integer(rbernoulli(N,prob)))) #N binomial samples with p = prob


P <- list()
lab=c("shale","sand")
for (i in (1:6)) {
  i=1
  seismic.df$lythology= realizations[,i]
  newplot = #should change the visual look of this
    ggplot(seismic.df, main="Text", aes(x=x, y=y)) + 
    coord_fixed(ratio = 1)+
    theme_minimal() +
    geom_tile(aes(fill=as.factor(lythology)))+
   scale_fill_manual(name = "Lythology",breaks = waiver(),labels=rev(lab),  values=c(0.1,1.1),  
                      na.value="black")
  i=i+1
  print(newplot)

  P=c(P,list(newplot))
}

library(gridExtra)
library(grid)
figure=mygrid <- grid.arrange(P[[1]],P[[2]],P[[3]],P[[4]],P[[5]],P[[6]],ncol=3,nrow=2) #arrange 4 plots

print(figure)

seismic.df$mean = prob #p
seismic.df$var = (c(rep(1,N))-prob)*prob #p(1-p)

meanplot =
  ggplot(seismic.df, aes(x=x, y=y, fill=mean)) +
  coord_fixed(ratio = 1)+ 
  theme_minimal() +
  geom_raster() +
  scale_fill_gradient2()
print(meanplot)

varplot =
  ggplot(seismic.df, aes(x=x, y=y, fill=var)) +
  coord_fixed(ratio = 1)+
  theme_minimal() +
  geom_raster() +
  scale_fill_gradient2()
print(varplot)