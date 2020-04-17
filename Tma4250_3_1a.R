seismic=read.delim("seismic.dat",sep=" ",header = FALSE)
complit=read.delim("complit.dat",sep=" ",header = FALSE)
library(ggplot2)

lithology=as.matrix(complit)
#the plotting works very slowly
for(i in 1:dim(lithology)[1])
  for(j in 1:dim(lithology)[2])
  { if(lithology[i,j]==1)
  {plot(i,j,pch=12,col='wheat4', bg='wheat2',asp=1,xlim=c(1,dim(lithology)[1]),ylim=c(1,dim(lithology)[2]),xlab="",ylab="")
    par(new=TRUE)
   }
    if(lithology[i,j]==0)
    {plot(i,j+0.2,pch='.',col='darkorange3',asp=1,xlim=c(1,dim(lithology)[1]),ylim=c(1,dim(lithology)[2]),xlab="",ylab="")
      par(new=TRUE)
      plot(i-0.2,j-0.1,pch='.',col='darkorange3',asp=1,xlim=c(1,dim(lithology)[1]),ylim=c(1,dim(lithology)[2]),xlab="",ylab="")
      par(new=TRUE)
      plot(i+0.2,j-0.1,pch='.',col='darkorange3',asp=1,xlim=c(1,dim(lithology)[1]),ylim=c(1,dim(lithology)[2]),xlab="",ylab="")
      par(new=TRUE)
      
    }
   
  }
par(new=FALSE)
legend("bottomright",legend=c("sand","shale"),
       col=c('darkorange3','wheat4'), pch=c('.',15))
