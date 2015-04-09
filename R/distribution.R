require(mvtnorm)
require(reshape2)
require(ggplot2)
mean=c(0,1)
sigma=cbind(c(1, 0),c(0, 1))

x<-seq(-2,2,by=0.01)
y<-seq(-2,2,by=0.01)
z<-outer(x,y,function(x,y)dmvnorm(cbind(x,y),mean=mean,sigma=sigma))

mz<-melt(z)

mz$Var1<- mz$Var1/100 - 2
mz$Var2<- mz$Var2/100 - 2

colnames(mz)<-c("x","y","z")

mz$intensity<-round(mz$z*1000)

g1 <- ggplot(mz, aes(x,y)) + 
  geom_contour(aes(z=z),colour="black") + 
  theme_bw() + xlab(expression(z[t1])) + ylab(expression(z[t2])) +
  ggtitle(expression(paste(b[1]==0, ",", b[2]==1))) +
  scale_x_continuous(limits = c(-2, 2)) +
  scale_y_continuous(limits = c(-2, 2)) +
  theme(axis.title = element_text(size=15), title = element_text(size=13))
plot(g1)