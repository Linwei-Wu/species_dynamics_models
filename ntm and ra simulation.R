
setwd("/Users/linwei/Dropbox/bioreactor/ratio model-new/1.3_ntm and neutral distribution")

C1<-read.csv("/Users/linwei/Dropbox/bioreactor/ratio model-new/1_model fitting/C1_fitting.csv",row.names = 1)
C2<-read.csv("/Users/linwei/Dropbox/bioreactor/ratio model-new/1_model fitting/C2_fitting.csv",row.names = 1)
C3<-read.csv("/Users/linwei/Dropbox/bioreactor/ratio model-new/1_model fitting/C3_fitting.csv",row.names = 1)

C1.NTm=(C1$x1.neutral/C1$SE.neutral/C1$SE.neutral+1)/C1$mean.ra
C2.NTm=(C2$x1.neutral/C2$SE.neutral/C2$SE.neutral+1)/C2$mean.ra
C3.NTm=(C3$x1.neutral/C3$SE.neutral/C3$SE.neutral+1)/C3$mean.ra



C1ntm<-data.frame(Ntm=C1.NTm, meanra=C1$mean.ra,modelP=C1$modP.neutral,reactor=rep("C1",nrow(C1)))
C2ntm<-data.frame(Ntm=C2.NTm, meanra=C2$mean.ra,modelP=C2$modP.neutral,reactor=rep("C2",nrow(C2)))
C3ntm<-data.frame(Ntm=C3.NTm, meanra=C3$mean.ra,modelP=C3$modP.neutral,reactor=rep("C3",nrow(C3)))
ntmdat<-rbind(C1ntm,C2ntm,C3ntm)
ntmdat<-ntmdat[!is.na(ntmdat$modelP),]

##1. plot ntm vs. relative abundance

library(ggplot2)

ggplot(ntmdat[ntmdat$modelP<0.05,], aes(x=meanra, y=Ntm))+
  geom_point(color="#1B1919FF",alpha=0.8)+
  geom_smooth(color="#BB0021FF")+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1))+
  scale_y_log10()+
  theme_bw()+
  facet_grid(~ reactor)

used<-ntmdat[ntmdat$modelP<0.05,]
cor.test(x=used$meanra[used$reactor=="C1"],y=used$Ntm[used$reactor=="C1"],method ="spearman")
###rho = -0.9192943, S = 937365594, p-value < 2.2e-16

cor.test(x=used$meanra[used$reactor=="C2"],y=used$Ntm[used$reactor=="C2"],method ="spearman")
###rho = -0.9431371,S = 669669014, p-value < 2.2e-16

cor.test(x=used$meanra[used$reactor=="C3"],y=used$Ntm[used$reactor=="C3"],method ="spearman")
###rho = -0.9475131,S = 1155702534, p-value < 2.2e-16

##2. individual species 
#read observed ra data##
otutab<-read.csv("/Users/linwei/Dropbox/bioreactor/ratio model-new/5_R code/otutab.csv",row.names = 1,check.names = F)
env<-read.csv("/Users/linwei/Dropbox/bioreactor/ratio model-new/5_R code/env.csv",row.names = 1)
otutab<-otutab[,match(row.names(env),colnames(otutab))]
comm<-t(otutab);
comm<-comm/rowSums(comm)
sum(row.names(comm)==row.names(env))

C1<-data.frame(NTm=C1.NTm,C1)
C2<-data.frame(NTm=C2.NTm,C2)
C3<-data.frame(NTm=C3.NTm,C3)

library(ggplot2)
##OTU2## #C1
C1["OTU_2",]
ntm=48.98771
pi=0.05744566
C1otu2<-comm[env$reactor=="C1","OTU_2"]

x<-seq(0,1,length=1000)
y<-dbeta(x,shape1=ntm*pi,shape2=ntm*(1-pi))

predict<-data.frame(x=x,y=y)
observed<-data.frame(OTU_2=C1otu2,reactors=rep("C1",length(C1otu2)))

ggplot(data=observed,aes(x=OTU_2))+
  geom_histogram(aes(y = ..density..),alpha=0.8,binwidth =0.01,fill="#1B1919FF")+  
  geom_area(data=predict,aes(x=x,y=y),fill="#3C5488FF",alpha=0.5)+
  geom_line(data=predict,aes(x=x,y=y,ymax=y),colour="#3C5488FF")+
  xlim(0,0.2)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



##OTU2## #C2
C2["OTU_2",]
ntm=136.511
pi=0.05156526
C1otu2<-comm[env$reactor=="C2","OTU_2"]

x<-seq(0,1,length=1000)
y<-dbeta(x,shape1=ntm*pi,shape2=ntm*(1-pi))

predict<-data.frame(x=x,y=y)
observed<-data.frame(OTU_2=C1otu2,reactors=rep("C1",length(C1otu2)))

ggplot(data=observed,aes(x=OTU_2))+
  geom_histogram(aes(y = ..density..),alpha=0.8,binwidth =0.01,fill="#1B1919FF")+  
  geom_area(data=predict,aes(x=x,y=y),fill="#3C5488FF",alpha=0.5)+
  geom_line(data=predict,aes(x=x,y=y,ymax=y),colour="#3C5488FF")+
  xlim(0,0.2)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

##OTU2## #C3
C3["OTU_2",]
ntm=65.47807
pi=0.04919329
C1otu2<-comm[env$reactor=="C3","OTU_2"]

x<-seq(0,1,length=1000)
y<-dbeta(x,shape1=ntm*pi,shape2=ntm*(1-pi))

predict<-data.frame(x=x,y=y)
observed<-data.frame(OTU_2=C1otu2,reactors=rep("C1",length(C1otu2)))

ggplot(data=observed,aes(x=OTU_2))+
  geom_histogram(aes(y = ..density..),alpha=0.8,binwidth =0.01,fill="#1B1919FF")+  
  geom_area(data=predict,aes(x=x,y=y),fill="#3C5488FF",alpha=0.5)+
  geom_line(data=predict,aes(x=x,y=y,ymax=y),colour="#3C5488FF")+
  xlim(0,0.2)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

###C1, OTU_102
C1["OTU_102",]
ntm=924.7942
pi=0.001677032
C1otu2<-comm[env$reactor=="C1","OTU_102"]

x<-seq(0,1,length=10000)
y<-dbeta(x,shape1=ntm*pi,shape2=ntm*(1-pi))

predict<-data.frame(x=x,y=y)
observed<-data.frame(OTU_2=C1otu2,reactors=rep("C2",length(C1otu2)))

ggplot(data=observed,aes(x=OTU_2))+
  geom_histogram(aes(y = ..density..),alpha=0.8,binwidth=0.0005,fill="#1B1919FF")+  
  geom_area(data=predict,aes(x=x,y=y),fill="#3C5488FF",alpha=0.5)+
  geom_line(data=predict,aes(x=x,y=y,ymax=y),colour="#3C5488FF")+
  xlim(0,0.0075)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

###C2, OTU_102
C2["OTU_102",]
ntm=1356.808
pi=0.001622726
C1otu2<-comm[env$reactor=="C2","OTU_102"]

x<-seq(0,1,length=10000)
y<-dbeta(x,shape1=ntm*pi,shape2=ntm*(1-pi))

predict<-data.frame(x=x,y=y)
observed<-data.frame(OTU_2=C1otu2,reactors=rep("C2",length(C1otu2)))

ggplot(data=observed,aes(x=OTU_2))+
  geom_histogram(aes(y = ..density..),alpha=0.8,binwidth=0.0005,fill="#1B1919FF")+  
  geom_area(data=predict,aes(x=x,y=y),fill="#3C5488FF",alpha=0.5)+
  geom_line(data=predict,aes(x=x,y=y,ymax=y),colour="#3C5488FF")+
  xlim(0,0.0075)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

###C3, OTU_102
C3["OTU_102",]
ntm=1011.033
pi=0.002037319
C1otu2<-comm[env$reactor=="C3","OTU_102"]

x<-seq(0,1,length=10000)
y<-dbeta(x,shape1=ntm*pi,shape2=ntm*(1-pi))

predict<-data.frame(x=x,y=y)
observed<-data.frame(OTU_2=C1otu2,reactors=rep("C2",length(C1otu2)))

ggplot(data=observed,aes(x=OTU_2))+
  geom_histogram(aes(y = ..density..),alpha=0.8,binwidth=0.0005,fill="#1B1919FF")+  
  geom_area(data=predict,aes(x=x,y=y),fill="#3C5488FF",alpha=0.5)+
  geom_line(data=predict,aes(x=x,y=y,ymax=y),colour="#3C5488FF")+
  xlim(0,0.0075)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


##C1 1002
C1["OTU_1003",]  ##OTU_1002 has too many missing values , thus use 1003 as the example instead
ntm=22918
pi=7.563576e-05
C1otu2<-comm[env$reactor=="C1","OTU_1003"]


x<-seq(0,1,length=200000)
y<-dbeta(x,shape1=ntm*pi,shape2=ntm*(1-pi))

predict<-data.frame(x=x,y=y)
observed<-data.frame(OTU_2=C1otu2,reactors=rep("C3",length(C1otu2)))

ggplot(data=observed,aes(x=OTU_2))+
  geom_histogram(aes(y = ..density..),alpha=0.8,binwidth=0.00002,fill="#1B1919FF")+  
  geom_area(data=predict,aes(x=x,y=y),fill="#3C5488FF",alpha=0.5)+
  geom_line(data=predict,aes(x=x,y=y,ymax=y),colour="#3C5488FF")+
  xlim(0,0.00025)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

##C2 1002
C2["OTU_1003",]  ##OTU_1002 has too many missing values , thus use 1003 as the example instead
ntm=22916
pi=8.714239e-05
C1otu2<-comm[env$reactor=="C2","OTU_1003"]


x<-seq(0,1,length=200000)
y<-dbeta(x,shape1=ntm*pi,shape2=ntm*(1-pi))

predict<-data.frame(x=x,y=y)
observed<-data.frame(OTU_2=C1otu2,reactors=rep("C3",length(C1otu2)))

ggplot(data=observed,aes(x=OTU_2))+
  geom_histogram(aes(y = ..density..),alpha=0.8,binwidth=0.00002,fill="#1B1919FF")+  
  geom_area(data=predict,aes(x=x,y=y),fill="#3C5488FF",alpha=0.5)+
  geom_line(data=predict,aes(x=x,y=y,ymax=y),colour="#3C5488FF")+
  xlim(0,0.00025)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


##C3 1002
C3["OTU_1003",]  ##OTU_1002 has too many missing values , thus use 1003 as the example instead
ntm=63036.67
pi=7.019597e-05
C1otu2<-comm[env$reactor=="C3","OTU_1003"]


x<-seq(0,1,length=200000)
y<-dbeta(x,shape1=ntm*pi,shape2=ntm*(1-pi))

predict<-data.frame(x=x,y=y)
observed<-data.frame(OTU_2=C1otu2,reactors=rep("C3",length(C1otu2)))

ggplot(data=observed,aes(x=OTU_2))+
  geom_histogram(aes(y = ..density..),alpha=0.8,binwidth=0.00002,fill="#1B1919FF")+  
  geom_area(data=predict,aes(x=x,y=y),fill="#3C5488FF",alpha=0.5)+
  geom_line(data=predict,aes(x=x,y=y,ymax=y),colour="#3C5488FF")+
  xlim(0,0.00025)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
