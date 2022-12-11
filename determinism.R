
setwd("/Users/linwei/Dropbox/bioreactor/neutral modeling/ratio model/2_determinism")

C1<-read.csv("../1_model fitting/C1_fitting.csv",header = T,row.names = 1)
C2<-read.csv("../1_model fitting/C2_fitting.csv",header = T,row.names = 1)
C3<-read.csv("../1_model fitting/C3_fitting.csv",header = T,row.names = 1)
D97<-read.csv("../1_model fitting/D_before100 ra fitting.csv",header = T,row.names = 1)

otutab<-read.csv("/Users/linwei/Dropbox/bioreactor/neutral modeling/ratio model/0_raw data/otutab.csv",header = T,row.names = 1,check.names = F)
env<-read.csv("/Users/linwei/Dropbox/bioreactor/neutral modeling/ratio model/0_raw data/env.csv",header = T,row.names = 1)

otutab<-otutab[,match(row.names(env),colnames(otutab))]
sum(row.names(env)==colnames(otutab))  #check
otutab<-t(otutab)
otutab<-otutab/rowSums(otutab)   #convert to relative abundance

### C1 result
otus<-otutab[env$reactor=="C1",]
envs<-env[env$reactor=="C1",]
sum(row.names(otus)==row.names(envs))

#determinism=mu/sigma
#d log x/xf /dt=b0+b1*(1/x)+b2*(1/xf)+b3*R + error  #the combined model

C1dtm.mean=sapply(row.names(C1), function(x){
  Cresult=C1[x,]
  a=(1/Cresult$SE.combine)^2
  b0=Cresult$X.Intercept..combine
  b1=Cresult$x1.combine
  b2=Cresult$x2.combine
  b3=Cresult$x3.combine
  xf=mean(otus[,"OTU_1"])
  x1=otus[,x]
  R=envs$VS;names(R)=row.names(envs)
  tpoints<-sapply(1:length(R), function(i){
    dtm=a*(b0+b1/x1[i]+b2/xf+b3*R[i])/sqrt(2/x1[i]+2/xf)
    abs(dtm) 
  })
  tpoints
})

C1dtm=sapply(row.names(C1), function(x){
  Cresult=C1[x,]
  a=(1/Cresult$SE.combine)^2
  b0=Cresult$X.Intercept..combine
  b1=Cresult$x1.combine
  b2=Cresult$x2.combine
  b3=Cresult$x3.combine
  xf=otus[,"OTU_1"]
  x1=otus[,x]
  R=envs$VS;names(R)=row.names(envs)
  tpoints<-sapply(1:length(R), function(i){
    dtm=a*(b0+b1/x1[i]+b2/xf[i]+b3*R[i])/sqrt(2/x1[i]+2/xf[i])
    abs(dtm) 
  })
  tpoints
})

### C2 result
otus<-otutab[env$reactor=="C2",]
envs<-env[env$reactor=="C2",]
sum(row.names(otus)==row.names(envs))

C2dtm.mean=sapply(row.names(C2), function(x){
  Cresult=C2[x,]
  a=(1/Cresult$SE.combine)^2
  b0=Cresult$X.Intercept..combine
  b1=Cresult$x1.combine
  b2=Cresult$x2.combine
  b3=Cresult$x3.combine
  xf=mean(otus[,"OTU_1"])
  x1=otus[,x]
  R=envs$VS;names(R)=row.names(envs)
  tpoints<-sapply(1:length(R), function(i){
    dtm=a*(b0+b1/x1[i]+b2/xf+b3*R[i])/sqrt(2/x1[i]+2/xf)
    abs(dtm) 
  })
  tpoints
})

C2dtm=sapply(row.names(C2), function(x){
  Cresult=C2[x,]
  a=(1/Cresult$SE.combine)^2
  b0=Cresult$X.Intercept..combine
  b1=Cresult$x1.combine
  b2=Cresult$x2.combine
  b3=Cresult$x3.combine
  xf=otus[,"OTU_1"]
  x1=otus[,x]
  R=envs$VS;names(R)=row.names(envs)
  tpoints<-sapply(1:length(R), function(i){
    dtm=a*(b0+b1/x1[i]+b2/xf[i]+b3*R[i])/sqrt(2/x1[i]+2/xf[i])
    abs(dtm) 
  })
  tpoints
})

### C3 result
otus<-otutab[env$reactor=="C3",]
envs<-env[env$reactor=="C3",]
sum(row.names(otus)==row.names(envs))

C3dtm.mean=sapply(row.names(C3), function(x){
  Cresult=C3[x,]
  a=(1/Cresult$SE.combine)^2
  b0=Cresult$X.Intercept..combine
  b1=Cresult$x1.combine
  b2=Cresult$x2.combine
  b3=Cresult$x3.combine
  xf=mean(otus[,"OTU_1"])
  x1=otus[,x]
  R=envs$VS;names(R)=row.names(envs)
  tpoints<-sapply(1:length(R), function(i){
    dtm=a*(b0+b1/x1[i]+b2/xf+b3*R[i])/sqrt(2/x1[i]+2/xf)
    abs(dtm) 
  })
  tpoints
})

C3dtm=sapply(row.names(C3), function(x){
  Cresult=C3[x,]
  a=(1/Cresult$SE.combine)^2
  b0=Cresult$X.Intercept..combine
  b1=Cresult$x1.combine
  b2=Cresult$x2.combine
  b3=Cresult$x3.combine
  xf=otus[,"OTU_1"]
  x1=otus[,x]
  R=envs$VS;names(R)=row.names(envs)
  tpoints<-sapply(1:length(R), function(i){
    dtm=a*(b0+b1/x1[i]+b2/xf[i]+b3*R[i])/sqrt(2/x1[i]+2/xf[i])
    abs(dtm) 
  })
  tpoints
})

### D97 result
otus<-otutab[env$reactor %in% c("D1","D2","D3"),]
envs<-env[env$reactor %in% c("D1","D2","D3"),]
sum(row.names(otus)==row.names(envs))

D97dtm.mean=sapply(row.names(D97), function(x){
  Cresult=D97[x,]
  a=(1/Cresult$SE.combine)^2
  b0=Cresult$X.Intercept..combine
  b1=Cresult$x1.combine
  b2=Cresult$x2.combine
  b3=Cresult$x3.combine
  xf=mean(otus[,"OTU_1"])
  x1=otus[,x]
  R=envs$VS;names(R)=row.names(envs)
  tpoints<-sapply(1:length(R), function(i){
    dtm=a*(b0+b1/x1[i]+b2/xf+b3*R[i])/sqrt(2/x1[i]+2/xf)
    abs(dtm) 
  })
  tpoints
})

D97dtm=sapply(row.names(D97), function(x){
  Cresult=D97[x,]
  a=(1/Cresult$SE.combine)^2
  b0=Cresult$X.Intercept..combine
  b1=Cresult$x1.combine
  b2=Cresult$x2.combine
  b3=Cresult$x3.combine
  xf=otus[,"OTU_1"]
  x1=otus[,x]
  R=envs$VS;names(R)=row.names(envs)
  tpoints<-sapply(1:length(R), function(i){
    dtm=a*(b0+b1/x1[i]+b2/xf[i]+b3*R[i])/sqrt(2/x1[i]+2/xf[i])
    abs(dtm) 
  })
  tpoints
})

determinism.meanxf=rbind(C1dtm.mean,C2dtm.mean,C3dtm.mean,D97dtm.mean)
determinism.dynamicxf=rbind(C1dtm,C2dtm,C3dtm,D97dtm)
write.csv(determinism.meanxf,"determinism.meanxf.species.csv")
write.csv(determinism.dynamicxf,"determinism.dynamicxf.species.csv")

##calculate community mean determinism
determinism.meanxf<-t(determinism.meanxf)
determinism.dynamicxf<-t(determinism.dynamicxf)

comm<-t(otutab[,-1])
sum(row.names(determinism.meanxf)==row.names(comm))
sum(row.names(determinism.dynamicxf)==row.names(comm))

dtmnsm.meanxf.unweighted=colMeans(determinism.meanxf,na.rm=T)
dtmnsm.dynamicxf.unweighted=colMeans(determinism.dynamicxf,na.rm=T)

dtmnsm.meanxf.weighted=colSums(determinism.meanxf*comm,na.rm=T)/colSums(comm,na.rm=T)
dtmnsm.dynamicxf.weighted=colSums(determinism.dynamicxf*comm,na.rm=T)/colSums(comm,na.rm=T)

community.determinism=data.frame(dtmnsm.meanxf.unweighted,dtmnsm.meanxf.weighted,dtmnsm.dynamicxf.unweighted,dtmnsm.dynamicxf.weighted)

sum(row.names(community.determinism)==row.names(env))
sum(env$reactor %in% c("C1","C2","C3"))
sum(env$reactor %in% c("D1","D2","D3"))

community.determinism<-data.frame(community.determinism,env,treat=c(rep("Control",159),rep("Treatment",33)))

write.csv(community.determinism,"community determinism.csv")

###rare species have higher determinism?
###figures
###1.  species determinism 
env<-read.csv("/Users/linwei/Dropbox/bioreactor/ratio model-new/5_R code/env.csv",row.names = 1)
Cft<-read.csv("/Users/linwei/Dropbox/bioreactor/ratio model-new/1_model fitting/C_before100 ra fitting.csv",row.names = 1)
Dft<-read.csv("/Users/linwei/Dropbox/bioreactor/ratio model-new/1_model fitting/D_before100 ra fitting.csv",row.names = 1)
dynamic.ref<-read.csv("/Users/linwei/Dropbox/bioreactor/ratio model-new/2_determinism/determinism.dynamicxf.species.csv",row.names = 1)
mean.ref<-read.csv("/Users/linwei/Dropbox/bioreactor/ratio model-new/2_determinism/determinism.meanxf.species.csv",row.names = 1)

#check if sample and species names matched
sum(row.names(dynamic.ref)==row.names(env))
sum(row.names(mean.ref)==row.names(env))
sum(colnames(dynamic.ref)==row.names(Cft))
sum(colnames(mean.ref)==row.names(Cft))
sum(colnames(mean.ref)==row.names(Dft))

library(reshape2)
dynamic_C<-dynamic.ref[env$reactor %in% c("C1","C2","C3") & env$time<100 & env$time>44 & env$time!=59,]
dynamic_C<-t(dynamic_C)
dynamic_C<-data.frame(mean.ra=Cft$mean.ra,dynamic_C)
rabins<-cut(Cft$mean.ra,c(0, 0.0001, 0.001,1))
dynamic_C<-data.frame(rabins=rabins,treat=rep("C",nrow(dynamic_C)),dynamic_C)
meltC<-melt(dynamic_C,id.vars = c("rabins","treat","mean.ra"))
meltC<-meltC[!is.na(meltC$value),]

dynamic_D<-dynamic.ref[env$reactor %in% c("D1","D2","D3") & env$time<100 & env$time>44 & env$time!=59,]
dynamic_D<-t(dynamic_D)
dynamic_D<-data.frame(mean.ra=Dft$mean.ra,dynamic_D)
rabins<-cut(Dft$mean.ra,c(0, 0.0001, 0.001,1))
dynamic_D<-data.frame(rabins=rabins,treat=rep("D",nrow(dynamic_D)),dynamic_D)
meltD<-melt(dynamic_D,id.vars = c("rabins","treat","mean.ra"))
meltD<-meltD[!is.na(meltD$value),]

allCD=rbind(meltC,meltD)
allCD<-allCD[!is.na(allCD$rabins),]

library(ggplot2)
ggplot(allCD, aes(x=treat, y=value))+
  geom_violin(aes(linetype=treat))+
  geom_dotplot(aes(fill=treat),binaxis='y', stackdir='center',position = "dodge",binwidth=0.55,dotsize = 1.5,colour=NA,binpositions = "all")+
  #geom_dotplot(aes(fill=treat),binaxis='y', stackdir='centerwhole',position = "dodge",binpositions = "all",binwidth = 0.045,colour=NA)+
  scale_fill_manual(values=c("C"="#999999", "D"="#E69F00"))+
  scale_linetype_manual(values=c("C"="dashed","D"= "solid"))+
  ylim(0,300)+
  xlab(NULL)+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  facet_grid(. ~ rabins)

cor.test(x=allCD$mean.ra[allCD$treat=="C"],y=allCD$value[allCD$treat=="C"],method ="spearman")
# S = 1.8531e+13, p-value < 2.2e-16
#       rho 
# -0.5299599
cor.test(x=allCD$mean.ra[allCD$treat=="D"],y=allCD$value[allCD$treat=="D"],method ="spearman")
#S = 1.6538e+13, p-value < 2.2e-16
#rho 
#-0.5542822 

wilcox.test(allCD$value[allCD$treat=="C" & allCD$rabins=="(0,0.0001]"], allCD$value[allCD$treat=="D" & allCD$rabins=="(0,0.0001]"])
t.test(allCD$value[allCD$treat=="C" & allCD$rabins=="(0,0.0001]"], allCD$value[allCD$treat=="D" & allCD$rabins=="(0,0.0001]"])

#W = 138452836, p-value = 9.942e-14
#t = -8.9344, df = 21845, p-value < 2.2e-16
#mean of x mean of y 
#152.1415  196.2012
wilcox.test(allCD$value[allCD$treat=="C" & allCD$rabins=="(0.0001,0.001]"], allCD$value[allCD$treat=="D" & allCD$rabins=="(0.0001,0.001]"])
t.test(allCD$value[allCD$treat=="C" & allCD$rabins=="(0.0001,0.001]"], allCD$value[allCD$treat=="D" & allCD$rabins=="(0.0001,0.001]"])

#W = 182687384, p-value = 0.002753
#t = -2.5259, df = 38485, p-value = 0.01154
#mean of x mean of y 
#53.97274  56.78406 
wilcox.test(allCD$value[allCD$treat=="C" & allCD$rabins=="(0.001,1]"], allCD$value[allCD$treat=="D" & allCD$rabins=="(0.001,1]"])
t.test(allCD$value[allCD$treat=="C" & allCD$rabins=="(0.001,1]"], allCD$value[allCD$treat=="D" & allCD$rabins=="(0.001,1]"])
#t = -8.7685, df = 8927.2, p-value < 2.2e-16
#W = 8436420, p-value < 2.2e-16
#mean of x mean of y 
#12.72425  15.86194 


ggplot(allCD[allCD$rabins=="(0,0.0001]",], aes(x=treat, y=value))+
  geom_violin(aes(linetype=treat))+
  geom_dotplot(aes(fill=treat),binaxis='y', stackdir='center',position = "dodge",binwidth=0.75,dotsize = 1.5,colour=NA,binpositions = "all")+
  #geom_dotplot(aes(fill=treat),binaxis='y', stackdir='centerwhole',position = "dodge",binpositions = "all",binwidth = 0.045,colour=NA)+
  scale_fill_manual(values=c("C"="#999999", "D"="#E69F00"))+
  scale_linetype_manual(values=c("C"="dashed","D"= "solid"))+
  ylim(0,300)+
  xlab(NULL)+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

##mean ref

mean_C<-mean.ref[env$reactor %in% c("C1","C2","C3") & env$time<100 & env$time>44 & env$time!=59,]
mean_C<-t(mean_C)
mean_C<-data.frame(mean.ra=Cft$mean.ra,mean_C)
rabins<-cut(Cft$mean.ra,c(0, 0.0001, 0.001,1))
mean_C<-data.frame(rabins=rabins,treat=rep("C",nrow(mean_C)),mean_C)
meltC<-melt(mean_C,id.vars = c("rabins","treat","mean.ra"))
meltC<-meltC[!is.na(meltC$value),]

mean_D<-mean.ref[env$reactor %in% c("D1","D2","D3") & env$time<100 & env$time>44 & env$time!=59,]
mean_D<-t(mean_D)
mean_D<-data.frame(mean.ra=Dft$mean.ra,mean_D)
rabins<-cut(Dft$mean.ra,c(0, 0.0001, 0.001,1))
mean_D<-data.frame(rabins=rabins,treat=rep("D",nrow(mean_D)),mean_D)
meltD<-melt(mean_D,id.vars = c("rabins","treat","mean.ra"))
meltD<-meltD[!is.na(meltD$value),]

allCD=rbind(meltC,meltD)
allCD<-allCD[!is.na(allCD$rabins),]

library(ggplot2)
ggplot(allCD, aes(x=treat, y=value))+
  geom_violin(aes(linetype=treat))+
  geom_dotplot(aes(fill=treat),binaxis='y', stackdir='center',position = "dodge",binwidth=0.55,dotsize = 1.5,colour=NA,binpositions = "all")+
  #geom_dotplot(aes(fill=treat),binaxis='y', stackdir='centerwhole',position = "dodge",binpositions = "all",binwidth = 0.045,colour=NA)+
  scale_fill_manual(values=c("C"="#999999", "D"="#E69F00"))+
  scale_linetype_manual(values=c("C"="dashed","D"= "solid"))+
  ylim(0,300)+
  xlab(NULL)+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  facet_grid(. ~ rabins)

cor.test(x=allCD$mean.ra[allCD$treat=="C"],y=allCD$value[allCD$treat=="C"],method ="spearman")
# S = 1.8531e+13, p-value < 2.2e-16
#       rho 
# -0.5559771 
cor.test(x=allCD$mean.ra[allCD$treat=="D"],y=allCD$value[allCD$treat=="D"],method ="spearman")
#S = 1.6455e+13, p-value < 2.2e-16
#rho 
#-0.5464876 

wilcox.test(allCD$value[allCD$treat=="C" & allCD$rabins=="(0,0.0001]"], allCD$value[allCD$treat=="D" & allCD$rabins=="(0,0.0001]"])
#W = 130583603, p-value < 2.2e-16
#mean of x mean of y 
#140.7962  188.2036  

wilcox.test(allCD$value[allCD$treat=="C" & allCD$rabins=="(0.0001,0.001]"], allCD$value[allCD$treat=="D" & allCD$rabins=="(0.0001,0.001]"])
#W = 168445765, p-value < 2.2e-16
#mean of x mean of y 
#47.28636  57.988666 
wilcox.test(allCD$value[allCD$treat=="C" & allCD$rabins=="(0.001,1]"], allCD$value[allCD$treat=="D" & allCD$rabins=="(0.001,1]"])
#W = 6641693, p-value < 2.2e-16
#mean of x mean of y 
#9.452811 19.106070 


ggplot(allCD[allCD$rabins=="(0,0.0001]",], aes(x=treat, y=value))+
  geom_violin(aes(linetype=treat))+
  geom_dotplot(aes(fill=treat),binaxis='y', stackdir='center',position = "dodge",binwidth=0.75,dotsize = 1.5,colour=NA,binpositions = "all")+
  #geom_dotplot(aes(fill=treat),binaxis='y', stackdir='centerwhole',position = "dodge",binpositions = "all",binwidth = 0.045,colour=NA)+
  scale_fill_manual(values=c("C"="#999999", "D"="#E69F00"))+
  scale_linetype_manual(values=c("C"="dashed","D"= "solid"))+
  ylim(0,300)+
  xlab(NULL)+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


####community determinism
comm.determinism<-read.csv("/Users/linwei/Dropbox/bioreactor/ratio model-new/2_determinism/community determinism.csv",row.names = 1)
env<-read.csv("/Users/linwei/Dropbox/bioreactor/ratio model-new/5_R code/env.csv",row.names = 1)
sum(row.names(comm.determinism)==row.names(env))

determ2<-comm.determinism[env$time<100 & env$time>44 & env$time!=59,] 
env2<-env[env$time<100 & env$time>44 & env$time!=59,]

library(ggpubr)

determ2$treat
ggline(determ2, x = "time", y = "dtmnsm.meanxf.unweighted", color = "treat",  
       add = c("mean_se", "jitter"),palette = c("#999999", "#E69F00"),add.params = list(group="treat"))
determ2$dtmnsm.meanxf.weighted
ggline(determ2, x = "time", y = "dtmnsm.meanxf.weighted", color = "treat",  
       add = c("mean_se", "jitter"),palette = c("#999999", "#E69F00"),add.params = list(group="treat"))
determ2$dtmnsm.dynamicxf.unweighted
ggline(determ2, x = "time", y = "dtmnsm.dynamicxf.unweighted", color = "treat",  
       add = c("mean_se", "jitter"),palette = c("#999999", "#E69F00"),add.params = list(group="treat"))

determ2$dtmnsm.dynamicxf.weighted
ggline(determ2, x = "time", y = "dtmnsm.dynamicxf.weighted", color = "treat",  
       add = c("mean_se", "jitter"),palette = c("#999999", "#E69F00"),add.params = list(group="treat"))

