
setwd("/Users/linwei/Dropbox/bioreactor/neutral modeling/new model/3_manuscript/20200115/R code")

otutab<-read.csv("otutab.csv",row.names = 1,check.names = F)
env<-read.csv("env.csv",row.names = 1)
sum(row.names(env)==colnames(otutab))  #check

comm<-t(otutab)
comm<-comm/rowSums(comm)      ##get the relative abundance

source("main modeling function/fitting and comparing models.R")

####1. modeling fitting on time series
##1.1 test on control reactors one by one
otus<-comm[env$reactor=="C1",] #C2, C3
Env<-env[env$reactor=="C1",]
tpoint<-Env$time

allresult<-one.t(otus=otus,tpoint = tpoint,Env = Env)
rownames(allresult)<-colnames(otus)
write.csv(allresult,"C1 ra fitting.csv")  ##the same for C2 and C3

##1.2 test on treatment reactors. Combine the three time series from D1, D2 and D3
D1otus<-comm[env$reactor=="D1",];D1Env<-env[env$reactor=="D1",]   
D1tpoint<-D1Env$time
D2otus<-comm[env$reactor=="D2",];D2Env<-env[env$reactor=="D2",]   
D2tpoint<-D2Env$time
D3otus<-comm[env$reactor=="D3",];D3Env<-env[env$reactor=="D3",]   
D3tpoint<-D3Env$time

allresult<-combine.t(D1otus,D1tpoint,D1Env,D2otus,D2tpoint,D2Env,D3otus,D3tpoint,D3Env)
rownames(allresult)<-colnames(D1otus)
write.csv(allresult,"D ra fitting_before 100.csv")  

##1.3 test on control reactors. Combine the three time series from C1, C2 
#and C3, which matched the time points in D reactors.
D1comm<-comm[env$reactor=="C1",]  #Although named as D1comm, it's data from the control reactor C1
D1env<-env[env$reactor=="C1",]
D2comm<-comm[env$reactor=="C2",]
D2env<-env[env$reactor=="C2",]
D3comm<-comm[env$reactor=="C3",]
D3env<-env[env$reactor=="C3",]
#C 45-97 days
D1id<-which(D1env$time<100 & D1env$time> 44 & (D1env$time !=59))  ##include all matched time points
D1otus<-D1comm[D1id,];D1Env<-D1env[D1id,]
D1tpoint<-D1Env$time

D2id<-which(D2env$time<100 & D2env$time> 44 & (D2env$time !=59))    ##include all matched time points
D2otus<-D2comm[D2id,];D2Env<-D2env[D2id,]
D2tpoint<-D2Env$time

D3id<-which(D3env$time<100 & D3env$time> 44 & (D3env$time !=59))   ##include all matched time points
D3otus<-D3comm[D3id,];D3Env<-D3env[D3id,]
D3tpoint<-D3Env$time

allresult<-combine.t(D1otus,D1tpoint,D1Env,D2otus,D2tpoint,D2Env,D3otus,D3tpoint,D3Env)
rownames(allresult)<-colnames(D1otus)
write.csv(allresult,"C ra fitting_before 100.csv")  

####2. plot the model fitting results
C1<-read.csv("C1 ra fitting.csv",header = T,row.names = 1)
C2<-read.csv("C2 ra fitting.csv",header = T,row.names = 1)
C3<-read.csv("C3 ra fitting.csv",header = T,row.names = 1)
D97<-read.csv("D ra fitting_before 100.csv",header = T,row.names = 1)
C97<-read.csv("C ra fitting_before 100.csv",header = T,row.names = 1)


aics<-data.frame(neutral=C1$aic.neutral,resource=C1$aic.resource,combine=C1$aic.combine)
models<-c("neutral","consumer-resource","combined")
C1model<-sapply(1:nrow(C1), function(i){
  result=ifelse(is.na(aics[i,1]), NA, models[which.min(aics[i,])])
  
})
table(C1model)

aics<-data.frame(neutral=C2$aic.neutral,resource=C2$aic.resource,combine=C2$aic.combine)
models<-c("neutral","consumer-resource","combined")
C2model<-sapply(1:nrow(C2), function(i){
  result=ifelse(is.na(aics[i,1]), NA, models[which.min(aics[i,])])
  
})
table(C2model)

aics<-data.frame(neutral=C3$aic.neutral,resource=C3$aic.resource,combine=C3$aic.combine)
models<-c("neutral","consumer-resource","combined")
C3model<-sapply(1:nrow(C3), function(i){
  result=ifelse(is.na(aics[i,1]), NA, models[which.min(aics[i,])])
  
})
table(C3model)

aics<-data.frame(neutral=D97$aic.neutral,resource=D97$aic.resource,combine=D97$aic.combine)
models<-c("neutral","consumer-resource","combined")
Dmodel<-sapply(1:nrow(D97), function(i){
  result=ifelse(is.na(aics[i,1]), NA, models[which.min(aics[i,])])
  
})
table(Dmodel)

aics<-data.frame(neutral=C97$aic.neutral,resource=C97$aic.resource,combine=C97$aic.combine)
models<-c("neutral","consumer-resource","combined")
Cmodel<-sapply(1:nrow(C97), function(i){
  result=ifelse(is.na(aics[i,1]), NA, models[which.min(aics[i,])])
  
})
table(Cmodel)

C1used<-data.frame(mean.ra=C1$mean.ra,r2.neutral=C1$r2psedo.neutral,r2.resource=C1$r2psedo.resource,r2.combine=C1$r2psedo.combine,modP.neutral=C1$modP.neutral,modP.resource=C1$modP.resource,modP.combine=C1$modP.combine,bestmodel=C1model,reactor=rep("C1",length(C1model)),treat=rep("C",length(C1model)))
C2used<-data.frame(mean.ra=C2$mean.ra,r2.neutral=C2$r2psedo.neutral,r2.resource=C2$r2psedo.resource,r2.combine=C2$r2psedo.combine,modP.neutral=C2$modP.neutral,modP.resource=C2$modP.resource,modP.combine=C2$modP.combine,bestmodel=C2model,reactor=rep("C2",length(C2model)),treat=rep("C",length(C2model)))
C3used<-data.frame(mean.ra=C3$mean.ra,r2.neutral=C3$r2psedo.neutral,r2.resource=C3$r2psedo.resource,r2.combine=C3$r2psedo.combine,modP.neutral=C3$modP.neutral,modP.resource=C3$modP.resource,modP.combine=C3$modP.combine,bestmodel=C3model,reactor=rep("C3",length(C3model)),treat=rep("C",length(C3model)))
Dused<-data.frame(mean.ra=D97$mean.ra,r2.neutral=D97$r2psedo.neutral,r2.resource=D97$r2psedo.resource,r2.combine=D97$r2psedo.combine,modP.neutral=D97$modP.neutral,modP.resource=D97$modP.resource,modP.combine=D97$modP.combine,bestmodel=Dmodel,reactor=rep("D",length(Dmodel)),treat=rep("D",length(Dmodel)))
Cused<-data.frame(mean.ra=C97$mean.ra,r2.neutral=C97$r2psedo.neutral,r2.resource=C97$r2psedo.resource,r2.combine=C97$r2psedo.combine,modP.neutral=C97$modP.neutral,modP.resource=C97$modP.resource,modP.combine=C97$modP.combine,bestmodel=Cmodel,reactor=rep("C",length(Dmodel)),treat=rep("C",length(Dmodel)))

alldata<-rbind(C1used,C2used,C3used,Dused)
alldata97<-rbind(Cused,Dused)  #match time points in C and D
alldata<-alldata[!is.na(alldata$mean.ra),]
alldata97<-alldata97[!is.na(alldata97$mean.ra),]

#####delete before uploading the R code####
rabins<-cut(alldata$mean.ra,c(0, 0.0001, 0.001,1)) 
alldatatest<-data.frame(alldata,rabins)
sum(alldatatest$reactor=="C1" & alldatatest$rabins=="(0.001,1]")
sum(alldatatest$reactor=="C1" & alldatatest$rabins=="(0.0001,0.001]")
sum(alldatatest$reactor=="C1" & alldatatest$rabins=="(0,0.0001]")

sum(alldatatest$reactor=="C2" & alldatatest$rabins=="(0.001,1]")
sum(alldatatest$reactor=="C2" & alldatatest$rabins=="(0.0001,0.001]")
sum(alldatatest$reactor=="C2" & alldatatest$rabins=="(0,0.0001]")

sum(alldatatest$reactor=="C3" & alldatatest$rabins=="(0.001,1]")
sum(alldatatest$reactor=="C3" & alldatatest$rabins=="(0.0001,0.001]")
sum(alldatatest$reactor=="C3" & alldatatest$rabins=="(0,0.0001]")

sum(alldatatest$reactor=="D" & alldatatest$rabins=="(0.001,1]")
sum(alldatatest$reactor=="D" & alldatatest$rabins=="(0.0001,0.001]")
sum(alldatatest$reactor=="D" & alldatatest$rabins=="(0,0.0001]")
#####delete the above code###



library(reshape2)
meltdata<-melt(alldata,id.vars = c("mean.ra","reactor","treat","bestmodel"))
measurement<-sapply(strsplit(as.character(meltdata$variable),"[.]"), function(x) x[1])
models<-sapply(strsplit(as.character(meltdata$variable),"[.]"), function(x) x[2])
meltdata<-data.frame(meltdata,measurement,models)

rabins<-cut(meltdata$mean.ra,c(0, 0.0001, 0.001,1))  ##classify to rare, moderate and abundant species based on relative abundance
meltdata<-data.frame(meltdata,rabins)
#write.csv(meltdata,"CD meltdata for plot.csv")

library(ggplot2)
meltdata$bestmodel=factor(meltdata$bestmodel,levels = c("neutral","consumer-resource","combined"))

##2.1 best model distribution (based on AIC values)
#simple stack plot so that inkscape can open it
meltdataC<-meltdata[meltdata$treat=="C",]
meltdataD<-meltdata[meltdata$treat=="D",]
tableC<-table(meltdataC$rabins,meltdataC$bestmodel)
tableC<-as.data.frame(tableC/rowSums(tableC))
tableD<-table(meltdataD$rabins,meltdataD$bestmodel)
tableD<-as.data.frame(tableD/rowSums(tableD))
tableC<-data.frame(tableC,treat=rep("C",9))
tableD<-data.frame(tableD,treat=rep("D",9))
tableCD<-rbind(tableC,tableD)
colnames(tableCD)<-c("rabins","bestmodel","Freq","treat")

tableCD$bestmodel=factor(tableCD$bestmodel,levels=c("neutral","consumer-resource","combined"))
ggplot(tableCD,aes(x = bestmodel, y = Freq,fill = bestmodel)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("#3C5488FF","#00A087FF","#E64B35FF"))+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())+
  facet_grid(treat ~ rabins)

## 2.2 P value distribution
#since p value is sensitive to sample numbers,
#we should use matched dataset (C and D have the same number of time points)
meltdata97<-melt(alldata97,id.vars = c("mean.ra","reactor","treat","bestmodel"))
measurement<-sapply(strsplit(as.character(meltdata97$variable),"[.]"), function(x) x[1])
models<-sapply(strsplit(as.character(meltdata97$variable),"[.]"), function(x) x[2])
meltdata97<-data.frame(meltdata97,measurement,models)
rabins<-cut(meltdata97$mean.ra,c(0, 0.0001, 0.001,1))  ##classify to rare, moderate and abundant species based on relative abundance
meltdata97<-data.frame(meltdata97,rabins)

pdata97<-meltdata97[meltdata97$measurement=="modP",]
pdata97$models<-factor(pdata97$models,levels = c("neutral","resource","combine"))

#rare species
ggplot(pdata97[pdata97$rabins=="(0,0.0001]",], aes(x=models, y=log10(value),alpha=treat))+
  geom_violin(aes(linetype=treat))+
  geom_dotplot(aes(fill=models),binaxis='y', stackdir='centerwhole',position = "dodge",binwidth = 0.033,dotsize = 2,colour=NA)+
  #geom_dotplot(aes(fill=models),binaxis='y', stackdir='center',position = "dodge",binpositions = "all",binwidth = 0.045,colour=NA)+
  geom_hline(yintercept = c(log10(0.05),log10(0.01)),linetype="dashed",color="grey30")+ 
  scale_fill_manual(values=c("#3C5488FF","#00A087FF","#E64B35FF"),labels=c("neutral","consumer-resource","combined"))+
  scale_alpha_manual(values=c("C"=0.3,"D"=1))+
  scale_linetype_manual(values=c("C"="dashed","D"= "solid"))+
  ylim(-6,0)+
  xlab(NULL)+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#moderate species
ggplot(pdata97[pdata97$rabins=="(0.0001,0.001]",], aes(x=models, y=log10(value),alpha=treat))+
  geom_violin(aes(linetype=treat))+
  geom_dotplot(aes(fill=models),binaxis='y', stackdir='centerwhole',position = "dodge",binwidth = 0.065,dotsize = 2,colour=NA)+
  #geom_dotplot(aes(fill=models),binaxis='y', stackdir='center',position = "dodge",binpositions = "all",binwidth = 0.045,colour=NA)+
  geom_hline(yintercept = c(log10(0.05),log10(0.01)),linetype="dashed",color="grey30")+ 
  scale_fill_manual(values=c("#3C5488FF","#00A087FF","#E64B35FF"),labels=c("neutral","consumer-resource","combined"))+
  scale_alpha_manual(values=c("C"=0.3,"D"=1))+
  scale_linetype_manual(values=c("C"="dashed","D"= "solid"))+
  ylim(-8,0)+
  xlab(NULL)+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#abundant species
ggplot(pdata97[pdata97$rabins=="(0.001,1]",], aes(x=models, y=log10(value),alpha=treat))+
  geom_violin(aes(linetype=treat))+
  geom_dotplot(aes(fill=models),binaxis='y', stackdir='centerwhole',position = "dodge",binwidth = 0.16,dotsize = 1.5,colour=NA)+
  #geom_dotplot(aes(fill=models),binaxis='y', stackdir='center',position = "dodge",binpositions = "all",binwidth = 0.045,colour=NA)+
  geom_hline(yintercept = c(log10(0.05),log10(0.01)),linetype="dashed",color="grey30")+ 
  scale_fill_manual(values=c("#3C5488FF","#00A087FF","#E64B35FF"),labels=c("neutral","consumer-resource","combined"))+
  scale_alpha_manual(values=c("C"=0.6,"D"=1))+
  scale_linetype_manual(values=c("C"="dashed","D"= "solid"))+
  ylim(-8,0)+
  xlab(NULL)+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

testdata<-pdata97[pdata97$rabins=="(0.001,1]",]
sum(testdata$value[testdata$models=="resource" & testdata$treat=="D"]<=0.05)/sum(testdata$models=="resource" & testdata$treat=="D")
sum(testdata$value[testdata$models=="resource" & testdata$treat=="C"]<=0.05)/sum(testdata$models=="resource" & testdata$treat=="C")

## 2.3 r-square plot
##r square values 
meltdata<-melt(alldata,id.vars = c("mean.ra","reactor","treat","bestmodel","modP.neutral","modP.resource","modP.combine"))
meltdata97<-melt(alldata97,id.vars = c("mean.ra","reactor","treat","bestmodel","modP.neutral","modP.resource","modP.combine"))

measurement<-sapply(strsplit(as.character(meltdata$variable),"[.]"), function(x) x[1])
models<-sapply(strsplit(as.character(meltdata$variable),"[.]"), function(x) x[2])
meltdata<-data.frame(meltdata,measurement,models)
rabins<-cut(meltdata$mean.ra,c(0, 0.0001, 0.001,1))
meltdata<-data.frame(meltdata,rabins)

measurement<-sapply(strsplit(as.character(meltdata97$variable),"[.]"), function(x) x[1])
models<-sapply(strsplit(as.character(meltdata97$variable),"[.]"), function(x) x[2])
meltdata97<-data.frame(meltdata97,measurement,models)
rabins<-cut(meltdata97$mean.ra,c(0, 0.0001, 0.001,1))
meltdata97<-data.frame(meltdata97,rabins)

meltdata$models<-factor(meltdata$models,levels = c("neutral","resource","combine"))
meltdata97$models<-factor(meltdata97$models,levels = c("neutral","resource","combine"))

mean(meltdata97$value[meltdata97$variable=="r2.resource" & meltdata97$rabins=="(0.001,1]" & meltdata97$treat=="D"])
mean(meltdata97$value[meltdata97$variable=="r2.resource" & meltdata97$rabins=="(0.001,1]" & meltdata97$treat=="C"])
sd(meltdata97$value[meltdata97$variable=="r2.resource" & meltdata97$rabins=="(0.001,1]" & meltdata97$treat=="D"])
sd(meltdata97$value[meltdata97$variable=="r2.resource" & meltdata97$rabins=="(0.001,1]" & meltdata97$treat=="C"])

mean(meltdata97$value[meltdata97$variable=="r2.combine" & meltdata97$rabins=="(0.001,1]" & meltdata97$treat=="D"])
mean(meltdata97$value[meltdata97$variable=="r2.combine" & meltdata97$rabins=="(0.001,1]" & meltdata97$treat=="C"])
sd(meltdata97$value[meltdata97$variable=="r2.combine" & meltdata97$rabins=="(0.001,1]" & meltdata97$treat=="D"])
sd(meltdata97$value[meltdata97$variable=="r2.combine" & meltdata97$rabins=="(0.001,1]" & meltdata97$treat=="C"])

t.test(x=meltdata97$value[meltdata97$variable=="r2.neutral" & meltdata97$rabins=="(0.001,1]" & meltdata97$treat=="D"],
       y=meltdata97$value[meltdata97$variable=="r2.neutral" & meltdata97$rabins=="(0.001,1]" & meltdata97$treat=="C"])

t.test(x=meltdata97$value[meltdata97$variable=="r2.resource" & meltdata97$rabins=="(0.001,1]" & meltdata97$treat=="D"],
       y=meltdata97$value[meltdata97$variable=="r2.resource" & meltdata97$rabins=="(0.001,1]" & meltdata97$treat=="C"])

t.test(x=meltdata97$value[meltdata97$variable=="r2.combine" & meltdata97$rabins=="(0.001,1]" & meltdata97$treat=="D"],
       y=meltdata97$value[meltdata97$variable=="r2.combine" & meltdata97$rabins=="(0.001,1]" & meltdata97$treat=="C"])

t.test(x=meltdata97$value[meltdata97$variable=="r2.neutral" & meltdata97$rabins=="(0.0001,0.001]" & meltdata97$treat=="D"],
       y=meltdata97$value[meltdata97$variable=="r2.neutral" & meltdata97$rabins=="(0.0001,0.001]" & meltdata97$treat=="C"])

t.test(x=meltdata97$value[meltdata97$variable=="r2.resource" & meltdata97$rabins=="(0.0001,0.001]" & meltdata97$treat=="D"],
       y=meltdata97$value[meltdata97$variable=="r2.resource" & meltdata97$rabins=="(0.0001,0.001]" & meltdata97$treat=="C"])

t.test(x=meltdata97$value[meltdata97$variable=="r2.combine" & meltdata97$rabins=="(0.0001,0.001]" & meltdata97$treat=="D"],
       y=meltdata97$value[meltdata97$variable=="r2.combine" & meltdata97$rabins=="(0.0001,0.001]" & meltdata97$treat=="C"])

t.test(x=meltdata97$value[meltdata97$variable=="r2.neutral" & meltdata97$rabins=="(0,0.0001]" & meltdata97$treat=="D"],
       y=meltdata97$value[meltdata97$variable=="r2.neutral" & meltdata97$rabins=="(0,0.0001]" & meltdata97$treat=="C"])

t.test(x=meltdata97$value[meltdata97$variable=="r2.resource" & meltdata97$rabins=="(0,0.0001]" & meltdata97$treat=="D"],
       y=meltdata97$value[meltdata97$variable=="r2.resource" & meltdata97$rabins=="(0,0.0001]" & meltdata97$treat=="C"])

t.test(x=meltdata97$value[meltdata97$variable=="r2.combine" & meltdata97$rabins=="(0,0.0001]" & meltdata97$treat=="D"],
       y=meltdata97$value[meltdata97$variable=="r2.combine" & meltdata97$rabins=="(0,0.0001]" & meltdata97$treat=="C"])

###using hist plot
ggplot(meltdata,aes(x =value,fill = models,linetype=treat,alpha=treat)) + 
  geom_density() +
  scale_fill_manual(values=c("#3C5488FF","#00A087FF","#E64B35FF"))+
  scale_color_manual(values=c("#3C5488FF","#00A087FF","#E64B35FF"))+
  scale_alpha_manual(values=c("C"=0.3,"D"=0.8))+
  scale_linetype_manual(values=c("C"="dashed","D"= "solid"))+
  xlim(-0.3,1.1)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())+
  facet_grid( models  ~  rabins)

ggplot(meltdata97,aes(x =value,fill = models,linetype=treat,alpha=treat)) + 
  geom_density() +
  scale_fill_manual(values=c("#3C5488FF","#00A087FF","#E64B35FF"))+
  scale_color_manual(values=c("#3C5488FF","#00A087FF","#E64B35FF"))+
  scale_alpha_manual(values=c("C"=0.3,"D"=0.8))+
  scale_linetype_manual(values=c("C"="dashed","D"= "solid"))+
  xlim(-0.5,1.2)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())+
  facet_grid( models  ~  rabins)

###using mean-error bar plot
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

meansd<-data_summary(meltdata97, varname="value", groupnames=c("models", "treat","rabins"))
meansd$rabins=factor(meansd$rabins,levels =c("(0.001,1]","(0.0001,0.001]","(0,0.0001]"))
ggplot(meansd, aes(x =value,y=treat,color = models,alpha=treat)) + 
  geom_point() +
  geom_errorbarh(aes(xmin=value-sd, xmax=value+sd), height=.2,)+
  scale_color_manual(values=c("#3C5488FF","#00A087FF","#E64B35FF"))+
  scale_alpha_manual(values=c("C"=0.3,"D"=1))+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())+
  facet_grid( models  ~  rabins)
  



## 2.4 Species relative abundance distribution
#compare with neutral model prediction (beta distribution)
#example on ESV1
##OTU1 in C1##
C1otu1<-comm[env$reactor=="C1","OTU_1"] #observed relative abundance

#NTm and pi predicted
id=which(row.names(C1)=="OTU_1")
ntm.all=(-C1$x1.neutral)/((C1$SE.neutral)^2)  
pi.all=C1$X.Intercept..neutral/(-C1$x1.neutral)
ntm=ntm.all[id]  #[1] 85.32622
pi=pi.all[id]    #[1] 0.0911959

x<-seq(0,1,length=1000)
y<-dbeta(x,shape1=7.781402,shape2=77.54482) ##predicted ra by beta-distribution
#shape1=ntm*pi; shape2=ntm*(1-pi)

predict<-data.frame(x=x,y=y)
observed<-data.frame(OTU_1=C1otu1,reactors=rep("C1",length(C1otu1)))

ggplot(data=observed,aes(x=OTU_1))+
  geom_histogram(aes(y = ..density..),alpha=0.8,binwidth =0.016,fill="#1B1919FF")+  
  geom_area(data=predict,aes(x=x,y=y),fill="#3C5488FF",alpha=0.5)+
  geom_line(data=predict,aes(x=x,y=y,ymax=y),colour="#3C5488FF")

ggplot(data=observed,aes(x=OTU_1))+
  geom_histogram(aes(y = ..density..),alpha=0.8,binwidth =0.016,fill="#1B1919FF")+  
  geom_area(data=predict,aes(x=x,y=y),fill="#3C5488FF",alpha=0.5)+
  geom_line(data=predict,aes(x=x,y=y,ymax=y),colour="#3C5488FF")+
  xlim(0,0.25)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



