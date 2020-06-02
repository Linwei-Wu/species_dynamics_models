
###1. read otu table and env table###
otutab<-read.csv("otutab.csv",row.names = 1,check.names = F)
env<-read.csv("env.csv",row.names = 1)
sum(row.names(env)==colnames(otutab))  #check
comm<-t(otutab)
comm<-comm/rowSums(comm)      ##get the relative abundance

###2. population and community stochasticity in control reactors
source("main modeling function/calculate community stochasticity.R")
#C1
otus<-comm[env$reactor=="C1",];Env<-env[env$reactor=="C1",]
tpoint<-Env$time
#whole community
cresult<-CSnormlz(otus=otus,Env=Env,tpoint = tpoint)
#different population
meanra=colMeans(otus)
rabins<-cut(meanra,c(0, 0.0001, 0.001,1))  ##classify to rare, moderate and abundant species 

presult<-sapply(split(data.frame(t(otus)),rabins),function(comms){
  CSnormlz(otus=t(comms),Env=Env,tpoint = tpoint)
})

Cresult.all<-data.frame(all=cresult,presult)
write.csv(Cresult.all,"C1_community stochasticity.csv")

#C2
otus<-comm[env$reactor=="C2",];Env<-env[env$reactor=="C2",]
tpoint<-Env$time
#whole community
cresult<-CSnormlz(otus=otus,Env=Env,tpoint = tpoint)
#different population
meanra=colMeans(otus)
rabins<-cut(meanra,c(0, 0.0001, 0.001,1))  ##classify to rare, moderate and abundant species 

presult<-sapply(split(data.frame(t(otus)),rabins),function(comms){
  CSnormlz(otus=t(comms),Env=Env,tpoint = tpoint)
})

Cresult.all<-data.frame(all=cresult,presult)
write.csv(Cresult.all,"C2_community stochasticity.csv")

#C3
otus<-comm[env$reactor=="C3",];Env<-env[env$reactor=="C3",]
tpoint<-Env$time
#whole community
cresult<-CSnormlz(otus=otus,Env=Env,tpoint = tpoint)
#different population
meanra=colMeans(otus)
rabins<-cut(meanra,c(0, 0.0001, 0.001,1))  ##classify to rare, moderate and abundant species 

presult<-sapply(split(data.frame(t(otus)),rabins),function(comms){
  CSnormlz(otus=t(comms),Env=Env,tpoint = tpoint)
})

Cresult.all<-data.frame(all=cresult,presult)
write.csv(Cresult.all,"C3_community stochasticity.csv")


###treatment reactors, combine the three time series together (more time points)
source("application on bioreactor data/calculate community stochasticity_combine reps.R")

D1otus<-comm[env$reactor=="D1",];D1Env<-env[env$reactor=="D1",]
D1tpoint<-D1Env$time
D2otus<-comm[env$reactor=="D2",];D2Env<-env[env$reactor=="D2",]
D2tpoint<-D2Env$time
D3otus<-comm[env$reactor=="D3",];D3Env<-env[env$reactor=="D3",]
D3tpoint<-D3Env$time

#whole community
dresult<-CSnormlz.reps(D1otus,D1tpoint,D1Env,D2otus,D2tpoint,D2Env,D3otus,D3tpoint,D3Env)

#different population
meanra=colMeans(rbind(D1otus,D2otus,D3otus))
rabins<-cut(meanra,c(0, 0.0001, 0.001,1))  ##classify to rare, moderate and abundant species 
names(rabins)<-names(meanra)

group=levels(rabins)
Comm<-comm[,!is.na(rabins)]
rabins=rabins[!is.na(rabins)]

dpresult<-sapply(1:length(group),function(i){
  message("i= ",i," in ",length(group))
  D1otus<-Comm[env$reactor=="D1",rabins==group[i]]
  D2otus<-Comm[env$reactor=="D2",rabins==group[i]]
  D3otus<-Comm[env$reactor=="D3",rabins==group[i]]
  
  D1Env<-env[env$reactor=="D1",]
  D2Env<-env[env$reactor=="D2",]
  D3Env<-env[env$reactor=="D3",]
  D1tpoint<-D1Env$time
  D2tpoint<-D2Env$time
  D3tpoint<-D3Env$time
  
  CSnormlz.reps(D1otus,D1tpoint,D1Env,D2otus,D2tpoint,D2Env,D3otus,D3tpoint,D3Env)

  })

colnames(dpresult)<-group
Dresult.all<-data.frame(all=dresult,dpresult)
write.csv(Dresult.all,"D_community stochasticity.csv")


###plot
C1cs<-read.csv("C1_community stochasticity.csv",row.names = 1)
C2cs<-read.csv("C2_community stochasticity.csv",row.names = 1)
C3cs<-read.csv("C3_community stochasticity.csv",row.names = 1)
Dcs<-read.csv("D_community stochasticity.csv",row.names = 1)
Ccs<-(C1cs+C2cs+C3cs)/3

colnames(Dcs)<-colnames(Ccs)<-c("all","rare","moderate","abundant")
library(reshape2)
Ccs.used<-data.frame(reactor=rep("C",4),process=rep(c("stochasticity","determinism"),each=2),
                     method=rep(c("unweighted","weighted"),2),rbind(Ccs[3:4,],1-Ccs[3:4,]))
Dcs.used<-data.frame(reactor=rep("D",4),process=rep(c("stochasticity","determinism"),each=2),
                     method=rep(c("unweighted","weighted"),2),rbind(Dcs[3:4,],1-Dcs[3:4,]))
Ccs.used<-melt(Ccs.used, id.vars = c("reactor","process","method"))
Dcs.used<-melt(Dcs.used, id.vars = c("reactor","process","method"))

df<-rbind(Ccs.used,Dcs.used)

library(ggplot2)
df$variable
df$variable=factor(df$variable,levels = rev(unique(df$variable)))

dat<-df[df$process=="stochasticity" & df$method=="weighted",]
ggplot(dat, aes(x=variable, y=value*100, fill=reactor)) +
  geom_bar(stat="identity", position=position_dodge(),width=0.7)+
  scale_fill_manual(values=c("#999999","#E69F00")) +
  ylab("Weighted community stochasticity (%)")+
  theme_bw()

dat<-df[df$process=="stochasticity" & df$method=="unweighted",]
ggplot(dat, aes(x=variable, y=value*100, fill=reactor)) +
  geom_bar(stat="identity", position=position_dodge(),width=0.7)+
  scale_fill_manual(values=c("#999999","#E69F00")) +
  ylab("Unweighted community stochasticity (%)")+
  theme_bw()
