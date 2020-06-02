
setwd("~/Dropbox/bioreactor/neutral modeling/new model/3_manuscript/20200115/R code")


otutab<-read.csv("otutab.csv",row.names = 1,check.names = F)
env<-read.csv("env.csv",row.names = 1)

sum(row.names(env)==colnames(otutab))  #check

comm<-t(otutab)
comm<-comm/rowSums(comm)      ##get the relative abundance


###1. pure neutral simulations. 
##ESV_3
##Ntm=140, p=0.02, a=7400, steps=1200##
#dx=Ntm*p/b*dt-Ntm*x/b*dt+sqt(2x(1-x)/b)dwt
#x0=0.022# 
#one time step: 2h# 
##2 hours one step, thus dt=2/24=1/12 and dwt ~ N(0, dt)
##then 'sample' the time series at exactly the same time points as the bioreactr project (11 time points for each reactor)
###D1 simulation
D1otus<-comm[env$reactor=="D1",];D1Env<-env[env$reactor=="D1",]
D1tpoint<-D1Env$time

t<-D1tpoint*12    

simuone<-function (j,Ntm,pi,a,x0){
  vec=numeric(1201)
  vec[1]=x0  ##initial relative abundance
  dwt<-rnorm(1200)/(12^0.5)
  for (i in 1:1200){
    dx=Ntm*pi/a/12-Ntm/a*vec[i]/12+((vec[i]*(1-vec[i])*2/a)^0.5)*dwt[i]
    vec[i+1]<-ifelse(vec[i]+dx>0,vec[i]+dx,0)
  }
  simu1<-c()
  for (i in 1:length(t)){
    simu1[i]<-vec[t[i]]
  }
  simu1
}

D1simu.neutral<-sapply(1:1000,simuone,Ntm=140,pi=0.02,a=7400,x0=0.021)
write.csv(D1simu.neutral,"fitting on simulated data/D1simu.neutral.csv")
D2simu.neutral<-sapply(1:1000,simuone,Ntm=140,pi=0.02,a=7400,x0=0.019)
write.csv(D2simu.neutral,"fitting on simulated data/D2simu.neutral.csv")
D3simu.neutral<-sapply(1:1000,simuone,Ntm=140,pi=0.02,a=7400,x0=0.020)
write.csv(D3simu.neutral,"fitting on simulated data/D3simu.neutral.csv")


###2. pure resource simulations

##D1 simulation
D1otus<-comm[env$reactor=="D1",];D1Env<-env[env$reactor=="D1",]
D1tpoint<-D1Env$time

steps<-(97-45)*12+1 

vssimu<-numeric(steps)  #simulate VS, assuming a step-wise change in VS with time.
for (i in 1:(length(D1tpoint)-1)){
  vssimu[((D1tpoint[i]-45)*12+1):((D1tpoint[i+1]-45)*12+1)] <- D1Env$VS[i]
}

##simulate: 

##ESV_3
##dx=b1*x*dt+b2*x*VS*dt+b3*x^2*dt
##b1=-0.7; b2= 0.06; b3=-15.0


simuonce<-function(x0,b1,b2,b3,vssimu){
  vec=numeric(steps)
  vec[1]=x0
  for (i in 1:(steps-1)){
    dx=b1*vec[i]/12+b2*vec[i]*vssimu[i]/12+b3*vec[i]*vec[i]/12
    vec[i+1]<-ifelse(vec[i]+dx>0,vec[i]+dx,0)
  }
  simu1<-c()
  for (i in 1:length(D1tpoint)){
    simu1[i]<-vec[((D1tpoint[i]-45)*12+1)]
  }
  simu1
}


D1simu.niche<-sapply(1:1000,function(i){
  simuonce(x0=0.009,b1=-0.7,b2=0.06,b3=-15.0,vssimu=vssimu)})

write.csv(D1simu.niche,"fitting on simulated data/D1simu.niche.csv")

##D2 simulation
D1otus<-comm[env$reactor=="D2",];D1Env<-env[env$reactor=="D2",]
D1tpoint<-D1Env$time


vssimu<-numeric(steps)  #simulate VS, assuming a step-wise change in VS with time.
for (i in 1:(length(D1tpoint)-1)){
  vssimu[((D1tpoint[i]-45)*12+1):((D1tpoint[i+1]-45)*12+1)] <- D1Env$VS[i]
}

##simulate: 
D2simu.niche<-sapply(1:1000,function(i){
  simuonce(x0=0.008,b1=-0.7,b2=0.06,b3=-15.0,vssimu=vssimu)})

write.csv(D2simu.niche,"fitting on simulated data/D2simu.niche.csv")

##D3 simulation
D1otus<-comm[env$reactor=="D3",];D1Env<-env[env$reactor=="D3",]
D1tpoint<-D1Env$time

vssimu<-numeric(steps)  #simulate VS, assuming a step-wise change in VS with time.
for (i in 1:(length(D1tpoint)-1)){
  vssimu[((D1tpoint[i]-45)*12+1):((D1tpoint[i+1]-45)*12+1)] <- D1Env$VS[i]
}
##simulate: 
D3simu.niche<-sapply(1:1000,function(i){
  simuonce(x0=0.01,b1=-0.7,b2=0.06,b3=-15.0,vssimu=vssimu)})

write.csv(D3simu.niche,"fitting on simulated data/D3simu.niche.csv")

###3. combined model simulations
##D3 simulation
D1otus<-comm[env$reactor=="D1",];D1Env<-env[env$reactor=="D1",]
D1tpoint<-D1Env$time

vssimu<-numeric(625)
for (i in 1:10){
  vssimu[((D1tpoint[i]-45)*12+1):((D1tpoint[i+1]-45)*12+1)] <- D1Env$VS[i]
}

##simulate: 

##ESV_3
##dx/dt=b0+b1*x+a1*(VS*x(1-x))+a2*x(1-x)+a3*(1-x)*(x^2)+error
##y=b0+b1*x1+b2*x2+b3*x3+b4*x4+error
##x1=x, x2=x(1-x), x3=vs*x(1-x),x4=x*x(1-x)
##b0=-0.0026, b1=2115.2; b2= -2114.93; b3=0.046; b4=-2205.5,timescale a=9273

simuonce.combine<-function(x0,b0,b1,b2,b3,b4,a,vssimu){
  vec=numeric(625)
  vec[1]=x0
  dwt<-rnorm(625)/(12^0.5)
  for (i in 1:624){
    dx=b0/12+b1*vec[i]/12+b2*vec[i]*(1-vec[i])/12+b3*vssimu[i]*vec[i]*(1-vec[i])/12+b4*vec[i]*vec[i]*(1-vec[i])/12+((vec[i]*(1-vec[i])*2/a)^0.5)*dwt[i]
    vec[i+1]<-ifelse(vec[i]+dx>0,vec[i]+dx,0)
  }
  simu1<-c()
  for (i in 1:11){
    simu1[i]<-vec[((D1tpoint[i]-45)*12+1)]
  }
  simu1
}

D1simu.combine<-sapply(1:1500,function(i){
  simuonce.combine(x0=0.005420811,b0=-0.002595874,b1=2115.200197,b2=-2114.933163,b3=0.045784234,b4=-2205.505617,a=9272.518406,vssimu=vssimu)
})
write.csv(D1simu.combine,"fitting on simulated data/D1simu.combine.csv")


##D2 simulation
D1otus<-comm[env$reactor=="D2",];D1Env<-env[env$reactor=="D2",]
D1tpoint<-D1Env$time

vssimu<-numeric(625)
for (i in 1:10){
  vssimu[((D1tpoint[i]-45)*12+1):((D1tpoint[i+1]-45)*12+1)] <- D1Env$VS[i]
}
D2simu.combine<-sapply(1:1500,function(i){
  simuonce.combine(x0=0.008667508,b0=-0.002595874,b1=2115.200197,b2=-2114.933163,b3=0.045784234,b4=-2205.505617,a=9272.518406,vssimu=vssimu)
})
write.csv(D2simu.combine,"fitting on simulated data/D2simu.combine.csv")

##D1 simulation
D1otus<-comm[env$reactor=="D3",];D1Env<-env[env$reactor=="D3",]
D1tpoint<-D1Env$time

vssimu<-numeric(625)
for (i in 1:10){
  vssimu[((D1tpoint[i]-45)*12+1):((D1tpoint[i+1]-45)*12+1)] <- D1Env$VS[i]
}

D3simu.combine<-sapply(1:1500,function(i){
  simuonce.combine(x0=0.0064,b0=-0.002595874,b1=2115.200197,b2=-2114.933163,b3=0.045784234,b4=-2205.505617,a=9272.518406,vssimu=vssimu)
})
write.csv(D3simu.combine,"fitting on simulated data/D3simu.combine.csv")


##############below is the parameter estimation and model fitting on
#############the simulated data
###4. parameter estimation and model fitting
 
#stochasticity=1   
sto=1  #change to 0, 0.1, 0.2,.., 0.9, 1 to test different species stochasticity level
D1=D1simu.neutral*sto+D1simu.niche*(1-sto)
#D1=D1simu.combine
#Noise=0
Noise<-sapply(1:1000,function(i){
  meanra=mean(D1[,i])
  noise<-rnorm(n=11,mean=0,sd=0.4*meanra)   #change sd=0.1*meanra, 0.2*meanra to indicate different noise level
  noise
})
D1<-D1+Noise

D2=D2simu.neutral*sto+D2simu.niche*(1-sto)
#D2=D2simu.combine
Noise<-sapply(1:1000,function(i){
  meanra=mean(D2[,i])
  noise<-rnorm(n=11,mean=0,sd=0.4*meanra)
  noise
})
D2<-D2+Noise

D3=D3simu.neutral*sto+D3simu.niche*(1-sto)
#D3=D3simu.combine
Noise<-sapply(1:1000,function(i){
  meanra=mean(D3[,i])
  noise<-rnorm(n=11,mean=0,sd=0.4*meanra)
  noise
})
D3<-D3+Noise


#prepare data to be combined

D1otus<-D1;D1Env<-env[env$reactor=="D1",]
D1tpoint<-D1Env$time
D2otus<-D2;D2Env<-env[env$reactor=="D2",]
D2tpoint<-D2Env$time
D3otus<-D3;D3Env<-env[env$reactor=="D3",]
D3tpoint<-D3Env$time

source("main modeling function/fitting and comparing models.R")
allresult<-combine.t(D1otus,D1tpoint,D1Env,D2otus,D2tpoint,D2Env,D3otus,D3tpoint,D3Env)
#write.csv(allresult,"fitting on simulated data/fitting result_stochastic0%_0%noise.csv")

##4.1 AIC-based best model distribution
C1<-data.frame(allresult)
aics<-data.frame(neutral=C1$aic.neutral,resource=C1$aic.resource,combine=C1$aic.combine)
models<-c("neutral","consumer-resource","combined")
C1model<-sapply(1:nrow(C1), function(i){
  result=ifelse(is.na(aics[i,1]), NA, models[which.min(aics[i,])])
  
})
table(C1model)

tableC<-data.frame(table(C1model)/1000)
# tableCD<-rbind(tableC, data.frame(C1model=c("neutral"),Freq=c(0)))
# tableCD<-rbind(tableC, data.frame(C1model=c("neutral","combined"),Freq=c(0,0)))
tableCD<-tableC
tableCD$C1model=factor(tableCD$C1model,levels=c("neutral","consumer-resource","combined"))
library(ggplot2)
ggplot(tableCD,aes(x = C1model, y = Freq,fill = C1model)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("#3C5488FF","#00A087FF","#E64B35FF"))+
  theme_bw()+
  #ylab("Probability of being the best model")+
  ylab(NULL)+
  xlab(NULL)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  ggtitle("stochasticity 0, 0% noise")



##4.2. plot for parameter estimation
#examples of neutral parameter
dat<-read.csv("fitting on simulated data/fitting result_stochastic100%_10%noise.csv",row.names = 1)
dat<-read.csv("/Users/linwei/Dropbox/bioreactor/neutral modeling/new model/2_stochasticity simulation/fitting on simulated data/simulation for parameter estimation/fitting result_stochastic100%_0%noise.csv",row.names = 1)


psimu<-dat$X.Intercept..neutral/(-dat$x1.neutral)
hist(psimu[psimu>0 ],breaks=40,main="noise0_p",xlim=c(0.01,0.04))
abline(v=c(0.02),col="red",lwd=2)

ntm<-(-dat$x1.neutral)/(dat$SE.neutral)^2
hist(ntm[ntm>0],breaks=30,main="noise0_ntm")
abline(v=c(140),col="red",lwd=2)


##example of niche paramters bici (competition strength)
#bici=0.06 for ESV3 that were simulated

dat2<-read.csv("fitting on simulated data/fitting result_stochastic0%_0%noise.csv",row.names = 1)

bc<- dat2$x5.resource
hist(bc,main="noise0_bici",xlim=c(0.045,0.065))
abline(v=c(0.06),col="red",lwd=2)


##example of combine model parameter bici
#bici=0.046 for ESV3 that were simulated

dat3<-read.csv("fitting on simulated data/fitting result_combined_0%noise.csv",row.names = 1)

bc<- dat3$x3.combine
hist(bc,breaks=40,main="noise0_bici",xlim=c(-0.1,0.3))
abline(v=c(0.046),col="red",lwd=2)

