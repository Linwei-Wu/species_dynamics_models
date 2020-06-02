

setwd("~/Dropbox/bioreactor/neutral modeling/new model/8_species simulation/combined model prediction/validation")


otutab<-read.csv("~/Dropbox/bioreactor/neutral modeling/new model/3_manuscript/20200115/R code/otutab.csv",row.names = 1,check.names = F)
env<-read.csv("~/Dropbox/bioreactor/neutral modeling/new model/3_manuscript/20200115/R code/env.csv",row.names = 1)

sum(row.names(env)==colnames(otutab))  #check

comm<-t(otutab)
comm<-comm/rowSums(comm)      ##get the relative abundance
comm[,"OTU_3"]
#0.006364657  0.008667508 0.005420811
###0. parameter estimation on one species, using data of two reactors 
source("/Users/linwei/Dropbox/bioreactor/neutral modeling/new model/3_manuscript/20200115/R code/main modeling function/fitting and comparing models.R")

oneSPtwoR<-function(D1otus,D1tpoint,D1Env,D2otus,D2tpoint,D2Env){
    D1x<-D1otus
    D1x0<-D1x[D1x>0];    #x1 initial value
    D1t0<-D1tpoint[D1x>0]
    D1x2<-D1x0*(1-D1x0)  #x2 initial value
    D1x3<-(D1Env$VS)*D1x*(1-D1x); D1x3<-D1x3[D1x>0] #x3 intitial value
    D1x4<-D1x0*D1x0*(1-D1x0)  #x4 initial value
    D1xvs<-(D1Env$VS)*D1x;D1x5<-D1xvs[D1x>0]  #x5 initial value
    D1x6<-D1x0^2   #x6 initial value
    D1dx<-D1x0[-1]-D1x0[-length(D1x0)];D1dt<-D1t0[-1]-D1t0[-length(D1x0)]
    D1y<-D1dx/D1dt  
    D1x1=D1x0[-length(D1x0)];   #x1
    D1x2=D1x2[-length(D1x0)];  #x2
    D1x3=D1x3[-length(D1x0)];   #x3
    D1x4=D1x4[-length(D1x0)];   #x4
    D1x5=D1x5[-length(D1x0)];   #x5
    D1x6=D1x6[-length(D1x0)];   #x6
    D1weight=D1dt/(2*D1x1*(1-D1x1))
    
    D2x<-D2otus
    D2x0<-D2x[D2x>0];    #x1 initial value
    D2t0<-D2tpoint[D2x>0]
    D2x2<-D2x0*(1-D2x0)  #x2 initial value
    D2x3<-(D2Env$VS)*D2x*(1-D2x); D2x3<-D2x3[D2x>0] #x3 intitial value
    D2x4<-D2x0*D2x0*(1-D2x0)  #x4 initial value
    D2xvs<-(D2Env$VS)*D2x;D2x5<-D2xvs[D2x>0]  #x5 initial value
    D2x6<-D2x0^2   #x6 initial value
    D2dx<-D2x0[-1]-D2x0[-length(D2x0)];D2dt<-D2t0[-1]-D2t0[-length(D2x0)]
    D2y<-D2dx/D2dt  
    D2x1=D2x0[-length(D2x0)];   #x1
    D2x2=D2x2[-length(D2x0)];  #x2
    D2x3=D2x3[-length(D2x0)];   #x3
    D2x4=D2x4[-length(D2x0)];   #x4
    D2x5=D2x5[-length(D2x0)];   #x5
    D2x6=D2x6[-length(D2x0)];   #x6
    D2weight=D2dt/(2*D2x1*(1-D2x1))
    
    y<-c(D1y,D2y)
    x1<-c(D1x1,D2x1)
    x2<-c(D1x2,D2x2)
    x3<-c(D1x3,D2x3)
    x4<-c(D1x4,D2x4)
    x5<-c(D1x5,D2x5)
    x6<-c(D1x6,D2x6)
    weight<-c(D1weight,D2weight)
    dt<-c(D1dt,D2dt)
    
    y[dt>20]<-NA   ##exclude y when time interval is more than 20 days
    if (sum(!is.na(y))<7 | length(unique(x1)) < 6){
      result<-rep(NA,42)
     } else{
      result<-comparemodel(y=y,x1=x1,x2=x2,x3=x3,x4=x4,x5=x5,x6=x6,Weight =weight )
      if (length(result)<42) result<-rep(NA,42)
    }
    result}

##0.a treatment reactors
D1otus<-comm[env$reactor=="D3","OTU_3"];D1Env<-env[env$reactor=="D3",]   
D1tpoint<-D1Env$time
D2otus<-comm[env$reactor=="D2","OTU_3"];D2Env<-env[env$reactor=="D2",]   
D2tpoint<-D2Env$time

pars<-oneSPtwoR(D1otus,D1tpoint,D1Env,D2otus,D2tpoint,D2Env)

##0.b control reactors
D1otus<-comm[env$reactor=="C1","OTU_3"];D1Env<-env[env$reactor=="C1",]   
D1tpoint<-D1Env$time
D2otus<-comm[env$reactor=="C3","OTU_3"];D2Env<-env[env$reactor=="C3",]   
D2tpoint<-D2Env$time

pars.C<-oneSPtwoR(D1otus,D1tpoint,D1Env,D2otus,D2tpoint,D2Env)


###1. pure neutral simulations. 
##ESV_3
ntm=(-pars["x1.neutral"])/((pars["SE.neutral"])^2)  
pi=-pars[1]/pars[2]
a=1/(pars["SE.neutral"])^2
ntm;pi;a
##Ntm=129, p=0.022, a=6113##
#dx=Ntm*p/b*dt-Ntm*x/b*dt+sqt(2x(1-x)/b)dwt
#one time step: 2h# 
##2 hours one step, thus dt=2/24=1/12 and dwt ~ N(0, dt)
##then 'sample' the time series at exactly the same time points as the bioreactr project (11 time points for each reactor)

###D1 prediction based on paramters estimated from D2 and D3
D1otus<-comm[env$reactor=="D1","OTU_3"];D1Env<-env[env$reactor=="D1",]
D1tpoint<-D1Env$time

t<-(D1tpoint-45)*12+1 ##this is the mapped t points in the simulated series

simuone<-function (j,Ntm,pi,a,x0){
  vec=numeric(1201)
  vec[1]=x0  ##initial relative abundance
  dwt<-rnorm(1200)/(12^0.5)
  for (i in 1:1200){
    dx=Ntm*pi/a/12-Ntm/a*vec[i]/12+((vec[i]*(1-vec[i])*2/a)^0.5)*dwt[i]
    vec[i+1]<-ifelse(vec[i]+dx>0,vec[i]+dx,0)
  }
  vec
}

D1simu.neutral<-sapply(1:1000,simuone,Ntm=129,pi=0.022,a=6113,x0=D1otus[1])
write.csv(D1simu.neutral,"D1prediction.neutral_ESV3.csv")


###2. pure resource simulations

##D1 simulation
D1Env<-env[env$reactor=="D1",]
D1tpoint<-D1Env$time

vssimu<-numeric(625)  #simulate VS, assuming a step-wise change in VS with time.
for (i in 1:10){
  vssimu[((D1tpoint[i]-45)*12+1):((D1tpoint[i+1]-45)*12+1)] <- D1Env$VS[i]
}

##simulate: 

##ESV_3
##dx=b1*x*dt+b2*x*VS*dt+b3*x^2*dt
##b1=-9.413135e-01; b2= 6.980018e-02; b3=-1.415508e+01

simuonce<-function(x0,b1,b2,b3,vssimu){
  vec=numeric(625)
  vec[1]=x0
  for (i in 1:624){
    dx=b1*vec[i]/12+b2*vec[i]*vssimu[i]/12+b3*vec[i]*vec[i]/12
    vec[i+1]<-ifelse(vec[i]+dx>0,vec[i]+dx,0)
  }
  vec
}


D1simu.niche<-sapply(1:1000,function(i){
  simuonce(x0=D1otus[1],b1=-9.413135e-01,b2=6.980018e-02,b3=-1.415508e+01,vssimu=vssimu)})

write.csv(D1simu.niche,"D1prediction.niche_ESV3.csv")


###3. combined model simulations
##D1 simulation
D1Env<-env[env$reactor=="D1",]
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
  vec
}

D1simu.combine<-sapply(1:1000,function(i){
  simuonce.combine(x0=D1otus[1],b0=-3.242007e-03,b1=2.801749e+03,b2=-2.801318e+03,b3=5.045712e-02,b4=-2.915794e+03,a=7133.569,vssimu=vssimu)
})
write.csv(D1simu.combine,"D1prediction.combine_ESV3.csv")


###compare predicted with observed
D1simu.neutral<-read.csv("D1prediction.neutral_ESV3.csv",row.names = 1)
D1simu.niche<-read.csv("D1prediction.niche_ESV3.csv",row.names = 1)
D1simu.combine<-read.csv("D1prediction.combine_ESV3.csv",row.names = 1)
#1. neutral model
# D1otus<-comm[env$reactor=="D1","OTU_3"]
time.d<-(1:nrow(D1simu.neutral)-1)/12+45
D1selected<-data.frame(Time=time.d,D1simu.neutral[,1:100])
D1selected<-D1selected[D1selected$Time<=97,]

means<-rowMeans(D1selected[,-1])
dat<-t(D1selected[,-1])
cf<-t(sapply(1:ncol(dat),function(i){
  quantile(dat[,i], probs = c(0.05, 0.5, 0.95))
}))

library(reshape2)
cfmean<-data.frame(Time=D1selected$Time,means=means,cf)  
head(cfmean)
cfmean<-cfmean[seq(1,625,by=12),]
cfmean<-melt(cfmean,id.vars ="Time")
write.csv(cfmean,"D1ESV3_neutral prediction_means and quantiles.csv")

D1otus<-comm[env$reactor=="D1","OTU_3"];D1Env<-env[env$reactor=="D1",]
D1tpoint<-D1Env$time
observed=data.frame(Time=D1tpoint,value=D1otus)


library(reshape2)
datamelt<-melt(D1selected,id.vars ="Time")
head(datamelt)
library(ggplot2)

ggplot()+
  geom_point(data=datamelt,aes(x=Time,y=value),color="grey74",alpha=0.1)+
  geom_line(data=cfmean[cfmean$variable=="X50.",],aes(x=Time,y=value),color="#3C5488FF")+
  geom_line(data=cfmean[cfmean$variable %in% c("X5.","X95."),],aes(x=Time,y=value,group=variable),linetype = "longdash",color="#3C5488FF")+
  geom_line(data=observed,aes(x=Time,y=value),linetype = "dotted")+
  geom_point(data=observed,aes(x=Time,y=value),color="orange")+
  theme_bw()+
  theme(legend.position = "none")+
  xlab("Time (days)")+
  ylab("Relative abunance")


#2. resource model
# D1otus<-comm[env$reactor=="D1","OTU_3"]
time.d<-(1:nrow(D1simu.niche)-1)/12+45
D1selected<-data.frame(Time=time.d,D1simu.niche[,1:100])
D1selected<-D1selected[D1selected$Time<=97,]

means<-rowMeans(D1selected[,-1])
dat<-t(D1selected[,-1])
cf<-t(sapply(1:ncol(dat),function(i){
  quantile(dat[,i], probs = c(0.05,0.5,  0.95))
}))

cfmean<-data.frame(Time=D1selected$Time,means=means,cf)  
head(cfmean)
cfmean<-cfmean[seq(1,625,by=12),]
cfmean<-melt(cfmean,id.vars ="Time")
write.csv(cfmean,"D1ESV3_niche prediction_means and quantiles.csv")


D1otus<-comm[env$reactor=="D1","OTU_3"];D1Env<-env[env$reactor=="D1",]
D1tpoint<-D1Env$time
observed=data.frame(Time=D1tpoint,value=D1otus)

datamelt<-melt(D1selected,id.vars ="Time")

ggplot()+
  #geom_point(data=datamelt,aes(x=Time,y=value),color="grey74",alpha=0.1)+
  geom_line(data=cfmean[cfmean$variable=="X50.",],aes(x=Time,y=value),color="#00A087FF")+
  #geom_line(data=cfmean[cfmean$variable!="means",],aes(x=Time,y=value,group=variable),linetype = "longdash",color="#00A087FF")+
  geom_line(data=observed,aes(x=Time,y=value),linetype = "dotted")+
  geom_point(data=observed,aes(x=Time,y=value),color="orange",size=2)+
  theme_bw()+
  theme(legend.position = "none")+
  xlab("Time (days)")+
  ylab("Relative abunance")+
  ylim(0,0.048)

#3. combined model
# D1otus<-comm[env$reactor=="D1","OTU_3"]
D1simu.combine[D1simu.combine==0]<-NA
time.d<-(1:nrow(D1simu.combine)-1)/12+45
D1selected<-data.frame(Time=time.d,D1simu.combine[,1:180])
D1selected<-D1selected[D1selected$Time<=97,]

means<-rowMeans(D1selected[,-1],na.rm = T)
dat<-t(D1selected[,-1])
cf<-t(sapply(1:ncol(dat),function(i){
  quantile(dat[,i], probs = c(0.05,0.5,  0.95),na.rm = T)
}))

cfmean<-data.frame(Time=D1selected$Time,means=means,cf)  
head(cfmean)
cfmean<-cfmean[seq(1,625,by=12),]
cfmean<-melt(cfmean,id.vars ="Time")
head(cfmean)
write.csv(cfmean,"D1ESV3_combined prediction_means and quantiles.csv")

D1otus<-comm[env$reactor=="D1","OTU_3"];D1Env<-env[env$reactor=="D1",]
D1tpoint<-D1Env$time
observed=data.frame(Time=D1tpoint,value=D1otus)
write.csv(observed,"D1ESV3_ra_observed.csv")

datamelt<-melt(D1selected,id.vars ="Time")

ggplot()+
  geom_point(data=datamelt,aes(x=Time,y=value),color="grey74",alpha=0.1)+
  geom_line(data=cfmean[cfmean$variable=="X50.",],aes(x=Time,y=value),color="#E64B35FF")+
  geom_line(data=cfmean[cfmean$variable %in% c("X5.","X95."),],aes(x=Time,y=value,group=variable),linetype = "longdash",color="#E64B35FF")+
  geom_line(data=observed,aes(x=Time,y=value),linetype = "dotted")+
  geom_point(data=observed,aes(x=Time,y=value),color="orange",size=2)+
  theme_bw()+
  theme(legend.position = "none")+
  xlab("Time (days)")+
  ylab("Relative abunance")+
  ylim(0,0.048)


######species prediction in control reactors
##Using C1 and C3 to predict C2
###1. pure neutral simulations. 
##ESV_3
ntm=(-pars.C["x1.neutral"])/((pars.C["SE.neutral"])^2)  
pi=-pars.C[1]/pars.C[2]
a=1/(pars.C["SE.neutral"])^2
ntm;pi;a
##Ntm=19.47077, p=0.04717401 , a=464.4287##
#dx=Ntm*p/b*dt-Ntm*x/b*dt+sqt(2x(1-x)/b)dwt
#one time step: 2h# 
##2 hours one step, thus dt=2/24=1/12 and dwt ~ N(0, dt)
##then 'sample' the time series at exactly the same time points as the bioreactr project (11 time points for each reactor)


###C2 prediction based on parameters estimated from C1 and C3
D1otus<-comm[env$reactor=="C2","OTU_3"];D1Env<-env[env$reactor=="C2",]
D1tpoint<-D1Env$time

steps<-(501-31)*12+1 ##this is the mapped t points in the simulated series
steps

simuone<-function (j,Ntm,pi,a,x0){
  vec=numeric(steps)
  vec[1]=x0  ##initial relative abundance
  dwt<-rnorm(steps)/(12^0.5)
  for (i in 1:steps){
    dx=Ntm*pi/a/12-Ntm/a*vec[i]/12+((vec[i]*(1-vec[i])*2/a)^0.5)*dwt[i]
    vec[i+1]<-ifelse(vec[i]+dx>0,vec[i]+dx,0)
  }
  vec
}

D1simu.neutral<-sapply(1:100,simuone,Ntm=19.47,pi=0.047,a=464.4,x0=D1otus[1])
write.csv(D1simu.neutral,"C2simu.neutral_ESV3.csv")


###2. pure resource simulations

##D1 simulation
D1otus<-comm[env$reactor=="C2","OTU_3"];D1Env<-env[env$reactor=="C2",]
D1tpoint<-D1Env$time

vssimu<-numeric(steps)  #simulate VS, assuming a step-wise change in VS with time.
for (i in 1:(length(D1tpoint)-1)){
  vssimu[((D1tpoint[i]-31)*12+1):((D1tpoint[i+1]-31)*12+1)] <- D1Env$VS[i]
}

##simulate: 

##ESV_3
##dx=b1*x*dt+b2*x*VS*dt+b3*x^2*dt
##b1=1.567850e-01; b2= -4.050253e-03 ; b3=-1.484848e+00

simuonce<-function(x0,b1,b2,b3,vssimu){
  vec=numeric(steps)
  vec[1]=x0
  for (i in 1:(steps-1)){
    dx=b1*vec[i]/12+b2*vec[i]*vssimu[i]/12+b3*vec[i]*vec[i]/12
    vec[i+1]<-ifelse(vec[i]+dx>0,vec[i]+dx,0)
  }
  vec
}

##b1=1.567850e-01; b2= -4.050253e-03 ; b3=-1.484848e+00
D1simu.niche<-sapply(1:100,function(i){
  simuonce(x0=D1otus[1],b1=1.567850e-01,b2=-4.050253e-03,b3=-1.484848e+00,vssimu=vssimu)})

write.csv(D1simu.niche,"C2simu.niche_ESV3.csv")


###3. combined model simulations
##D1 simulation
D1otus<-comm[env$reactor=="C2","OTU_3"];D1Env<-env[env$reactor=="C2",]
D1tpoint<-D1Env$time

vssimu<-numeric(steps)
for (i in 1:(length(D1tpoint)-1)){
  vssimu[((D1tpoint[i]-31)*12+1):((D1tpoint[i+1]-31)*12+1)] <- D1Env$VS[i]
}

##simulate: 

##ESV_3
##dx/dt=b0+b1*x+a1*(VS*x(1-x))+a2*x(1-x)+a3*(1-x)*(x^2)+error
##y=b0+b1*x1+b2*x2+b3*x3+b4*x4+error
##x1=x, x2=x(1-x), x3=vs*x(1-x),x4=x*x(1-x)
##b0=-0.00100405, b1=28.48534975; b2= -28.21675614; b3=0.002237907; b4=-35.24420597,timescale a=406

simuonce.combine<-function(x0,b0,b1,b2,b3,b4,a,vssimu){
  vec=numeric(steps)
  vec[1]=x0
  dwt<-rnorm(steps)/(12^0.5)
  for (i in 1:(steps-1)){
    dx=b0/12+b1*vec[i]/12+b2*vec[i]*(1-vec[i])/12+b3*vssimu[i]*vec[i]*(1-vec[i])/12+b4*vec[i]*vec[i]*(1-vec[i])/12+((vec[i]*(1-vec[i])*2/a)^0.5)*dwt[i]
    vec[i+1]<-ifelse(vec[i]+dx>0,vec[i]+dx,0)
  }
  vec
}

D1simu.combine<-sapply(1:300,function(i){
  simuonce.combine(x0=D1otus[1],b0=-1.313223e-03,b1=3.159177e+01,b2=-3.137698e+01 ,b3=6.100720e-03,b4=-3.873245e+01,a=486.7,vssimu=vssimu)
})
write.csv(D1simu.combine,"C2simu.combine_ESV3.csv")

###compare predicted with observed
D1simu.neutral<-read.csv("C1simu.neutral_ESV3.csv",row.names = 1)
D1simu.niche<-read.csv("C1simu.niche_ESV3.csv",row.names = 1)
D1simu.combine<-read.csv("C1simu.combine_ESV3.csv",row.names = 1)

#1. neutral model
# D1otus<-comm[env$reactor=="D1","OTU_3"]
time.d<-(1:nrow(D1simu.neutral)-1)/12+31
D1selected<-data.frame(Time=time.d,D1simu.neutral[,1:100])

means<-rowMeans(D1selected[,-1])
dat<-t(D1selected[,-1])
cf<-t(sapply(1:ncol(dat),function(i){
  quantile(dat[,i], probs = c(0.05, 0.5, 0.95))
}))

cfmean<-data.frame(Time=D1selected$Time,means=means,cf)  
head(cfmean)
cfmean<-cfmean[seq(1,5641,by=12),]
cfmean<-melt(cfmean,id.vars ="Time")
head(cfmean)
write.csv(cfmean,"C2ESV3_neutral prediction_means and quantiles.csv")


D1otus<-comm[env$reactor=="C2","OTU_3"];D1Env<-env[env$reactor=="C2",]
D1tpoint<-D1Env$time
observed=data.frame(Time=D1tpoint,value=D1otus)
write.csv(observed,"C2ESV3_ra_observed.csv")


D1selected<-D1selected[seq(1,5641,by=12),]
datamelt<-melt(D1selected,id.vars ="Time")

ggplot()+
  geom_point(data=datamelt,aes(x=Time,y=value),color="grey74",alpha=0.1)+
  geom_line(data=cfmean[cfmean$variable=="X50.",],aes(x=Time,y=value),color="#3C5488FF")+
  geom_line(data=cfmean[cfmean$variable %in% c("X5.","X95."),],aes(x=Time,y=value,group=variable),linetype = "longdash",color="#3C5488FF")+
  geom_line(data=observed,aes(x=Time,y=value),linetype = "dotted")+
  geom_point(data=observed,aes(x=Time,y=value),color="orange")+
  theme_bw()+
  theme(legend.position = "none")+
  xlab("Time (days)")+
  ylab("Relative abunance")


#2. resource model
# D1otus<-comm[env$reactor=="D1","OTU_3"]
time.d<-(1:nrow(D1simu.niche)-1)/12+31
D1selected<-data.frame(Time=time.d,D1simu.niche[,1:100])

means<-rowMeans(D1selected[,-1])
dat<-t(D1selected[,-1])
cf<-t(sapply(1:ncol(dat),function(i){
  quantile(dat[,i], probs = c(0.05,0.5,  0.95))
}))

cfmean<-data.frame(Time=D1selected$Time,means=means,cf)  
head(cfmean)
cfmean<-cfmean[seq(1,5641,by=12),]
cfmean<-melt(cfmean,id.vars ="Time")
write.csv(cfmean,"C2ESV3_niche prediction_means and quantiles.csv")


D1otus<-comm[env$reactor=="C2","OTU_3"];D1Env<-env[env$reactor=="C2",]
D1tpoint<-D1Env$time
observed=data.frame(Time=D1tpoint,value=D1otus)

D1selected<-D1selected[seq(1,5641,by=12),]
datamelt<-melt(D1selected,id.vars ="Time")

ggplot()+
  #geom_point(data=datamelt,aes(x=Time,y=value),color="grey74",alpha=0.1)+
  geom_line(data=cfmean[cfmean$variable=="X50.",],aes(x=Time,y=value),color="#00A087FF",size=1)+
  #geom_line(data=cfmean[cfmean$variable!="means",],aes(x=Time,y=value,group=variable),linetype = "longdash",color="#00A087FF")+
  geom_line(data=observed,aes(x=Time,y=value),linetype = "dotted")+
  geom_point(data=observed,aes(x=Time,y=value),color="orange",size=1)+
  theme_bw()+
  theme(legend.position = "none")+
  xlab("Time (days)")+
  ylab("Relative abunance")+
  ylim(0,0.40)

#3. combined model
# D1otus<-comm[env$reactor=="D1","OTU_3"]
#D1simu.combine[D1simu.combine==0]<-NA
time.d<-(1:nrow(D1simu.combine)-1)/12+31
D1selected<-data.frame(Time=time.d,D1simu.combine)

means<-rowMeans(D1selected[,-1],na.rm = T)
dat<-t(D1selected[,-1])
cf<-t(sapply(1:ncol(dat),function(i){
  quantile(dat[,i], probs = c(0.05,0.5,  0.95),na.rm = T)
}))

cfmean<-data.frame(Time=D1selected$Time,means=means,cf)  
head(cfmean)
cfmean<-cfmean[seq(1,5641,by=12),]
cfmean<-melt(cfmean,id.vars ="Time")
write.csv(cfmean,"C2ESV3_combine prediction_means and quantiles.csv")

D1otus<-comm[env$reactor=="C2","OTU_3"];D1Env<-env[env$reactor=="C2",]
D1tpoint<-D1Env$time
observed=data.frame(Time=D1tpoint,value=D1otus)

D1selected<-D1selected[seq(1,5641,by=12),]
datamelt<-melt(D1selected,id.vars ="Time")

ggplot()+
  geom_point(data=datamelt,aes(x=Time,y=value),color="grey74",alpha=0.1)+
  geom_line(data=cfmean[cfmean$variable=="X50.",],aes(x=Time,y=value),color="#E64B35FF",size=1)+
  geom_line(data=cfmean[cfmean$variable %in% c("X5.","X95."),],aes(x=Time,y=value,group=variable),linetype = "longdash",color="#E64B35FF")+
  geom_line(data=observed,aes(x=Time,y=value),linetype = "dotted")+
  geom_point(data=observed,aes(x=Time,y=value),color="orange",size=1)+
  theme_bw()+
  theme(legend.position = "none")+
  xlab("Time (days)")+
  ylab("Relative abunance")+
  ylim(0,0.40)


