


###1. generate the initial community###
otutab<-read.csv("otutab.csv",row.names = 1,check.names = F)
env<-read.csv("env.csv",row.names = 1)
sum(row.names(env)==colnames(otutab))  #check
comm<-t(otutab)
comm<-comm/rowSums(comm)      ##get the relative abundance

#D1 initial
D1otus<-comm[env$reactor=="D1",];D1Env<-env[env$reactor=="D1",]
D1tpoint<-D1Env$time
D1intl<-colMeans(D1otus) #initial value of relative abundance

#D2 initial
D2otus<-comm[env$reactor=="D2",];D2Env<-env[env$reactor=="D2",]
D2tpoint<-D2Env$time
D2intl<-colMeans(D2otus)  #initial value

#D3 initial
D3otus<-comm[env$reactor=="D3",];D3Env<-env[env$reactor=="D3",]
D3tpoint<-D3Env$time
D3intl<-colMeans(D3otus)  #initial value

meanra<-rowMeans(cbind(D1intl,D2intl,D2intl))  ##the initial relative abundace for each species

###2. simulate niche relative abundance 
#get niche parameters for each species

fitting.result<-read.csv("~/Dropbox/bioreactor/neutral modeling/new model/1_model fitting/D_bf100 ra fitting_r2.csv",row.names = 1)
OTUid<-row.names(fitting.result)[!is.na(fitting.result$x1.neutral)]  ##OTU id used

fitting.result<-fitting.result[OTUid,]
meanra<-meanra[OTUid]
sum(names(meanra)==row.names(fitting.result))  #check

niche.p<-data.frame(meanra=meanra,fitting.result[,c("x1.resource","x5.resource","x6.resource")])
#niche.p is the niche parameter used to generate the ra simulations

#D1 niche simulation
vssimu<-numeric(625)  ##simulate VS (volatile solids), 2 hours per step
for (i in 1:10){
  vssimu[((D1tpoint[i]-45)*12+1):((D1tpoint[i+1]-45)*12+1)] <- D1Env$VS[i]
}

simuonce<-function(x0,b1,b2,b3,vssimu){
  vec=numeric(625)
  vec[1]=x0
  for (i in 1:624){
    dx=b1*vec[i]/12+b2*vec[i]*vssimu[i]/12+b3*vec[i]*vec[i]/12
    vec[i+1]<-ifelse(vec[i]+dx>0,vec[i]+dx,0)
  }
  simu1<-c()
  for (i in 1:11){
    simu1[i]<-vec[((D1tpoint[i]-45)*12+1)]
  }
  simu1
}

D1simu.niche<-sapply(1:nrow(niche.p),function(i){
  simuonce(x0=niche.p[i,1],b1=niche.p[i,2],b2=niche.p[i,3],b3=niche.p[i,4],vssimu=vssimu)})
colnames(D1simu.niche)<-row.names(niche.p)
D1simu.niche[is.na(D1simu.niche)]<-0
#write.csv(D1simu.niche,"D1simu.niche.csv")

#D2 niche simulation
vssimu<-numeric(625)
for (i in 1:10){
  vssimu[((D2tpoint[i]-45)*12+1):((D2tpoint[i+1]-45)*12+1)] <- D2Env$VS[i]
}

D2simu.niche<-sapply(1:nrow(niche.p),function(i){
  simuonce(x0=niche.p[i,1],b1=niche.p[i,2],b2=niche.p[i,3],b3=niche.p[i,4],vssimu=vssimu)})
colnames(D2simu.niche)<-row.names(niche.p)
D2simu.niche[is.na(D2simu.niche)]<-0
#write.csv(D2simu.niche,"D2simu.niche.csv")

#D3 niche simulation
vssimu<-numeric(625)
for (i in 1:10){
  vssimu[((D3tpoint[i]-45)*12+1):((D3tpoint[i+1]-45)*12+1)] <- D3Env$VS[i]
}

D3simu.niche<-sapply(1:nrow(niche.p),function(i){
  simuonce(x0=niche.p[i,1],b1=niche.p[i,2],b2=niche.p[i,3],b3=niche.p[i,4],vssimu=vssimu)})
colnames(D3simu.niche)<-row.names(niche.p)
D3simu.niche[is.na(D3simu.niche)]<-0
#write.csv(D3simu.niche,"D3simu.niche.csv")

###3. neutral simulation##
neutral.p<-data.frame(meanra=meanra,Ntm=(-fitting.result$x1.neutral)/(fitting.result$SE.neutral)^2,
                      pi=fitting.result$X.Intercept..neutral/(-fitting.result$x1.neutral),a=1/(fitting.result$SE.neutral)^2)
#neutral.p is the neutral parameters

#D1 simulation
t<-D1tpoint*12  ##2-h per step
simuone<-function (Ntm,pi,a,x0){
  vec=numeric(1201)
  vec[1]=x0  ##initial relative abundance
  dwt<-rnorm(1200)/(12^0.5)
  for (i in 1:1200){
    dx=Ntm*pi/a/12-Ntm/a*vec[i]/12+((vec[i]*(1-vec[i])*2/a)^0.5)*dwt[i]
    vec[i+1]<-ifelse(vec[i]+dx>0,vec[i]+dx,0)
  }
  simu1<-c()
  for (i in 1:11){
    simu1[i]<-vec[t[i]]
  }
  simu1
}

D1simu.neutral<-sapply(1:nrow(neutral.p),function(i){
  simuone(Ntm=neutral.p[i,2],pi=neutral.p[i,3],a=neutral.p[i,4],x0=neutral.p[i,1])
})
colnames(D1simu.neutral)<-row.names(neutral.p)
D1simu.neutral[is.na(D1simu.neutral)]<-0
#write.csv(D1simu.neutral,"D1simu.neutral.csv")

#D2 simulation
D2simu.neutral<-sapply(1:nrow(neutral.p),function(i){
  simuone(Ntm=neutral.p[i,2],pi=neutral.p[i,3],a=neutral.p[i,4],x0=neutral.p[i,1])
})
colnames(D2simu.neutral)<-row.names(neutral.p)
D2simu.neutral[is.na(D2simu.neutral)]<-0
#write.csv(D2simu.neutral,"D2simu.neutral.csv")

#D3 simulation
D3simu.neutral<-sapply(1:nrow(neutral.p),function(i){
  simuone(Ntm=neutral.p[i,2],pi=neutral.p[i,3],a=neutral.p[i,4],x0=neutral.p[i,1])
})
colnames(D3simu.neutral)<-row.names(neutral.p)
D3simu.neutral[is.na(D3simu.neutral)]<-0
#write.csv(D3simu.neutral,"D3simu.neutral.csv")


###4. generate community RA under a certain community stochasticity (CS) 
### and fit the models
#add 0, 10%, 20%, 30% and 40% noise
#permute 1000 times

#CS is the community stochasticity level
source("main modeling function/fitting and comparing models.R")
##4.1 test unweighted CS
CSsimu.one<-function(k,CS){
  message("k=",k)
  sspecies.n<-round(ncol(D1simu.neutral)*CS)
  sspecies.id<-sample(1:ncol(D1simu.neutral), sspecies.n)
  
  D1<-matrix(0,nrow=nrow(D1simu.neutral),ncol=ncol(D1simu.neutral))
  D2<-matrix(0,nrow=nrow(D2simu.neutral),ncol=ncol(D2simu.neutral))
  D3<-matrix(0,nrow=nrow(D3simu.neutral),ncol=ncol(D3simu.neutral))
  
  for(i in 1:ncol(D1simu.neutral)){
    sto<-ifelse((i %in% sspecies.id),runif(1,min = 0.55,max=1),runif(1,min = 0,max=0.45))
    D1[,i]=D1simu.neutral[,i]*sto+D1simu.niche[,i]*(1-sto)
    D2[,i]=D2simu.neutral[,i]*sto+D2simu.niche[,i]*(1-sto)
    D3[,i]=D3simu.neutral[,i]*sto+D3simu.niche[,i]*(1-sto)
  }
  
  predicted.CS<-sapply(c(0,0.1,0.2,0.3,0.4),function(noise){
    message("noise=",noise)
    if (noise==0){
      Noise=0
    } else {
      Noise<-sapply(1:ncol(D1), function(i){
        meanra=mean(D1[,i])
        Ns<-rnorm(n=11,mean=0,sd=noise*meanra)
        Ns
      })}
    D1noise=D1+Noise
    D1noise[D1noise<0]<-0
    D1noise[D1noise>1]<-0
    
    if (noise==0){
      Noise=0
    } else {
      Noise<-sapply(1:ncol(D2), function(i){
        meanra=mean(D2[,i])
        Ns<-rnorm(n=11,mean=0,sd=noise*meanra)
        Ns
      })}
    D2noise=D2+Noise
    D2noise[D2noise<0]<-0
    D2noise[D2noise>1]<-0
    
    if (noise==0){
      Noise=0
    } else {
      Noise<-sapply(1:ncol(D3), function(i){
        meanra=mean(D3[,i])
        Ns<-rnorm(n=11,mean=0,sd=noise*meanra)
        Ns
      })}
    D3noise=D3+Noise
    D3noise[D3noise<0]<-0
    D3noise[D3noise>1]<-0
    
    D1otus<-D1noise;D1Env<-env[env$reactor=="D1",]
    D1tpoint<-D1Env$time
    D2otus<-D2noise;D2Env<-env[env$reactor=="D2",]
    D2tpoint<-D2Env$time
    D3otus<-D3noise;D3Env<-env[env$reactor=="D3",]
    D3tpoint<-D3Env$time
    
    allresult<-data.frame(combine.t(D1otus,D1tpoint,D1Env,D2otus,D2tpoint,D2Env,D3otus,D3tpoint,D3Env))
    aics<-data.frame(neutral=allresult$aic.neutral,resource=allresult$aic.resource,combine=allresult$aic.combine)
    models<-c("neutral","consumer-resource","combined")
    allmodel<-sapply(1:nrow(allresult), function(i){
      result=ifelse(is.na(aics[i,1]), NA, models[which.min(aics[i,])])
    })
    testtable<-table(allmodel)
    cs=testtable[3]/sum(testtable)
    cs
  })
  names(predicted.CS)<-c("no noise","noise10%","noise20%","noise30%","noise40%")
  predicted.CS
}

CS0<-sapply(1:1000,CSsimu.one,CS=0);write.csv(CS0,"CS0_update.csv")
CS0.1<-sapply(1:1000,CSsimu.one,CS=0.1);write.csv(CS0.1,"CS0.1_update.csv")
CS0.2<-sapply(1:1000,CSsimu.one,CS=0.2);write.csv(CS0.2,"CS0.2_update.csv")
CS0.3<-sapply(1:1000,CSsimu.one,CS=0.3);write.csv(CS0.3,"CS0.3_update.csv")
CS0.4<-sapply(1:1000,CSsimu.one,CS=0.4);write.csv(CS0.4,"CS0.4_update.csv")
CS0.5<-sapply(1:1000,CSsimu.one,CS=0.5);write.csv(CS0.5,"CS0.5_update.csv")
CS0.6<-sapply(1:1000,CSsimu.one,CS=0.6);write.csv(CS0.6,"CS0.6_update.csv")
CS0.7<-sapply(1:1000,CSsimu.one,CS=0.7);write.csv(CS0.7,"CS0.7_update.csv")
CS0.8<-sapply(1:1000,CSsimu.one,CS=0.8);write.csv(CS0.8,"CS0.8_update.csv")
CS0.9<-sapply(1:1000,CSsimu.one,CS=0.9);write.csv(CS0.9,"CS0.9_update.csv")
CS1<-sapply(1:1000,CSsimu.one,CS=1);write.csv(CS1,"CS1_update.csv")

##4.2 test weighted CS (weighted by species abundance)
ra0<-meanra;sum(names(ra0)==colnames(D1simu.niche))
CSsimuw.one<-function(j,CS){
  message("j=",j)
  sspecies.n<-round(ncol(D1simu.neutral)*CS)
  sspecies.id<-sample(1:ncol(D1simu.neutral), sspecies.n)
  truecs=sum(ra0[sspecies.id])/sum(ra0)
  
  D1<-matrix(0,nrow=nrow(D1simu.neutral),ncol=ncol(D1simu.neutral))
  D2<-matrix(0,nrow=nrow(D2simu.neutral),ncol=ncol(D2simu.neutral))
  D3<-matrix(0,nrow=nrow(D3simu.neutral),ncol=ncol(D3simu.neutral))
  
  for(i in 1:ncol(D1simu.neutral)){
    sto<-ifelse((i %in% sspecies.id),runif(1,min = 0.55,max=1),runif(1,min = 0,max=0.45))
    D1[,i]=D1simu.neutral[,i]*sto+D1simu.niche[,i]*(1-sto)
    D2[,i]=D2simu.neutral[,i]*sto+D2simu.niche[,i]*(1-sto)
    D3[,i]=D3simu.neutral[,i]*sto+D3simu.niche[,i]*(1-sto)
  }
  predicted.CS<-sapply(c(0,0.1,0.2,0.3,0.4),function(noise){
    message("noise=",noise)
    if (noise==0){
      Noise=0
    } else {
      Noise<-sapply(1:ncol(D1), function(i){
        meanra=mean(D1[,i])
        Ns<-rnorm(n=11,mean=0,sd=noise*meanra)
        Ns
      })}
    D1noise=D1+Noise
    D1noise[D1noise<0]<-0
    D1noise[D1noise>1]<-0
    
    if (noise==0){
      Noise=0
    } else {
      Noise<-sapply(1:ncol(D2), function(i){
        meanra=mean(D2[,i])
        Ns<-rnorm(n=11,mean=0,sd=noise*meanra)
        Ns
      })}
    D2noise=D2+Noise
    D2noise[D2noise<0]<-0
    D2noise[D2noise>1]<-0
    
    if (noise==0){
      Noise=0
    } else {
      Noise<-sapply(1:ncol(D3), function(i){
        meanra=mean(D3[,i])
        Ns<-rnorm(n=11,mean=0,sd=noise*meanra)
        Ns
      })}
    D3noise=D3+Noise
    D3noise[D3noise<0]<-0
    D3noise[D3noise>1]<-0
     
    D1otus<-D1noise;D1Env<-env[env$reactor=="D1",]
    D1tpoint<-D1Env$time
    D2otus<-D2noise;D2Env<-env[env$reactor=="D2",]
    D2tpoint<-D2Env$time
    D3otus<-D3noise;D3Env<-env[env$reactor=="D3",]
    D3tpoint<-D3Env$time
    
    allresult<-data.frame(combine.t(D1otus,D1tpoint,D1Env,D2otus,D2tpoint,D2Env,D3otus,D3tpoint,D3Env))
    aics<-data.frame(neutral=allresult$aic.neutral,resource=allresult$aic.resource,combine=allresult$aic.combine)
    models<-c("neutral","consumer-resource","combined")
    allmodel<-sapply(1:nrow(allresult), function(i){
      result=ifelse(is.na(aics[i,1]), NA, models[which.min(aics[i,])])
      
    })
    cs= sum((allmodel=="neutral")*ra0,na.rm = T)/sum(ra0)
    cs
  })
  names(predicted.CS)<-c("no noise","noise10%","noise20%","noise30%","noise40%")
  result<-c(truecs,predicted.CS)
  names(result)<-c("trueCS","no noise","noise10%","noise20%","noise30%","noise40%")
  result
}

CS0<-sapply(1:1000,CSsimuw.one,CS=0);write.csv(CS0,"CS0.w_update.csv")
CS0.1<-sapply(1:1000,CSsimuw.one,CS=0.1);write.csv(CS0.1,"CS0.1.w_update.csv")
CS0.2<-sapply(1:1000,CSsimuw.one,CS=0.2);write.csv(CS0.2,"CS0.2.w_update.csv")
CS0.3<-sapply(1:1000,CSsimuw.one,CS=0.3);write.csv(CS0.3,"CS0.3.w_update.csv")
CS0.4<-sapply(1:1000,CSsimuw.one,CS=0.4);write.csv(CS0.4,"CS0.4.w_update.csv")
CS0.5<-sapply(1:1000,CSsimuw.one,CS=0.5);write.csv(CS0.5,"CS0.5.w_update.csv")
CS0.6<-sapply(1:1000,CSsimuw.one,CS=0.6);write.csv(CS0.6,"CS0.6.w_update.csv")
CS0.7<-sapply(1:1000,CSsimuw.one,CS=0.7);write.csv(CS0.7,"CS0.7.w_update.csv")
CS0.8<-sapply(1:1000,CSsimuw.one,CS=0.8);write.csv(CS0.8,"CS0.8.w_update.csv")
CS0.9<-sapply(1:1000,CSsimuw.one,CS=0.9);write.csv(CS0.9,"CS0.9.w_update.csv")
CS1<-sapply(1:1000,CSsimuw.one,CS=1);write.csv(CS1,"CS1.w_update.csv")

###5. compare predicted CS and true CS
##5.1 unweighted stochasticity###
CS0<-read.csv("CS0_update.csv",row.names = 1)
CS0<-data.frame(True.stochasticity=rep(0,100),t(CS0))
CS0.1<-read.csv("CS0.1_update.csv",row.names = 1)
CS0.1<-data.frame(True.stochasticity=rep(0.1,100),t(CS0.1))
CS0.2<-read.csv("CS0.2_update.csv",row.names = 1)
CS0.2<-data.frame(True.stochasticity=rep(0.2,100),t(CS0.2))
CS0.3<-read.csv("CS0.3_update.csv",row.names = 1)
CS0.3<-data.frame(True.stochasticity=rep(0.3,100),t(CS0.3))
CS0.4<-read.csv("CS0.4_update.csv",row.names = 1)
CS0.4<-data.frame(True.stochasticity=rep(0.4,100),t(CS0.4))
CS0.5<-read.csv("CS0.5_update.csv",row.names = 1)
CS0.5<-data.frame(True.stochasticity=rep(0.5,100),t(CS0.5))
CS0.6<-read.csv("CS0.6_update.csv",row.names = 1)
CS0.6<-data.frame(True.stochasticity=rep(0.6,100),t(CS0.6))
CS0.7<-read.csv("CS0.7_update.csv",row.names = 1)
CS0.7<-data.frame(True.stochasticity=rep(0.7,100),t(CS0.7))
CS0.8<-read.csv("CS0.8_update.csv",row.names = 1)
CS0.8<-data.frame(True.stochasticity=rep(0.8,100),t(CS0.8))
CS0.9<-read.csv("CS0.9_update.csv",row.names = 1)
CS0.9<-data.frame(True.stochasticity=rep(0.9,100),t(CS0.9))
CS1<-read.csv("CS1_update.csv",row.names = 1)
CS1<-data.frame(True.stochasticity=rep(1,100),t(CS1))
alldat<-data.frame(rbind(CS0,CS0.1,CS0.2,CS0.3,CS0.4,CS0.5,CS0.6,
                         CS0.7,CS0.8,CS0.9,CS1))

#normalizaion based on STmin and STmax from the no-noise data simulated
#ST.n=(ST-STmin)/(STmax-STmin)
STmin=mean(alldat$no.noise[alldat$True.stochasticity==0]);
STmax=mean(alldat$no.noise[alldat$True.stochasticity==1])
all.normalized=(alldat[,-1]-STmin)/(STmax-STmin)
all.n<-data.frame(True.stochasticity=alldat$True.stochasticity,all.normalized)

library(ieggr)
plot(x=all.n$True.stochasticity,y=all.n$noise40.,col="blue",ylim=c(0, 1), xlab="True ST.uw",ylab="Predicted ST.uw")
abline(lm(all.n$noise40.~all.n$True.stochasticity), col="red", lwd=2) 
ccc(x=all.n$True.stochasticity,y=all.n$noise40.) ###get the precision and accuracy data

##5.2  weighted community stochasticity
CS0<-read.csv("CS0.w_update.csv",row.names = 1)
CS0<-data.frame(True.stochasticity=rep(0,100),t(CS0))
CS0.1<-read.csv("CS0.1.w_update.csv",row.names = 1)
CS0.1<-data.frame(True.stochasticity=rep(0.1,100),t(CS0.1))
CS0.2<-read.csv("CS0.2.w_update.csv",row.names = 1)
CS0.2<-data.frame(True.stochasticity=rep(0.2,100),t(CS0.2))
CS0.3<-read.csv("CS0.3.w_update.csv",row.names = 1)
CS0.3<-data.frame(True.stochasticity=rep(0.3,100),t(CS0.3))
CS0.4<-read.csv("CS0.4.w_update.csv",row.names = 1)
CS0.4<-data.frame(True.stochasticity=rep(0.4,100),t(CS0.4))
CS0.5<-read.csv("CS0.5.w_update.csv",row.names = 1)
CS0.5<-data.frame(True.stochasticity=rep(0.5,100),t(CS0.5))
CS0.6<-read.csv("CS0.6.w_update.csv",row.names = 1)
CS0.6<-data.frame(True.stochasticity=rep(0.6,100),t(CS0.6))
CS0.7<-read.csv("CS0.7.w_update.csv",row.names = 1)
CS0.7<-data.frame(True.stochasticity=rep(0.7,100),t(CS0.7))
CS0.8<-read.csv("CS0.8.w_update.csv",row.names = 1)
CS0.8<-data.frame(True.stochasticity=rep(0.8,100),t(CS0.8))
CS0.9<-read.csv("CS0.9.w_update.csv",row.names = 1)
CS0.9<-data.frame(True.stochasticity=rep(0.9,100),t(CS0.9))
CS1<-read.csv("CS1.w_update.csv",row.names = 1)
CS1<-data.frame(True.stochasticity=rep(1,100),t(CS1))
alldat<-data.frame(rbind(CS0,CS0.1,CS0.2,CS0.3,CS0.4,CS0.5,CS0.6,
                         CS0.7,CS0.8,CS0.9,CS1))

#normalizaion based on STmin and STmax from the no-noise data simulated
#ST.n=(ST-STmin)/(STmax-STmin)
STmin=mean(alldat$no.noise[alldat$True.stochasticity==0]);
STmax=mean(alldat$no.noise[alldat$True.stochasticity==1])
all.normalized=(alldat[,-(1:2)]-STmin)/(STmax-STmin)
all.n<-data.frame(True.stochasticity=alldat$trueCS,all.normalized)

library(ieggr)
plot(x=all.n$True.stochasticity,y=all.n$noise40.,col="blue",ylim=c(0, 1), xlab="True ST.w",ylab="Predicted ST.w")
abline(lm(all.n$noise40.~all.n$True.stochasticity), col="red", lwd=2) 
ccc(x=all.n$True.stochasticity,y=all.n$noise40.)  ###get the precision and accuracy data

##if you cannot install ieggr
##then load the ccc function as below  
ccc<-function (x, y, out.vector = FALSE) {
  d <- y - x
  m1 <- mean(y)
  m2 <- mean(x)
  v1 <- var(y)
  v2 <- var(x)
  n <- length(d)
  e2 <- sum(d^2)/n
  if (length(y) != length(x)) 
    stop("x and y should have the same length!")
  mu_d <- m1 - m2
  d2 <- mu_d^2
  s12 <- v1 * (n - 1)/n
  s22 <- v2 * (n - 1)/n
  U <- mu_d/sqrt(sqrt(s12 * s22))
  V <- sqrt(s12/s22)
  Ca <- 2/(V + 1/V + U^2)
  rc <- 1 - e2/(d2 + s12 + s22)
  r <- (rc/Ca)
  out = list(accuracy = Ca, precision = r, ccc = rc)
  if (out.vector) {
    out = unlist(out)
  }
  out
}

