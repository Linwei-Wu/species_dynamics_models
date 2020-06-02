
####this is the function to calculate STmin and STmax for the three replicated time-series in treatment reactors,
####then calculated the normalized community stochasticity 

##input otus (species relative abundance data), Env that contains the resource 
##variable, tpoint for time point (days)
##otus and Env should have matched rownames (samples)

CSnormlz.reps<-function(D1otus,D1tpoint,D1Env,D2otus,D2tpoint,D2Env,D3otus,D3tpoint,D3Env){
  ###1. generate the initial community###
  meanra<-rowMeans(cbind(colMeans(D1otus),colMeans(D2otus),colMeans(D3otus)))  #initial value of relative abundance
  
  ###2. simulate niche relative abundance 
  #get niche parameters for each species
  fitting.result<-data.frame(combine.t(D1otus,D1tpoint,D1Env,D2otus,D2tpoint,D2Env,D3otus,D3tpoint,D3Env))
  rownames(fitting.result)<-colnames(D1otus)
  OTUid<-row.names(fitting.result)[!is.na(fitting.result$x1.neutral)]  ##OTU id used
  fitting.result<-fitting.result[OTUid,]
  meanra<-meanra[OTUid]
  niche.p<-data.frame(meanra=meanra,fitting.result[,c("x1.resource","x5.resource","x6.resource")])
  #niche.p is the niche parameter used to generate the ra simulations
  #niche simulation
  message("now simulate species ra under pure niche model")
  #D1
  nstep=(max(D1tpoint)-min(D1tpoint))*12+1
  vssimu<-numeric(nstep)  ##simulate VS (volatile solids), 2 hours per step
  for (i in 1:(length(D1tpoint)-1)){
    vssimu[((D1tpoint[i]-min(D1tpoint))*12+1):((D1tpoint[i+1]-min(D1tpoint))*12+1)] <- D1Env$VS[i]
  }
  
  D1simu.niche<-sapply(1:nrow(niche.p),function(i){
    simuonce(x0=niche.p[i,1],b1=niche.p[i,2],b2=niche.p[i,3],b3=niche.p[i,4],vssimu=vssimu,tpoint=D1tpoint)})
  colnames(D1simu.niche)<-row.names(niche.p)
  D1simu.niche[is.na(D1simu.niche)]<-0
  #D2
  nstep=(max(D2tpoint)-min(D2tpoint))*12+1
  vssimu<-numeric(nstep)  ##simulate VS (volatile solids), 2 hours per step
  for (i in 1:(length(D2tpoint)-1)){
    vssimu[((D2tpoint[i]-min(D2tpoint))*12+1):((D2tpoint[i+1]-min(D2tpoint))*12+1)] <- D2Env$VS[i]
  }
  
  D2simu.niche<-sapply(1:nrow(niche.p),function(i){
    simuonce(x0=niche.p[i,1],b1=niche.p[i,2],b2=niche.p[i,3],b3=niche.p[i,4],vssimu=vssimu,tpoint=D2tpoint)})
  colnames(D2simu.niche)<-row.names(niche.p)
  D2simu.niche[is.na(D2simu.niche)]<-0
  #D3
  nstep=(max(D3tpoint)-min(D3tpoint))*12+1
  vssimu<-numeric(nstep)  ##simulate VS (volatile solids), 2 hours per step
  for (i in 1:(length(D3tpoint)-1)){
    vssimu[((D3tpoint[i]-min(D3tpoint))*12+1):((D3tpoint[i+1]-min(D3tpoint))*12+1)] <- D3Env$VS[i]
  }
  
  D3simu.niche<-sapply(1:nrow(niche.p),function(i){
    simuonce(x0=niche.p[i,1],b1=niche.p[i,2],b2=niche.p[i,3],b3=niche.p[i,4],vssimu=vssimu,tpoint=D3tpoint)})
  colnames(D3simu.niche)<-row.names(niche.p)
  D3simu.niche[is.na(D3simu.niche)]<-0
   ###3. neutral simulation##
  neutral.p<-data.frame(meanra=meanra,Ntm=(-fitting.result$x1.neutral)/(fitting.result$SE.neutral)^2,
                        pi=fitting.result$X.Intercept..neutral/(-fitting.result$x1.neutral),a=1/(fitting.result$SE.neutral)^2)
  #neutral.p is the neutral parameters
  message("now simulate species ra under pure neutral model")
  D1simu.neutral<-sapply(1:nrow(neutral.p),function(i){
    simuone(Ntm=neutral.p[i,2],pi=neutral.p[i,3],a=neutral.p[i,4],x0=neutral.p[i,1],tpoint=D1tpoint)
  })
  colnames(D1simu.neutral)<-row.names(neutral.p)
  D1simu.neutral[is.na(D1simu.neutral)]<-0
  
  D2simu.neutral<-sapply(1:nrow(neutral.p),function(i){
    simuone(Ntm=neutral.p[i,2],pi=neutral.p[i,3],a=neutral.p[i,4],x0=neutral.p[i,1],tpoint=D2tpoint)
  })
  colnames(D2simu.neutral)<-row.names(neutral.p)
  D2simu.neutral[is.na(D2simu.neutral)]<-0
  
  D3simu.neutral<-sapply(1:nrow(neutral.p),function(i){
    simuone(Ntm=neutral.p[i,2],pi=neutral.p[i,3],a=neutral.p[i,4],x0=neutral.p[i,1],tpoint=D2tpoint)
  })
  colnames(D3simu.neutral)<-row.names(neutral.p)
  D3simu.neutral[is.na(D3simu.neutral)]<-0
  ###4. get STmin and STmax under CS=0 and CS=1
  ##4.1 unweighted CS
  message("now doing simulation analyses to get STmin and STmax")
  message("simulation on unweighted CS")
  CS0.uw<-sapply(1:100,function(k){
    CSsimu.one.reps(CS=0,D1simu.neutral,D1simu.niche,D1Env,D1tpoint,
                    D2simu.neutral,D2simu.niche,D2Env,D2tpoint,
                    D3simu.neutral,D3simu.niche,D3Env,D3tpoint)})
  CS1.uw<-sapply(1:100,function(m){
    CSsimu.one.reps(CS=1,D1simu.neutral,D1simu.niche,D1Env,D1tpoint,
                    D2simu.neutral,D2simu.niche,D2Env,D2tpoint,
                    D3simu.neutral,D3simu.niche,D3Env,D3tpoint)})
  STmin.uw=mean(CS0.uw,na.rm = T)
  STmax.uw=mean(CS1.uw,na.rm = T)
  ##4.2 weighted CS
  message("simulation on weighted CS")
  CS0.w<-sapply(1:100,function(k){
    CSsimuw.one.reps(CS=0,ra0=meanra,D1simu.neutral,D1simu.niche,D1Env,D1tpoint,
                D2simu.neutral,D2simu.niche,D2Env,D2tpoint,
                D3simu.neutral,D3simu.niche,D3Env,D3tpoint)})
  CS1.w<-sapply(1:100,function(m){
    CSsimuw.one.reps(CS=1,ra0=meanra,D1simu.neutral,D1simu.niche,D1Env,D1tpoint,
                D2simu.neutral,D2simu.niche,D2Env,D2tpoint,
                D3simu.neutral,D3simu.niche,D3Env,D3tpoint)})
  STmin.w=mean(CS0.w,na.rm = T)
  STmax.w=mean(CS1.w,na.rm = T)
  
  ###5. normalize observed CS
  message("now normalize observed CS using STmin and STmax")
  aics<-data.frame(neutral=fitting.result$aic.neutral,resource=fitting.result$aic.resource,combine=fitting.result$aic.combine)
  models<-c("neutral","consumer-resource","combined")
  allmodel<-sapply(1:nrow(fitting.result), function(i){
    result=ifelse(is.na(aics[i,1]), NA, models[which.min(aics[i,])])
  })
  cs.w= sum((allmodel=="neutral")*meanra,na.rm = T)/sum(meanra)
  cs.uw= sum((allmodel=="neutral"),na.rm = T)/sum(!is.na(allmodel))
  res<-c(cs.uw,cs.w,(cs.uw-STmin.uw)/(STmax.uw-STmin.uw),(cs.w-STmin.w)/(STmax.w-STmin.w),
            STmin.uw,STmax.uw,STmin.w,STmax.w,
            sd(CS0.uw,na.rm = T),sd(CS1.uw,na.rm = T),sd(CS0.w,na.rm = T),sd(CS1.w,na.rm = T))
  names(res)<-c("cs.uw","cs.w","CS.uw.n","CS.w.n","STmin.uw","STmax.uw","STmin.w","STmax.w",
                "STmin.uw.sd","STmax.uw.sd","STmin.w.sd","STmax.w.sd")
  res
}


#function for niche simulation
simuonce<-function(x0,b1,b2,b3,vssimu,tpoint){
  vec=numeric(length(vssimu))
  vec[1]=x0
  for (i in 1:(length(vssimu)-1)){
    dx=b1*vec[i]/12+b2*vec[i]*vssimu[i]/12+b3*vec[i]*vec[i]/12
    vec[i+1]<-ifelse(vec[i]+dx>0,vec[i]+dx,0)
  }
  simu1<-c()
  for (i in 1:length(tpoint)){
    simu1[i]<-vec[((tpoint[i]-tpoint[1])*12+1)]
  }
  simu1
}

#function for neutral simulation
simuone<-function (Ntm,pi,a,x0,tpoint){
  vec=numeric(max(tpoint)*12+1)
  vec[1]=x0  ##initial relative abundance
  dwt<-rnorm(max(tpoint)*12)/(12^0.5)
  for (i in 1:(max(tpoint)*12)){
    dx=Ntm*pi/a/12-Ntm/a*vec[i]/12+((vec[i]*(1-vec[i])*2/a)^0.5)*dwt[i]
    vec[i+1]<-ifelse(vec[i]+dx>0,vec[i]+dx,0)
  }
  simu1<-c()
  for (i in 1:length(tpoint)){
    simu1[i]<-vec[tpoint[i]*12]
  }
  simu1
}

#function for model fitting
source("main modeling function/fitting and comparing models.R")

##function to get predicted CS (community stochasticity), unweighted
CSsimu.one.reps<-function(CS,D1simu.neutral,D1simu.niche,D1Env,D1tpoint,
                     D2simu.neutral,D2simu.niche,D2Env,D2tpoint,
                     D3simu.neutral,D3simu.niche,D3Env,D3tpoint){
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
  D1[D1<0]<-0;D2[D2<0]<-0;D3[D3<0]<-0
  D1[D1>1]<-0;D2[D2>1]<-0;D3[D3>1]<-0
  allresult<-data.frame(combine.t(D1,D1tpoint,D1Env,D2,D2tpoint,D2Env,D3,D3tpoint,D3Env))
  aics<-data.frame(neutral=allresult$aic.neutral,resource=allresult$aic.resource,combine=allresult$aic.combine)
  models<-c("neutral","consumer-resource","combined")
  allmodel<-sapply(1:nrow(allresult), function(i){
    result=ifelse(is.na(aics[i,1]), NA, models[which.min(aics[i,])])
  })
  testtable<-table(allmodel)
  cs=testtable[3]/sum(testtable)
  cs
}

##function to get predicted CS, weighted
CSsimuw.one.reps<-function(CS,ra0,D1simu.neutral,D1simu.niche,D1Env,D1tpoint,
                           D2simu.neutral,D2simu.niche,D2Env,D2tpoint,
                           D3simu.neutral,D3simu.niche,D3Env,D3tpoint){
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
  D1[D1<0]<-0;D2[D2<0]<-0;D3[D3<0]<-0
  D1[D1>1]<-0;D2[D2>1]<-0;D3[D3>1]<-0
  
  allresult<-data.frame(combine.t(D1,D1tpoint,D1Env,D2,D2tpoint,D2Env,D3,D3tpoint,D3Env))
  aics<-data.frame(neutral=allresult$aic.neutral,resource=allresult$aic.resource,combine=allresult$aic.combine)
  models<-c("neutral","consumer-resource","combined")
  allmodel<-sapply(1:nrow(allresult), function(i){
    result=ifelse(is.na(aics[i,1]), NA, models[which.min(aics[i,])])
  })
  cs= sum((allmodel=="neutral")*ra0,na.rm = T)/sum(ra0)
  cs
}

