

setwd("/Users/linwei/Dropbox/bioreactor/neutral modeling/ratio model/1_model fitting")
source("modelcompare.R")

#read otu table
otutab<-read.csv("/Users/linwei/Dropbox/bioreactor/neutral modeling/ratio model/0_raw data/otutab.csv",header = T,row.names = 1,check.names = F)
env<-read.csv("/Users/linwei/Dropbox/bioreactor/neutral modeling/ratio model/0_raw data/env.csv",header = T,row.names = 1)
midasclassifier<-read.csv("~/Dropbox/bioreactor/neutral modeling/0_raw data/Midas.rdp_assigned_taxonomy/otus_resampled_tax_assignments_Midas.txt",sep="\t",header = F,row.names = 1)


rdpclassifier<-otutab$taxonomy
names(rdpclassifier)<-row.names(otutab)
midasclassifier<-midasclassifier[match(names(rdpclassifier),row.names(midasclassifier)),]

otutab<-otutab[,match(row.names(env),colnames(otutab))]
sum(row.names(env)==colnames(otutab))  #check


otutab<-t(otutab)
otutab<-otutab/rowSums(otutab)   #convert to relative abundance

sum(row.names(env)==row.names(otutab))
###select the reference otu： OTU_1, which is detected across all samples


##test on control reactors one by one
otus<-otutab[env$reactor=="C3",]
Env<-env[env$reactor=="C3",]
tpoint<-Env$time

# y = delta(log x/xf)/delta t
# x1= 1/x; x2=1/xf; x3=vs;
# neutral model: y=b0+b1*x1+b2*x2+error
# consumer model: y=b0 + b3*x3
# combined model: y=b0+b1*x1+b2*x2+ b3*x3 + error
#weight=delta t/(2*x1+2*x2)

allresult<-sapply(2:ncol(otus),function(i){
  message("i =", i)
  xf<-otus[,"OTU_1"]
  x<-otus[,i]
  x0<-x[x>0]; #x initial value
  x1<-1/x0;  #x1 initial value
  xf<-xf[x>0];   #xf initial value
  x2<-1/xf;    #x2 initial value
  t0<-tpoint[x>0]
  
  x3<-(Env$VS); x3<-x3[x>0] #x3 intitial value
  lratio=log(x0/xf)  ##log-ratio
  dlrati<-lratio[-1]-lratio[-length(x0)];dt<-t0[-1]-t0[-length(x0)]
  y<-dlrati/dt  
  x1=x1[-length(x0)];   #x1
  x2=x2[-length(x0)];  #x2
  x3=x3[-length(x0)];   #x3
  weight=dt/(2*x1+2*x2)
  y[dt>20]<-NA   ##exclude y when time interval is more than 20 days
  if (sum(!is.na(y))<7 | length(unique(x1)) < 6){
    result<-rep(NA,40)
  } else{
    result<-comparemodel(y=y,x1=x1,x2=x2,x3=x3,Weight =weight )
    if (length(result)<40) result<-rep(NA,40)
  }
  result
})

colnames(allresult)<-colnames(otus)[-1]
allresult<-t(allresult)
sum(row.names(allresult)==names(rdpclassifier)[-1])
allresult<-data.frame(allresult,rdpclassifier=rdpclassifier[-1],midasclassifier=midasclassifier[-1,1])
write.csv(allresult,"C3_fitting.csv")



######  D reactors   ##############################


D1comm<-otutab[env$reactor=="D1",]
D1env<-env[env$reactor=="D1",]
D2comm<-otutab[env$reactor=="D2",]
D2env<-env[env$reactor=="D2",]
D3comm<-otutab[env$reactor=="D3",]
D3env<-env[env$reactor=="D3",]


#D 45-100 days
D1id<-which(D1env$time<100)  ##include all time points
D1otus<-D1comm[D1id,];D1Env<-D1env[D1id,]
D1tpoint<-D1Env$time

D2id<-which(D2env$time<100)  ##include all time points
D2otus<-D2comm[D2id,];D2Env<-D2env[D2id,]
D2tpoint<-D2Env$time

D3id<-which(D3env$time<100)  ##include all time points
D3otus<-D3comm[D3id,];D3Env<-D3env[D3id,]
D3tpoint<-D3Env$time

##For comparison, select C reactors of the same time points as D
##remove '#' for the below code to run for C samples
#D1comm<-otutab[env$reactor=="C1",]
#D1env<-env[env$reactor=="C1",]
#D2comm<-otutab[env$reactor=="C2",]
#D2env<-env[env$reactor=="C2",]
#D3comm<-otutab[env$reactor=="C3",]
#D3env<-env[env$reactor=="C3",]

#D1id<-which(D1env$time<100 & D1env$time>44 & D1env$time!=59)  ##include all time points
#D1otus<-D1comm[D1id,];D1Env<-D1env[D1id,]
#D1tpoint<-D1Env$time

#D2id<-which(D2env$time<100 & D2env$time>44 & D2env$time!=59)  ##include all time points
#D2otus<-D2comm[D2id,];D2Env<-D2env[D2id,]
#D2tpoint<-D2Env$time

#D3id<-which(D3env$time<100 & D3env$time>44 & D3env$time!=59)  ##include all time points
#D3otus<-D3comm[D3id,];D3Env<-D3env[D3id,]
#D3tpoint<-D3Env$time

allresult<-sapply(2:ncol(D1otus),function(i){
  message("i =", i)
  
  D1xf<-D1otus[,"OTU_1"]
  D1x<-D1otus[,i]
  D1x0<-D1x[D1x>0]; #x initial value
  D1x1<-1/D1x0;  #x1 initial value
  D1xf<-D1xf[D1x>0];   #xf initial value
  D1x2<-1/D1xf;    #x2 initial value
  D1t0<-D1tpoint[D1x>0]
  D1x3<-(D1Env$VS); D1x3<-D1x3[D1x>0] #x3 intitial value
  D1lratio=log(D1x0/D1xf)  ##log-ratio
  D1dlrati<-D1lratio[-1]-D1lratio[-length(D1x0)];D1dt<-D1t0[-1]-D1t0[-length(D1x0)]
  D1y<-D1dlrati/D1dt  
  D1x1=D1x1[-length(D1x0)];   #x1
  D1x2=D1x2[-length(D1x0)];  #x2
  D1x3=D1x3[-length(D1x0)];   #x3
  D1weight=D1dt/(2*D1x1+2*D1x2)
  
  D2xf<-D2otus[,"OTU_1"]
  D2x<-D2otus[,i]
  D2x0<-D2x[D2x>0]; #x initial value
  D2x1<-1/D2x0;  #x1 initial value
  D2xf<-D2xf[D2x>0];   #xf initial value
  D2x2<-1/D2xf;    #x2 initial value
  D2t0<-D2tpoint[D2x>0]
  D2x3<-(D2Env$VS); D2x3<-D2x3[D2x>0] #x3 intitial value
  D2lratio=log(D2x0/D2xf)  ##log-ratio
  D2dlrati<-D2lratio[-1]-D2lratio[-length(D2x0)];D2dt<-D2t0[-1]-D2t0[-length(D2x0)]
  D2y<-D2dlrati/D2dt  
  D2x1=D2x1[-length(D2x0)];   #x1
  D2x2=D2x2[-length(D2x0)];  #x2
  D2x3=D2x3[-length(D2x0)];   #x3
  D2weight=D2dt/(2*D2x1+2*D2x2)
  
  D3xf<-D3otus[,"OTU_1"]
  D3x<-D3otus[,i]
  D3x0<-D3x[D3x>0]; #x initial value
  D3x1<-1/D3x0;  #x1 initial value
  D3xf<-D3xf[D3x>0];   #xf initial value
  D3x2<-1/D3xf;    #x2 initial value
  D3t0<-D3tpoint[D3x>0]
  D3x3<-(D3Env$VS); D3x3<-D3x3[D3x>0] #x3 intitial value
  D3lratio=log(D3x0/D3xf)  ##log-ratio
  D3dlrati<-D3lratio[-1]-D3lratio[-length(D3x0)];D3dt<-D3t0[-1]-D3t0[-length(D3x0)]
  D3y<-D3dlrati/D3dt  
  D3x1=D3x1[-length(D3x0)];   #x1
  D3x2=D3x2[-length(D3x0)];  #x2
  D3x3=D3x3[-length(D3x0)];   #x3
  D3weight=D3dt/(2*D3x1+2*D3x2)
  
  y<-c(D1y,D2y,D3y)
  x1<-c(D1x1,D2x1,D3x1)
  x2<-c(D1x2,D2x2,D3x2)
  x3<-c(D1x3,D2x3,D3x3)
  weight<-c(D1weight,D2weight,D3weight)
  dt<-c(D1dt,D2dt,D3dt)
  
  y[dt>20]<-NA   ##exclude y when time interval is more than 20 days
  if (sum(!is.na(y))<7 | length(unique(x1)) < 6){
    result<-rep(NA,40)
  } else{
    result<-comparemodel(y=y,x1=x1,x2=x2,x3=x3,Weight =weight )
    if (length(result)<40) result<-rep(NA,40)
  }
  result
})
colnames(allresult)<-colnames(D1otus)[-1]
allresult<-t(allresult)
sum(row.names(allresult)==names(rdpclassifier)[-1])
allresult<-data.frame(allresult,rdpclassifier=rdpclassifier[-1],midasclassifier=midasclassifier[-1,1])
write.csv(allresult,"D_before100 ra fitting.csv")
#write.csv(allresult,"C_before100 ra fitting.csv")

