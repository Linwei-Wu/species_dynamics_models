##the function 'comparemodel' is to compare AIC of different model, and return the  parameters of each model##

##dx/dt=b0+b1*x+b2*x(1-x)+b3*(VS*x(1-x))+b4*(1-x)*(x^2)+error  #the combined model
#VS is volatile solids, representing the resource level in the bioreators
#if you have more than one resource variable, you should modify the function accordingly (i.e., adding more predictors).

##x1=x, x2=x(1-x), x3=vs*x(1-x),x4=x*x(1-x)
##weight=dt/(2*x1*(1-x1))
#x5=vs*x,x6=x^2 for the consumer-resource model

comparemodel<-function(y,x1,x2,x3,x4,x5,x6,Weight){
  neutral.model<-lm(y~x1,weights = Weight)   ##pure neutral model
  resource.model<-lm(y~x1+x5+x6+0)       ##resource-consumer model
  combine.model<-lm(y~x1+x2+x3+x4,weights = Weight)   ##combined model
  
  neutral.summary<-summary(neutral.model)
  resource.summary<-summary(resource.model)
  combine.summary<-summary(combine.model)
  
  #neutral model
  coef.neutral<-coef(neutral.summary)[ , "Estimate"] 
  coefP.neutral<-coef(neutral.summary)[ , "Pr(>|t|)"] 
  names(coef.neutral)<-paste0(names(coef.neutral),".neutral")
  names(coefP.neutral)<-paste0(names(coefP.neutral),".neutral.P")
  r2.neutral<- neutral.summary$r.squared
  SE.neutral<- neutral.summary$sigma
  aic.neutral<-AIC(neutral.model)
  bic.neutral<-BIC(neutral.model)
  f.neutral <-  neutral.summary$fstatistic
  modP.neutral <- pf(f.neutral[1],f.neutral[2],f.neutral[3],lower.tail=F)
  attributes( modP.neutral) <- NULL
  df.neutral<-neutral.model$df.residual
  other.neutral<-c(r2.neutral,SE.neutral,aic.neutral,bic.neutral,modP.neutral,df.neutral)
  names(other.neutral)<-c("r2.neutral","SE.neutral","aic.neutral","bic.neutral","modP.neutral","df.neutral")
  
  #resource model
  coef.resource<-coef(resource.summary)[ , "Estimate"] 
  coefP.resource<-coef(resource.summary)[ , "Pr(>|t|)"] 
  names(coef.resource)<-paste0(names(coef.resource),".resource")
  names(coefP.resource)<-paste0(names(coefP.resource),".resource.P")
  r2.resource<- resource.summary$r.squared
  SE.resource<- resource.summary$sigma
  aic.resource<-AIC(resource.model)
  bic.resource<-BIC(resource.model)
  f.resource <-  resource.summary$fstatistic
  modP.resource <- pf(f.resource[1],f.resource[2],f.resource[3],lower.tail=F)
  attributes( modP.resource) <- NULL
  df.resource<-resource.model$df.residual
  other.resource<-c(r2.resource,SE.resource,aic.resource,bic.resource,modP.resource,df.resource)
  names(other.resource)<-c("r2.resource","SE.resource","aic.resource","bic.resource","modP.resource","df.resource")
  
  #combined model
  coef.combine<-coef(combine.summary)[ , "Estimate"] 
  coefP.combine<-coef(combine.summary)[ , "Pr(>|t|)"] 
  names(coef.combine)<-paste0(names(coef.combine),".combine")
  names(coefP.combine)<-paste0(names(coefP.combine),".combine.P")
  r2.combine<- combine.summary$r.squared
  SE.combine<- combine.summary$sigma
  aic.combine<-AIC(combine.model)
  bic.combine<-BIC(combine.model)
  f.combine <-  combine.summary$fstatistic
  modP.combine <- pf(f.combine[1],f.combine[2],f.combine[3],lower.tail=F)
  attributes( modP.combine) <- NULL
  df.combine<-combine.model$df.residual
  other.combine<-c(r2.combine,SE.combine,aic.combine,bic.combine,modP.combine,df.combine)
  names(other.combine)<-c("r2.combine","SE.combine","aic.combine","bic.combine","modP.combine","df.combine")
  
  r2.psedo<-c(r2(neutral.model),r2(resource.model),r2(combine.model))
  names(r2.psedo)<-c("r2psedo.neutral","r2psedo.resource","r2psedo.combine")
  
  mean.ra<-mean(x1,na.rm = T)
  names(mean.ra)<-"mean.ra"
  
  
  result<-c(coef.neutral,coefP.neutral,other.neutral,coef.resource,coefP.resource,other.resource,
            coef.combine,coefP.combine,other.combine,r2.psedo,mean.ra)
  result
}

r2 <- function(x){  
  SSe <- sum((x$residuals)^2);  
  observed <- x$residuals+x$fitted.values;  
  SSt <- sum((observed-mean(observed))^2);  
  value <- 1-SSe/SSt;  
  return(value);  
} 

r2ww <- function(x){
  SSe <- sum(x$weights*(x$residuals)^2); #the residual sum of squares is weighted
  observed <- x$residuals+x$fitted.values;
  SSt <- sum(x$weights*(observed-weighted.mean(observed,x$weights))^2); #the total sum of squares is weighted      
  value <- 1-SSe/SSt;
  return(value);
}

#####below is the function of model fitting on a single time series dataset
#otus is the community matrix, rows are samples ordered by time, columns are species, entries are relative abundances.
#Env contains the variable on resource level. In this study, we only have one variable
#that is VS (volatile solids).

one.t<-function(otus,tpoint,Env){
  t(sapply(1:ncol(otus),function(i){
    #message("i =", i)
    x<-otus[,i]
    x0<-x[x>0];    #x1 initial value
    t0<-tpoint[x>0]
    x2<-x0*(1-x0)  #x2 initial value
    x3<-(Env$VS)*x*(1-x); x3<-x3[x>0] #x3 intitial value
    x4<-x0*x0*(1-x0)  #x4 initial value
    xvs<-(Env$VS)*x;x5<-xvs[x>0]  #x5 initial value
    x6<-x0^2   #x6 initial value
    dx<-x0[-1]-x0[-length(x0)];dt<-t0[-1]-t0[-length(x0)]
    y<-dx/dt  
    x1=x0[-length(x0)];   #x1
    x2=x2[-length(x0)];  #x2
    x3=x3[-length(x0)];   #x3
    x4=x4[-length(x0)];   #x4
    x5=x5[-length(x0)];   #x5
    x6=x6[-length(x0)];   #x6
    weight=dt/(2*x1*(1-x1))
    y[dt>20]<-NA   ##exclude y when time interval is more than 20 days
    if (sum(!is.na(y))<7 | length(unique(x1)) < 6){
      result<-rep(NA,42)
      names(result)<-c("X.Intercept..neutral","x1.neutral","X.Intercept..neutral.P","x1.neutral.P","r2.neutral","SE.neutral",            
                        "aic.neutral", "bic.neutral","modP.neutral","df.neutral", "x1.resource", "x5.resource", "x6.resource",           
                       "x1.resource.P" ,"x5.resource.P", "x6.resource.P", "r2.resource" , "SE.resource" ,"aic.resource",           
                       "bic.resource","modP.resource", "df.resource", "X.Intercept..combine","x1.combine", "x2.combine",       
                       "x3.combine","x4.combine" , "X.Intercept..combine.P" ,"x1.combine.P" ,"x2.combine.P" ,          
                       "x3.combine.P" ,"x4.combine.P" , "r2.combine","SE.combine","aic.combine","bic.combine",             
                       "modP.combine" , "df.combine" , "r2psedo.neutral", "r2psedo.resource" ,"r2psedo.combine" ,"mean.ra" )      

    } else{
      result<-comparemodel(y=y,x1=x1,x2=x2,x3=x3,x4=x4,x5=x5,x6=x6,Weight =weight )
      if (length(result)<42) result<-rep(NA,42)
    }
    result
  }))
}


#####below is the function to combine the replicated time series from treatment reactors, then fit the models##

combine.t<-function(D1otus,D1tpoint,D1Env,D2otus,D2tpoint,D2Env,D3otus,D3tpoint,D3Env){
  t(sapply(1:ncol(D1otus),function(i){
    D1x<-D1otus[,i]
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
    
    D2x<-D2otus[,i]
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
    
    D3x<-D3otus[,i]
    D3x0<-D3x[D3x>0];    #x1 initial value
    D3t0<-D3tpoint[D3x>0]
    D3x2<-D3x0*(1-D3x0)  #x2 initial value
    D3x3<-(D3Env$VS)*D3x*(1-D3x); D3x3<-D3x3[D3x>0] #x3 intitial value
    D3x4<-D3x0*D3x0*(1-D3x0)  #x4 initial value
    D3xvs<-(D3Env$VS)*D3x;D3x5<-D3xvs[D3x>0]  #x5 initial value
    D3x6<-D3x0^2   #x6 initial value
    D3dx<-D3x0[-1]-D3x0[-length(D3x0)];D3dt<-D3t0[-1]-D3t0[-length(D3x0)]
    D3y<-D3dx/D3dt  
    D3x1=D3x0[-length(D3x0)];   #x1
    D3x2=D3x2[-length(D3x0)];  #x2
    D3x3=D3x3[-length(D3x0)];   #x3
    D3x4=D3x4[-length(D3x0)];   #x4
    D3x5=D3x5[-length(D3x0)];   #x5
    D3x6=D3x6[-length(D3x0)];   #x6
    D3weight=D3dt/(2*D3x1*(1-D3x1))
    
    y<-c(D1y,D2y,D3y)
    x1<-c(D1x1,D2x1,D3x1)
    x2<-c(D1x2,D2x2,D3x2)
    x3<-c(D1x3,D2x3,D3x3)
    x4<-c(D1x4,D2x4,D3x4)
    x5<-c(D1x5,D2x5,D3x5)
    x6<-c(D1x6,D2x6,D3x6)
    weight<-c(D1weight,D2weight,D3weight)
    dt<-c(D1dt,D2dt,D3dt)
    
    y[dt>20]<-NA   ##exclude y when time interval is more than 20 days
    if (sum(!is.na(y))<7 | length(unique(x1)) < 6){
      result<-rep(NA,42)
      names(result)<-c("X.Intercept..neutral","x1.neutral","X.Intercept..neutral.P","x1.neutral.P","r2.neutral","SE.neutral",            
                       "aic.neutral", "bic.neutral","modP.neutral","df.neutral", "x1.resource", "x5.resource", "x6.resource",           
                       "x1.resource.P" ,"x5.resource.P", "x6.resource.P", "r2.resource" , "SE.resource" ,"aic.resource",           
                       "bic.resource","modP.resource", "df.resource", "X.Intercept..combine","x1.combine", "x2.combine",       
                       "x3.combine","x4.combine" , "X.Intercept..combine.P" ,"x1.combine.P" ,"x2.combine.P" ,          
                       "x3.combine.P" ,"x4.combine.P" , "r2.combine","SE.combine","aic.combine","bic.combine",             
                       "modP.combine" , "df.combine" , "r2psedo.neutral", "r2psedo.resource" ,"r2psedo.combine" ,"mean.ra" )      
    } else{
      result<-comparemodel(y=y,x1=x1,x2=x2,x3=x3,x4=x4,x5=x5,x6=x6,Weight =weight )
      if (length(result)<42) result<-rep(NA,42)
    }
    result
  }))
  
}




