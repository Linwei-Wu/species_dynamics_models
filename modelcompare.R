
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


# x focal species
# xf reference species
#d log x/xf /dt=b0+b1*(1/x)+b2*(1/xf)+b3*R + error  #the combined model
#weight = dt/（2/x+2/xf)

# y = delta(log x/xf)/delta t
# x1= 1/x; x2=1/xf; x3=vs;

# neutral model: y=b0+b1*x1+b2*x2+error
# consumer model: y=b0 + b3*x3
# combined model: y=b0+b1*x1+b2*x2+ b3*x3 + error
#weight=delta t/(2*x1+2*x2)

##compare AIC of different model, and return parameters of each model##
comparemodel<-function(y,x1,x2,x3,Weight){
  neutral.model<-lm(y~x1+x2,weights = Weight)   ##pure neutral model
  resource.model<-lm(y~x3)       ##resource-consumer model
  combine.model<-lm(y~x1+x2+x3,weights = Weight)   ##combined model
  
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
  
  mean.ra<-mean(1/x1,na.rm = T)
  names(mean.ra)<-"mean.ra"
  
  
  result<-c(coef.neutral,coefP.neutral,other.neutral,coef.resource,coefP.resource,other.resource,
            coef.combine,coefP.combine,other.combine,r2.psedo,mean.ra)
  result
}

