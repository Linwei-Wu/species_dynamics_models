
###model fitting for C and D, all data points
setwd("/Users/linwei/Dropbox/bioreactor/ratio model-new/1_model fitting")

C1<-read.csv("C1_fitting.csv",header = T,row.names = 1)
C2<-read.csv("C2_fitting.csv",header = T,row.names = 1)
C3<-read.csv("C3_fitting.csv",header = T,row.names = 1)
D97<-read.csv("D_before100 ra fitting.csv",header = T,row.names = 1)
C97<-read.csv("C_before100 ra fitting.csv",header = T,row.names = 1)


###############three models####

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



C1used<-data.frame(mean.ra=C1$mean.ra,r2.neutral=C1$r2psedo.neutral,r2.resource=C1$r2psedo.resource,r2.combine=C1$r2psedo.combine,modP.neutral=C1$modP.neutral,modP.resource=C1$modP.resource,modP.combine=C1$modP.combine,bestmodel=C1model,reactor=rep("C1",length(C1model)),treat=rep("C",length(C1model)),class=C1$rdpclassifier)
C2used<-data.frame(mean.ra=C2$mean.ra,r2.neutral=C2$r2psedo.neutral,r2.resource=C2$r2psedo.resource,r2.combine=C2$r2psedo.combine,modP.neutral=C2$modP.neutral,modP.resource=C2$modP.resource,modP.combine=C2$modP.combine,bestmodel=C2model,reactor=rep("C2",length(C2model)),treat=rep("C",length(C2model)),class=C2$rdpclassifier)
C3used<-data.frame(mean.ra=C3$mean.ra,r2.neutral=C3$r2psedo.neutral,r2.resource=C3$r2psedo.resource,r2.combine=C3$r2psedo.combine,modP.neutral=C3$modP.neutral,modP.resource=C3$modP.resource,modP.combine=C3$modP.combine,bestmodel=C3model,reactor=rep("C3",length(C3model)),treat=rep("C",length(C3model)),class=C3$rdpclassifier)
Dused<-data.frame(mean.ra=D97$mean.ra,r2.neutral=D97$r2psedo.neutral,r2.resource=D97$r2psedo.resource,r2.combine=D97$r2psedo.combine,modP.neutral=D97$modP.neutral,modP.resource=D97$modP.resource,modP.combine=D97$modP.combine,bestmodel=Dmodel,reactor=rep("D",length(Dmodel)),treat=rep("D",length(Dmodel)),class=C97$rdpclassifier)
Cused<-data.frame(mean.ra=C97$mean.ra,r2.neutral=C97$r2psedo.neutral,r2.resource=C97$r2psedo.resource,r2.combine=C97$r2psedo.combine,modP.neutral=C97$modP.neutral,modP.resource=C97$modP.resource,modP.combine=C97$modP.combine,bestmodel=Cmodel,reactor=rep("C",length(Dmodel)),treat=rep("C",length(Dmodel)),class=D97$rdpclassifier)


alldata<-rbind(C1used,C2used,C3used,Dused)
alldata<-alldata[!is.na(alldata$mean.ra),]
alldata97<-rbind(Cused,Dused)
alldata97<-alldata97[!is.na(alldata97$mean.ra),]
#alldata<-alldata97   #match time points in C and D, for p value distribution figure


library(reshape2)
meltdata<-melt(alldata,id.vars = c("mean.ra","reactor","treat","bestmodel","class"))
measurement<-sapply(strsplit(as.character(meltdata$variable),"[.]"), function(x) x[1])
models<-sapply(strsplit(as.character(meltdata$variable),"[.]"), function(x) x[2])
phylum<-sapply(strsplit(as.character(meltdata$class),"; "), function(x) x[2])

meltdata<-data.frame(meltdata,measurement,models,phylum)

### 1. plot percentage of the best models
rabins<-cut(meltdata$mean.ra,c(0, 0.0001, 0.001,1))
meltdata<-data.frame(meltdata,rabins)

#write.csv(meltdata,"CD100 meltdata for plot.csv")  ##data for p value plot (control the sample number)
write.csv(meltdata,"CD meltdata for plot.csv")


library(ggplot2)
meltdata$bestmodel=factor(meltdata$bestmodel,levels = c("neutral","consumer-resource","combined"))

#simple stack plot so that inkscape can open it
meltdataC<-meltdata[meltdata$treat=="C",]
meltdataD<-meltdata[meltdata$treat=="D",]

tableC<-table(meltdataC$rabins,meltdataC$bestmodel)
tableC<-tableC/rowSums(tableC)

tableD<-table(meltdataD$rabins,meltdataD$bestmodel)
tableD<-tableD/rowSums(tableD)

tableC<-as.data.frame(tableC)
tableD<-as.data.frame(tableD)

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



###2. r square plot
##r square values 
alldata<-alldata[,-11];
meltdata<-melt(alldata,id.vars = c("mean.ra","reactor","treat","bestmodel","modP.neutral","modP.resource","modP.combine"))

measurement<-sapply(strsplit(as.character(meltdata$variable),"[.]"), function(x) x[1])
models<-sapply(strsplit(as.character(meltdata$variable),"[.]"), function(x) x[2])
meltdata<-data.frame(meltdata,measurement,models)
rabins<-cut(meltdata$mean.ra,c(0, 0.0001, 0.001,1))
meltdata<-data.frame(meltdata,rabins)


meltdata$models<-factor(meltdata$models,levels = c("neutral","resource","combine"))

#using hist plot
ggplot(meltdata,aes(x =value,fill = models,linetype=treat,alpha=treat)) + 
  geom_density() +
  scale_fill_manual(values=c("#3C5488FF","#00A087FF","#E64B35FF"))+
  scale_color_manual(values=c("#3C5488FF","#00A087FF","#E64B35FF"))+
  scale_alpha_manual(values=c("C"=0.3,"D"=0.8))+
  scale_linetype_manual(values=c("C"="dashed","D"= "solid"))+
  xlim(-0.2,1.1)+
  #ylim(0,5)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())+
  facet_grid( models  ~  rabins, scale="free_y")


#get the statistics
quantile(meltdata$value[meltdata$models=="neutral"], probs = c(0.25, 0.5,0.75))
#25%       50%       75% 
#0.1878248 0.2887918 0.4117356 
quantile(meltdata$value[meltdata$models=="resource"], probs = c(0.25, 0.5,0.75))
#25%         50%         75% 
#0.002543855 0.013331008 0.044597147 
quantile(meltdata$value[meltdata$models=="combine"], probs = c(0.25, 0.5,0.75))
#25%       50%       75% 
#0.2273200 0.3268498 0.4600736 

mean(meltdata$value[meltdata$models=="neutral"]);
sd(meltdata$value[meltdata$models=="neutral"]);
#[1] 0.3093077
#[1] 0.1867827
mean(meltdata$value[meltdata$models=="resource"]);
sd(meltdata$value[meltdata$models=="resource"]);
#[1] 0.04183539
#[1] 0.07931047
mean(meltdata$value[meltdata$models=="combine"]);
sd(meltdata$value[meltdata$models=="combine"]);
#[1] 0.3555939
#[1] 0.1968022

##combined model, control vs treatment 
#rare
wilcox.test(meltdata$value[meltdata$models=="combine" & meltdata$treat=="C" & meltdata$rabins=="(0,0.0001]"], meltdata$value[meltdata$models=="combine" & meltdata$treat=="D" & meltdata$rabins=="(0,0.0001]"])
t.test(meltdata$value[meltdata$models=="combine" & meltdata$treat=="C" & meltdata$rabins=="(0,0.0001]"], meltdata$value[meltdata$models=="combine" & meltdata$treat=="D" & meltdata$rabins=="(0,0.0001]"])
#W = 1821113, p-value = 0.0004726
#t = -3.7552, df = 1415.8, p-value = 0.0001802
#mean of x mean of y 
#0.4033995 0.4308452 

#moderate
wilcox.test(meltdata$value[meltdata$models=="combine" & meltdata$treat=="C" & meltdata$rabins=="(0.0001,0.001]"], meltdata$value[meltdata$models=="combine" & meltdata$treat=="D" & meltdata$rabins=="(0.0001,0.001]"])
t.test(meltdata$value[meltdata$models=="combine" & meltdata$treat=="C" & meltdata$rabins=="(0.0001,0.001]"], meltdata$value[meltdata$models=="combine" & meltdata$treat=="D" & meltdata$rabins=="(0.0001,0.001]"])
#W = 589967, p-value = 2.868e-13
#t = -7.5098, df = 1109.5, p-value = 1.214e-13
# mean of x mean of y 
#0.2672344 0.3165298


#abundant
wilcox.test(meltdata$value[meltdata$models=="combine" & meltdata$treat=="C" & meltdata$rabins=="(0.001,1]"], meltdata$value[meltdata$models=="combine" & meltdata$treat=="D" & meltdata$rabins=="(0.001,1]"])
t.test(meltdata$value[meltdata$models=="combine" & meltdata$treat=="C" & meltdata$rabins=="(0.001,1]"], meltdata$value[meltdata$models=="combine" & meltdata$treat=="D" & meltdata$rabins=="(0.001,1]"])
# W = 25535, p-value = 0.1908
# t = -2.2867, df = 306.54, p-value = 0.02289
#mean of x mean of y 
#0.2347288 0.2617985

##neutral model, control vs treatment 
wilcox.test(meltdata$value[meltdata$models=="neutral" & meltdata$treat=="C" & meltdata$rabins=="(0,0.0001]"], meltdata$value[meltdata$models=="neutral" & meltdata$treat=="D" & meltdata$rabins=="(0,0.0001]"])
t.test(meltdata$value[meltdata$models=="neutral" & meltdata$treat=="C" & meltdata$rabins=="(0,0.0001]"], meltdata$value[meltdata$models=="neutral" & meltdata$treat=="D" & meltdata$rabins=="(0,0.0001]"])
# W = 2204105, p-value = 5.492e-09
# t = 5.2981, df = 1399.8, p-value = 1.358e-07
# mean of x mean of y 
# 0.3560671 0.3187186 


wilcox.test(meltdata$value[meltdata$models=="neutral" & meltdata$treat=="C" & meltdata$rabins=="(0.0001,0.001]"], meltdata$value[meltdata$models=="neutral" & meltdata$treat=="D" & meltdata$rabins=="(0.0001,0.001]"])
t.test(meltdata$value[meltdata$models=="neutral" & meltdata$treat=="C" & meltdata$rabins=="(0.0001,0.001]"], meltdata$value[meltdata$models=="neutral" & meltdata$treat=="D" & meltdata$rabins=="(0.0001,0.001]"])
# W = 727157, p-value = 0.8409
# t = -0.81416, df = 1071.4, p-value = 0.4157
# mean of x mean of y 
# 0.254026  0.259513

wilcox.test(meltdata$value[meltdata$models=="neutral" & meltdata$treat=="C" & meltdata$rabins=="(0.001,1]"], meltdata$value[meltdata$models=="neutral" & meltdata$treat=="D" & meltdata$rabins=="(0.001,1]"])
t.test(meltdata$value[meltdata$models=="neutral" & meltdata$treat=="C" & meltdata$rabins=="(0.001,1]"], meltdata$value[meltdata$models=="neutral" & meltdata$treat=="D" & meltdata$rabins=="(0.001,1]"])
#W = 38104, p-value = 2.241e-11
# t = 4.511, df = 283.04, p-value = 9.467e-06
#mean of x mean of y 
#0.2199217 0.1585845 

##consumer-resource model, control vs treatment 
wilcox.test(meltdata$value[meltdata$models=="resource" & meltdata$treat=="C" & meltdata$rabins=="(0,0.0001]"], meltdata$value[meltdata$models=="resource" & meltdata$treat=="D" & meltdata$rabins=="(0,0.0001]"])
t.test(meltdata$value[meltdata$models=="resource" & meltdata$treat=="C" & meltdata$rabins=="(0,0.0001]"], meltdata$value[meltdata$models=="resource" & meltdata$treat=="D" & meltdata$rabins=="(0,0.0001]"])
# W = 1490517, p-value < 2.2e-16
# t = -4.7906, df = 1462.9, p-value = 1.831e-06
# mean of x  mean of y 
#0.05370573 0.06928013 

wilcox.test(meltdata$value[meltdata$models=="resource" & meltdata$treat=="C" & meltdata$rabins=="(0.0001,0.001]"], meltdata$value[meltdata$models=="resource" & meltdata$treat=="D" & meltdata$rabins=="(0.0001,0.001]"])
t.test(meltdata$value[meltdata$models=="resource" & meltdata$treat=="C" & meltdata$rabins=="(0.0001,0.001]"], meltdata$value[meltdata$models=="resource" & meltdata$treat=="D" & meltdata$rabins=="(0.0001,0.001]"])
# W = 252274, p-value < 2.2e-16
# t = -15.138, df = 1014.1, p-value < 2.2e-16
#  mean of x  mean of y 
# 0.01324735 0.03830038 


wilcox.test(meltdata$value[meltdata$models=="resource" & meltdata$treat=="C" & meltdata$rabins=="(0.001,1]"], meltdata$value[meltdata$models=="resource" & meltdata$treat=="D" & meltdata$rabins=="(0.001,1]"])
t.test(meltdata$value[meltdata$models=="resource" & meltdata$treat=="C" & meltdata$rabins=="(0.001,1]"], meltdata$value[meltdata$models=="resource" & meltdata$treat=="D" & meltdata$rabins=="(0.001,1]"])
# W = 1911, p-value < 2.2e-16
# t = -18.945, df = 154.66, p-value < 2.2e-16
# mean of x   mean of y 
# 0.006114477 0.044497282


###3. Best model (AIC-based) across phylums
#simple stack plot so that inkscape can open it
meltdata<-read.csv("CD meltdata for plot.csv",row.names = 1)

meltdataD<-meltdata[meltdata$treat=="D",]
phylumcount<-table(meltdataD$phylum)
top10<-names(sort(phylumcount,decreasing = T))[1:11]
top10<-top10[-c(6,11)] #remove unclassified and archae group

##rare species
meltdataC<-meltdata[meltdata$treat=="C" & meltdata$rabins=="(0,0.0001]",]   ##rare species
meltdataD<-meltdata[meltdata$treat=="D" & meltdata$rabins=="(0,0.0001]",]

tableC<-table(meltdataC$phylum,meltdataC$bestmodel)
tableC<-tableC/rowSums(tableC)

tableD<-table(meltdataD$phylum,meltdataD$bestmodel)
tableD<-tableD/rowSums(tableD)

tableC<-as.data.frame(tableC)
tableD<-as.data.frame(tableD)

tableC<-data.frame(tableC,treat=rep("C",63))
tableD<-data.frame(tableD,treat=rep("D",66))
tableCD<-rbind(tableC,tableD)
colnames(tableCD)<-c("rabins","bestmodel","Freq","treat")

tableCDused<-tableCD[tableCD$rabins %in% top10,]
tableCDused$bestmodel=factor(tableCDused$bestmodel,levels=c("neutral","consumer-resource","combined"))
tableCDused$rabins=factor(tableCDused$rabins,levels=top10)

ggplot(tableCDused,aes(x = bestmodel, y = Freq,fill = bestmodel)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("#3C5488FF","#00A087FF","#E64B35FF"))+
  theme_bw()+
  facet_grid(treat ~ rabins)

##moderate species
meltdataC<-meltdata[meltdata$treat=="C" & meltdata$rabins=="(0.0001,0.001]",]   
meltdataD<-meltdata[meltdata$treat=="D" & meltdata$rabins=="(0.0001,0.001]",]

tableC<-table(meltdataC$phylum,meltdataC$bestmodel)
tableC<-tableC/rowSums(tableC)

tableD<-table(meltdataD$phylum,meltdataD$bestmodel)
tableD<-tableD/rowSums(tableD)

tableC<-as.data.frame(tableC)
tableD<-as.data.frame(tableD)

tableC<-data.frame(tableC,treat=rep("C",72))
tableD<-data.frame(tableD,treat=rep("D",60))
tableCD<-rbind(tableC,tableD)
colnames(tableCD)<-c("rabins","bestmodel","Freq","treat")

tableCDused<-tableCD[tableCD$rabins %in% top10,]
tableCDused$bestmodel=factor(tableCDused$bestmodel,levels=c("neutral","consumer-resource","combined"))
tableCDused$rabins=factor(tableCDused$rabins,levels=top10)

ggplot(tableCDused,aes(x = bestmodel, y = Freq,fill = bestmodel)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("#3C5488FF","#00A087FF","#E64B35FF"))+
  theme_bw()+
  facet_grid(treat ~ rabins)

##abundant species
meltdataC<-meltdata[meltdata$treat=="C" & meltdata$rabins=="(0.001,1]",]   
meltdataD<-meltdata[meltdata$treat=="D" & meltdata$rabins=="(0.001,1]",]

tableC<-table(meltdataC$phylum,meltdataC$bestmodel)
tableC<-tableC/rowSums(tableC)

tableD<-table(meltdataD$phylum,meltdataD$bestmodel)
tableD<-tableD/rowSums(tableD)

tableC<-as.data.frame(tableC)
tableD<-as.data.frame(tableD)

tableC<-data.frame(tableC,treat=rep("C",36))
tableD<-data.frame(tableD,treat=rep("D",30))
tableCD<-rbind(tableC,tableD)
colnames(tableCD)<-c("rabins","bestmodel","Freq","treat")

tableCDused<-tableCD[tableCD$rabins %in% top10,]
tableCDused$bestmodel=factor(tableCDused$bestmodel,levels=c("neutral","consumer-resource","combined"))
tableCDused$rabins=factor(tableCDused$rabins,levels=top10)

ggplot(tableCDused,aes(x = bestmodel, y = Freq,fill = bestmodel)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("#3C5488FF","#00A087FF","#E64B35FF"))+
  theme_bw()+
  facet_grid(treat ~ rabins)






