

setwd("/Users/linwei/Dropbox/bioreactor/ratio model-new/1.2_competition and seq diff")

#read model fitting result
D97<-read.csv("/Users/linwei/Dropbox/bioreactor/ratio model-new/1_model fitting/D_before100 ra fitting.csv",header = T,row.names = 1)
C97<-read.csv("/Users/linwei/Dropbox/bioreactor/ratio model-new/1_model fitting/C_before100 ra fitting.csv",header = T,row.names = 1)
midasclassifier<-read.csv("/Users/linwei/Dropbox/bioreactor/ratio model-new/0_raw data/otus_resampled_tax_assignments_Midas.txt",sep="\t",header = F,row.names = 1)

###sequence similarity
seqdis<-read.csv("/Users/linwei/Dropbox/bioreactor/ratio model-new/0_raw data/Rep_clean_sequence.distmat",sep = "",header = F,row.names = 1,skip = 1)
colnames(seqdis)<-row.names(seqdis)
seqdis<-100-seqdis  ###convert similarity to dissimilarity


######D reactor
selected<-D97[D97$modP.combine <= 0.05 & (!is.na(D97$modP.combine)),]
plot(selected$x3.combine, selected$x3.resource)

seqdis.selected<-seqdis[match(row.names(selected),row.names(seqdis)), match(row.names(selected),colnames(seqdis))]

library(vegan)
b.dis<-vegdist(selected$x3.combine, "euclidean")
b.dis<-as.matrix(b.dis); row.names(b.dis)<-colnames(b.dis)<-row.names(selected)

mean(b.dis); sd(b.dis)
#[1] 0.1593356
#[1] 0.1366137
##transfer b distance matrix to list##
xy <- t(combn(colnames(b.dis), 2))
bdislist<-data.frame(xy, dist=b.dis[xy]);
min(bdislist[,3]);max(bdislist[,3])

##transfer sequence dissimilarity to list##
xy2 <- t(combn(colnames(seqdis.selected), 2));sum(xy[,1]==xy2[,1]);sum(xy[,2]==xy2[,2]);
seqdislist<-data.frame(xy2, dist=seqdis.selected[xy2])
min(seqdislist[,3]);max(seqdislist[,3])

compareD<-data.frame(seqdislist,bdislist,otu1.class=midasclassifier[match(seqdislist$X1,row.names(midasclassifier)),1],
                     otu2.class=midasclassifier[match(seqdislist$X2,row.names(midasclassifier)),1])

colnames(compareD)<-c("otu1","otu2","seqdis","otu1","otu2","bdis","otu1.class","otu2.class")
# plot(compareD$seqdis,compareD$bdis)

#####C97 reactor##

selected<-C97[C97$modP.combine <= 0.05 & (!is.na(C97$modP.combine)),]
plot(selected$x3.combine, selected$x3.resource)

seqdis.selected<-seqdis[match(row.names(selected),row.names(seqdis)), match(row.names(selected),colnames(seqdis))]


b.dis<-vegdist(selected$x3.combine, "euclidean")
b.dis<-as.matrix(b.dis); row.names(b.dis)<-colnames(b.dis)<-row.names(selected)
mean(b.dis); sd(b.dis)
##transfer b distance matrix to list##
xy <- t(combn(colnames(b.dis), 2))
bdislist<-data.frame(xy, dist=b.dis[xy]);
min(bdislist[,3]);max(bdislist[,3])

##transfer sequence dissimilarity to list##
xy2 <- t(combn(colnames(seqdis.selected), 2));sum(xy[,1]==xy2[,1]);sum(xy[,2]==xy2[,2]);
seqdislist<-data.frame(xy2, dist=seqdis.selected[xy2])
min(seqdislist[,3]);max(seqdislist[,3])

compareC<-data.frame(seqdislist,bdislist,otu1.class=midasclassifier[match(seqdislist$X1,row.names(midasclassifier)),1],
                     otu2.class=midasclassifier[match(seqdislist$X2,row.names(midasclassifier)),1])

colnames(compareC)<-c("otu1","otu2","seqdis","otu1","otu2","bdis","otu1.class","otu2.class")
#plot(compareC$seqdis,compareC$bdis)

###ggplot
compareC<-data.frame(compareC, treat=rep("C",nrow(compareC)))
compareD<-data.frame(compareD, treat=rep("D",nrow(compareD)))

compareCD<-rbind(compareC,compareD)

library(ggplot2)


compareCD$bdis
ggplot(compareCD, aes(x=seqdis, y=bdis,linetype=treat,colour=treat,fill=treat))+
  geom_smooth()+ 
  scale_color_manual(values=c("#999999", "#E69F00"))+
  scale_fill_manual(values=c("#999999", "#E69F00"))+
  scale_linetype_manual(values=c("C"="dashed","D"= "solid"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot(compareCD[compareCD$seqdis<53,], aes(x=seqdis, y=bdis,linetype=treat,colour=treat,fill=treat))+
  geom_smooth()+ 
  scale_color_manual(values=c("#999999", "#E69F00"))+
  scale_fill_manual(values=c("#999999", "#E69F00"))+
  scale_linetype_manual(values=c("C"="dashed","D"= "solid"))+
  scale_x_continuous(breaks=c(0,10,20,30,40,50))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


ggplot(compareCD[compareCD$seqdis<30,], aes(x=seqdis, y=bdis,linetype=treat,colour=treat,fill=treat))+
  geom_smooth()+ 
  scale_color_manual(values=c("#999999", "#E69F00"))+
  scale_fill_manual(values=c("#999999", "#E69F00"))+
  scale_linetype_manual(values=c("C"="dashed","D"= "solid"))+
  scale_x_continuous(breaks=c(0,10,20,30,40,50))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

