library("ggplot2")
library("reshape2")
library(easyGgplot2)
library(scales)
main_theme = theme(
  panel.grid=element_blank(),
  axis.line.x=element_line(size=0.5, colour="black"),
  axis.line.y=element_line(size=0.5, colour="black"),
  axis.ticks=element_line(color="black"),
  axis.text=element_text(color="black", size=14),
  legend.position="right",
  legend.background=element_blank(),
  legend.key=element_blank(),
  legend.text= element_text(size=12),
  text=element_text(family="sans", size=12))

otu<-read.csv("otu_taxa_table.csv",head=T,row.names=1)
otu<-as.data.frame(t(otu[,-dim(otu)[2]]))
for (i in 1:length(rownames(otu))){
  otu$group[i]<-paste(unlist(strsplit(rownames(otu)[i],split='\\.'))[1:2],collapse=".")
}

otu.mean = aggregate(otu[,-dim(otu)[2]], by=list(otu$group), FUN=mean)
rownames(otu.mean)<-otu.mean[,1]
otu.mean<-t(otu.mean[,-1])
write.csv(otu.mean,"otu.mean.csv")
otu.mean<-read.csv("otu.mean.csv",header=T)
colnames(otu.mean)[1]="Name"

for(i in 1:dim(otu.mean)[1]){
otu.mean$Control.37[i]<-mean(c(otu.mean$Control.3[i] , otu.mean$Control.7[i]))
otu.mean$Control.2130[i]<-mean(c(otu.mean$Control.21[i] , otu.mean$Control.30[i]))
otu.mean$PGP5.37[i]<-mean(c(otu.mean$PGP5.3[i] , otu.mean$PGP5.7[i]))
otu.mean$PGP5.2130[i]<-mean(c(otu.mean$PGP5.21[i] , otu.mean$PGP5.30[i]))
otu.mean$PGP41.37[i]<-mean(c(otu.mean$PGP41.3[i] , otu.mean$PGP41.7[i]))
otu.mean$PGP41.2130[i]<-mean(c(otu.mean$PGP41.21[i] , otu.mean$PGP41.30[i]))
}

ck3<-read.table("otu_ck3_node.txt",header=T,sep="\t")
pgp53<-read.table("otu_pgp53_node.txt",header=T,sep="\t")
pgp413<-read.table("otu_pgp413_node.txt",header=T,sep="\t")
ck30<-read.table("otu_ck30_node.txt",header=T,sep="\t")
pgp530<-read.table("otu_pgp530_node.txt",header=T,sep="\t")
pgp4130<-read.table("otu_pgp4130_node.txt",header=T,sep="\t")

ck3$abun<-otu.mean[which(otu.mean$Name %in% ck3$Name), "Control.15"]
pgp53$abun<-otu.mean[which(otu.mean$Name %in% pgp53$Name), "PGP5.15"]
pgp413$abun<-otu.mean[which(otu.mean$Name %in% pgp413$Name), "PGP41.15"]
ck30$abun<-otu.mean[which(otu.mean$Name %in% ck30$Name), "Control.15"]
pgp530$abun<-otu.mean[which(otu.mean$Name %in% pgp530$Name), "PGP5.15"]
pgp4130$abun<-otu.mean[which(otu.mean$Name %in% pgp4130$Name), "PGP41.15"]

overlape.ck.pgp53<-merge(ck3,pgp53,by="Name",all=FALSE)
overlape.ck.pgp413<-merge(ck3,pgp413,by="Name",all=FALSE)
overlape.ck.pgp530<-merge(ck30,pgp530,by="Name",all=FALSE)
overlape.ck.pgp4130<-merge(ck30,pgp4130,by="Name",all=FALSE)

ck$abun<-otu.mean[which(otu.mean$Name %in% ck3$Name), "Control.15"]
pgp53$abun<-otu.mean[which(otu.mean$Name %in% pgp53$Name), "PGP5.15"]
pgp413$abun<-otu.mean[which(otu.mean$Name %in% pgp413$Name), "PGP41.15"]
ck30$abun<-otu.mean[which(otu.mean$Name %in% ck30$Name), "Control.15"]
pgp530$abun<-otu.mean[which(otu.mean$Name %in% pgp530$Name), "PGP5.15"]
pgp4130$abun<-otu.mean[which(otu.mean$Name %in% pgp4130$Name), "PGP41.15"]

wilcox.test(ck3$node.degree, pgp53$node.degree, exact=FALSE)
wilcox.test(ck3$node.betw, pgp53$node.betw, exact=FALSE)
wilcox.test(ck3$node.stress, pgp53$node.stress, exact=FALSE)
wilcox.test(ck3$node.evcent, pgp53$node.evcent, exact=FALSE)
wilcox.test(ck3$Clustering.Coefficient, pgp53$Clustering.Coefficient, exact=FALSE)
wilcox.test(ck3$abun, pgp53$abun, exact=FALSE)

wilcox.test(ck3$node.degree, pgp413$node.degree, exact=FALSE)
wilcox.test(ck3$node.betw, pgp413$node.betw, exact=FALSE)
wilcox.test(ck3$node.stress, pgp413$node.stress, exact=FALSE)
wilcox.test(ck3$node.evcent, pgp413$node.evcent, exact=FALSE)
wilcox.test(ck3$Clustering.Coefficient, pgp413$Clustering.Coefficient, exact=FALSE)
wilcox.test(ck3$abun, pgp413$abun, exact=FALSE)

wilcox.test(ck30$node.degree, pgp530$node.degree, exact=FALSE)
wilcox.test(ck30$node.betw, pgp530$node.betw, exact=FALSE)
wilcox.test(ck30$node.stress, pgp530$node.stress, exact=FALSE)
wilcox.test(ck30$node.evcent, pgp530$node.evcent, exact=FALSE)
wilcox.test(ck30$Clustering.Coefficient, pgp530$Clustering.Coefficient, exact=FALSE)
wilcox.test(ck30$abun, pgp530$abun, exact=FALSE)

wilcox.test(ck30$node.degree, pgp4130$node.degree, exact=FALSE)
wilcox.test(ck30$node.betw, pgp4130$node.betw, exact=FALSE)
wilcox.test(ck30$node.stress, pgp4130$node.stress, exact=FALSE)
wilcox.test(ck30$node.evcent, pgp4130$node.evcent, exact=FALSE)
wilcox.test(ck30$Clustering.Coefficient, pgp4130$Clustering.Coefficient, exact=FALSE)
wilcox.test(ck30$abun, pgp4130$abun, exact=FALSE)

degree.ck.pgp53<-rbind(ck3[,c("Name","node.degree")], pgp53[,c("Name","node.degree")])
degree.ck.pgp53$type<-c(rep("CK",dim(ck3)[1]),rep("PGP53",dim(pgp53)[1]))
stress.ck.pgp53<-rbind(ck3[,c("Name","node.stress")], pgp53[,c("Name","node.stress")])
stress.ck.pgp53$type<-c(rep("CK",dim(ck3)[1]),rep("PGP53",dim(pgp53)[1]))
evcent.ck.pgp53<-rbind(ck3[,c("Name","node.evcent")], pgp53[,c("Name","node.evcent")])
evcent.ck.pgp53$type<-c(rep("CK",dim(ck3)[1]),rep("PGP53",dim(pgp53)[1]))

stress.ck.pgp413<-rbind(ck3[,c("Name","node.stress")], pgp413[,c("Name","node.stress")])
stress.ck.pgp413$type<-c(rep("CK",dim(ck3)[1]),rep("PGP413",dim(pgp413)[1]))
Clustering.Coefficient.ck.pgp413<-rbind(ck3[,c("Name","Clustering.Coefficient")], pgp413[,c("Name","Clustering.Coefficient")])
Clustering.Coefficient.ck.pgp413$type<-c(rep("CK",dim(ck3)[1]),rep("PGP413",dim(pgp413)[1]))

p1<-ggplot(degree.ck.pgp53, aes(type,node.degree,fill=type))+
  geom_boxplot(outlier.size=0.2)+
  main_theme+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.border=element_rect(fill='transparent', color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.position='none')+
  labs(y = "")
p2<-ggplot(stress.ck.pgp53, aes(type,node.stress,fill=type))+
  geom_boxplot(outlier.size=0.2)+
  main_theme+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.border=element_rect(fill='transparent', color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.position='none')+
        labs(y = "")+
        scale_y_continuous(limits=c(0,1e+6),labels =scientific)
p3<-ggplot(evcent.ck.pgp53, aes(type,node.evcent,fill=type))+
  geom_boxplot(outlier.size=0.2)+
  main_theme+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.border=element_rect(fill='transparent', color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.position='none')+
        labs(y = "")

p4<-ggplot(stress.ck.pgp413, aes(type,node.stress,fill=type))+
  geom_boxplot(outlier.size=0.2)+
  main_theme+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.border=element_rect(fill='transparent', color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.position='none')+
  labs(y = "")+
  scale_y_continuous(labels =scientific)
p5<-ggplot(Clustering.Coefficient.ck.pgp413, aes(type,Clustering.Coefficient,fill=type))+
  geom_boxplot(outlier.size=0.2)+
  main_theme+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.border=element_rect(fill='transparent', color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.position='none')+
  labs(y = "")
  
ggsave("1.pdf", p1, width = 1.5, height = 2.5)
ggsave("2.pdf", p2, width = 2, height = 2.5)
ggsave("3.pdf", p3, width = 1.66, height = 2.5)
ggsave("4.pdf", p4, width = 2, height = 2.5)
ggsave("5.pdf", p5, width = 1.66, height = 2.5)
pdf(file="network.node.property.pdf", width=2, height=12)
ggplot2.multiplot(p1, p2, p3, p4,p5,cols=1)
dev.off()






#基于共有otu的检验和作图

wilcox.test(overlape.ck.pgp53$node.degree.x, overlape.ck.pgp53$node.degree.y, exact=FALSE)#p-value = 0.007166
wilcox.test(overlape.ck.pgp53$node.betw.x, overlape.ck.pgp53$node.betw.y, exact=FALSE)#p-value = 0.009354
wilcox.test(overlape.ck.pgp53$node.stress.x, overlape.ck.pgp53$node.stress.y, exact=FALSE)#p-value = 2.362e-08
wilcox.test(overlape.ck.pgp53$node.evcent.x, overlape.ck.pgp53$node.evcent.y, exact=FALSE)#p-value = 0.7934
wilcox.test(overlape.ck.pgp53$Clustering.Coefficient.x, overlape.ck.pgp53$Clustering.Coefficient.y, exact=FALSE)#p-value = 0.3041

wilcox.test(overlape.ck.pgp413$node.degree.x, overlape.ck.pgp413$node.degree.y, exact=FALSE)#p-value = 0.001137
wilcox.test(overlape.ck.pgp413$node.betw.x, overlape.ck.pgp413$node.betw.y, exact=FALSE)#p-value = 0.2648
wilcox.test(overlape.ck.pgp413$node.stress.x, overlape.ck.pgp413$node.stress.y, exact=FALSE)#p-value = 0.2295
wilcox.test(overlape.ck.pgp413$node.evcent.x, overlape.ck.pgp413$node.evcent.y, exact=FALSE)#p-value = 0.704
wilcox.test(overlape.ck.pgp413$Clustering.Coefficient.x, overlape.ck.pgp413$Clustering.Coefficient.y, exact=FALSE)#p-value = 2.319e-06

wilcox.test(overlape.ck.pgp530$node.degree.x, overlape.ck.pgp530$node.degree.y, exact=FALSE)#p-value = 0.08243
wilcox.test(overlape.ck.pgp530$node.betw.x, overlape.ck.pgp530$node.betw.y, exact=FALSE)#p-value = 0.8775
wilcox.test(overlape.ck.pgp530$node.stress.x, overlape.ck.pgp530$node.stress.y, exact=FALSE)#p-value = 0.3876
wilcox.test(overlape.ck.pgp530$node.evcent.x, overlape.ck.pgp530$node.evcent.y, exact=FALSE)#p-value = 4.827e-05
wilcox.test(overlape.ck.pgp530$Clustering.Coefficient.x, overlape.ck.pgp530$Clustering.Coefficient.y, exact=FALSE)#p-value = 0.2918

wilcox.test(overlape.ck.pgp4130$node.degree.x, overlape.ck.pgp4130$node.degree.y, exact=FALSE)#p-value = 0.1033
wilcox.test(overlape.ck.pgp4130$node.betw.x, overlape.ck.pgp4130$node.betw.y, exact=FALSE)#p-value = 0.0007026
wilcox.test(overlape.ck.pgp4130$node.stress.x, overlape.ck.pgp4130$node.stress.y, exact=FALSE)#p-value = 0.001754
wilcox.test(overlape.ck.pgp4130$node.evcent.x, overlape.ck.pgp4130$node.evcent.y, exact=FALSE)#p-value = 0.4901
wilcox.test(overlape.ck.pgp4130$Clustering.Coefficient.x, overlape.ck.pgp4130$Clustering.Coefficient.y, exact=FALSE)#p-value = 0.823


ck.pgp53.abun<-otu.mean[which(otu.mean$Name %in% overlape.ck.pgp53$Name), "Control.21"]
pgp53.abun<-otu.mean[which(otu.mean$Name %in% overlape.ck.pgp53$Name), "PGP5.21"]
pgp413.abun<-otu.mean[which(otu.mean$Name %in% overlape.ck.pgp413$Name), "PGP41.3"]
ck.pgp413.abun<-otu.mean[which(otu.mean$Name %in% overlape.ck.pgp413$Name), "Control.3"]

ck.pgp530.abun<-otu.mean[which(otu.mean$Name %in% overlape.ck.pgp530$Name), "Control.3"]
pgp530.abun<-otu.mean[which(otu.mean$Name %in% overlape.ck.pgp530$Name), "PGP5.3"]
pgp4130.abun<-otu.mean[which(otu.mean$Name %in% overlape.ck.pgp4130$Name), "PGP41.3"]
ck.pgp4130<-otu.mean[which(otu.mean$Name %in% overlape.ck.pgp4130$Name), "Control.3"]

wilcox.test(ck.pgp53.abun, pgp53.abun, exact=FALSE)#p-value = 0.3575
wilcox.test(ck.pgp413.abun, pgp413.abun, exact=FALSE)#p-value = 0.3575

#作图

pgp53.net<-overlape.ck.pgp53[,c("Name","node.degree.x", "node.degree.y")]
pgp53.net.degree<-melt(pgp53.net,id.vars=c("Name"),variable.name = "OTU",value.name = "node.degree")
pgp53.net.degree$type<-c(rep("CK",dim(overlape.ck.pgp53)[1]),rep("PGP53",dim(overlape.ck.pgp53)[1]))   

pgp53.net<-overlape.ck.pgp53[,c("Name","node.betw.x", "node.betw.y")]
pgp53.net.betw<-melt(pgp53.net,id.vars=c("Name"),variable.name = "OTU",value.name = "node.betw")
pgp53.net.betw$type<-c(rep("CK",dim(overlape.ck.pgp53)[1]),rep("PGP53",dim(overlape.ck.pgp53)[1])) 

pgp53.net<-overlape.ck.pgp53[,c("Name","node.stress.x", "node.stress.y")]
pgp53.net.stress<-melt(pgp53.net,id.vars=c("Name"),variable.name = "OTU",value.name = "node.stress")
pgp53.net.stress$type<-c(rep("CK",dim(overlape.ck.pgp53)[1]),rep("PGP53",dim(overlape.ck.pgp53)[1]))


pgp413.net<-overlape.ck.pgp413[,c("Name","node.degree.x", "node.degree.y")]
pgp413.net.degree<-melt(pgp413.net,id.vars=c("Name"),variable.name = "OTU",value.name = "node.degree")
pgp413.net.degree$type<-c(rep("CK",dim(overlape.ck.pgp413)[1]),rep("PGP413",dim(overlape.ck.pgp413)[1])) 

pgp413.net<-overlape.ck.pgp413[,c("Name","Clustering.Coefficient.x", "Clustering.Coefficient.y")]
pgp413.net.Clustering<-melt(pgp413.net,id.vars=c("Name"),variable.name = "OTU",value.name = "Clustering.Coefficient")
pgp413.net.Clustering$type<-c(rep("CK",dim(overlape.ck.pgp413)[1]),rep("PGP413",dim(overlape.ck.pgp413)[1])) 


p11<-ggplot(pgp53.net.degree, aes(type,node.degree,fill=type))+
  geom_boxplot()+
#  scale_y_continuous(limits=c(0,3.5))+
  main_theme+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.border=element_rect(fill='transparent', color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.position='none')+
        labs(y = "")

p12<-ggplot(pgp53.net.betw, aes(type,node.betw,fill=type))+
  geom_boxplot()+
    scale_y_continuous(limits=c(0,5000))+
  main_theme+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.border=element_rect(fill='transparent', color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.position='none')+
  labs(y = "")

p13<-ggplot(pgp53.net.stress, aes(type,node.stress,fill=type))+
  geom_boxplot()+
  scale_y_continuous(limits=c(0,50000))+
  main_theme+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.border=element_rect(fill='transparent', color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.position='none')+
        labs(y = "")

p14<-ggplot(pgp413.net.degree, aes(type,node.degree,fill=type))+
  geom_boxplot()+
  #  scale_y_continuous(limits=c(0,3.5))+
  main_theme+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.border=element_rect(fill='transparent', color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.position='none')+
  labs(y = "")

p15<-ggplot(pgp413.net.Clustering, aes(type,Clustering.Coefficient,fill=type))+
  geom_boxplot()+
  #  scale_y_continuous(limits=c(0,3.5))+
  main_theme+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.border=element_rect(fill='transparent', color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.position='none')+
  labs(y = "")
