main_theme = theme(
  panel.grid=element_blank(),
  axis.line.x=element_line(size=0.5, colour="black"),
  axis.line.y=element_line(size=0.5, colour="black"),
  axis.ticks=element_line(color="black"),
  axis.text=element_text(color="black", size=12),
  legend.position="right",
  legend.background=element_blank(),
  legend.key=element_blank(),
  legend.text= element_text(size=12),
  text=element_text(family="sans", size=12))
                    
                  
dmr.pgp53.chh<-read.csv("CK-4-R_A5-2-R.CHH.dep4.win200.dmc.anno.csv",header=T)
dmr.pgp413.chh<-read.csv("CK-4-R_S41-5-R.CHH.dep4.win200.dmc.anno.csv",header=T)

dmr.pgp530.chh<-read.csv("CK-30_A5-30.CHH.dmr.anno.csv",header=T)
dmr.pgp4130.chh<-read.csv("CK-30_S41-30.CHH.dmr.anno.csv",header=T)

gene.exp<-read.csv("gene.expression.mean.csv",header=T)
head(gene.exp)
gene.exp$Geneid<-gene.exp$X
gene.exp$log2pgp53<-abs(log2(gene.exp$PGP5/gene.exp$CK3))
gene.exp$log2pgp413<-abs(log2(gene.exp$PGP413/gene.exp$CK3))
gene.exp$log2pgp530<-abs(log2(gene.exp$PGP530/gene.exp$CK30))
gene.exp$log2pgp4130<-abs(log2(gene.exp$PGP4130/gene.exp$CK30))

overlape.dmr.pgp5<-merge(dmr.pgp53.chh,dmr.pgp530.chh,by="Geneid",all=FALSE)
overlape.dmr.pgp41<-merge(dmr.pgp413.chh,dmr.pgp4130.chh,by="Geneid",all=FALSE)
nonoverlape.dmr.pgp53<-dmr.pgp53.chh[!(dmr.pgp53.chh$Geneid %in% overlape.dmr.pgp5$Geneid),]
nonoverlape.dmr.pgp530<-dmr.pgp530.chh[!(dmr.pgp530.chh$Geneid %in% overlape.dmr.pgp5$Geneid),]
nonoverlape.dmr.pgp413<-dmr.pgp413.chh[!(dmr.pgp413.chh$Geneid %in% overlape.dmr.pgp41$Geneid),]
nonoverlape.dmr.pgp4130<-dmr.pgp4130.chh[!(dmr.pgp4130.chh$Geneid %in% overlape.dmr.pgp41$Geneid),]

overlape.dmr.pgp5.exp<-merge(overlape.dmr.pgp5,gene.exp,by="Geneid",all=FALSE)
overlape.dmr.pgp41.exp<-merge(overlape.dmr.pgp41,gene.exp,by="Geneid",all=FALSE)
nonoverlape.dmr.pgp53.exp<-merge(nonoverlape.dmr.pgp53,gene.exp,by="Geneid",all=FALSE)
nonoverlape.dmr.pgp530.exp<-merge(nonoverlape.dmr.pgp530,gene.exp,by="Geneid",all=FALSE)
nonoverlape.dmr.pgp413.exp<-merge(nonoverlape.dmr.pgp413,gene.exp,by="Geneid",all=FALSE)
nonoverlape.dmr.pgp4130.exp<-merge(nonoverlape.dmr.pgp4130,gene.exp,by="Geneid",all=FALSE)


pgp53.m<-read.csv("A5-2-R.gene.methy.csv",header=T,row.names=1)
ck3.m<-read.csv("CK-4-R.gene.methy.csv",header=T,row.names=1)
pgp413.m<-read.csv("S41-5-R.gene.methy.csv",header=T,row.names=1)
ck3.m$Geneid<-rownames(ck3.m)
pgp53.m$Geneid<-rownames(pgp53.m)
pgp413.m$Geneid<-rownames(pgp413.m)
over.ck.pgp53<-merge(ck3.m,pgp53.m,by="Geneid",all=FALSE)
over.ck.pgp413<-merge(ck3.m,pgp413.m,by="Geneid",all=FALSE)

pgp530.m<-read.csv("A5-30.gene.methy.csv",header=T,row.names=1)
ck30.m<-read.csv("CK-30.gene.methy.csv",header=T,row.names=1)
pgp4130.m<-read.csv("S41-30.gene.methy.csv",header=T,row.names=1)
ck30.m$Geneid<-rownames(ck30.m)
pgp530.m$Geneid<-rownames(pgp530.m)
pgp4130.m$Geneid<-rownames(pgp4130.m)
over.ck.pgp530<-merge(ck30.m,pgp530.m,by="Geneid",all=FALSE)
over.ck.pgp4130<-merge(ck30.m,pgp4130.m,by="Geneid",all=FALSE)

over.ck.pgp53.exp<-merge(over.ck.pgp53,gene.exp,by="Geneid",all=FALSE)
over.ck.pgp413.exp<-merge(over.ck.pgp413,gene.exp,by="Geneid",all=FALSE)
over.ck.pgp530.exp<-merge(over.ck.pgp530,gene.exp,by="Geneid",all=FALSE)
over.ck.pgp4130.exp<-merge(over.ck.pgp4130,gene.exp,by="Geneid",all=FALSE)

wilcox.test(overlape.dmr.pgp5.exp$log2pgp53,over.ck.pgp53.exp$log2pgp53,exact=FALSE)#p-value = 0.4693
wilcox.test(nonoverlape.dmr.pgp53.exp$log2pgp53,over.ck.pgp53.exp$log2pgp53,exact=FALSE)#p-value = 0.02499
wilcox.test(overlape.dmr.pgp5.exp$log2pgp530,over.ck.pgp530.exp$log2pgp530,exact=FALSE)#p-value = 0.9949
wilcox.test(nonoverlape.dmr.pgp530.exp$log2pgp530,over.ck.pgp530.exp$log2pgp530,exact=FALSE)#p-value = 0.1051

overlape.dmr.pgp53.hypo.exp<-overlape.dmr.pgp5.exp[overlape.dmr.pgp5.exp$type.x=="hypo",]
overlape.dmr.pgp53.hyper.exp<-overlape.dmr.pgp5.exp[overlape.dmr.pgp5.exp$type.x=="hyper",]
overlape.dmr.pgp530.hypo.exp<-overlape.dmr.pgp5.exp[overlape.dmr.pgp5.exp$type.y=="hypo",]
overlape.dmr.pgp530.hyper.exp<-overlape.dmr.pgp5.exp[overlape.dmr.pgp5.exp$type.y=="hyper",]

wilcox.test(overlape.dmr.pgp53.hypo.exp$log2pgp53,over.ck.pgp53.exp$log2pgp53,exact=FALSE)#p-value = 0.6575
wilcox.test(overlape.dmr.pgp53.hyper.exp$log2pgp53,over.ck.pgp53.exp$log2pgp53,exact=FALSE)#p-value = 0.5654
wilcox.test(overlape.dmr.pgp530.hypo.exp$log2pgp530,over.ck.pgp530.exp$log2pgp530,exact=FALSE)#p-value = 0.6575
wilcox.test(overlape.dmr.pgp530.hyper.exp$log2pgp530,over.ck.pgp530.exp$log2pgp530,exact=FALSE)#p-value = 0.5654


wilcox.test(overlape.dmr.pgp5.exp$log2pgp53,nonoverlape.dmr.pgp53.exp$log2pgp53,exact=FALSE)#p-value = 0.4693
wilcox.test(overlape.dmr.pgp5.exp$log2pgp530,nonoverlape.dmr.pgp530.exp$log2pgp530,exact=FALSE)#p-value = 0.6128

#
dmr.pgp53.exp<-merge(dmr.pgp53.chh,gene.exp,by="Geneid",all=FALSE)
dmr.pgp413.exp<-merge(dmr.pgp413.chh,gene.exp,by="Geneid",all=FALSE)
dmr.pgp530.exp<-merge(dmr.pgp530.chh,gene.exp,by="Geneid",all=FALSE)
dmr.pgp4130.exp<-merge(dmr.pgp4130.chh,gene.exp,by="Geneid",all=FALSE)


dmr.pgp53.exp.hypo<-dmr.pgp53.exp[dmr.pgp53.exp$type=="hypo",]
dmr.pgp53.exp.hyper<-dmr.pgp53.exp[dmr.pgp53.exp$type=="hyper",]
dmr.pgp530.exp.hypo<-dmr.pgp530.exp[dmr.pgp530.exp$type=="hypo",]
dmr.pgp530.exp.hyper<-dmr.pgp530.exp[dmr.pgp530.exp$type=="hyper",]

dmr.pgp413.exp.hypo<-dmr.pgp413.exp[dmr.pgp413.exp$type=="hypo",]
dmr.pgp413.exp.hyper<-dmr.pgp413.exp[dmr.pgp413.exp$type=="hyper",]
dmr.pgp4130.exp.hypo<-dmr.pgp4130.exp[dmr.pgp4130.exp$type=="hypo",]
dmr.pgp4130.exp.hyper<-dmr.pgp4130.exp[dmr.pgp4130.exp$type=="hyper",]

wilcox.test(dmr.pgp53.exp.hypo$log2pgp53,over.ck.pgp53.exp$log2pgp53,exact=FALSE)#p-value = 0.6518
wilcox.test(dmr.pgp53.exp.hyper$log2pgp53,over.ck.pgp53.exp$log2pgp53,exact=FALSE)#p-value = 0.008441
wilcox.test(dmr.pgp530.exp.hypo$log2pgp530,over.ck.pgp530.exp$log2pgp530,exact=FALSE)#p-value = 0.04869
wilcox.test(dmr.pgp530.exp.hyper$log2pgp530,over.ck.pgp530.exp$log2pgp530,exact=FALSE)#p-value = 0.926

wilcox.test(dmr.pgp413.exp.hypo$log2pgp413,over.ck.pgp413.exp$log2pgp413,exact=FALSE)#p-value = 0.02436
wilcox.test(dmr.pgp413.exp.hyper$log2pgp413,over.ck.pgp413.exp$log2pgp413,exact=FALSE)#p-value = 0.01684
wilcox.test(dmr.pgp4130.exp.hypo$log2pgp4130,over.ck.pgp4130.exp$log2pgp4130,exact=FALSE)#p-value = 0.01202
wilcox.test(dmr.pgp4130.exp.hyper$log2pgp4130,over.ck.pgp4130.exp$log2pgp4130,exact=FALSE)#p-value = 0.2172

library("ggplot2")
library("reshape2")
library(easyGgplot2)
pgp53.chh<-rbind(dmr.pgp53.exp.hypo[,c("Geneid","log2pgp53")],
                dmr.pgp53.exp.hyper[,c("Geneid","log2pgp53")],
                over.ck.pgp53.exp[,c("Geneid","log2pgp53")])
pgp53.chh$type<-c(rep("hypo",dim(dmr.pgp53.exp.hypo)[1]),
                 rep("hyper",dim(dmr.pgp53.exp.hyper)[1]),
                 rep("all",dim(over.ck.pgp53.exp)[1]))

p1<-ggplot(pgp53.chh, aes(type,log2pgp53,fill=type))+
  geom_boxplot()+
  scale_y_continuous(limits=c(0,3.5))+
  main_theme+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
      panel.border=element_rect(fill='transparent', color="black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.background=element_blank(),
      legend.key=element_blank(),
      legend.position='none')+
  labs(y = "")
  
#ggsave("pgp5.chh.pdf", p1, width = 3, height = 3)

pgp530.chh<-rbind(dmr.pgp530.exp.hypo[,c("Geneid","log2pgp530")],
                dmr.pgp530.exp.hyper[,c("Geneid","log2pgp530")],
                over.ck.pgp530.exp[,c("Geneid","log2pgp530")])
pgp530.chh$type<-c(rep("hypo",dim(dmr.pgp530.exp.hypo)[1]),
                 rep("hyper",dim(dmr.pgp530.exp.hyper)[1]),
                 rep("all",dim(over.ck.pgp530.exp)[1]))

p2<-ggplot(pgp530.chh, aes(type,log2pgp530,fill=type))+
  geom_boxplot()+
  scale_y_continuous(limits=c(0,1.5))+
  main_theme+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.border=element_rect(fill='transparent', color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.position='none')+
  labs(y = "")

pgp413.chh<-rbind(dmr.pgp413.exp.hypo[,c("Geneid","log2pgp413")],
                  dmr.pgp413.exp.hyper[,c("Geneid","log2pgp413")],
                  over.ck.pgp413.exp[,c("Geneid","log2pgp413")])
pgp413.chh$type<-c(rep("hypo",dim(dmr.pgp413.exp.hypo)[1]),
                   rep("hyper",dim(dmr.pgp413.exp.hyper)[1]),
                   rep("all",dim(over.ck.pgp413.exp)[1]))

p3<-ggplot(pgp413.chh, aes(type,log2pgp413,fill=type))+
  geom_boxplot()+
  scale_y_continuous(limits=c(0,3.5))+
  main_theme+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.border=element_rect(fill='transparent', color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.position='none')+
  labs(y = "")
#ggsave("pgp41.chh.pdf", p1, width = 3, height = 3)

pgp4130.chh<-rbind(dmr.pgp4130.exp.hypo[,c("Geneid","log2pgp4130")],
                   dmr.pgp4130.exp.hyper[,c("Geneid","log2pgp4130")],
                   over.ck.pgp4130.exp[,c("Geneid","log2pgp4130")])
pgp4130.chh$type<-c(rep("hypo",dim(dmr.pgp4130.exp.hypo)[1]),
                    rep("hyper",dim(dmr.pgp4130.exp.hyper)[1]),
                    rep("all",dim(over.ck.pgp4130.exp)[1]))

p4<-ggplot(pgp4130.chh, aes(type,log2pgp4130,fill=type))+
  geom_boxplot()+
  scale_y_continuous(limits=c(0,1.5))+
  main_theme+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.border=element_rect(fill='transparent', color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.position='none')+
  labs(y = "")

pdf(file="deg.dmr.pdf", width=6, height=8)
ggplot2.multiplot(p1, p3, p2, p4, cols=2)
dev.off()

count<-data.frame(Hyper=c("612","506","680","577"),Hypo=c("398","498","1381","820"))
rownames(count)<-c("PGP5-DAI3","PGP41-DAI3","PGP5-DAI30","PGP41-DAI30")
count$Sample<-rownames(count)
count.p<-melt(count,id.vars=c("Sample"),variable.name = "Type",value.name = "Count")
count.p$Count<-as.numeric(count.p$Count)
count.p$Sample<-factor(count$Sample,levels<-c("PGP5-DAI3","PGP41-DAI3","PGP5-DAI30","PGP41-DAI30"))

p5<-ggplot(count.p,aes(Sample,Count))+
  geom_bar(stat="identity",position="dodge",aes(fill=factor(Type)))+
  main_theme+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.border=element_rect(fill='transparent', color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background=element_blank(),
        legend.key=element_blank())

count.c<-read.table("cgchgchh.txt",sep="\t",header=T)
count.c$Sample<-factor(count.c$Sample,levels<-c(rep(c("PGP53","PGP413"),1),rep(c("PGP530","PGP4130"),1)))

p6<-ggplot(count.c,aes(Sample,Count))+
  geom_bar(stat="identity",position="dodge",aes(fill=factor(Type)))+
  main_theme+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.border=element_rect(fill='transparent', color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background=element_blank(),
        legend.key=element_blank())

pdf(file="deg.dmr5.pdf", width=6, height=3)
ggplot2.multiplot(p6,cols=1)
dev.off()
      
pdf(file="deg.dmr2.pdf", width=6, height=3)
ggplot2.multiplot(p1, p3, p2, p4,cols=4)
dev.off()
pdf(file="deg.dmr3.pdf", width=4, height=4)
ggplot2.multiplot(p6,p8,p7,p9, cols=2)
dev.off()
pdf(file="deg.dmr4.pdf", width=5, height=3)
ggplot2.multiplot(p5,cols=1)
dev.off()

methy<-read.csv("methy.GO.csv",header=T)
main_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   axis.line.x=element_line(size=0.5, colour="black"),
                   axis.line.y=element_line(size=0.5, colour="black"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=10),
                   legend.position="right",
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=10),
                   text=element_text(family="sans", size=10))

p=ggplot(data = methy, mapping = aes(x=methy$description,y=methy$p.value,fill=Treatments)) + 
  geom_bar(stat="identity")+coord_flip()+main_theme+
  labs(y = "-log10(p-value)")+
  scale_fill_manual(values=c("grey","grey50"))
ggsave("methy.go.pdf", p, width = 8, height = 3.7)	

