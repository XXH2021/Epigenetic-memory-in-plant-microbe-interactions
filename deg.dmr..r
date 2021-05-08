library("ggplot2")
library(devtools)
library(easyGgplot2)
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

dmr.pgp5.chh<-read.csv("CK-4-R_A5-2-R.CHH.dep4.win200.dmc.anno.csv",header=T)
dmr.pgp41.chh<-read.csv("CK-4-R_S41-5-R.CHH.dep4.win200.dmc.anno.csv",header=T)
deg.pgp5<-read.csv("PGP5-3_vs_CK-3.DEG.csv",header=T)
deg.pgp41<-read.csv("PGP41-3_vs_CK-3.DEG.csv",header=T)
colnames(dmr.pgp5.chh)[1]<-colnames(deg.pgp5)[1]
colnames(dmr.pgp41.chh)[1]<-colnames(deg.pgp41)[1]

pgp5.deg.dmr<-merge(deg.pgp5,dmr.pgp5.chh,by="GeneID",all=FALSE)
pgp41.deg.dmr<-merge(deg.pgp41,dmr.pgp41.chh,by="GeneID",all=FALSE)

#write.csv(pgp5.deg.dmr,"pgp5.deg.dmr.csv")
#write.csv(pgp41.deg.dmr,"pgp41.deg.dmr.csv")

for(i in 1:length(pgp5.deg.dmr$type.y)){
  if(pgp5.deg.dmr$type.y[i]=="hypo"){
    pgp5.deg.dmr$fc[i]<- 1/pgp5.deg.dmr$fc[i]
  }else{
  pgp5.deg.dmr$fc[i]<- pgp5.deg.dmr$fc[i]  
  }  
}
pgp5.deg.dmr$logFCm<-log2(pgp5.deg.dmr$fc)

for(i in 1:length(pgp41.deg.dmr$type.y)){
  if(pgp41.deg.dmr$type.y[i]=="hypo"){
    pgp41.deg.dmr$fc[i]<- 1/pgp41.deg.dmr$fc[i]
  }else{
    pgp41.deg.dmr$fc[i]<- pgp41.deg.dmr$fc[i]  
  }  
}
pgp41.deg.dmr$logFCm<-log2(pgp41.deg.dmr$fc)

#write.csv(pgp5.deg.dmr,"pgp5.deg.dmr.csv")
#write.csv(pgp41.deg.dmr,"pgp41.deg.dmr.csv")

pgp5.deg.dmr.p<-pgp5.deg.dmr[(pgp5.deg.dmr$type.x=="up" & pgp5.deg.dmr$type.y=="hyper") | 
                             (pgp5.deg.dmr$type.x=="down" & pgp5.deg.dmr$type.y=="hypo"),]
pgp5.deg.dmr.n<-pgp5.deg.dmr[(pgp5.deg.dmr$type.x=="up" & pgp5.deg.dmr$type.y=="hypo") | 
                               (pgp5.deg.dmr$type.x=="down" & pgp5.deg.dmr$type.y=="hyper"),]
pgp41.deg.dmr.p<-pgp41.deg.dmr[(pgp41.deg.dmr$type.x=="up" & pgp41.deg.dmr$type.y=="hyper") | 
                               (pgp41.deg.dmr$type.x=="down" & pgp41.deg.dmr$type.y=="hypo"),]
pgp41.deg.dmr.n<-pgp41.deg.dmr[(pgp41.deg.dmr$type.x=="up" & pgp41.deg.dmr$type.y=="hypo") | 
                               (pgp41.deg.dmr$type.x=="down" & pgp41.deg.dmr$type.y=="hyper"),]

#write.csv(pgp5.deg.dmr.p,"pgp5.deg.dmr.p.csv")
#write.csv(pgp5.deg.dmr.n,"pgp5.deg.dmr.n.csv")
#write.csv(pgp41.deg.dmr.p,"pgp41.deg.dmr.p.csv")
#write.csv(pgp41.deg.dmr.n,"pgp41.deg.dmr.n.csv")

pgp5.p<-read.csv("pgp5.deg.dmr.p.csv",header=T)
pgp5.n<-read.csv("pgp5.deg.dmr.n.csv",header=T)
pgp41.p<-read.csv("pgp41.deg.dmr.p.csv",header=T)
pgp41.n<-read.csv("pgp41.deg.dmr.n.csv",header=T)

colnames(pgp5.n)
pgp5.p<-na.omit(pgp5.p)
pgp5.n<-na.omit(pgp5.n)
pgp41.p<-na.omit(pgp41.p)
pgp41.n<-na.omit(pgp41.n)

r.pgp5.n<-cor(pgp5.n$logFC,pgp5.n$logFCm,method="pearson")
p.pgp5.n<-cor.test(pgp5.n$logFC,pgp5.n$logFCm,method="pearson")
p.pgp5.n<-p.pgp5.n[[3]]

r.pgp5.p<-cor(pgp5.p$logFC,pgp5.p$logFCm,method="pearson")
p.pgp5.p<-cor.test(pgp5.p$logFC,pgp5.p$logFCm,method="pearson")
p.pgp5.p<-p.pgp5.p[[3]]

r.pgp41.n<-cor(pgp41.n$logFC,pgp41.n$logFCm,method="pearson")
p.pgp41.n<-cor.test(pgp41.n$logFC,pgp41.n$logFCm,method="pearson")
p.pgp41.n<-p.pgp41.n[[3]]

r.pgp41.p<-cor(pgp41.p$logFC,pgp41.p$logFCm,method="pearson")
p.pgp41.p<-cor.test(pgp41.p$logFC,pgp41.p$logFCm,method="pearson")
p.pgp41.p<-p.pgp41.p[[3]]

p6<-ggplot(pgp5.n,aes(logFC,logFCm))+
  geom_point(stat="identity",size=0.7,color="red")+
  geom_smooth(method="lm",color="black")+
  geom_text(aes(x=2,y=4.5,label=paste('R','=',signif(r.pgp5.n,3),seq="")))+
  geom_text(aes(x=2,y=3.5,label=paste('P','=',signif(p.pgp5.n,3),seq="")))+
  theme_bw()+
  main_theme+
  labs(y = "")

p7<-ggplot(pgp5.p,aes(logFC,logFCm))+
  geom_point(stat="identity",size=0.7,color="blue")+
  geom_smooth(method="lm",color="black")+
  geom_text(aes(x=-2,y=6.5,label=paste('R','=',signif(r.pgp5.p,3),seq="")))+
  geom_text(aes(x=-2,y=5.5,label=paste('P','=',signif(p.pgp5.p,3),seq="")))+
  theme_bw()+
  main_theme+
  labs(y = "")

p8<-ggplot(pgp41.n,aes(logFC,logFCm))+
  geom_point(stat="identity",size=0.7,color="red")+
  geom_smooth(method="lm",color="black")+
  geom_text(aes(x=1,y=4,label=paste('R','=',signif(r.pgp41.n,3),seq="")))+
  geom_text(aes(x=1,y=3,label=paste('P','=',signif(p.pgp41.n,3),seq="")))+
  theme_bw()+
  main_theme+
  labs(y = "")

p9<-ggplot(pgp41.p,aes(logFC,logFCm))+
  geom_point(stat="identity",size=0.7,color="blue")+
  geom_smooth(method="lm",color="black")+
  geom_text(aes(x=-2,y=8,label=paste('R','=',signif(r.pgp41.p,3),seq="")))+
  geom_text(aes(x=-2,y=7,label=paste('P','=',signif(p.pgp41.p,3),seq="")))+
  theme_bw()+
  main_theme+
  labs(y = "")

pdf(file="deg.dmr.day3.pdf", width=7, height=5)
ggplot2.multiplot(p6, p7, p8, p9, cols=2)
dev.off()

dmr.pgp5.chh<-read.csv("CK-4-R_A5-2-R.CHH.dep4.win200.dmc.anno.csv",header=T)
dmr.pgp41.chh<-read.csv("CK-4-R_S41-5-R.CHH.dep4.win200.dmc.anno.csv",header=T)
gene.exp<-read.csv("gene.expression.mean.csv",header=T)
head(gene.exp)

overlape.dmr.chh<-merge(dmr.pgp5.chh,dmr.pgp41.chh,by="Geneid",all=FALSE)

pgp5.m<-read.csv("A5-2-R.gene.methy.csv",header=T,row.names=1)
ck.m<-read.csv("CK-4-R.gene.methy.csv",header=T,row.names=1)
pgp41.m<-read.csv("S41-5-R.gene.methy.csv",header=T,row.names=1)

ck.m$GeneID<-rownames(ck.m)
pgp5.m$GeneID<-rownames(pgp5.m)
pgp41.m$GeneID<-rownames(pgp41.m)

over.ck.pgp5<-merge(ck.m,pgp5.m,by="GeneID",all=FALSE)
over.ck.pgp41<-merge(ck.m,pgp41.m,by="GeneID",all=FALSE)

gene.exp$GeneID<-gene.exp$X
over.pgp5.expre<-merge(gene.exp,over.ck.pgp5,by="GeneID",all=FALSE)
over.pgp5.expre<-over.pgp5.expre[,-1]
colnames(over.pgp5.expre)[1]<-"Geneid"

over.pgp41.expre<-merge(gene.exp,over.ck.pgp41,by="GeneID",all=FALSE)
over.pgp41.expre<-over.pgp41.expre[,-1]
colnames(over.pgp41.expre)[1]<-"Geneid"

over.pgp5.expre$log2pgp5<-abs(log2(over.pgp5.expre$PGP5/over.pgp5.expre$CK3))
over.pgp41.expre$log2pgp41<-abs(log2(over.pgp41.expre$PGP413/over.pgp41.expre$CK3))
colnames(over.pgp41.expre)

dmr.pgp5.chh.exp<-merge(dmr.pgp5.chh,over.pgp5.expre,by="Geneid",all=FALSE)
dmr.pgp5.chh.exp.hypo<-dmr.pgp5.chh.exp[dmr.pgp5.chh.exp$type=="hypo",]
dmr.pgp5.chh.exp.hyper<-dmr.pgp5.chh.exp[dmr.pgp5.chh.exp$type=="hyper",]

dmr.pgp41.chh.exp<-merge(dmr.pgp41.chh,over.pgp41.expre,by="Geneid",all=FALSE)
dmr.pgp41.chh.exp.hypo<-dmr.pgp41.chh.exp[dmr.pgp41.chh.exp$type=="hypo",]
dmr.pgp41.chh.exp.hyper<-dmr.pgp41.chh.exp[dmr.pgp41.chh.exp$type=="hyper",]

wilcox.test(dmr.pgp5.chh.exp.hypo$log2pgp5,over.pgp5.expre$log2pgp5,exact=FALSE)#p-value = 0.6518
wilcox.test(dmr.pgp5.chh.exp.hyper$log2pgp5,over.pgp5.expre$log2pgp5,exact=FALSE)#p-value = 0.008441

wilcox.test(dmr.pgp41.chh.exp.hypo$log2pgp41,over.pgp41.expre$log2pgp41,exact=FALSE)#p-value = 0.02436
wilcox.test(dmr.pgp41.chh.exp.hyper$log2pgp41,over.pgp41.expre$log2pgp41,exact=FALSE)#p-value = 0.01684



