library("ggplot2") # load related packages
library("grid")
library("scales")
library("vegan")
library("agricolae")
library("ggbiplot")
library("dplyr")
library("ggrepel")
main_theme = theme(panel.background=element_blank(),
                    panel.grid=element_blank(),
                    axis.line.x=element_line(size=.5, colour="black"),
                    axis.line.y=element_line(size=.5, colour="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(color="black", size=14),
                    legend.position="right",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    legend.text= element_text(size=14),
                    text=element_text(family="sans", size=14))
otu_table = read.delim("otu.txt", row.names= 1,  header=T, sep="\t")
design = read.delim("design.txt", row.names= 1,  header=T, sep="\t")
if (TRUE){
    design = subset(design,groupID %in% c("d0_CK","d0_PGP41","d0_PGP5","d3_CK","d3_PGP41","d3_PGP5","d7_CK","d7_PGP41","d7_PGP5","d15_CK","d15_PGP41","d15_PGP5","d21_CK","d21_PGP41","d21_PGP5","d30_CK","d30_PGP41","d30_PGP5") ) # select group1
}
calc_dis = function(subset){
	idx = grepl(subset, design$groupID)
	sub_design = design[idx,]
	print(paste("Number of group: ",length(unique(sub_design$group)),sep="")) # show group numbers

	idx=rownames(sub_design) %in% colnames(otu_table)
	sub_design=sub_design[idx,]
	count = otu_table[, rownames(sub_design)]
	norm = t(t(count)/colSums(count,na=T)) * 100 # normalization to total 100
	norm=as.data.frame(norm)

	final=rownames(sub_design[sub_design$day2 %in% 30,])
	ck=norm[,final] 
	ck_mean= rowMeans(ck)
	norm$ck_mean=ck_mean 

	bray_curtis = vegdist(t(norm), method = "bray")
	bray_curtis= as.matrix(bray_curtis)

	dat = t(as.data.frame(c("sampleA","sampleB","0.15","group","genosite")))
	colnames(dat) = c("sampleA","sampleB","distance","group","type")
	rownames(dat) = c("test")

	for (i in sort(unique(sub_design$day2))){
		group = rownames(sub_design[sub_design$day2 %in% i,])
		for (m in 1:(length(group)-1)) {
			x = c(group[m],"ck_mean",bray_curtis[group[m],"ck_mean"],i,subset)
			dat=rbind(dat,x)
		}
	}

	dat = as.data.frame(dat[-1,], stringsAsFactors=F) 
	dat$distance=round(as.numeric(as.character(dat$distance)), digits=3)
	dat$group = factor(dat$group, levels=unique(dat$group))
	return(dat)
}
dat = calc_dis("CK")
all = dat
for (i in c("PGP41", "PGP5")){
  dat = calc_dis(i)
  all = rbind(all, dat)
}
write.csv(all,file="distance.csv")
all<-read.csv("distance.csv",row.names= 1,  header=T)

library(ggpubr)
ck<-all[all$type=="CK",]
ck_mean = aggregate(ck[,3], by=list(ck$group), FUN=mean)
ck_sd = aggregate(ck[,3], by=list(ck$group), FUN=sd)

pgp5<-all[all$type=="PGP5",]
pgp5_mean = aggregate(pgp5[,3], by=list(pgp5$group), FUN=mean)
pgp5_sd = aggregate(pgp5[,3], by=list(pgp5$group), FUN=sd)

pgp41<-all[all$type=="PGP41",]
pgp41_mean = aggregate(pgp41[,3], by=list(pgp41$group), FUN=mean)
pgp41_sd = aggregate(pgp41[,3], by=list(pgp41$group), FUN=sd)

ck<-cbind(ck_mean,ck_sd)
pgp5<-cbind(pgp5_mean,pgp5_sd)
pgp41<-cbind(pgp41_mean,pgp41_sd)

div<-rbind(ck,pgp5,pgp41)
colnames(div)<-c("Days","means","Day2","sds")
div$Treatments<-c(rep("CK",6),rep("PGP5",6),rep("PGP41",6))
p4<-ggplot(div, aes(x=Days, y=means, color=Treatments)) + 
      geom_line() + 
	  geom_point() +
	  geom_errorbar(aes(ymin=means - sds , ymax=means + sds ), width=.5) +
	  scale_color_manual(values=c("red","green","blue")) +
	  labs(x="Days",y = "Bray-Curtis distances",colour="Treatments") +
	  scale_x_continuous(breaks=c(0,3,7,15,21,30))+
	  main_theme 
ggsave("Bray-Curtis distances.pdf", p4, width = 5, height = 4)