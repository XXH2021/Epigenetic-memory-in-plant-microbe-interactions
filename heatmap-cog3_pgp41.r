library("ggplot2")
library(pheatmap)
data = read.csv("diff_cog_CK_PGP413_anno.csv", header=TRUE, row.names = 2)
data<-t(data[,8:25])
treat<-c(rep("CK-30",3),rep("PGP41-30",3),rep("CK-3",3),rep("PGP41-3",3),rep("PGP5-30",3),rep("PGP5-3",3))
na<-rownames(data)
tr<-data.frame(na,treat)
count.data<-data.frame(tr,data)[,-1]
mat_mean = aggregate(count.data[,-1], by=list(treat=count.data$treat), FUN=mean) # mean
rownames(mat_mean)=mat_mean[,1]
mat_mean<-t(mat_mean[,-1])

fun<-read.table("NOG.funccat.txt",sep="\t")
fun<-fun[1:4873,]
colnames(fun)<-c("cog","cato")
fun$cog<-as.vector(fun$cog)
fun$cato<-as.vector(fun$cato)
rownames(fun)<-fun$cog
fun<-fun[rownames(mat_mean),]
write.csv(fun,"cog3.pgp41.fun.csv")

annotation_row = data.frame(COGClass=factor(fun$cato))
rownames(annotation_row) = rownames(fun)
fun.count<-table(annotation_row$COGClass)
write.csv(fun.count,"cog3_pgp41.csv")

cato.5<-levels(annotation_row$COGClass)[table(annotation_row$COGClass)<6]
cato.6<-levels(annotation_row$COGClass)[table(annotation_row$COGClass)>5]
annotation_row[which(annotation_row$COGClass %in% cato.5),1]<-"Others"
annotation_row$COGClass<-as.vector(annotation_row$COGClass)
annotation_row[which(is.na(annotation_row$COGClass)),1]<-"Others"

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999","#000000")
GeneClass<-matrix(ncol=(length(cato.6)+1),nrow=1)
colnames(GeneClass)<-c(cato.6,"Others")
GeneClass[1,]<-c(cbPalette[1:(length(cato.6))],"#999999")
da<-list()
for(i in 1:nrow(GeneClass)) da[[i]]<-GeneClass[i,]
annotation_colors<-list(COGClass=da[[1]]) 

pheatmap(mat_mean, scale = "row", annotation_row = annotation_row,annotation_colors=annotation_colors,border=FALSE,
    fontsize = 9,treeheight_row=25,treeheight_col =15,show_rownames=F,
	width = 4, height = 4,filename = "heatmap_diff_cog_3_pgp41_.pdf")

