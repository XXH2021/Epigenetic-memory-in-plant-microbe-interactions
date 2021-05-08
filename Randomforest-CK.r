library(ggplot2)
main_theme = theme(panel.background=element_blank(),
                    panel.grid=element_blank(),
                    axis.line.x=element_line(size=.5, colour="black"),
                    axis.line.y=element_line(size=.5, colour="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(color="black", size=7),
                    legend.position="right",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    legend.text= element_text(size=7),
                    text=element_text(family="sans", size=7)
					)

tc_map =read.table("design.txt",header = T, row.names = 1)
otu_table = read.csv("otu_taxa_table.csv", row.names= 1,  header=T)
colnames(otu_table)<-gsub("\\.","",colnames(otu_table))
sub_map = tc_map[tc_map$Treatment %in% c("CK"),]
idx = rownames(sub_map) %in% colnames(otu_table)
sub_map = sub_map[idx,]
sub_otu = otu_table[, rownames(sub_map)]
sub_otu <- sub_otu[which(rowSums(sub_otu) > 30), ]
library(randomForest)
set.seed(315)
rf = randomForest(t(sub_otu), sub_map$day, importance=TRUE, proximity=TRUE, ntree = 1000)
print(rf)

set.seed(315) 
result = rfcv(Treatment~., data=sub_otu_pgp5, cv.fold=10)
result$error.cv
pdf(file="rfcv_ck.pdf", height = 4, width = 4)
with(result, plot(n.var, error.cv, log="x", type="l", lwd=1))
dev.off()

imp= as.data.frame(rf$importance)
imp = imp[order(imp[,1],decreasing = T),]
head(imp)
write.table(imp,file = "importance_class_ck.txt",quote = F,sep = '\t', row.names = T, col.names = T)
imp = read.table("importance_class_ck.txt", header=T, row.names= 1, sep="\t") 
imp = head(imp, n=22)
imp=imp[order(1:22,decreasing = T),]
imp$temp = rownames(imp)

name_otu=otu_table[rownames(otu_table) %in% rownames(imp),]
imp <- cbind(rownames(imp), imp)
name_otu<-cbind(rownames(name_otu), name_otu)
imp<-merge(imp,name_otu,by.x="rownames(imp)",by.y="rownames(name_otu)")
write.table(imp,file = "name_otu_ck.txt",quote = F,sep = '\t', row.names = T, col.names = T)
imp$taxonomy = gsub("^.*p__","",imp$taxonomy,perl=TRUE) 
imp$phylum = gsub(";.*","",imp$taxonomy,perl=TRUE) 
imp$temp=factor(imp$temp)
p=ggplot(data = imp, mapping = aes(x=reorder(imp$temp,imp$X.IncMSE,sum),y=X.IncMSE,fill=phylum)) + 
  geom_bar(stat="identity")+coord_flip()+main_theme
p
ggsave(paste("rf_imp_feature_ck_22",".pdf", sep=""), p, width = 4, height =3)

library(pheatmap)
rownames(imp)<-imp$temp
sub_abu = sub_otu[rownames(imp),]
pheatmap(sub_abu, scale = "row")

sampFile = as.data.frame(sub_map$day,row.names = row.names(sub_map))
colnames(sampFile)[1] = "group"
mat_t = t(sub_abu)
mat_t2 = merge(sampFile, mat_t, by="row.names")
mat_t2 = mat_t2[,-1]
mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
otu_norm_group = do.call(rbind, mat_mean)[-1,]
colnames(otu_norm_group) = mat_mean$group
rownames(otu_norm_group)=gsub("^.*\\;f","f",rownames(otu_norm_group),perl=TRUE) 
myheatmap<-pheatmap(otu_norm_group, scale = "row", cluster_cols = F, cluster_rows = T, border=FALSE,gaps_col=c(3),cutree_row = 2,filename = "heatmap_groups_ck_13.pdf", width = 3.15, height = 3)

row_order<-myheatmap$tree_row$order
otu.max.ck<-otu_norm_group[row_order,]
otu.max.ck.pgp5<-t(otu_table[rownames(otu.max.ck), 37:54])
mat.ck.pgp5 = cbind(sampFile, otu.max.ck.pgp5)
mat_mean.ck.pgp5 = aggregate(mat.ck.pgp5[,-1], by=mat.ck.pgp5[1], FUN=mean)
rownames(mat_mean.ck.pgp5) = mat_mean.ck.pgp5$group
pheatmap(t(mat_mean.ck.pgp5[,-1]), scale = "row", cluster_cols = F, cluster_rows = F, border=FALSE,gaps_col=c(3),gaps_row=c(1),filename = "heatmap_groups_ck_13_pgp5.pdf", width = 2.5, height = 3)

otu.max.ck.pgp41<-t(otu_table[rownames(otu.max.ck), 19:36])
mat.ck.pgp41 = cbind(sampFile, otu.max.ck.pgp41)
mat_mean.ck.pgp41 = aggregate(mat.ck.pgp41[,-1], by=mat.ck.pgp41[1], FUN=mean)
rownames(mat_mean.ck.pgp41) = mat_mean.ck.pgp41$group
pheatmap(t(mat_mean.ck.pgp41[,-1]), scale = "row", cluster_cols = F, cluster_rows = F, border=FALSE,gaps_col=c(3),gaps_row=c(1),filename = "heatmap_groups_ck_13_pgp41.pdf", width = 2.5, height = 3)
