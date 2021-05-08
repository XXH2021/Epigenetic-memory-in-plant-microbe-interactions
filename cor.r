rm(list=ls()) # clean enviroment object
library("corrplot")
library("pheatmap")
library(ggcorrplot)
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
                    text=element_text(family="sans", size=7))
otu_table = read.delim("otu2.txt", row.names= 1,  header=T, sep="\t")
design = read.delim("design.txt", row.names= 1,  header=T, sep="\t")
# Set group order
if ("TRUE" == "TRUE") {
    design$Group  = factor(design$Group, levels=c("d0_CK","d0_PGP41","d0_PGP5","d3_CK","d3_PGP41","d3_PGP5","d7_CK","d7_PGP41","d7_PGP5","d15_CK","d15_PGP41","d15_PGP5","d21_CK","d21_PGP41","d21_PGP5","d30_CK","d30_PGP41","d30_PGP5"))
    }else{design$Group  = as.factor(design$Group)}
print(paste("Number of group: ",length(unique(design$Group)),sep="")) # show group numbers
count = otu_table[,-1]
norm = count/rowSums(count,na=T) * 100 # normalization to total 100
sampFile = as.data.frame(design$Group,row.names = row.names(design))
colnames(sampFile)[1] = "Group"
mat = merge(sampFile, norm, by="row.names")
mat_mean = aggregate(mat[,-1], by=mat[2], FUN=mean)
mat_mean_final = mat_mean[,-2]
rownames(mat_mean_final)=mat_mean_final[,1]
mat_mean_final=t(mat_mean_final[,-1])
sim=cor(mat_mean_final,method="pearson")
sim=round(sim,3)
col1 <- colorRampPalette(c("green", "green", "red"))
pdf(file="corplot_pie_pearson.pdf", height = 12, width = 12)
corrplot(sim, method="circle", type="lower", col=col1(100),addrect = 4) # , diag=F , na.label = "1"
dev.off()
time1 = c(0,0,0,3,3,3,7,7,7,15,,15,15,21,21,21,30,30,30)
time2 = c(0,3,7,15,21,30)
time=data.frame(time1,time2)
pheatmap(time, cluster_rows = F,  cluster_cols = F)
pheatmap(time, cluster_rows = F,  cluster_cols = F, filename = "corplot_pie_legend_time.pdf" ,width=2, height=4)