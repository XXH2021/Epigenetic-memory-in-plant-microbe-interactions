library("DESeq2")
data=read.csv("eggNOG.profile.cog.csv",header=T,row.names=1)
countData<-data[,1:18]
head(countData,2)

countData_PGP53<-cbind(countData[,7:9],countData[,16:18])
countData_PGP53<-floor(countData_PGP53)
condition <- factor(c(rep("CK3",3),rep("PGP53",3)), levels = c("CK3","PGP53"))
coldata <- data.frame(condition,row.names=colnames(countData_PGP53))
all(rownames(coldata)==colnames(countData_PGP53))

countData_PGP53<-floor(countData_PGP53)
ds<-DESeqDataSetFromMatrix(countData=countData_PGP53,colData=coldata,design = ~ condition)

keep <- rowSums(counts(ds)) >= 1
ds <- ds[keep,]
dds <- DESeq(ds)
res = results(dds, contrast=c("condition", "CK3", "PGP53"))
resOrdered <- res[order(res$pvalue),]
res = res[order(res$pvalue),]
table(res$padj<0.05)
write.csv(res,file="cog_CK_PGP53_results.csv")

countData_PGP530<-cbind(countData[,1:3],countData[,13:15])
countData_PGP530<-floor(countData_PGP530)
condition <- factor(c(rep("CK30",3),rep("PGP530",3)), levels = c("CK30","PGP530"))
coldata <- data.frame(condition,row.names=colnames(countData_PGP530))
all(rownames(coldata)==colnames(countData_PGP530))
countData_PGP530<-floor(countData_PGP530)
ds<-DESeqDataSetFromMatrix(countData=countData_PGP530,colData=coldata,design = ~ condition)
keep <- rowSums(counts(ds)) >= 1
ds <- ds[keep,]
dds <- DESeq(ds)
res = results(dds, contrast=c("condition", "CK30", "PGP530"))
resOrdered <- res[order(res$pvalue),]
res = res[order(res$pvalue),]
table(res$padj<0.05)
write.csv(res,file="cog_CK_PGP530_results.csv")

countData_PGP413<-cbind(countData[,7:9],countData[,10:12])
countData_PGP413<-floor(countData_PGP413)
condition <- factor(c(rep("CK3",3),rep("PGP413",3)), levels = c("CK3","PGP413"))
coldata <- data.frame(condition,row.names=colnames(countData_PGP413))
all(rownames(coldata)==colnames(countData_PGP413))
countData_PGP413<-floor(countData_PGP413)
ds<-DESeqDataSetFromMatrix(countData=countData_PGP413,colData=coldata,design = ~ condition)
keep <- rowSums(counts(ds)) >= 1
ds <- ds[keep,]
dds <- DESeq(ds)
res = results(dds, contrast=c("condition", "CK3", "PGP413"))
resOrdered <- res[order(res$pvalue),]
res = res[order(res$pvalue),]
table(res$padj<0.05)
write.csv(res,file="cog_CK_PGP413_results.csv")

countData_PGP4130<-cbind(countData[,1:3],countData[,4:6])
countData_PGP4130<-floor(countData_PGP4130)
condition <- factor(c(rep("CK30",3),rep("PGP4130",3)), levels = c("CK30","PGP4130"))
coldata <- data.frame(condition,row.names=colnames(countData_PGP4130))
all(rownames(coldata)==colnames(countData_PGP4130))
countData_PGP4130<-floor(countData_PGP4130)
ds<-DESeqDataSetFromMatrix(countData=countData_PGP4130,colData=coldata,design = ~ condition)
keep <- rowSums(counts(ds)) >= 1
ds <- ds[keep,]
dds <- DESeq(ds)
res = results(dds, contrast=c("condition", "CK30", "PGP4130"))
resOrdered <- res[order(res$pvalue),]
res = res[order(res$pvalue),]
table(res$padj<0.05)
write.csv(res,file="cog_CK_PGP4130_results.csv")

countData_PGP4130<-cbind(countData[,1:6],countData[,13:15],countData[,7:12],countData[,16:18])
countData_PGP4130<-floor(countData_PGP4130)
condition <- factor(c(rep("CK30",9),rep("PGP4130",9)), levels = c("CK30","PGP4130"))
coldata <- data.frame(condition,row.names=colnames(countData_PGP4130))
all(rownames(coldata)==colnames(countData_PGP4130))
countData_PGP4130<-floor(countData_PGP4130)
ds<-DESeqDataSetFromMatrix(countData=countData_PGP4130,colData=coldata,design = ~ condition)
keep <- rowSums(counts(ds)) >= 1
ds <- ds[keep,]
dds <- DESeq(ds)
res = results(dds, contrast=c("condition", "CK30", "PGP4130"))
resOrdered <- res[order(res$pvalue),]
res = res[order(res$pvalue),]
table(res$padj<0.05)
write.csv(res,file="cog_all_day3_day30.csv")

cazy.anno<-read.csv("eggNOG.profile.cog.csv",header=T)
colnames(cazy.anno)[1]<-"ID"
dif.cazy<-read.csv("cog_all_day3_day30.csv",header=T)
colnames(dif.cazy)[1]<-"ID"
cazy.a<-dif.cazy[which(dif.cazy$padj<0.05),]
cazy.b<-cazy.a[which(cazy.a[,3]<=-0.176 | cazy.a[,3]>=0.585),]
write.csv(cazy.b,file="dif_cog_all.csv")
dif.cazy.anno<-merge(cazy.b,cazy.anno,by="ID",all=FALSE)
write.csv(dif.cazy.anno,file="diff_cog_all_anno.csv")

cazy.anno<-read.csv("eggNOG.profile.cog.csv",header=T)
colnames(cazy.anno)[1]<-"ID"
dif.cazy<-read.csv("cog_CK_PGP53_results.csv",header=T)
colnames(dif.cazy)[1]<-"ID"
cazy.a<-dif.cazy[which(dif.cazy$padj<0.05),]
cazy.b<-cazy.a[which(cazy.a[,3]<=-0.176 | cazy.a[,3]>=0.585),]
write.csv(cazy.b,file="dif_cog_CK_PGP53.csv")
dif.cazy.anno<-merge(cazy.b,cazy.anno,by="ID",all=FALSE)
write.csv(dif.cazy.anno,file="diff_cog_CK_PGP53_anno.csv")

cazy.anno<-read.csv("eggNOG.profile.cog.csv",header=T)
colnames(cazy.anno)[1]<-"ID"
dif.cazy<-read.csv("cog_CK_PGP413_results.csv",header=T)
colnames(dif.cazy)[1]<-"ID"
cazy.a<-dif.cazy[which(dif.cazy$padj<0.05),]
cazy.b<-cazy.a[which(cazy.a[,3]<=-0.176 | cazy.a[,3]>=0.585),]
write.csv(cazy.b,file="dif_cog_CK_PGP413.csv")
dif.cazy.anno<-merge(cazy.b,cazy.anno,by="ID",all=FALSE)
write.csv(dif.cazy.anno,file="diff_cog_CK_PGP413_anno.csv")



