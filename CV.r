cv<-read.csv("kegg.profile.pathway.csv",header=T,row.names=1)

cv <- cv[which(rowSums(cv[,1:18]) > 18), ]
cv[cv<1]<-NA

calcv<-function(x){sd(x)/mean(x)*100}

cv_ck3<-apply(cv[,7:9],1,calcv)
cv_pgp413<-apply(cv[,10:12],1,calcv)
cv_pgp53<-apply(cv[,16:18],1,calcv)
cv_ck30<-apply(cv[,1:3],1,calcv)
cv_pgp4130<-apply(cv[,4:6],1,calcv)
cv_pgp530<-apply(cv[,13:15],1,calcv)

cv_cv1<-cbind(cv_ck3,cv_pgp413,cv_pgp53,cv_ck30,cv_pgp4130,cv_pgp530)
boxplot(cv_cv1,names=c("CK-d3","PGP41-d3","PGP5-d3","CK-d30","PGP41-d30","PGP5-d30"),ylab="Coefficient of variation (%)", notch=T,col=c("steelblue","mediumturquoise","sandybrown","hotpink","mediumpurple","grey"),ylim=c(0,60),cex=0.5,cex.axis=0.75)

day3<-cbind(cv[,7:9],cv[,10:12],cv[,16:18])
day30<-cbind(cv[,7:9],cv[,4:6],cv[,13:15])
ck<-cbind(cv[,1:3],cv[,7:9])
pgp41<-cbind(cv[,4:6],cv[,10:12])
pgp5<-cbind(cv[,16:18],cv[,13:15])

cv_day3<-apply(day3,1,calcv)
cv_day30<-apply(day30,1,calcv)
cv_ck<-apply(ck,1,calcv)
cv_pgp41<-apply(pgp41,1,calcv)
cv_pgp5<-apply(pgp5,1,calcv)

cv_cv2<-cbind(cv_day3,cv_day30,cv_ck,cv_pgp41,cv_pgp5)
boxplot(cv_cv2,names=c("Day3","Day30","CK","PGP41","PGP5"),ylab="Coefficient of variation (%)", notch=T,col=c("steelblue","mediumturquoise","sandybrown","hotpink","mediumpurple"),ylim=c(0,30),cex=0.5,cex.axis=0.75,las=2)

kegg<-read.csv("kegg_mata.csv",header=T,row.names=1)
mergedcv<-merge(kegg,cv_cv2,by="row.names")
row.names(mergedcv)=mergedcv$keggid
write.csv(mergedcv,file="combined.csv")

cv_cv3<-cv_cv1[rownames(mergedcv),]
cv_cv4<-cv_cv2[rownames(mergedcv),]

pdf("combined_single3.pdf", width=2, height=4.3)
boxplot(cv_cv4[,1:2],names=c("Day3","Day30"),ylab="Coefficient of variation (%)", notch=T,col=c("steelblue","mediumturquoise"),ylim=c(0,30),cex=0.5,cex.axis=0.75,las=2)
dev.off()

pdf("sigle_combined_cv4.pdf",width=10, height=8)
layout(matrix(c(1:12),2,6))
boxplot(cv_cv4[,1:2],names=c("Day3","Day30"), notch=T,main="All",cex.main=1,col=c("steelblue","mediumturquoise"),ylim=c(0,30),cex=0.5,cex.axis=0.75,las=2)
for (i in (2:13)){
	cv_cv<-subset(mergedcv,mergedcv$group==i)
	data=matrix(as.numeric(as.matrix(cv_cv[,5:6])),nrow=nrow(cv_cv))
boxplot(data, names=c("Day3","Day30"), notch=T,main=as.character(cv_cv[1,2]),cex.main=1, col=c("steelblue","mediumturquoise"),las=2,ylim=c(0,30),cex=0.5,cex.axis=0.75)
}
dev.off()










