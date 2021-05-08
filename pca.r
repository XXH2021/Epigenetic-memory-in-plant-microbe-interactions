library(ggplot2)
library(easyGgplot2)
main_theme = theme(panel.background=element_blank(),
                    panel.grid=element_blank(),
                    axis.line.x=element_line(size=.5, colour="black"),
                    axis.line.y=element_line(size=.5, colour="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(color="black", size=7),
                    legend.position="right",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
					)

kegg_table = read.csv("eggNOG.profile.cog.csv", row.names= 1,  header=T)
kegg_table <-t(kegg_table[,-dim(kegg_table)[2]])

Day=c("30","30","30","30","30","30","3","3","3","3","3","3","30","30","30","3","3","3")
Treatment=c(rep("CK",3),rep("PGP41",3),rep("CK",3),rep("PGP41",3),rep("PGP5",6))
kegg_table<-kegg_table[,colSums(kegg_table)>0]
pca.ck<-prcomp(kegg_table,scale=T,center=T)
pca.ck$x
points.ck<-as.data.frame(pca.ck$x)

write.csv(points.ck, "points.ck.csv")

p.ck = ggplot(points.ck, aes(x=PC1, y=PC2,color=Day,shape = Treatment)) + 
  geom_point(alpha=.7, size=2) +
  labs(x=paste("PCA 1 (", summary(pca.ck)$importance[2,1]*100, "%)", sep=""),
       y=paste("PCA 2 (", summary(pca.ck)$importance[2,2]*100, "%)", sep="")
	   )+
	scale_color_manual(values=c("red","blue"))

p<-p.ck + main_theme	
ggsave("kegg_PCA.pdf",p, width = 4, height = 4)

	
