library(ggplot2)
library(devtools)
library(easyGgplot2)
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
                    text=element_text(family="sans", size=14)
					)
alpha_div =read.table("alpha_div.txt",header = T, row.names = 1)
div_mean = aggregate(alpha_div[,1:6], by=list(alpha_div$groupID), FUN=mean)
colnames(div_mean)=paste(colnames(div_mean),"mean",sep="")
div_std = aggregate(alpha_div[,1:6], by=list(alpha_div$groupID), FUN=sd)
colnames(div_std)=paste(colnames(div_std),"std",sep="")
div_comb=cbind(div_mean,div_std)
div_comb$Treatment<-c("CK","PGP41","PGP5","CK","PGP41","PGP5","CK","PGP41","PGP5","CK","PGP41","PGP5","CK","PGP41","PGP5","CK","PGP41","PGP5")
div_comb$Day<-c(rep(0,3),rep(15,3),rep(21,3),rep(3,3),rep(30,3),rep(7,3))
p1<-ggplot(div_comb, aes(x=Day, y=OTUmean, color=Treatment)) + 
      geom_line() + 
	  geom_point() +
	  geom_errorbar(aes(ymin=OTUmean - OTUstd, ymax=OTUmean + OTUstd), width=.5) +
	  labs(y = "OTU")
p2<-ggplot(div_comb, aes(x=Day, y=chaomean, color=Treatment)) + 
      geom_line() + 
	  geom_point() +
	  geom_errorbar(aes(ymin=chaomean - chaostd, ymax=chaomean + chaostd), width=.5) +
	  scale_color_manual(values=c("red","green","blue")) +
	  labs(x="Days",y = "Chao1",colour="Treatments")+
	  scale_x_continuous(breaks=c(0,3,7,15,21,30))+
	  main_theme 
p3<-ggplot(div_comb, aes(x=Day, y=shannonmean, color=Treatment)) + 
      geom_line() + 
	  geom_point() +
	  geom_errorbar(aes(ymin=shannonmean - shannonstd , ymax=shannonmean + shannonstd ), width=.5) +
	  scale_color_manual(values=c("red","green","blue")) +
	  labs(x="Days",y = "Shannon",colour="Treatments") +
	  scale_x_continuous(breaks=c(0,3,7,15,21,30))+
	  main_theme 
p4<-ggplot(div_comb, aes(x=Day, y=simpsonmean, color=Treatment)) + 
      geom_line() + 
	  geom_point() +
	  geom_errorbar(aes(ymin=simpsonmean - simpsonstd , ymax=simpsonmean + simpsonstd ), width=.5) +
	  scale_color_manual(values=c("red","green","blue")) +
	  labs(x="Days",y = "Simpson",colour="Treatments") +
	  scale_x_continuous(breaks=c(0,3,7,15,21,30))+
	  main_theme  
pdf(file="alpha_dev.pdf", width=8, height=5)	  
ggplot2.multiplot(p1, p2, p3, p4, cols=2)
dev.off()
ggsave("Chao1.pdf", p2, width = 5, height = 4)
ggsave("Shannon.pdf", p3, width = 5, height = 4)
ggsave("Simpson.pdf", p4, width = 5, height = 4)  