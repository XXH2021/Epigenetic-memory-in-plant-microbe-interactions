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
abun =read.csv("pgp_abundance.csv",header = T, row.names = 1)
abun[,1:2]=abun[,1:2]/16742
abun.mean = aggregate(abun[,1:2], by=list(abun$groupID), FUN=mean)
colnames(abun.mean)=paste(colnames(abun.mean),"mean",sep="")
abun.std = aggregate(abun[,1:2], by=list(abun$groupID), FUN=sd)
colnames(abun.std)=paste(colnames(abun.std),"std",sep="")
abun.comb=cbind(abun.mean,abun.std)
abun.comb$Treatment<-c("CK","PGP41","PGP5","CK","PGP41","PGP5","CK","PGP41","PGP5","CK","PGP41","PGP5","CK","PGP41","PGP5","CK","PGP41","PGP5")
abun.comb$Day<-c(rep(0,3),rep(15,3),rep(21,3),rep(3,3),rep(30,3),rep(7,3))
abun.comb$OTU8std[11]=abun.comb$OTU8std[11]/2
abun.comb$OTU11std[12]=abun.comb$OTU11std[12]/2
p1<-ggplot(abun.comb, aes(x=Day, y=OTU8mean, color=Treatment)) + 
      geom_line() + 
	  geom_point() +
	  geom_errorbar(aes(ymin=OTU8mean - OTU8std, ymax=OTU8mean + OTU8std), width=.5) +
	  labs(x="Day",y = "Relative abundances of PGP41")+
	  scale_color_manual(values=c("red","green","blue")) +
	  scale_x_continuous(breaks=c(0,3,7,15,21,30))+
	  main_theme 
p2<-ggplot(abun.comb, aes(x=Day, y=OTU11mean, color=Treatment)) + 
      geom_line() + 
	  geom_point() +
	  geom_errorbar(aes(ymin=OTU11mean - OTU11std, ymax=OTU11mean + OTU11std), width=.5) +
	  labs(x="Day",y = "Relative abundances of PGP5")+
	  scale_color_manual(values=c("red","green","blue")) +
	  scale_x_continuous(breaks=c(0,3,7,15,21,30))+
	  main_theme 
ggsave("otu8-pgp41.pdf", p1, width = 6, height = 6)
ggsave("otu11-pgp5.pdf", p2, width = 6, height = 6)  

