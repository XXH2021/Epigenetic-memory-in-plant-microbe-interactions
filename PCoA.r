# clean enviroment object
rm(list=ls()) 
# load related packages
library("ggplot2") 
library("vegan")
# Set ggplot2 drawing parameter, such as axis line and text size, lengend and title size, and so on.
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
points = read.table("points.txt", header=T, row.names= 1, sep="\t") 
# plot PCo 1 and 2
p = ggplot(points, aes(x=PC1, y=PC2,color=Day,shape = Treatment))
p = p + geom_point(alpha=.7, size=2) +
  labs(x=paste("PCoA 1 (", 64.46, "%)", sep=""),
       y=paste("PCoA 2 (", 14.29, "%)", sep="")
	   ) + main_theme
q1= p + scale_color_manual(values=rainbow(6))
p = ggplot(points, aes(x=PC1, y=PC3, color=Day, shape = Treatment))
p = p + geom_point(alpha=.7, size=2) +
  labs(x=paste("PCoA 1 (", 64.46, "%)", sep=""),
       y=paste("PCoA 3 (", 3.36, "%)", sep="")
	   ) + main_theme
q2= p + scale_color_manual(values=rainbow(6))
p = ggplot(points, aes(x=PC2, y=PC3, color=Day, shape = Treatment))
p = p + geom_point(alpha=.7, size=2) +
  labs(x=paste("PCoA 2 (", 14.29, "%)", sep=""),
       y=paste("PCoA 3 (", 3.36, "%)", sep="")
	   ) + main_theme
q3= p + scale_color_manual(values=rainbow(6))

ggsave("beta_pcoa_day_bray_curtis1.pdf", q1, width = 4, height = 4)
ggsave("beta_pcoa_day_bray_curtis2.pdf", q2, width = 4, height = 4)
ggsave("beta_pcoa_day_bray_curtis3.pdf", q3, width = 4, height = 4)

