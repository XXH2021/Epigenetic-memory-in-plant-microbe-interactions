deg.pgp53<-read.csv("PGP5-3_vs_CK-3.DEG.csv",header=T,row.names=1)
deg.pgp530<-read.csv("PGP5-30_vs_CK-30.DEG.csv",header=T,row.names=1)
deg.pgp413<-read.csv("PGP41-3_vs_CK-3.DEG.csv",header=T,row.names=1)
deg.pgp4130<-read.csv("PGP41-30_vs_CK-30.DEG.csv",header=T,row.names=1)

dmr.pgp5.c<-read.csv("CK-30_A5-30.C.dmr.anno.csv",header=T,row.names=1)
dmr.pgp5.cg<-read.csv("CK-30_A5-30.CG.dmr.anno.csv",header=T,row.names=1)
dmr.pgp5.chg<-read.csv("CK-30_A5-30.CHG.dmr.anno.csv",header=T,row.names=1)
dmr.pgp5.chh<-read.csv("CK-30_A5-30.CHH.dmr.anno.csv",header=T,row.names=1)

dmr.pgp41.c<-read.csv("CK-30_S41-30.C.dmr.anno.csv",header=T,row.names=1)
dmr.pgp41.cg<-read.csv("CK-30_S41-30.CG.dmr.anno.csv",header=T,row.names=1)
dmr.pgp41.chg<-read.csv("CK-30_S41-30.CHG.dmr.anno.csv",header=T,row.names=1)
dmr.pgp41.chh<-read.csv("CK-30_S41-30.CHH.dmr.anno.csv",header=T,row.names=1)

library(VennDiagram)

venn.diagram(list(deg.PGP5.3=rownames(deg.pgp53),deg.PGP5.30=rownames(deg.pgp530),deg.PGP41.3=rownames(deg.pgp413),deg.PGP41.30=rownames(deg.pgp4130)), 
    resolution = 300, imagetype = "tiff", alpha=c(0.5,0.5,0.5,0.5), cex = 2,cat.cex = 2,
    fill=c("red","yellow","blue","green"), cat.fontface=4,fontfamily=3,
    main="Shared DEG by different groups",
    main.cex = 2, main.fontface = 2, main.fontfamily = 3,
    filename = "deg.VennDiagram.tif")


venn.diagram(list(dmr.pgp5.c=rownames(dmr.pgp5.c),dmr.pgp5.cg=rownames(dmr.pgp5.cg),dmr.pgp5.chg=rownames(dmr.pgp5.chg),dmr.pgp5.chh=rownames(dmr.pgp5.chh)), 
    resolution = 300, imagetype = "tiff", alpha=c(0.5,0.5,0.5,0.5), cex = 2,cat.cex = 2,
    fill=c("red","yellow","blue","green"), cat.fontface=4,fontfamily=3,
    main="Shared DMR by different groups",
    main.cex = 2, main.fontface = 2, main.fontfamily = 3,
    filename = "dmr.pgp5.VennDiagram.tif")
	
venn.diagram(list(dmr.pgp5.cg=rownames(dmr.pgp5.cg),dmr.pgp5.chg=rownames(dmr.pgp5.chg),dmr.pgp5.chh=rownames(dmr.pgp5.chh)), 
    resolution = 300, imagetype = "tiff", alpha=c(0.5,0.5,0.5), cex = 2,cat.cex = 2,
    fill=c("red","yellow","blue"), cat.fontface=4,fontfamily=3,
    main="Shared DMR by different groups",
    main.cex = 2, main.fontface = 2, main.fontfamily = 3,
    filename = "dmr.pgp5.VennDiagram2.tif")

	
venn.diagram(list(dmr.pgp41.c=rownames(dmr.pgp41.c),dmr.pgp41.cg=rownames(dmr.pgp41.cg),dmr.pgp41.chg=rownames(dmr.pgp41.chg),dmr.pgp41.chh=rownames(dmr.pgp41.chh)), 
    resolution = 300, imagetype = "tiff", alpha=c(0.5,0.5,0.5,0.5), cex = 2,cat.cex = 2,
    fill=c("red","yellow","blue","green"), cat.fontface=4,fontfamily=3,
    main="Shared DMR by different groups",
    main.cex = 2, main.fontface = 2, main.fontfamily = 3,
    filename = "deg.pgp41.VennDiagram.tif")
	
venn.diagram(list(dmr.pgp41.cg=rownames(dmr.pgp41.cg),dmr.pgp41.chg=rownames(dmr.pgp41.chg),dmr.pgp41.chh=rownames(dmr.pgp41.chh)), 
    resolution = 300, imagetype = "tiff", alpha=c(0.5,0.5,0.5), cex = 2,cat.cex = 2,
    fill=c("red","yellow","blue"), cat.fontface=4,fontfamily=3,
    main="Shared DMR by different groups",
    main.cex = 2, main.fontface = 2, main.fontfamily = 3,
    filename = "deg.pgp41.VennDiagram2.tif")
	

venn.diagram(list(deg.PGP5.3=rownames(deg.pgp53),dmr.pgp5.c=rownames(dmr.pgp5.c)), 
    resolution = 300, imagetype = "tiff", alpha=c(0.5,0.5), cex = 2,cat.cex = 2,
    fill=c("red","blue"), cat.fontface=4,fontfamily=3,
    main="deg.dmr.pgp53.c",
    main.cex = 2, main.fontface = 2, main.fontfamily = 3,
    filename = "deg.dmr.pgp53.c.VennDiagram.tif")	

