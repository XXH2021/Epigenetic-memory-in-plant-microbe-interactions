library(ggplot2)
library(ggraph)
library(igraph)
library(tidyverse)
library(viridis)

lesfe =read.table("Galaxy71-[B)_LDA_Effect_Size_(LEfSe)_on_data_70].lefse_internal_res", sep="\t",row.names=1)
abundance =read.table("18.txt",sep="\t",row.names=1)
for (i in 1:dim(abundance)[1]){
	a<-rownames(abundance)[i]
	a<-gsub("\\|", "\\.",a)
	a<-gsub(" ", "",a)
	a<-gsub("[\\(\\)]","_",a,perl=TRUE)
	a<-gsub("[\\[\\]]","_",a,perl=TRUE)
	a<-gsub("-","_",a,perl=TRUE)
	rownames(abundance)[i]<-a	
	}
colnames(abundance)<-c("ck31","ck32","ck33","pgp531","pgp532","pgp533")
abundance<-t(abundance[-1,])
abundance<-data.frame(abundance,stringsAsFactors = F)
abundance<-as.data.frame(lapply(abundance,as.numeric))

abundance$mean=c(rep("CK3",3),rep("PGP53",3))
abundance.mean<-aggregate(abundance[,-dim(abundance)[2]],by=list(abundance$mean),FUN=mean)
abundance.mean<-t(abundance.mean)
colnames(abundance.mean)<-abundance.mean[1,]
abundance.mean<-as.data.frame(abundance.mean[-1,])
abundance.mean$taxon<-rownames(abundance.mean)
lesfe$taxon<-rownames(lesfe)
comb<-merge(abundance.mean,lesfe,by="taxon")

for (i in 1:dim(comb)[1]){
  spln<-unlist(strsplit(comb$taxon[i],"\\."))
  if (i >1){
    spll<-unlist(strsplit(comb$taxon[i-1],"\\."))
    for (j in 1:length(spln)){
      if (!setequal(grep("norank",spln[j]),integer(0))){
        if (!setequal(grep("norank",spll[j]),integer(0))){
          spln[j] = spll[j]
        }
        else if (!setequal(grep("norank",spln[j-1]),integer(0))){
          spln[j] = spln[j-1]
        }
        else {
          spln[j] = paste(spln[j],i,sep = "_")
        }
      }
    }
    comb$taxon[i] = paste(spln,collapse=".")
  }
}

for (i in 1:dim(comb)[1]){
  spln<-unlist(strsplit(comb$taxon[i],"\\."))
  if (i >1){
    spll<-unlist(strsplit(comb$taxon[i-1],"\\."))
    for (j in 1:length(spln)){
      if (!setequal(grep("Other",spln[j]),integer(0))){
        if (!setequal(grep("Other",spll[j]),integer(0))){
          spln[j] = spll[j]
        }
        else if (!setequal(grep("Other",spln[j-1]),integer(0))){
          spln[j] = spln[j-1]
        }
        else {
          spln[j] = paste(spln[j],i,sep = "_")
        }
      }
    }
    comb$taxon[i] = paste(spln,collapse=".")
  }
}

for (i in 1:dim(comb)[1]){
  spln<-unlist(strsplit(comb$taxon[i],"\\."))
  if (i >1){
    spll<-unlist(strsplit(comb$taxon[i-1],"\\."))
    for (j in 1:length(spln)){
      if (!setequal(grep("uncultured",spln[j]),integer(0))){
        if (!setequal(grep("uncultured",spll[j]),integer(0))){
          spln[j] = spll[j]
        }
        else if (!setequal(grep("uncultured",spln[j-1]),integer(0))){
          spln[j] = spln[j-1]
        }
        else {
          spln[j] = paste(spln[j],i,sep = "_")
        }
      }
    }
    comb$taxon[i] = paste(spln,collapse=".")
  }
}
write.csv(comb,"combmean.csv")

comb.edges<-comb
for (i in 1:dim(comb.edges)[1]){
	b<-comb.edges$taxon[i]
	comb.edges$to[i]<-comb.edges$taxon[i]
	spl<-unlist(strsplit(b,"\\."))
	comb.edges$from[i]<-paste(spl[-length(spl)],collapse=".")
}
rownames(comb.edges)<-c(1:dim(comb.edges)[1])
edges<-comb.edges[,c("from","to")]
edges<-edges[which(edges$from !=""),]
rownames(edges)<-c(1:dim(edges)[1])
write.csv(edges,"combmean.edges.csv")

comb.vertices<-comb
comb.vertices$shortName<-gsub("^.*\\.","",comb.vertices$taxon)
colnames(comb.vertices)[5]<-c("enrich")
comb.vertices$enrich<-as.vector(comb.vertices$enrich)

for (i in 1:dim(comb.vertices)[1]){
	spl<-unlist(strsplit(comb.vertices$taxon[i],"\\."))
	if (comb.vertices$enrich[i]==""){
		if (length(spl) >= 6){
		comb.vertices$enrich[i]<-c("NO")
		}		
	}
	if (length(spl) >= 6){
	  comb.vertices$shortName[i]<-""
	}
	else{
		if (comb.vertices$enrich[i] %in% c("CK-3","CK-30")){
		  comb.vertices$enrich[i]<-""
		}
	}
}
write.csv(comb.vertices,"comb.vertices.combmean.vertices.csv")

k = 1
genera = data.frame()
for (i in 1:dim(comb.vertices)[1]){
  spl<-unlist(strsplit(comb.vertices$taxon[i],"\\."))
  if (length(spl)==6){
   genera[k,1] = comb.vertices$taxon[i]
   genera[k,2] = spl[1]
   genera[k,3] = spl[2]
   genera[k,4] = spl[3]
   genera[k,5] = spl[4]
   genera[k,6] = spl[5]
   genera[k,7] = spl[6]
   genera[k,8] = comb.vertices$CK3[i]
   genera[k,9] = comb.vertices$PGP53[i]
   k = k+1
  }
}  
colnames(genera)=c('taxon','level1','level2','level3','level4','level5','level6','CK3','PGP5')  
comb.vertices$CK3 = as.numeric(as.character(comb.vertices$CK3))
comb.vertices$PGP53 = as.numeric(as.character(comb.vertices$PGP53))
write.csv(genera,"genera.comb.vertices.combmean.vertices.csv")

for (i in 1:dim(comb.vertices)[1]){
  spl<-unlist(strsplit(comb.vertices$taxon[i],"\\."))
  if (length(spl)==1){
    tempdatack3 = genera$CK3[which(genera$level1 %in% spl[1])]
    tempsumck3 = sum(as.numeric(as.character(tempdatack3)))
    tempdatapgp5 = genera$PGP5[which(genera$level1 %in% spl[1])]
    tempsumpgp5 = sum(as.numeric(as.character(tempdatapgp5)))   
    comb.vertices$CK3[i] = tempsumck3
    comb.vertices$PGP53[i] = tempsumpgp5
  }else if (length(spl)==2){
    tempdatack3 = genera$CK3[which(genera$level2 %in% spl[2])]
    tempsumck3 = sum(as.numeric(as.character(tempdatack3)))
    tempdatapgp5 = genera$PGP5[which(genera$level2 %in% spl[2])]
    tempsumpgp5 = sum(as.numeric(as.character(tempdatapgp5)))   
    comb.vertices$CK3[i] = tempsumck3
    comb.vertices$PGP53[i] = tempsumpgp5
  }
  else if (length(spl)==3){
  tempdatack3 = genera$CK3[which(genera$level3 %in% spl[3])]
  tempsumck3 = sum(as.numeric(as.character(tempdatack3)))
  tempdatapgp5 = genera$PGP5[which(genera$level3 %in% spl[3])]
  tempsumpgp5 = sum(as.numeric(as.character(tempdatapgp5)))   
  comb.vertices$CK3[i] = tempsumck3
  comb.vertices$PGP53[i] = tempsumpgp5
  }  
  else if (length(spl)==4){
    tempdatack3 = genera$CK3[which(genera$level4 %in% spl[4])]
    tempsumck3 = sum(as.numeric(as.character(tempdatack3)))
    tempdatapgp5 = genera$PGP5[which(genera$level4 %in% spl[4])]
    tempsumpgp5 = sum(as.numeric(as.character(tempdatapgp5)))   
    comb.vertices$CK3[i] = tempsumck3
    comb.vertices$PGP53[i] = tempsumpgp5
  }
  else if (length(spl)==5){
    tempdatack3 = genera$CK3[which(genera$level5 %in% spl[5])]
    tempsumck3 = sum(as.numeric(as.character(tempdatack3)))
    tempdatapgp5 = genera$PGP5[which(genera$level5 %in% spl[5])]
    tempsumpgp5 = sum(as.numeric(as.character(tempdatapgp5)))   
    comb.vertices$CK3[i] = tempsumck3
    comb.vertices$PGP53[i] = tempsumpgp5
  }  
}


rownames(comb.vertices)<-c(1:dim(comb.vertices)[1])
vertices<-comb.vertices[,c("taxon","PGP53","shortName","enrich")]
colnames(vertices)<-c("taxon","size","shortName","enrich")
vertices$size<-as.vector(vertices$size)
vertices$size<-as.numeric(vertices$size)

phylum<-c("Archaea","Bacteria","Bacteria.Proteobacteria","Bacteria.Actinobacteria","Bacteria.Chloroflexi","Bacteria.Acidobacteria","Bacteria.Gemmatimonadetes","Bacteria.Bacteroidetes","Bacteria.Firmicutes")
vertices[which(!(vertices$taxon %in% phylum)),c("shortName")]<-NA
write.csv(vertices,"combmean.vertices.csv")

edges<-read.csv("combmean.edges.csv",header = T)
edges<-edges[,-1]
vertices<-read.csv("combmean.vertices.csv",header = T)
vertices<-vertices[,-1]

vertices<-vertices[which(vertices$size>0),]
edges<-edges[which(edges$to %in% vertices$taxon),]

del.from<-which(!(edges$from %in% vertices$taxon))

if (!setequal(del.from,integer(0))){
  edges<-edges[-del.from, ]
}

vertices<-vertices[order(vertices$size,decreasing=T),]
write.csv(vertices,"1vertices.csv")
write.csv(edges,"1edges.csv")

mygraph <- graph_from_data_frame(edges, vertices=vertices )
 
main_theme = theme(panel.background=element_blank(),
                    panel.grid=element_blank(),
                    legend.position="right",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    legend.text= element_text(size=7),
                    text=element_text(family="sans", size=7)
					) 
color<-c("white","blue","yellow","green")
p <- ggraph(mygraph, layout = 'circlepack', weight=size) + 
  geom_node_circle(aes(fill = as.factor(enrich))) +
  scale_fill_manual(values=color)+
  geom_node_label(aes(label=shortName, size=size))+
  main_theme +
  theme(legend.position="FALSE")
ggsave("pgpCK32mean.pdf", p, width = 5.6, height = 5.6) 


  