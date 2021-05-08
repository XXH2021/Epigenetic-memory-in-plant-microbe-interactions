library(WGCNA)
library(reshape2)
library(stringr)
enableWGCNAThreads()
ALLOW_WGCNA_THREADS=4
type = "unsigned"
corType = "pearson"
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
maxPOutliers = ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)

exprMat <- "nrep_expr.csv"
dataExpr <- read.csv(exprMat,header=T,row.names=1)
dim(dataExpr)

dataExprVar <- dataExpr
zero.count1 <- function(x){
			a<-c()
			b<-c()
			out<-c()
			x<-as.numeric(x)
			for (i in 1:6){
				a[i]<-table(x[((i-1)*3+1):((i-1)*3+3)])[[1]]
				b[i]<-mean(x[((i-1)*3+1):((i-1)*3+3)])			
				out[i]=ifelse ((a[i] > 1 | b[i]< 1),0,1)		
				}
			sum(out)
			}
m.zero <- apply(dataExprVar,1,zero.count1)
dataExprVar2 <- dataExprVar[which(m.zero>1),]
dim(dataExprVar2)

dataExpr <- as.data.frame(t(dataExprVar2))
gsg = goodSamplesGenes(dataExpr, verbose = 3)
gsg$allOK
if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}
dim(dataExpr)

sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dataExpr<-dataExpr[c(-11,-18),]
dim(dataExpr)
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
rownames(dataExpr)
dataExpr<-dataExpr[c(-2,-6),]
rownames(dataExpr)
dim(dataExpr)
write.csv(dataExpr,"dataExpr12.09.csv")
sampleTree = hclust(dist(dataExpr), method = "average")
pdf("Sample clustering to detect outliers2.pdf")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
powers = c(c(1:10), seq(from = 12, to=30, by=1))
sft = pickSoftThreshold(dataExpr, powerVector=powers, networkType=type, verbose=5)
pdf("powerEstimate.pdf")
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")
dev.off()
power = sft$powerEstimate
power
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType, 
                       maxPOutliers=maxPOutliers, loadTOM =TRUE,
                       saveTOMFileBase = paste0(exprMat, ".tom"),
                       verbose = 3)

table(net$colors)

moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
pdf("dendrograms.pdf")
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
pdf("Eigengene adjacency heatmap.pdf")
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()

trait <- "dataTraits.txt"
if(trait != "") {
  traitData <- read.table(file=trait, sep='\t', header=T, row.names=1)
  traitData<-traitData[,c(1,8,9,10)]
  traitData <-data.frame(traitData)
  row.names(traitData)<-gsub("-",".",row.names(traitData))
  sampleName = rownames(dataExpr)
  traitData = traitData[match(sampleName, rownames(traitData)), ]
}

if (corType=="pearsoon") {
  modTraitCor = cor(MEs_col, traitData, use = "p")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
  modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
  modTraitCor = modTraitCorP$bicor
  modTraitP   = modTraitCorP$p
}

textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
pdf("gene_modules_trait2.pdf",width=4,height=8)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
               yLabels = colnames(MEs_col), 
               cex.lab = 0.5, 
               ySymbols = colnames(MEs_col), colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 0.5, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

if (corType=="pearsoon") {
  geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
             as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}

if (corType=="pearsoon") {
  geneTraitCor = as.data.frame(cor(dataExpr, traitData, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
             as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(dataExpr, traitData, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}

module = "lightcyan"
pheno = "taxonPC3"
modNames = substring(colnames(MEs_col), 3)

module_column = match(module, modNames)
pheno_column = match(pheno,colnames(traitData))

moduleGenes = moduleColors == module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))

pdf("magenta.pdf")
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

load(net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)
dissTOM = 1-TOM
plotTOM = dissTOM^power
diag(plotTOM) = NA

probes = colnames(dataExpr)
dimnames(TOM) <- list(probes, probes)
cyt = exportNetworkToCytoscape(TOM,
             edgeFile = paste(exprMat, ".edges.txt", sep=""),
             nodeFile = paste(exprMat, ".nodes.txt", sep=""),
             weighted = TRUE, threshold = 0,
             nodeNames = probes, nodeAttr = moduleColors)

probes = colnames(dataExpr)
colorlevels=unique(moduleColors)
modProbes<-matrix(nrow=max(table(net$colors)),ncol=length(colorlevels))
modProbes<-as.data.frame(modProbes)
colnames(modProbes)=colorlevels
for (i in colorlevels){     
     inModule = (moduleColors==i)
     modProbes[1:length(probes[inModule]),colnames(modProbes)==i] = probes[inModule]
}
count<-function(x){length(na.omit(x))}
apply(modProbes,2,count)
write.csv(modProbes,file="modProbes.csv")

pval<-as.data.frame(modTraitP)
write.csv(pval,file="modProbes_psig.csv")
groD3<-pval[pval$Group<0.05,]
psigD3<-modProbes[,sub('^..','',rownames(groD3))]
write.csv(psigD3,file="modProbes_psig_Group.csv")

ptaxonPC2<-pval[pval$taxonPC2<0.05,]
taxonPC2<-modProbes[,sub('^..','',rownames(ptaxonPC2))]
write.csv(taxonPC2,file="modProbes_psig_taxonPC2.csv")

ptaxonPC3<-pval[pval$taxonPC3<0.05,]
taxonPC3<-modProbes[,sub('^..','',rownames(ptaxonPC3))]
write.csv(taxonPC3,file="modProbes_psig_taxonPC3.csv")

par(mfrow=c(3,as.integer(1+ncol(psigD3)/4)))
pdf("Module.membership.gene.significance.Group.pdf")
pheno = "Group"
for (i in c(1:ncol(psigD3))){
	module=colnames(psigD3)[i]
	modNames = substring(colnames(MEs_col), 3)
	module_column = match(module, modNames)
	pheno_column = match(pheno,colnames(traitData))
	moduleGenes = moduleColors == module
	verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1, cex.lab = 1, cex.axis = 1, col = module,abline=TRUE)
}
dev.off()

par(mfrow=c(3,as.integer(1+ncol(taxonPC2)/4)))
pdf("Module.membership.gene.significance.taxonPC2.pdf")
pheno = "taxonPC2"
for (i in c(1:ncol(taxonPC2))){
	module=colnames(taxonPC2)[i]
	modNames = substring(colnames(MEs_col), 3)
	module_column = match(module, modNames)
	pheno_column = match(pheno,colnames(traitData))
	moduleGenes = moduleColors == module
	verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1, cex.lab = 1, cex.axis = 1, col = module,abline=TRUE)
}
dev.off()

par(mfrow=c(3,as.integer(1+ncol(taxonPC3)/4)))
pdf("Module.membership.gene.significance.taxonPC3.pdf")
pheno = "taxonPC3"
for (i in c(1:ncol(taxonPC3))){
	module=colnames(taxonPC3)[i]
	modNames = substring(colnames(MEs_col), 3)
	module_column = match(module, modNames)
	pheno_column = match(pheno,colnames(traitData))
	moduleGenes = moduleColors == module
	verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1, cex.lab = 1, cex.axis = 1, col = module,abline=TRUE)
}
dev.off()

connet=abs(cor(dataExpr,use="p"))^power
Alldegrees1=intramodularConnectivity(connet, moduleColors)
head(Alldegrees1)
EB= as.data.frame(traitData[,3]);# change specific 
names(EB) = "EB"
GS1 = as.numeric(cor(EB,dataExpr, use="p"))
GeneSignificance=abs(GS1)
colorlevels=unique(moduleColors)
sizeGrWindow(9,6)
par(mfrow=c(4,as.integer(1+length(colorlevels)/4)))
for (i in c(1:length(colorlevels)))
{
  whichmodule=colorlevels[[i]];
  restrict1 = (moduleColors==whichmodule);
  verboseScatterplot(Alldegrees1$kWithin[restrict1],
                     GeneSignificance[restrict1], col=moduleColors[restrict1],
                     main=whichmodule,
                     xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
}

datKME=signedKME(dataExpr, MEs_col, outputColumnName="MM.")

head(datKME)
EB= as.data.frame(traitData[,3]);# change specific 
names(EB) = "EB"
GS1 = as.numeric(cor(EB,dataExpr, use="p"))

FilterGenes= abs(GS1)> .5 & abs(datKME$MM.purple)>.8
table(FilterGenes)
trait_hubGenes_spe <- colnames(dataExpr)[FilterGenes]
write.csv(trait_hubGenes_spe, "trait_hubGenes_spe_purple_pc2_0.5_0.8.csv") 

FilterGenes= abs(GS1)> .7 & abs(datKME$MM.tan)>.8
table(FilterGenes)
trait_hubGenes_spe <- colnames(dataExpr)[FilterGenes]
write.csv(trait_hubGenes_spe, "trait_hubGenes_spe_tan_PC2_0.7_0.8.csv") 

FilterGenes= abs(GS1)> .9 & abs(datKME$MM.salmon)>.8
table(FilterGenes)
trait_hubGenes_spe <- colnames(dataExpr)[FilterGenes]
write.csv(trait_hubGenes_spe, "trait_hubGenes_spe_salmon_PC2_0.9_0.8.csv") 

FilterGenes= abs(GS1)> .9 & abs(datKME$MM.turquoise)>.8
table(FilterGenes)
trait_hubGenes_spe <- colnames(dataExpr)[FilterGenes]
write.csv(trait_hubGenes_spe, "trait_hubGenes_spe_turquoise_PC2_0.9_0.8.csv")

FilterGenes= abs(GS1)> .8 & abs(datKME$MM.brown)>.9
table(FilterGenes)
trait_hubGenes_spe <- colnames(dataExpr)[FilterGenes]
write.csv(trait_hubGenes_spe, "trait_hubGenes_spe_brown_PC2_0.8_0.9.csv") 

FilterGenes= abs(GS1)> .6 & abs(datKME$MM.black)>.8
table(FilterGenes)
trait_hubGenes_spe <- colnames(dataExpr)[FilterGenes]
write.csv(trait_hubGenes_spe, "trait_hubGenes_spe_black_PC2_0.6_0.8.csv")

EB= as.data.frame(traitData[,4]);# change specific 
names(EB) = "EB"
GS2 = as.numeric(cor(EB,dataExpr, use="p"))

FilterGenes= abs(GS2)> .5 & abs(datKME$MM.pink)>.8
table(FilterGenes)
trait_hubGenes_spe <- colnames(dataExpr)[FilterGenes]
write.csv(trait_hubGenes_spe, "trait_hubGenes_spe_pink_PC3_0.5_0.8.csv")

FilterGenes= abs(GS2)> .8 & abs(datKME$MM.magenta)>.9
table(FilterGenes)
trait_hubGenes_spe <- colnames(dataExpr)[FilterGenes]
write.csv(trait_hubGenes_spe, "trait_hubGenes_spe_magenta_PC3_0.8_0.9.csv")

FilterGenes= abs(GS2)> .6 & abs(datKME$MM.salmon)>.8
table(FilterGenes)
trait_hubGenes_spe <- colnames(dataExpr)[FilterGenes]
write.csv(trait_hubGenes_spe, "trait_hubGenes_spe_salmon_PC3_0.6_0.8.csv")

FilterGenes= abs(GS2)> .8 & abs(datKME$MM.lightcyan)>.8
table(FilterGenes)
trait_hubGenes_spe <- colnames(dataExpr)[FilterGenes]
write.csv(trait_hubGenes_spe, "trait_hubGenes_spe_lightcyan_PC3_0.8_0.8.csv")


EB= as.data.frame(traitData[,1]);# change specific 
names(EB) = "EB"
GS3 = as.numeric(cor(EB,dataExpr, use="p"))

FilterGenes= abs(GS3)> .8 & abs(datKME$MM.purple)>.8
table(FilterGenes)
trait_hubGenes_spe <- colnames(dataExpr)[FilterGenes]
write.csv(trait_hubGenes_spe, "trait_hubGenes_spe_purple_group.csv") 

load(net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)

module = "turquoise";
probes = colnames(dataExpr) 
inModule = (moduleColors==module);
modProbes = probes[inModule]; 
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
vis = exportNetworkToVisANT(modTOM,
                            file = paste("VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0)
cyt = exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
  weighted = FALSE,
  threshold = 0.05,
  nodeNames = modProbes, 
)




