#libraries and functions====
library(Biobase)
library(bseqsc)
library(data.table)
library(DESeq2)
library(ggplot2)
library(Matrix)
library(plyr)
library(reshape2)
library(xbioc)
require(MuSiC)

#1: create an ExpressionSet to hold our tissue level data----
	AtrialRNAData <- readRDS("Data/DESeq2.OvAr31+ERCC.2FilesRM.batchCutRINOX3Group.maxit2000.7nonconverge.RDS")
	AtrialRNACounts = AtrialRNACounts = readRDS("Data/AtrialRNACountMatrix.rds")
	#AtrialRNACounts = AtrialRNACounts[which(rowSums(AtrialRNACounts)!=0),]
	metadata = colData(AtrialRNAData)
	metadata=metadata[,which(colnames(metadata) %in% c("sample", "batch", "Group"))]
	metadata=data.frame(metadata)
	metadataLabels = data.frame(row.names = c("sample", "batch", "Group"),
													 labelDescription=c("Unique identifier of each RNA sample", "batch ID (taken from day of extraction)",
													 									 "Group membership (timepoint)"))
	metadata=AnnotatedDataFrame(metadata, metadataLabels)

#2: preparing single-cell RNA seq data from Skelly et al 's data====
	SkellyClusterIDs = fread("Data/cluster_labels.tsv.txt")
	SkellyClusterIDs.DF = as.data.frame(SkellyClusterIDs)
	rownames(SkellyClusterIDs.DF)=SkellyClusterIDs$cell
	readMatrix = fread("Data/E-MTAB-6173.processed.1/full_count_matrix.txt")
	readMatrix = as.matrix(readMatrix, rownames="V1")
	#Skelly appear to have filtered out some columns, which are not represented in the cluster data; we too have to remove them
	readMatrix = readMatrix[, which(colnames(readMatrix) %in% rownames(SkellyClusterIDs.DF))]

	#skelly's data uses mouse gene symbols; our internal data uses ENSEMBL IDs, ergo:
	#We require a 'translation'; a simple toUpper() won't work due to differences in naming schemes between mouse and sheep (thanks, -science-)
	#We get the necessary 'rosetta stone' from BIOMART:
	MouseHumanSheepIDs = fread("Data/181205_SheepMouseHuman_geneIDs_geneNames_OrthologTypes.txt")
	MouseHumanSheepIDs = unique(MouseHumanSheepIDs[`Mouse homology type`=="ortholog_one2one",
																						.(`Gene stable ID`,`Gene name`, `Mouse gene stable ID`, `Mouse gene name`, `Mouse homology type`)])
		#We filter our the certain, 1:1 relationships
	MouseHumanSheepIDs[,length(unique(`Gene stable ID`))]

	#Now we filter the atrial data to only those that have a 1:1 mouse ortholog
	AtrialRNACounts = AtrialRNACounts[which(rownames(AtrialRNACounts) %in% MouseHumanSheepIDs[,`Gene stable ID`]),]
	#dim(AtrialRNACounts) #Matches the length of the filtered gene list (15956) exactly if you do NOT filter out rowSums()==0
	#Now we need to replace the row names
	SheepRowNames = data.table(`Gene stable ID`=rownames(AtrialRNACounts))
	SheepRowNames = merge(SheepRowNames, MouseHumanSheepIDs[,.(`Gene stable ID`,`Mouse gene name`)], by="Gene stable ID", all.x=T)
	rownames(AtrialRNACounts) = SheepRowNames[,`Mouse gene name`]

	#now we only retain those genes we know we have a sheep 1:1 ortholog of
	readMatrix = readMatrix[which(rownames(readMatrix) %in% MouseHumanSheepIDs[,`Mouse gene name`]),]
	SkellyNonCMNonEndoData <- ExpressionSet(assayData = readMatrix, phenoData = AnnotatedDataFrame(SkellyClusterIDs.DF))
	SkellyCellTypeProportions <- SkellyClusterIDs[,.N,.(cluster)][,.(cluster, "Number of Cells"=N, "Percent of Dataset"=round((N/sum(N))*100, 2))]

	AtrialTissueExpressionData <- ExpressionSet(assayData=AtrialRNACounts, phenoData=metadata)
	saveRDS(SkellyNonCMNonEndoData, "Data/Skelly-scRNAData.rds")
	saveRDS(AtrialTissueExpressionData, "Data/AtrialTissue-whole-RNA-seq-Data.rds")
	rm(readMatrix, SkellyClusterIDs, SkellyClusterIDs.DF)

#3: TODO preparing single-cell RNA seq data from DeLaughter====
#TODO

#4: TODO merging the two murine scRNA-datasets???====

#5: MuSiC====

	CardiacSingleCellData.Expressed <- rownames(SkellyNonCMNonEndoData)[rowMeans(exprs(SkellyNonCMNonEndoData)) != 0] #Filters for non-zero expression?
	sc_basis <- music_basis(x = SkellyNonCMNonEndoData, non.zero = TRUE, markers = CardiacSingleCellData.Expressed, clusters = "cluster", samples = "cell", verbose = T)
	saveRDS(sc_basis, "Data/Skelly-Sc-Basis.rds")

singleCell_sigma <- sc_basis$Sigma


	# Estimate cell type proportions of artificial bulk data
	props <- music_prop(AtrialTissueExpressionData, SkellyNonCMNonEndoData, clusters="cluster", samples="cell")
	saveRDS(props, "Data/AtrialTissue-CellProportionsMkI.rds")

# Fri Dec 07 16:18:25 2018 works

#props$r.squared.full

stop()
mb <- music_basis(sc_es, clusters="cluster", samples="cell")
#mb$Disgn.mtx = gene x cell_type Design matrix
#mb$S = subject x cell_type Library sizes
#mb$M.S = Cell type average library sizes
#mb$M.theta = Gene x cell_type matrix of average relative abundance
#mb$Sigma = Gene x cell_type matrix of cross-sample variation

d <- dist(t(log(mb$Disgn.mtx + 1e-6)), method = "euclidean")
# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "complete" )
# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1, main = 'Cluster log(Design Matrix)')
phc <- as.phylo(hc1)
(ggtobj <- ggtree(phc,size=1.5) + xlim(c(NA,600)) + geom_tiplab(size=3,offset = 5, lineheight = 0.8) + geom_treescale())
ggsave("tree-basic.pdf",plot=ggtobj,width=4,height=9)

d <- dist(t(log(mb$M.theta + 1e-8)), method = "euclidean")
# Hierarchical clustering using Complete Linkage
# hc2 <- hclust(d, method = "complete" )
hc2 <- hclust(d, method = "complete")
# Plot the obtained dendrogram
plot(hc2, cex = 0.6, hang = -1, main = 'Cluster log(Mean of RA)')

pr <- readRDS("music-props.rds")
str(pr,max.level = 2)

mp <- data.table(melt(pr[["Weight.gene"]],varnames = c("ensembl_gene_id","cluster_id"),value.name = "prop"))

mp <- merge(mp,refined_clusters[,.(cluster_id,cluster_group,cluster_colour)],by="cluster_id")

mp <- mp[,.(prop=sum(prop)),.(cluster_group, sample_id, cluster_colour)]
mp[,sample_type:=ifelse(sample_id %like% "POZT|N","Biopsy","Nephrectomy")]
mp[,study:=ifelse(sample_id %like% "-","TCGA","MTG")]
mp[,old_prop:=prop]
mp[,prop:=prop/sum(prop),sample_id]
mp[,sum(prop),sample_id]

mat <- acast(mp,cluster_group ~ sample_id, value.var = "prop")
hc <- hclust(dist(mat))
plot(hc)
ordered_samples <- hc$labels[hc$order]
