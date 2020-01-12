#libraries and functions====
library(Biobase); library(bseqsc); library(data.table); library(DESeq2); library(ggplot2)
library(Matrix); library(plyr); library(reshape2); library(xbioc); require(MuSiC)
source("Code/Desktop/functions.r")
rm(AskYN, give.n, se, GenerateGeneSummaryPlots, completeDT, DoHeatmap, LeafCutterContainerToResultsObj)

#skelly's data uses mouse gene symbols; our internal data uses ENSEMBL IDs, ergo:
#We require a 'translation'; a simple toUpper() won't work due to differences in naming schemes between mouse and sheep (thanks, -science-)
#We get the necessary 'rosetta stone' from BIOMART:
tmsg("Filtering for 1:1 orthologs...")
MouseHumanSheepIDs = fread("Data/biomart/181205_SheepMouseHuman_geneIDs_geneNames_OrthologTypes.txt")
MouseHumanSheepIDs = unique(MouseHumanSheepIDs[`Mouse homology type`=="ortholog_one2one",
                                               .(`Gene stable ID`,`Gene name`, `Mouse gene stable ID`, `Mouse gene name`, `Mouse homology type`)])
#ortholog == "ortholog_one2one" leaves only certain, 1:1 relationships

#1: preparing single-cell RNA seq data from Skelly et al 's data====
tmsg("Loading Skelly's murine whole heart data...")
SkellyClusterIDs = fread("Data/MuSiC/Skelly/cluster_labels.tsv.txt")
SkellyClusterIDs.DF = as.data.frame(SkellyClusterIDs)
rownames(SkellyClusterIDs.DF)=SkellyClusterIDs$cell
readMatrix = fread("Data/MuSiC/Skelly/full_count_matrix.txt")
#readMatrix[,V1:=toupper(V1)] #PISD
readMatrix = as.matrix(readMatrix, rownames="V1")
readMatrix = readMatrix[, which(colnames(readMatrix) %in% rownames(SkellyClusterIDs.DF))]

tmsg("Filtering for 1:1 orthologs...")
#now we only retain those genes we know we have a sheep 1:1 ortholog of
readMatrix = readMatrix[which(rownames(readMatrix) %in% MouseHumanSheepIDs[,`Mouse gene name`]),]
#SkellyNonCMNonEndoData <- ExpressionSet(assayData = readMatrix, phenoData = AnnotatedDataFrame(SkellyClusterIDs.DF))
#SkellyCellTypeProportions <- SkellyClusterIDs[,.N,.(cluster)][,.(cluster, "Number of Cells"=N, "Percent of Dataset"=round((N/sum(N))*100, 2))]

#2: preparing DeLaughter dataset (See also: Auxiliary/3-CleanDeLaughterData.R)====
#Using the biomart export, which delivers, presumably, the 'best' entry for each ID <-> Gene Symbol pair...
DeLaughterMatrix = readRDS("Data/MuSiC/DeLaughter/DevCell_LKB-STAR_p21-CM-FB-EC.rds")
DeLaughter.Meta = fread("Data/MuSiC/DeLaughter/190128-110504_DeLaughter.Metadata.tsv")
DeLaughter.Meta = DeLaughter.Meta[,.(cell=Library_Name, cluster=Cell_Type)]
DeLaughter.Meta = as.data.frame(DeLaughter.Meta)
rownames(DeLaughter.Meta)=DeLaughter.Meta$cell

CMRowNames = merge(rownames(DeLaughterMatrix), MouseHumanSheepIDs[,.(`Mouse gene stable ID`,`Mouse gene name`)], by.x="x", by.y="Mouse gene stable ID", all.x=T)
setDT(CMRowNames)
CMRowNames[is.na(`Mouse gene name`), `Mouse gene name`:=x]
CMRowNames=CMRowNames[,`Mouse gene name`]
rownames(DeLaughterMatrix)=unname(unlist2(CMRowNames))

#Science fiction, double filter: we can only take those genes forward that are in Skelly's data AND our own data
dim(DeLaughterMatrix) #31053
DeLaughterMatrix = DeLaughterMatrix[which(rownames(DeLaughterMatrix) %in% MouseHumanSheepIDs[,`Mouse gene name`]),]
dim(DeLaughterMatrix) #15136

#Now this needs to be merged with the other matrix via cbind()...
FullCardiacCellCounts=merge(DeLaughterMatrix, readMatrix, by="row.names")
rownames(FullCardiacCellCounts)=FullCardiacCellCounts$Row.names
FullCardiacCellCounts$Row.names=NULL
FullCardiacCellCounts = as.matrix(FullCardiacCellCounts)

FullCardiacCellCounts.Meta = rbind(DeLaughter.Meta, SkellyClusterIDs.DF)

FullCardiacCellSet = ExpressionSet(FullCardiacCellCounts, AnnotatedDataFrame(FullCardiacCellCounts.Meta))

#Run Music_Theta()
#FullCardiacCellSet.Theta = music_Theta(FullCardiacCellSet, non.zero = TRUE, clusters="cluster", samples="cell")
FullCardiacCellSet.Basis = music_basis(FullCardiacCellSet, non.zero = TRUE, clusters="cluster", samples="cell")

#Export Result====
#tmsg("Saving expression set...")
#saveRDS(FullCardiacCellSet, "Data/MuSiC/Skelly+DeLaughter-scRNAData-DESeqNorm.rds")
