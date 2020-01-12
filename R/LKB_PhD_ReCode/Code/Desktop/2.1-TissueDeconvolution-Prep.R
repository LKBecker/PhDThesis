#libraries and functions====
library(Biobase); library(bseqsc); library(data.table); library(DESeq2); library(ggplot2)
library(Matrix); library(plyr); library(reshape2); library(xbioc)
require(MuSiC)
source("Code/Desktop/functions.r")
rm(AskYN, give.n, se, GenerateGeneSummaryPlots, completeDT, DoHeatmap, LeafCutterContainerToResultsObj)
library(DESeq2)

#1: create an ExpressionSet to hold our tissue level data----
tmsg("Loading in-house sheep atrium data...")
AtrialRNAData <- readRDS("Data/RNAseq/Output/SplineModel.ERCC.ImputeB4Exclude.2excl.maxit5000.batchCutRINOXnsImputedAgeDaysdf3.LRT.RDS")
AtrialRNACounts  <- DESeq2::counts(AtrialRNAData, normalized=FALSE)

metadata = colData(AtrialRNAData)
metadata$sample = rownames(metadata)
metadata=metadata[,which(colnames(metadata) %in% c("sample", "batch", "Timepoint"))]
metadata=as.data.frame(metadata)
metadataLabels = data.frame(row.names = c("sample", "batch", "Group"),
                            labelDescription=c("Unique identifier of each RNA sample", "batch ID (taken from day of extraction)",
                                               "Group membership (timepoint)"))
metadata=AnnotatedDataFrame(metadata, metadataLabels)

#2: preparing single-cell RNA seq data from Skelly et al 's data====
tmsg("Loading Skelly's murine whole heart data...")
SkellyClusterIDs = fread("Data/MuSiC/Skelly/cluster_labels.tsv.txt")
SkellyClusterIDs.DF = as.data.frame(SkellyClusterIDs)
rownames(SkellyClusterIDs.DF)=SkellyClusterIDs$cell
readMatrix = fread("Data/MuSiC/Skelly/full_count_matrix.txt")
#readMatrix[,V1:=toupper(V1)] #PISD
readMatrix = as.matrix(readMatrix, rownames="V1")
readMatrix = readMatrix[, which(colnames(readMatrix) %in% rownames(SkellyClusterIDs.DF))]

#skelly's data uses mouse gene symbols; our internal data uses ENSEMBL IDs, ergo:
#We require a 'translation'; a simple toUpper() won't work due to differences in naming schemes between mouse and sheep (thanks, -science-)
#We get the necessary 'rosetta stone' from BIOMART:
tmsg("Filtering for 1:1 orthologs...")
MouseHumanSheepIDs = fread("Data/biomart/181205_SheepMouseHuman_geneIDs_geneNames_OrthologTypes.txt")
MouseHumanSheepIDs = unique(MouseHumanSheepIDs[`Mouse homology type`=="ortholog_one2one",
                                               .(`Gene stable ID`,`Gene name`, `Mouse gene stable ID`, `Mouse gene name`, `Mouse homology type`)])
#ortholog == "ortholog_one2one" leaves only certain, 1:1 relationships
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
#SkellyNonCMNonEndoData <- ExpressionSet(assayData = readMatrix, phenoData = AnnotatedDataFrame(SkellyClusterIDs.DF))
#SkellyCellTypeProportions <- SkellyClusterIDs[,.N,.(cluster)][,.(cluster, "Number of Cells"=N, "Percent of Dataset"=round((N/sum(N))*100, 2))]

#3: preparing DeLaughter dataset====
#See also: Auxiliary/3-CleanDeLaughterData.R

#To nobody's surprise, the gene IDs here are an archaic and badly-documented format, because why should I have nice things
#I know what will easily and quickly solve my problem easier: 
#SQL!
#library(org.Mm.eg.db)
#DeLaughterIDMap = select(org.Mm.eg.db, keys = rownames(DeLaughterMatrix), keytype = "ENSEMBL", columns = c("SYMBOL", "ENSEMBL"))

#The above gives one:many relationships so let's just use the biomart export, which delivers, presumably, the 'best' entry for each ID <-> Gene Symbol pair...
DeLaughterMatrix = readRDS("Data/MuSiC/DeLaughter/DevCell_LKB-STAR_p21CM-Only.rds")
DeLaughter.Meta = fread("Data/MuSiC/DeLaughter/190115-195902_DeLaughter.Metadata.tsv")
DeLaughter.Meta = DeLaughter.Meta[,.(cell=Library_Name, cluster=Cell_Type)]
DeLaughter.Meta = as.data.frame(DeLaughter.Meta)
rownames(DeLaughter.Meta)=DeLaughter.Meta$cell
DeLaughterDESeq = DESeqDataSetFromMatrix(DeLaughterMatrix, DeLaughter.Meta, design=~1)
DeLaughterDESeq = estimateSizeFactors(DeLaughterDESeq)
DeLaughterDESeq.Norm = DESeq2::counts(DeLaughterDESeq, normalized=TRUE)

CMRowNames = rownames(DeLaughterDESeq.Norm)
CMRowNames = merge(rownames(DeLaughterDESeq.Norm), MouseHumanSheepIDs[,.(`Mouse gene stable ID`,`Mouse gene name`)], by.x="x", by.y="Mouse gene stable ID", all.x=T)
setDT(CMRowNames)
CMRowNames[is.na(`Mouse gene name`), `Mouse gene name`:=x]
CMRowNames=CMRowNames[,`Mouse gene name`]
rownames(DeLaughterDESeq.Norm)=unname(unlist2(CMRowNames))

#Science fiction, double filter: we can only take those genes forward that are in Skelly's data AND our own data
dim(DeLaughterDESeq.Norm) #31053
DeLaughterDESeq.Norm = DeLaughterDESeq.Norm[which(rownames(DeLaughterDESeq.Norm) %in% MouseHumanSheepIDs[,`Mouse gene name`]),]
dim(DeLaughterDESeq.Norm) #15136

#Now this needs to be merged with the other matrix via cbind()...
FullCardiacCellCounts=merge(DeLaughterDESeq.Norm, readMatrix, by="row.names")
rownames(FullCardiacCellCounts)=FullCardiacCellCounts$Row.names
FullCardiacCellCounts$Row.names=NULL
FullCardiacCellCounts = as.matrix(FullCardiacCellCounts)

FullCardiacCellCounts.Meta = rbind(DeLaughter.Meta, SkellyClusterIDs.DF)

FullCardiacCellSet = ExpressionSet(FullCardiacCellCounts, AnnotatedDataFrame(FullCardiacCellCounts.Meta))
    
#Export Result====
AtrialTissueExpressionData <- ExpressionSet(assayData=AtrialRNACounts, phenoData=metadata)
tmsg("Saving expression sets...")
#saveRDS(SkellyNonCMNonEndoData, "Data/MuSic/Skelly-scRNAData.rds")
saveRDS(FullCardiacCellSet, "Data/MuSiC/Skelly+DeLaughter-scRNAData-DESeqNorm.rds")
saveRDS(AtrialTissueExpressionData, "Data/MuSic/AtrialTissue-whole-RNA-seq-Data.rds")
tmsg("Complete.")