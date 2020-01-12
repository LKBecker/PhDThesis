#libraries====
library(data.table)
library(DESeq2)
library(ggplot2)
library(splines)
library(stringr)
#functions====
source("Code/Desktop/functions.r")
source("Code/Desktop/functions.DESeq2.r")
MakeRNAModelOutputname <- function(ModelName, TestUsed){
    if (! ModelName %in% c("Spline", "Group")) { stop("Model must be either 'Group' or 'Spline'")}
    Filename = "Data/RNAseq/Output/"
    NExcluded = paste0(length(CONF.RNAseq$FILE_EXCLUDE_LIST),"excl")
    MaxIter   = paste0("maxit", CONF.RNAseq$MAX_WALD_ITERATIONS)
    if(ModelName=="Spline"){
        Model = as.character(CONF.RNAseq$SPLINE_FORMULA)[2]    
    }else{
        Model = as.character(CONF.RNAseq$GROUP_FORMULA)[2]    
    }
    Model = str_remove_all(Model, "\\+")
    Model = str_remove(Model, "\\.")
    Model = str_remove(Model, "\\(")
    Model = str_remove(Model, "\\)")
    Model = str_remove(Model, "\\,")
    Model = str_remove(Model, "=")
    Model = str_remove_all(Model, " ")
    RNAModelOutputName = paste(paste0(Filename, ModelName, "GWPM"), NExcluded, 
                               MaxIter, Model, TestUsed, "RDS", sep = ".") 
    RNAModelOutputName
}
#Config====
CONF.RNAseq = list()
CONF.RNAseq[["VERSION"]]		        = "0.0-r1"
CONF.RNAseq[["FILE_EXCLUDE_LIST"]]	    = c("12C", "13R")
CONF.RNAseq[["MIN_FDR"]]				= 0.05
CONF.RNAseq[["MIN_LFC"]]				= 0.0
CONF.RNAseq[["MAX_WALD_ITERATIONS"]]    = 5000 #do not set to 999999, it takes four+ hours and still doesn't converge
CONF.RNAseq[["COUNTS_FOLDER"]]          = "Data/RNAseq/Input/CountFiles.GWPM/"
CONF.RNAseq[["GROUP_FORMULA"]]          = formula(~Batch+Group)
#Loading covariates====
tmsg(paste0("Preparing to generate GWPM RNAseq models, script version ", CONF.RNAseq$VERSION, "..."))
tmsg("Loading saved covariates...")
GWPMBatchCovariate <- fread("Data/RNAseq/Input/181025-160817_GWPM.Metadata.tsv", na.strings = "#N/A") 
#two missing samples (12C, 13R) are excluded by PCA; metadatafile must've been written after

#Listing .counts files====
FileNames.GWPM <- data.table(fileName = dir(pattern = "*.counts", path = CONF.RNAseq$COUNTS_FOLDER, full.names = T))
FileNames.GWPM[,Sample:=factor(str_match(fileName, ".*/((\\d{1,2})(C|R|HF))")[,2])]
FileNames.GWPM = merge(FileNames.GWPM, GWPMBatchCovariate, by="Sample", all=T)

FileNames.GWPM[,Batch:=factor(Batch)]
FileNames.GWPM[,Group:=factor(Group)]
setcolorder(FileNames.GWPM, unique(c("Sample", "fileName", colnames(FileNames.GWPM))))
tmsg(paste0("Scanned .counts folder '", getwd(), .Platform$file.sep, CONF.RNAseq$COUNTS_FOLDER, "'. ", nrow(FileNames.GWPM), " files found."))

#PCA====
DESeqData.Raw = DESeqDataSetFromHTSeqCount(sampleTable = FileNames.GWPM,  design = ~1) #No design needed for PCA
DESeqData.Raw.VST   <- vst(DESeqData.Raw) #Variance-Stabilising Transformation

rv <- rowVars(assay(DESeqData.Raw.VST))
select <- order(rv, decreasing = TRUE)[seq_len(min(20000, length(rv)))]
pca <- prcomp(t(assay(DESeqData.Raw.VST)[select, ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar <- percentVar[1:2]
rm(select, rv, pca)

PCA = plotPCA(DESeqData.Raw.VST, intgroup="Group", ntop=20000, returnData=TRUE)
PCA = ggplot(PCA, aes(x=PC1, y=PC2, color=Group))+geom_label(aes(label=name))+theme(legend.position="none", plot.margin = margin(5, 5, 10, 5)) 
#ntop=20000 is excessive (5000 gives same result) but why not
PCA
rm(DESeqData.Raw, DESeqData.Raw.VST)
#This recommends excluding 13R, 12C

#Excluding according to PCA====
FileNames.GWPM.PostPCA = FileNames.GWPM[!( Sample %in% CONF.RNAseq$FILE_EXCLUDE_LIST )]
tmsg(paste0("Excluding files [", paste0(CONF.RNAseq$FILE_EXCLUDE_LIST, collapse = ", "), 
            "] according to PCA results, remainder [", nrow(FileNames.GWPM.PostPCA), "] files."))

#Group Model====
GroupModel <- DESeqDataSetFromHTSeqCount(sampleTable = FileNames.GWPM.PostPCA,  design = CONF.RNAseq$GROUP_FORMULA)
GroupModel <- GroupModel[ rowSums(counts(GroupModel)) > 1, ]
GroupModel <- estimateSizeFactors(GroupModel)
GroupModel <- estimateDispersions(GroupModel)
tmsg("Running group model...")
GroupModel <- nbinomWaldTest(GroupModel, maxit = CONF.RNAseq$MAX_WALD_ITERATIONS)

#Export to RDS====
tmsg("Exporting to RDS...")
saveRDS(GroupModel, MakeRNAModelOutputname("Group", "Wald"))
