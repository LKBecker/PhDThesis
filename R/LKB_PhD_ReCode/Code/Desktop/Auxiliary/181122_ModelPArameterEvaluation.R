rm(list = ls()) #clear environment
source("Code/libraries.r"); source("Code/functions.R"); source("Code/config.R"); source("Code/Graphs.ThemeKat.R")
library(DESeq2)
LKBRNACovariates <- fread("Data/Metadata/20170911-133126_RNASeq-LKB-Covariates.csv", na.strings = "#N/A")

FileNames <- data.table(fileName = dir(pattern = "*.counts", path = "Data/Count Files/", full.names = T))
tmsg(paste0(nrow(FileNames), " files (NO ERCC!) were found."))
FileNames[,Timepoint:=str_match(fileName, "Sample-\\d{6}-\\d{2}_(\\w{1,2})_")[,2]]
FileNames[,Timepoint:=factor(Timepoint, levels = c("1W", "1M", "3M", "AD"))]
FileNames[,sample:=str_match(fileName, "Sample-(\\d{6}-\\d{2})")[,2]]
FileNames[,shortName:= paste(sample, Timepoint, sep="_")]
FileNames[,batch:=factor(str_match(sample, "\\d{4}(\\d{2})-\\d{2}")[,2])]

Metadata <- merge(LKBRNACovariates, FileNames, by.x="ID_OX", by.y="shortName")
setcolorder(Metadata, unique(c("ID_OX", "fileName","sample", "batch", colnames(Metadata)))) 
set(Metadata, j = Metadata[, which((colnames(Metadata) %in%
                                      c("Sex","Location","HasEM", "HasPFA", "HasInVivo" ,"RNA_ID", "Vol", "260/230",
                                        "260/280", "ng/ul", "Subject")))], value=NULL)
Metadata <- Metadata[!(ID_OX %in% CONF$FILE_EXCLUDE_LIST)] #based on PCA,which indeed holds up with ERCC
tmsg(paste("Removing files [", paste0(CONF$FILE_EXCLUDE_LIST, collapse = ", "), " ] as outliers..."))
tmsg(paste0("Final count: ", nrow(Metadata), " files going into fitting and analysis."))
source("Code/AgeImputation.R")
Metadata[, CutRINOX:=cut(RIN_OX, 3) ] #/!\
rm(FileNames, LKBRNACovariates)

DESeqData.LKB <- DESeqDataSetFromHTSeqCount(sampleTable = Metadata,  design = ~batch+CutRINOX+ns(ImputedAgeDays)) 

#DESeqData.LKB <- DESeqData.LKB[ rowSums(counts(DESeqData.LKB)) > 1, ] # this changes stat outcomes 
DESeqData.LKB	<- estimateSizeFactors(DESeqData.LKB)
DESeqData.LKB	<- estimateDispersions(DESeqData.LKB)
DESeqData.LKB	<- nbinomLRT(DESeqData.LKB, reduced = ~batch+CutRINOX, maxit = 5000)

Spline.LRT.ERCC = GrabResult(DESeqData.LKB, "ns.ImputedAgeDays.", comparisonStr = "LRT", 0.05, TRUE)
message(Spline.LRT.ERCC[,length(unique(geneID))])
