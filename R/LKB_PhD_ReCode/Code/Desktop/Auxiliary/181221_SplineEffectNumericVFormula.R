#RNA-seq models for further analysis
#libraries and functions====
library(data.table)
library(DESeq2)
library(ggplot2)
library(splines)
library(stringr)
source("Code/Desktop/functions.r")
MakeRNAModelOutputname <- function(ModelName, TestUsed){
    if (! ModelName %in% c("Spline", "Group")) { stop("Model must be either 'Group' or 'Spline'")}
    Filename = "Data/RNAseq/Output"
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
    RNAModelOutputName = paste(paste0("Data/", ModelName, "Model.ERCC.ImputeB4Exclude"), NExcluded, 
                               MaxIter, Model, TestUsed, "RDS", sep = ".") 
    RNAModelOutputName
}
#Config====
CONF.RNAseq = list()
CONF.RNAseq[["VERSION"]]		        = "0.1-INIT"
CONF.RNAseq[["FILE_EXCLUDE_LIST"]]	    = c("161206-07_3M", "161207-05_3M")
CONF.RNAseq[["MIN_FDR"]]				= 0.05
CONF.RNAseq[["MIN_LFC"]]				= 0.0
CONF.RNAseq[["MAX_WALD_ITERATIONS"]]    = 5000 #do not set to 999999, it takes four hours and still doesn't converge
CONF.RNAseq[["COUNTS_FOLDER"]]          = "Data/RNAseq/Input/CountFiles.ERCC/"
CONF.RNAseq[["EXPERIMENTAL_FORMULA"]]   = formula(~batch+CutRINOX+Spline1+Spline2+Spline3)
CONF.RNAseq[["REDUCED_FORMULA"]]        = formula(~batch+CutRINOX)
#Loading covariates====
tmsg(paste0("Preparing to run 2-EffectOfSplineNumericVPart.R, version ", CONF.RNAseq$VERSION, "..."))
tmsg("Loading saved covariates...")
LKBRNACovariates <- fread("Data/RNAseq/Input/20170911-133126_RNASeq-LKB-Covariates.csv", na.strings = "#N/A")
#Listing .counts files====
FileNames.ERCC <- data.table(fileName = dir(pattern = "*.counts", path = CONF.RNAseq$COUNTS_FOLDER, full.names = T))
FileNames.ERCC[,Timepoint:=str_match(fileName, "Sample\\d{6}-\\d{2}_(\\w{1,2})_")[,2]]
FileNames.ERCC[,sample:=str_match(fileName, "Sample(\\d{6}-\\d{2})")[,2]]
FileNames.ERCC[,Timepoint:=factor(Timepoint, levels = c("1W", "1M", "3M", "AD"))]
FileNames.ERCC[,shortName:= paste(sample, Timepoint, sep="_")]
FileNames.ERCC[,batch:=factor(str_match(sample, "\\d{4}(\\d{2})-\\d{2}")[,2])]
tmsg(paste0("Scanned .counts folder '", getwd(), .Platform$file.sep, CONF.RNAseq$COUNTS_FOLDER, "'. ", nrow(FileNames.ERCC), " files found."))
#Creating Metadata by merging counts file list with table of corresponding covariates====
Metadata <- merge(LKBRNACovariates, FileNames.ERCC, by.x="ID_OX", by.y="shortName")
if ( !(nrow(Metadata) == nrow(FileNames.ERCC)) ) { 
    stop("Critical error: Unequal length of file list pre- and post covariate merge. Please check covariates are sufficient.") 
}
#Editing metadata, imputing subject age ====
Metadata[, CutRINOX:=cut(RIN_OX, 3) ] #/!\
setcolorder(Metadata, unique(c("ID_OX", "fileName","sample", "batch", colnames(Metadata))))
set(Metadata, j = Metadata[, which((colnames(Metadata) %in%
                                        c("Sex","Location","HasEM", "HasPFA", "HasInVivo" ,"RNA_ID", "Vol", "260/230", "260/280", "ng/ul", "Subject", "Group")))], 
    value=NULL)
Metadata[,S1Date:=as.Date.character(S1Date, "%d/%m/%Y")]
Metadata[,Born:=as.Date.character(Born, "%d/%m/%Y")]
#Imputing birthday 
Metadata[!is.na(Born),ImputationBDay := Born]
Metadata[is.na(Born) & Timepoint == "AD" & !is.na(S1Date),
         ImputationBDay := as.Date.character(paste0("01/05/", year(S1Date)-1), "%d/%m/%Y")]
#Assumption: All lambs are born on the first of May (middle of lambing season); all AD samples are at least one year old
#e.g. computing age in days as distance between S1date and first of May
Metadata[!is.na(DaysOld), ImputedAgeDays:=DaysOld] #If an age in days is known, we use the certain value of course
Metadata[is.na(DaysOld) & !is.na(S1Date) & !is.na(ImputationBDay), ImputedAgeDays:= as.integer(S1Date - ImputationBDay)]
Metadata[,GroupMean:=mean(DaysOld, na.rm = T), Timepoint] #We compute a group mean from the known and imputed values
Metadata[is.na(ImputedAgeDays), ImputedAgeDays:=as.integer(GroupMean)] #Those without a birtsday gain an approximate age
set(Metadata, j = Metadata[, which((colnames(Metadata) %in% c("S1Date", "ImputedBDay", "GroupMean", "BodyWght", "HeartWght",
                                                              "ImputationBDay", "Born", "DaysOld", "sample")))], value=NULL)
Metadata.PostPCA = Metadata[!( ID_OX %in% CONF.RNAseq$FILE_EXCLUDE_LIST )] #If removing BEFORE imputation, diff result!
tmsg(paste0("Excluding files \"", paste0(CONF.RNAseq$FILE_EXCLUDE_LIST, collapse = ", "), 
            " according to PCR results, remainder ", nrow(Metadata.PostPCA), " files."))
#EXPERIMENTAL====
SplineDF = as.data.table(Metadata.PostPCA[,ns(ImputedAgeDays, df = 3)]) #To compare to december's spline mode, df must be 3!
colnames(SplineDF)<- c("Spline1", "Spline2", "Spline3")
Metadata.nsTable = cbind(Metadata.PostPCA, SplineDF)
#Experimental model====
ExperimentalModel <- DESeqDataSetFromHTSeqCount(sampleTable = Metadata.nsTable,  design = CONF.RNAseq$EXPERIMENTAL_FORMULA)
ExperimentalModel <- ExperimentalModel[ rowSums(counts(ExperimentalModel)) > 1, ]
ExperimentalModel <- estimateSizeFactors(ExperimentalModel)
ExperimentalModel <- estimateDispersions(ExperimentalModel)
tmsg(paste0("Running Experimental model with formula", as.character(CONF.RNAseq$EXPERIMENTAL_FORMULA)[2]))
ExperimentalModel <- nbinomLRT(ExperimentalModel, reduced = CONF.RNAseq$REDUCED_FORMULA, maxit = CONF.RNAseq$MAX_WALD_ITERATIONS)

SplineModel = readRDS("Data/RNAseq/Output/SplineModel.ERCC.ImputeB4Exclude.2excl.maxit5000.batchCutRINOXnsImputedAgeDaysdf3.LRT.RDS")
# SplineResults1=GrabResult(SplineModel, "ns.ImputedAgeDays..df...3.3", comparisonStr = "LRT", Alpha = 0.05, Filter = TRUE)
# SplineResults2=GrabResult(SplineModel, "ns.ImputedAgeDays..df...3.2", comparisonStr = "LRT", Alpha = 0.05, Filter = TRUE)
# SplineResults3=results(SplineModel)
#TESTING====
#all(SplineResults1==SplineResults2)                     #Different, filtered and unfiltered. Same nrow() tho. 
#all(SplineResults1[,geneID]==SplineResults2[,geneID])   #TRUE = pvalues are constant
#Ergo, the p-value calculation does not change, but the log2FC does due to who is the reference value...

#DEEP HURTING====
#resutls() uses DEseq2:::getCoef(), which needs name; default is DESeq2:::lastCoefName(object), e.g. last formula term is default results level...
#So calls to results(ExperimentalModel) and results(SplineModel) should produce identical results IF spline being internal does not matter
SplineResults=results(SplineModel)
ExperimentalResults=results(ExperimentalModel)
all(SplineResults==ExperimentalResults) #TRUE TRUE TRUE TRUE TRUE NA

#SOME GRAPHS====