#RNA-seq models for further analysis
#libraries and functions====
library(cowplot)
library(data.table)
library(DESeq2)
library(extrafont)
library(ggplot2)
library(ggsignif)
library(RColorBrewer)
library(splines)
library(stringr)

source("Code/Desktop/functions.r")
source("Code/Desktop/functions.DESeq2.r")
#Theme====
TextSize = 10 #10 for thesis, 20+ for presentation
geom.text.size =  TextSize #/ ggplot2::.pt #geom_text uses "mm", ggplot theme uses point... yikes.
geom.text.size.small = geom.text.size / 2
theme_set(theme_cowplot(line_size = 1))
theme_update( text = element_text(family = windowsFont("Arial"), size = geom.text.size),  line = element_line(size=unit(0.7, "cm")), 
              legend.text.align = 1, plot.margin=margin(10,2,-15,0, "pt")
)
theme_update(
    axis.ticks.length=unit(-0.2, "cm"), strip.background = element_rect(fill = NA, colour = "#000000", linetype = "solid"),
    axis.text.x = element_text(family = "Arial", size = geom.text.size, margin=margin(8,0,0,0,"pt")), 
    axis.text.y = element_text(family = "Arial", size = geom.text.size, margin=margin(0,10,0,0,"pt")), 
    strip.text = element_text(family = "Arial", size = geom.text.size, margin = margin(2,2,2,2, "pt"))
)
#script----
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
    RNAModelOutputName = paste(paste0(Filename, ModelName, "Model.ERCC.ImputeB4Exclude"), NExcluded, 
                               MaxIter, Model, TestUsed, "RDS", sep = ".") 
    RNAModelOutputName
}
#Config====
CONF.RNAseq = list()
CONF.RNAseq[["VERSION"]]		        = "1.0-Recode"
CONF.RNAseq[["FILE_EXCLUDE_LIST"]]	    = c("161206-07_3M", "161207-05_3M")
CONF.RNAseq[["MIN_FDR"]]				= 0.05
CONF.RNAseq[["MIN_LFC"]]				= 0.0
CONF.RNAseq[["MAX_WALD_ITERATIONS"]]    = 5000 #do not set to 999999, it takes four hours and still doesn't converge
CONF.RNAseq[["COUNTS_FOLDER"]]          = "Data/RNAseq/Input/CountFiles.ERCC/"
CONF.RNAseq[["SPLINE_FORMULA"]]         = formula(~batch+CutRINOX+ns(ImputedAgeDays, df=3))
CONF.RNAseq[["GROUP_FORMULA"]]          = formula(~batch+CutRINOX+Timepoint)
CONF.RNAseq[["REDUCED_FORMULA"]]        = formula(~batch+CutRINOX)
#Loading covariates====
tmsg(paste0("Preparing to run 1-GenerateRNASeqModels.R, version ", CONF.RNAseq$VERSION, "..."))
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
#Imputing various covariates from the group mean: #Currently unused
# Metadata[,GroupMean:=mean(HeartWght, na.rm = T), Timepoint]
# Metadata[is.na(HeartWght), HeartWght:=GroupMean]
# Metadata[,GroupMean:=mean(BodyWght, na.rm = T), Timepoint]
# Metadata[is.na(BodyWght), BodyWght:=GroupMean]
#Imputing birthday 
Metadata[!is.na(Born),ImputationBDay := Born]
Metadata[is.na(Born) & Timepoint == "AD" & !is.na(S1Date),
         ImputationBDay := as.Date.character(paste0("01/05/", year(S1Date)-1), "%d/%m/%Y")]
#Assumption: All lambs are born on the first of May (middle of lambing season); all AD samples are at least one year old
#e.g. computing age in days as distance between S1date and first of May
Metadata[!is.na(DaysOld), ImputedAgeDays:=DaysOld] #If an age in days is known, we use the certain value of course
Metadata[is.na(DaysOld) & !is.na(S1Date) & !is.na(ImputationBDay), ImputedAgeDays:= as.integer(S1Date - ImputationBDay)]
Metadata[,GroupMean:=mean(DaysOld, na.rm = T), Timepoint] #We compute a group mean from the known and imputed values
Metadata[is.na(ImputedAgeDays), ImputedAgeDays:=as.integer(GroupMean)] #Those without a birthday gain an approximate age

#Plot for Methods section====
BasePlot = ggplot(Metadata, aes(x=Timepoint, y=ImputedAgeDays))+geom_boxplot(outlier.shape = NA)+geom_jitter(aes(color=Timepoint), data = Metadata[!is.na(DaysOld)], width=0.1)+
    geom_jitter(data = Metadata[is.na(DaysOld)], width=0.1, color="black", shape=2)
Plot1 = BasePlot+coord_cartesian(xlim=c(1, 4), ylim = c(0, 120))+theme(legend.position="none", plot.margin = margin(10, 0, -5, 5))+ylab("")+xlab("")+
    scale_y_continuous(breaks = seq(0, 120, length.out = 4), expand = c(0,0))
Plot2 = BasePlot+coord_cartesian(xlim=c(1, 4), ylim = c(370,470))+
    theme(legend.position="none", axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
          axis.title.y = element_blank(), plot.margin = margin(5, 0, 10, 22))+scale_y_continuous(breaks = seq(370, 470, length.out = 5), expand = c(0,0))
Fig2.4 = plot_grid(Plot2, Plot1, ncol = 1, axis = "l")+draw_label("Age (days)", angle = 90, 0.02, 0.5, size = geom.text.size)

set(Metadata, j = Metadata[, which((colnames(Metadata) %in% c("S1Date", "ImputedBDay", "GroupMean", "BodyWght", "HeartWght",
                                                              "ImputationBDay", "Born", "DaysOld", "sample")))], value=NULL)

#PCA====
DESeqData.Raw = DESeqDataSetFromHTSeqCount(sampleTable = Metadata,  design = ~1) #No design needed for PCA
DESeqData.Raw.VST   <- vst(DESeqData.Raw) #Variance-Stabilising Transformation
PCAData = plotPCA(DESeqData.Raw.VST, intgroup="Timepoint", ntop=20000, returnData=T)
setDT(PCAData, keep.rownames = "RNAsample")
PCAData[inrange(PC2, mean(PC2)-sd(PC2)*1.5, mean(PC2)+sd(PC2)*1.5), name:=NA]

Fig2.5 = ggplot(PCAData, aes(x=PC1, y=PC2, color=Timepoint))+geom_point()+geom_label(aes(label=name), nudge_x=5, nudge_y=5)+
    theme(legend.position = "bottom", plot.margin = margin(5,5,10,5))+coord_cartesian(xlim = c(-60, 60)) 
Fig2.5
rm(DESeqData.Raw, DESeqData.Raw.VST)

#Excluding according to PCA====
Metadata.PostPCA = Metadata[!( ID_OX %in% CONF.RNAseq$FILE_EXCLUDE_LIST )] #If removing BEFORE imputation, diff result!
tmsg(paste0("Excluding files \"", paste0(CONF.RNAseq$FILE_EXCLUDE_LIST, collapse = ", "), 
            " according to PCR results, remainder ", nrow(Metadata.PostPCA), " files."))
#Final Metadata====
#ffwrite(Metadata.PostPCA[,.(ID_OX, Timepoint, RIN_OX, CutRINOX, ImputedAgeDays, batch)], "clipboard") 
#EXPERIMENTAL====
SplineDF = as.data.table(Metadata.PostPCA[,ns(ImputedAgeDays, df = 3)])
colnames(SplineDF)<- c("Spline1", "Spline2", "Spline3")
Metadata.nsTable = cbind(Metadata.PostPCA, SplineDF)

#Spline Models====
SplineModel <- DESeqDataSetFromHTSeqCount(sampleTable = Metadata.PostPCA,  design = CONF.RNAseq$SPLINE_FORMULA)
SplineModel <- SplineModel[ rowSums(counts(SplineModel)) > 1, ]
SplineModel	<- estimateSizeFactors(SplineModel)
SplineModel	<- estimateDispersions(SplineModel)
tmsg(paste0("Running spline model with formula ", as.character(CONF.RNAseq$SPLINE_FORMULA)[2]))
SplineModel	<- nbinomLRT(SplineModel, reduced = CONF.RNAseq$REDUCED_FORMULA, maxit = CONF.RNAseq$MAX_WALD_ITERATIONS)

#Group Models====
GroupModel <- DESeqDataSetFromHTSeqCount(sampleTable = Metadata.PostPCA,  design = CONF.RNAseq$GROUP_FORMULA)
GroupModel <- GroupModel[ rowSums(counts(GroupModel)) > 1, ]
GroupModel <- estimateSizeFactors(GroupModel)
GroupModel <- estimateDispersions(GroupModel)
tmsg("Running group model...")
GroupModel <- nbinomWaldTest(GroupModel, maxit = CONF.RNAseq$MAX_WALD_ITERATIONS)

#Export to RDS====
tmsg("Exporting to RDS...")
saveRDS(SplineModel, MakeRNAModelOutputname("Spline", "LRT"))
saveRDS(GroupModel, MakeRNAModelOutputname("Group", "Wald"))