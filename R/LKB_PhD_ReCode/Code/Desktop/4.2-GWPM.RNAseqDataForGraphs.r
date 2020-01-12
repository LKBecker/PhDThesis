#CONFIG====
CONF.GWPMRNA4Graphs = list()
CONF.GWPMRNA4Graphs[["VERSION"]]		        = "0.1-Init"
CONF.GWPMRNA4Graphs[["EXPRESSED_MIN_COUNTS"]]	= 6
CONF.GWPMRNA4Graphs[["EXPRESSED_PCT_LIMIT"]]	= 0.5
#Libraries and functions====
library(data.table); library(DESeq2); library(stringr); 
require(GenomicRanges); require(rtracklayer)
source("Code/Desktop/functions.DESeq2.r")
source("Code/Desktop/functions.r")
#Load====
GWPMGroupModel = readRDS("Data/RNAseq/Output/GroupGWPM.2excl.maxit5000.BatchGroup.Wald.RDS") 
DEGs.GWPM = GetAllResultsForVariable(Dataset = GWPMGroupModel, item = "Group", Alpha = 0.05, filterNonPass = TRUE)

#Generate objects neccesary for graphing====
#raw counts, allowing graphing per sample
Counts.GWPM    = DESeq2::counts(GWPMGroupModel, normalized=TRUE)
Counts.GWPM.M  = data.table(melt(Counts.GWPM))
colnames(Counts.GWPM.M) <- c("gene", "sample", "value")
Counts.GWPM.M[,Timepoint:=str_match(sample, "\\d{1,2}(C|HF|R)")[,2]]
Counts.GWPM.M[,Is_Expr:=(value>=CONF.GWPMRNA4Graphs$EXPRESSED_MIN_COUNTS)] 
#'legacy' algo; since we do not have TPM values, we can't fully emulate GTEX' methodology

SamplesPerGroup = unique(Counts.GWPM.M[,sample, Timepoint])[,.(MaxSamples=.N),Timepoint]

#In how many samples is any gene expressed at any timepoint?====
ExpressPct.GWPM=Counts.GWPM.M[,.N,.(gene, Timepoint, Is_Expr)]
ExpressPct.GWPM=dcast(ExpressPct.GWPM, gene+Timepoint~Is_Expr, value.var = "N")
colnames(ExpressPct.GWPM) = c("geneID", "Timepoint", "samplesNoExpression", "samplesExpressed")
ExpressPct.GWPM[is.na(samplesExpressed), samplesExpressed:=0]
ExpressPct.GWPM[,samplesNoExpression:=NULL]
ExpressPct.GWPM=merge(ExpressPct.GWPM, SamplesPerGroup, by="Timepoint", all.x=T)
ExpressPct.GWPM[,ExpressedRatio:=samplesExpressed/MaxSamples]
ExpressPct.GWPM[,IsExpressed:=ExpressedRatio>CONF.GWPMRNA4Graphs$EXPRESSED_PCT_LIMIT]

#Which patterns does this gene expression have?====
GeneExpression.S = ExpressPct.GWPM[IsExpressed==T,.N,Timepoint]
GeneExpression.S[,Timepoint:=factor(Timepoint, levels=c("C", "HF", "R"), ordered = T)]

#Scaling counts for a heatmap====
ScaledCounts.GWPM = MakeNormalisedCountTable(GWPMGroupModel, ExpressPct.GWPM, isGWPM = TRUE)
ScaledCounts.GWPM = ScaledCounts.GWPM[,.(geneID, geneName, Sample, value, Timepoint)]
setnames(ScaledCounts.GWPM, "Timepoint", "group")

#TPM calculation====
#Load FastQC Data and calculate average Read length
BHB.TPMCountTable = DESeq2::counts(GWPMGroupModel, normalized=FALSE)

SeqQCFiles <- list.files("Data/PICARD/FragmentLengths/BHB_RNA/", pattern = "*alignment.txt", full.names = T)
SeqLenDistribution = data.table()
for(file in SeqQCFiles){
    FileName = str_match(file, "(\\d{1,2}(C|R|HF)).*\\.txt")[,2]
    tryCatch(expr = { QCData = fread(file)}, error = function (e) { QCData = data.table()  })
    QCData[,Sample:=FileName]
    SeqLenDistribution = rbind(SeqLenDistribution, QCData, fill=T)
} #for file in SeqQCFiles
rm(file, FileName, QCData, SeqQCFiles)
SeqLenDistribution = SeqLenDistribution[CATEGORY=="PAIR", .(Sample, MEAN_READ_LENGTH)]
setkey(SeqLenDistribution, "Sample")

#Load PICARD data on mean insert length
InsertLengths = data.table()
PICARDFiles <- list.files("Data/PICARD/FragmentLengths/BHB_RNA/", pattern = "*Inserts.txt", full.names = T)
for(file in PICARDFiles){
    FileName = str_match(file, "(\\d{1,2}(C|R|HF)).*\\.txt")[,2]
    tryCatch(expr = { PICARDData = fread(file, skip = 6, nrows = 1)}, error = function (e) { PICARDData = data.table()  })
    PICARDData[,Sample:=FileName]
    InsertLengths = rbind(InsertLengths, PICARDData, fill=T)
} #for file in InsertLengths
rm(file, FileName, PICARDFiles, PICARDData)
InsertLengths = InsertLengths[,.(Sample, MEAN_INSERT_SIZE)]
setkey(InsertLengths, "Sample")

InsertLengths = InsertLengths[SeqLenDistribution]
InsertLengths[,MeanFragmentLength := MEAN_INSERT_SIZE + 2 * MEAN_READ_LENGTH]

#Get FEATURE (sum of union of exons) lengths
#Counts.ERCC[which(!(rownames(Counts.ERCC) %in% GeneLengthData[,`Gene stable ID`])),]
OvisAriesGTF = rtracklayer::import.gff("Data/biomart/Ovis_aries.Oar_v3.1.95.sorted.gtf", format="gtf", feature.type="exon")
OvisAriesGTF = GenomicRanges::reduce(split(OvisAriesGTF, elementMetadata(OvisAriesGTF)$gene_id))
OvisAriesGTF.UL = unlist(OvisAriesGTF, use.names=T)
elementMetadata(OvisAriesGTF.UL)$gene_id <- rep(names(OvisAriesGTF), elementNROWS(OvisAriesGTF))
OvisAriesGTF.UL =as.data.table(OvisAriesGTF.UL)
GeneLengthData = OvisAriesGTF.UL[,.(featureLength=sum(width)),gene_id]
rm(OvisAriesGTF, OvisAriesGTF.UL)

GeneLengthData = GeneLengthData[gene_id %in% rownames(BHB.TPMCountTable)]
stopifnot(length(rownames(BHB.TPMCountTable)) == length(GeneLengthData[,gene_id]))
stopifnot(all(rownames(BHB.TPMCountTable) == GeneLengthData[,gene_id]))

InsertLengths = InsertLengths[Sample %in% colnames(BHB.TPMCountTable)]
stopifnot(length(colnames(BHB.TPMCountTable)) == length(InsertLengths[,Sample])) #Can they be equal?
stopifnot(all(colnames(BHB.TPMCountTable) == InsertLengths[,Sample]))            #ARE they equal in order?

source("Code/counts_to_tpm.R")
BHB.TPMCountTable = counts_to_tpm(BHB.TPMCountTable, GeneLengthData[,featureLength], InsertLengths[,MeanFragmentLength])
BHB.TPMCountTable = data.table(BHB.TPMCountTable, keep.rownames="geneID")

#Export====
rm(Dictionary, CONF.GWPMRNA4Graphs, DESeqResultToDataTable, GrabAllPossibleResults, GrabResult,
   MakeNormalisedCountTable, AskYN, completeDT, DoHeatmap, GenerateGeneSummaryPlots, give.n, LeafCutterContainerToResultsObj, se, 
   SplitDataTableWithMultiRows, tmsg, tstamp, ffwrite, GetAllResultsForVariable, GeneExpression.S, SamplesPerGroup, Counts.GWPM,
   counts_to_tpm, InsertLengths, SeqLenDistribution, GroupNonOverlappingIntervals, GWPMGroupModel, GeneLengthData)
save(list = ls(), file = "Data/RNAseq/Output/GWPM.GeneSummaryObjects.rda") #Export to a single object
