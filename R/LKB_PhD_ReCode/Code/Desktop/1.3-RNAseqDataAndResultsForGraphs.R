#Before running: check filepaths, e.g. for SplineModel import and PANTHER data to combine with the DEGs. 
#Manually insert headers into any PANTHER output files.
#CONFIG====
CONF.RNA4Graphs = list()
CONF.RNA4Graphs[["VERSION"]]		        = "0.1-Init"
CONF.RNA4Graphs[["EXPRESSED_MIN_COUNTS"]]	= 6
CONF.RNA4Graphs[["EXPRESSED_PCT_LIMIT"]]	= 0.5
#Libraries and functions====
library(data.table); library(DESeq2); library(stringr); 
require(rtracklayer); require(GenomicRanges)
source("Code/Desktop/functions.DESeq2.r")
source("Code/Desktop/functions.r")
#Dictionary <- fread("Data/biomart/BIOMART-OvAr3.1-GeneIDTxIDNameDescrLocation-FINAL.txt")
#Load in a model (which does not matter as counts are independent of model/question)====
SplineModel = readRDS("Data/RNAseq/Output/SplineModel.ERCC.ImputeB4Exclude.2excl.maxit5000.batchCutRINOXnsImputedAgeDaysdf3.LRT.RDS")
GroupModel  = readRDS("Data/RNAseq/Output/GroupModel.ERCC.ImputeB4Exclude.2excl.maxit5000.batchCutRINOXTimepoint.Wald.RDS")
#Generated in file 1.1

#Basic exports====
#The log2fc values of the above will differ depending on the df...x.x used, but pValues and pAdj are constant. hm.

DEGs.ERCC.SplineDF3 = GrabResult(DESeqObject = SplineModel, ContrastOrItem = "ns.ImputedAgeDays..df...3.3", 
                                 comparisonStr = "All", Alpha = 0.05, Filter = TRUE)
# DEGs.ERCC.SplineDF3.Top10Counts = DEGs.ERCC.SplineDF3[order(baseMean, decreasing = T)][1:10]
# DEGs.ERCC.SplineDF3.Top10Log2FC = DEGs.ERCC.SplineDF3[baseMean>500][order(Log2FCAbs, decreasing = T)][1:10]


DEGs.ERCC.1WAD      = GrabResult(DESeqObject = GroupModel, ContrastOrItem = c("Timepoint", "AD", "1W"), 
                                 comparisonStr = "1WAD", Alpha = 0.05, Filter = TRUE)
# DEGs.ERCC.1WAD.Top10Counts = DEGs.ERCC.1WAD[order(baseMean, decreasing = T)][1:10]
# DEGs.ERCC.1WAD.Top10Log2FC = DEGs.ERCC.1WAD[baseMean>500][order(Log2FCAbs, decreasing = T)][1:10]

#Generate objects neccesary for graphing====
#All this data is based on counts, so the model is irrelevant
#raw counts, allowing graphing per sample
Counts.ERCC    = DESeq2::counts(SplineModel, normalized=TRUE)
Counts.ERCC.M  = data.table(melt(Counts.ERCC))
colnames(Counts.ERCC.M) <- c("gene", "sample", "value")
Counts.ERCC.M[,Timepoint:=str_match(sample, "\\d{6}-\\d{2}_(1W|1M|3M|AD)")[,2]]

#Calculations 1====
SamplesPerGroup = unique(Counts.ERCC.M[,sample, Timepoint])[,.(MaxSamples=.N),Timepoint]

Counts.ERCC.M[,Is_Expr:=(value>=CONF.RNA4Graphs$EXPRESSED_MIN_COUNTS)] 
#'legacy' algo; since we do not have TPM values, we can't fully emulate GTEX' methodology

#In how many samples is any gene expressed at any timepoint?
ExpressPct=Counts.ERCC.M[,.N,.(gene, Timepoint, Is_Expr)]
ExpressPct=dcast(ExpressPct, gene+Timepoint~Is_Expr, value.var = "N")
colnames(ExpressPct) = c("geneID", "Timepoint", "samplesNoExpression", "samplesExpressed")
ExpressPct[is.na(samplesExpressed), samplesExpressed:=0]
ExpressPct[,samplesNoExpression:=NULL]
ExpressPct=merge(ExpressPct, SamplesPerGroup, by="Timepoint", all.x=T)
ExpressPct[,ExpressedRatio:=samplesExpressed/MaxSamples]
ExpressPct[,IsExpressed:=ExpressedRatio>CONF.RNA4Graphs$EXPRESSED_PCT_LIMIT]

#Which patterns does this gene expression have?
GeneExpression.S = ExpressPct[IsExpressed==T,.N,Timepoint]
GeneExpression.S[,Timepoint:=factor(Timepoint, levels=c("1W", "1M", "3M", "AD"), ordered = T)]

#Scaling counts for a heatmap
ScaledCounts = MakeNormalisedCountTable(SplineModel, ExpressPct)
ScaledCounts = ScaledCounts[,.(geneID, geneName, Sample, value, Timepoint)]
setnames(ScaledCounts, "Timepoint", "group")

#TPM calculation using PICARD average Insert length + 2x FastQC Read Length for Fragment length====
LKB.TPMCountTable = DESeq2::counts(SplineModel, normalized=FALSE)
LKB.TPMCountTable = LKB.TPMCountTable[!grepl("^ERCC", rownames(LKB.TPMCountTable)),]

#Load FastQC Data and calculate average Read length
SeqQCFiles <- list.files("Data/FastQC/RawData/", pattern = "*Sequence Length Distribution.txt", full.names = T)
SeqLenDistribution = data.table()
for(file in SeqQCFiles){
    FileName = str_match(file, "(WTCHG_\\d+_\\d+)_\\d_fastqc_.*\\.txt")[,2]
    tryCatch(expr = { QCData = fread(file)}, error = function (e) { QCData = data.table()  })
    QCData[,Sample:=FileName]
    SeqLenDistribution = rbind(SeqLenDistribution, QCData, fill=T)
} #for file in SeqQCFiles
rm(file, FileName, QCData, SeqQCFiles)

FilesToSamples = fread("Data/RNAseq/Input/2017-04-25_RNAFiles2.txt")
setkey(FilesToSamples, Filename)
setkey(SeqLenDistribution, Sample)
SeqLenDistribution = SeqLenDistribution[FilesToSamples]
#all(SeqLenDistribution[,Count==i.Count])
SeqLenDistribution[,i.Count:=NULL]
SeqLenDistribution = dcast.data.table(SeqLenDistribution, SampleID~., value.var = "Length", fun.aggregate = mean)
setnames(SeqLenDistribution, ".", "AverageReadLength")
setkey(SeqLenDistribution, "SampleID")
#It's 75 bp as ordered. But we checked anyway.

#Load PICARD data on mean insert length
InsertLengths = data.table()
PICARDFiles <- list.files("Data/PICARD/FragmentLengths/LKB_RNA/", pattern = "*.txt", full.names = T)
for(file in PICARDFiles){
    FileName = str_match(file, "^.*Sample(\\d{6}-\\d{2}_..).*\\.txt")[,2]
    tryCatch(expr = { PICARDData = fread(file, skip = 6, nrows = 1)}, error = function (e) { PICARDData = data.table()  })
    PICARDData[,Sample:=FileName]
    InsertLengths = rbind(InsertLengths, PICARDData, fill=T)
} #for file in InsertLengths
rm(file, FileName, PICARDFiles, PICARDData)
InsertLengths = InsertLengths[,.(Sample, MEAN_INSERT_SIZE)]
setkey(InsertLengths, "Sample")

InsertLengths = InsertLengths[SeqLenDistribution]
InsertLengths[,MeanFragmentLength := MEAN_INSERT_SIZE + 2 * AverageReadLength]

#Get gene lengths
#We need to remove ERCCs, see:
#Counts.ERCC[which(!(rownames(Counts.ERCC) %in% GeneLengthData[,`Gene stable ID`])),]

#from https://www.biostars.org/p/91218/#91233 - validate using 'in house' code
OvisAriesGTF = rtracklayer::import.gff("Data/biomart/Ovis_aries.Oar_v3.1.95.sorted.gtf", format="gtf", feature.type="exon")
OvisAriesGTF = GenomicRanges::reduce(split(OvisAriesGTF, elementMetadata(OvisAriesGTF)$gene_id))
OvisAriesGTF.UL = unlist(OvisAriesGTF, use.names=T)
elementMetadata(OvisAriesGTF.UL)$gene_id <- rep(names(OvisAriesGTF), elementNROWS(OvisAriesGTF))
OvisAriesGTF.UL =as.data.table(OvisAriesGTF.UL)
GeneLengthData = OvisAriesGTF.UL[,.(featureLength=sum(width)),gene_id]
rm(OvisAriesGTF, OvisAriesGTF.UL)

#GeneLengthData = fread("Data/biomart/parsedOvAr3.1GTF.tsv")[feature=="exon"]
#GeneLengthData[,source:=NULL]
#setorder(GeneLengthData, gene_id, exon_number)
# Inter range transformations (e.g. reduce()) transform all the ranges together as a set to produce a new set of ranges. 
# They return an object that is generally NOT parallel to the input object. 
# Those transformations are described in the inter-range-methods man page (see ?`inter-range-methods`).
#[...]
#reduce
#reduce first orders the ranges in x from left to right, then merges the overlapping or adjacent ones.
#order by gene_id, start, end
#for each gene_id:
# call unique() to exclude perfect overlaps
# then check the list of start, end pairs X, from i .. n, whether X[i].start <= X[i+1].end AND X[i].end >= X[i+1].start
## if TRUE, then make start = max(X[i].start, X[i+1].start), same for end
# then repeat with the NEW list until all are processed
#similar to leafcutter clustering algo
#...
#... but because there's no inherent value in suffering and GenomicRanges::reduce() already does this, why reinvent the wheel?
GeneLengthData = GeneLengthData[gene_id %in% rownames(LKB.TPMCountTable)]

stopifnot(length(rownames(LKB.TPMCountTable)) == length(GeneLengthData[,gene_id]))
stopifnot(all(rownames(LKB.TPMCountTable) == GeneLengthData[,gene_id]))

InsertLengths = InsertLengths[Sample %in% colnames(LKB.TPMCountTable)]
stopifnot(length(colnames(LKB.TPMCountTable)) == length(InsertLengths[,Sample])) #Can they be equal?
stopifnot(all(colnames(LKB.TPMCountTable) == InsertLengths[,Sample]))            #ARE they equal in order?

source("Code/counts_to_tpm.R") #includes code which removes genes with fragment length =< gene length
LKB.TPMCountTable = counts_to_tpm(LKB.TPMCountTable, GeneLengthData[,featureLength], InsertLengths[,MeanFragmentLength])
LKB.TPMCountTable = data.table(LKB.TPMCountTable, keep.rownames="geneID")

rm(Dictionary, CONF.RNA4Graphs, Counts.ERCC, DESeqResultToDataTable, GrabAllPossibleResults, GrabResult, GroupModel,
   MakeNormalisedCountTable, AskYN, completeDT, DoHeatmap, GenerateGeneSummaryPlots, give.n, LeafCutterContainerToResultsObj, se, 
   SplineModel, SplitDataTableWithMultiRows, tmsg, tstamp, ffwrite, GetAllResultsForVariable, InsertLengths, GeneLengthData, counts_to_tpm, 
   FilesToSamples, GroupNonOverlappingIntervals, SeqLenDistribution)
save(list = ls(), file = "Data/RNAseq/Output/LKB.ERCC.GeneSummaryObjects.rda") #Export to a single object

#Final size should be approx 19987 genes (which is nrow() of a matrix with rowSums() > 0).
#TPM table has less, likely requiring a greater RowSum?