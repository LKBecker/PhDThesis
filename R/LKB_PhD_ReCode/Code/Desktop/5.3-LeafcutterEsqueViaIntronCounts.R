DATASET = "LKB"
VERSION = "v2.0"

library(data.table)
library(limma)
library(splines)
source("Code/Desktop/functions.R")
rm(AskYN, completeDT, DoHeatmap, GenerateGeneSummaryPlots, give.n, LeafCutterContainerToResultsObj, se, SplitDataTableWithMultiRows)
source("Code/Desktop/functions.DESeq2.r")
tmsg(sprintf("Preparing to run 5.3-LeafcutterEsqueViaIntronCounts.R, version %s, with dataset %s...", VERSION, DATASET))
#Loading in the gene dictionary====
Dictionary = fread("Data/biomart/BIOMART-OvAr3.1-GeneIDTxIDNameDescrLocation-FINAL.txt")
Dictionary[,`Transcript stable ID`:=NULL]
Dictionary = unique(Dictionary)
setkey(Dictionary, "Chromosome/scaffold name", "Gene start (bp)", "Gene end (bp)")

#Loading intron count matrix====
tmsg("Loading intron count matrix, created from .junc files using a custom C# leafcutter implementation...")
switch(DATASET,
    "LKB"= { Introns = fread("Data/leafcutter/juncs_perIntronCounts-Raw.out") },
    "BHB"= { Introns = fread("Data/leafcutter/BHB/juncs-BHB_perIntronCounts-Raw.out") }
)

colnames(Introns) = gsub("_Aligned.out.bam", "", colnames(Introns))
Introns = foverlaps(Introns, Dictionary[,.(`Chromosome/scaffold name`, `Gene start (bp)`, `Gene end (bp)`, `Gene stable ID`)], 
                          by.x = c("Intron.Chr", "Intron.Start", "Intron.Stop"), 
                          by.y = c("Chromosome/scaffold name", "Gene start (bp)", "Gene end (bp)"))
Introns[, geneID:= sprintf("%s-%d", `Gene stable ID`, seq_len(.N)), by = `Gene stable ID`]
IntronMap = Introns[,.(Intron.Chr, Intron.Start, Intron.Stop, `Gene stable ID`, geneID)]
IntronIDs = Introns[,geneID]
set(Introns, j=Introns[, which((colnames(Introns) %in% 
    c("Intron.Chr", "Gene start (bp)", "Gene end (bp)", "Intron.Start", "Intron.Stop", "geneID", "Gene stable ID")))], value=NULL)
Introns = as.matrix(Introns)
rownames(Introns) = IntronIDs
saveRDS(Introns, sprintf("Data/leafcutter/%s-Raw_Intron_Matrix.rds", DATASET))
Introns = Introns[which(rowSums(Introns)>50),]
#Introns = Introns[which(substr(rownames(Introns),1,2)!="NA"),]

#Loading covariates====
tmsg("Loading saved covariates...")
switch(DATASET, "LKB"={
        Metadata <- fread("Data/RNAseq/Input/20170911-133126_RNASeq-LKB-Covariates.csv", na.strings = "#N/A")
        Metadata[,batch:=factor(str_match(ID_OX, "\\d{4}(\\d{2})-\\d{2}")[,2])]
        Metadata[,SampleName:=sprintf("Sample-%s", ID_OX)]
        Metadata[,Timepoint:=factor(Group, levels=c("1W", "1M", "3M", "AD"))]
        Metadata = Metadata[,.(ID_OX, SampleName, Timepoint, batch, S1Date, Born, DaysOld, RIN_OX)]
        #Editing metadata, imputing subject age ====
        Metadata[, CutRINOX:=cut(RIN_OX, 3) ] #/!\
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
        set(Metadata, j = Metadata[, which((colnames(Metadata) %in% 
            c("S1Date", "ImputedBDay", "GroupMean",  "ImputationBDay", "Born", "DaysOld", "sample")))], value=NULL)
        setorder(Metadata, "ID_OX")
        stopifnot(all(Metadata[,SampleName]==colnames(Introns))) #
        limmaDesign = model.matrix(~batch+CutRINOX+Timepoint, data = Metadata)
        DESeqData.Introns.Group = DESeqDataSetFromMatrix(Introns, Metadata, design = ~batch+CutRINOX+Timepoint)
    }, 
    "BHB"= {
         Metadata <- fread("Data/RNAseq/Input/181025-160817_GWPM.Metadata.tsv")   
         Metadata[,Sample:=sprintf("Sample-%s", Sample)]
         Metadata[,Batch:=factor(Batch)]
         setorder(Metadata, Sample)
         Introns = Introns[,which( !(colnames(Introns) %in% c("Sample-12C", "Sample-13R")))]
         stopifnot(all(Metadata[,Sample]==colnames(Introns))) # DESeq2 does NOT check if matrix columns and metadata columns are equal... for some reason...
         #But they must be or you're gonna create wank.
         limmaDesign = model.matrix(~Batch+Group, data = Metadata)
         #Samples 12 and 13 were removed according to the PCA and thus there is no data available
         DESeqData.Introns.Group = DESeqDataSetFromMatrix(Introns, Metadata, design = ~Batch+Group)
    }
)
#Running RNA-seq on intron data====
tmsg("Running DESeq with 5000 iterations, Wald test...")
DESeqData.Introns.Group = DESeq(DESeqData.Introns.Group, test = "Wald")
DESeqData.Introns.Group = nbinomWaldTest(DESeqData.Introns.Group, maxit = 5000) #fair warning, maxit=5000 takes a long long long time
#default - 118, 5000 - 98
tmsg("Extracting results...")
switch(DATASET, 
       "LKB"={ IntronsResults.Group = GetAllResultsForVariable(DESeqData.Introns.Group, "Timepoint", 0.05, FALSE) },
       "BHB"={ IntronsResults.Group = GetAllResultsForVariable(DESeqData.Introns.Group, "Group", 0.05, FALSE) }
)
tmsg("Output to file...")
saveRDS(DESeqData.Introns.Group, sprintf("Data/leafcutter/leafcutterIntrons-%s-GroupModel-Unfiltered.rds", DATASET))

rm(IntronIDs, tmsg, tstamp, GrabResult, GetAllResultsForVariable, GrabAllPossibleResults)
save(list = ls(), file = sprintf("Data/leafcutter/%s-DESeq2wIntrons-v1-Unfiltered.rda", DATASET))
message("Complete.")