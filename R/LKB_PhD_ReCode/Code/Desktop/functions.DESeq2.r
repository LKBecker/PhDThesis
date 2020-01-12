require(data.table)
require(DESeq2)
require(aroma.light)
require(GenABEL); 

#Almost the same as GetAllResults for Variable but a) with 'translation', b) one defined contrast only
Dictionary <- fread("Data/biomart/BIOMART-OvAr3.1-GeneIDNameDescrLocation-FINAL.txt")

DESeqResultToDataTable<-function(DESeqResultObject, FDRLimit=CONF$MIN_FDR, Filter=TRUE){
    if (! "DESeqResults" %in% class(DESeqResultObject)) { stop("Error in DESeqResultToDataTable(): Object is not of class DESeqResults.") }
    if (Filter==TRUE) { DESeqResultObject 	<- DESeqResultObject[which(DESeqResultObject$padj < FDRLimit),] }
    ResultsDT		<- data.table(geneID=rownames(DESeqResultObject), pValue=DESeqResultObject$pvalue, pValAdjust= DESeqResultObject$padj,
                             baseMean=DESeqResultObject$baseMean, log2FC=DESeqResultObject$log2FoldChange, Log2FCStdErr=DESeqResultObject$lfcSE,
                             t_stat=DESeqResultObject$stat)
    ResultsDT[,Log2FCAbs:=abs(log2FC)]
    ResultsDT[is.na(pValAdjust),pValAdjust:=1]
    ResultsDT[,foldChange:=2^log2FC]
    #ResultsDT[log2FC>=0,foldChange:=2^log2FC]
    #ResultsDT[log2FC<0,foldChange:=-1/(2^log2FC)]
    return(ResultsDT)
}

#Different results for 
#CountData <- readRDS("Data/RDS_DESeq2Data_batchCutRINOX-Group_2FileRemoved-7NonConverge.RDS")
#a=GetAllResultsForVariable(CountData, "Group", 0.05, TRUE)[Comparison == "Group|1M|1W"] #1W1M = 1002
#b=DoDESeqResult(c("Group", "1M", "1W"), DESeqObject = CountData, 0.05, TRUE)#[,.N] #1109, same as used for older graphs
GetAllResultsForVariable <- function(Dataset, item, Alpha=CONF$MIN_FDR, filterNonPass=TRUE){
    if(!(item %in% colnames(colData(Dataset)))) { stop(paste0("Item has to be one of ", paste(colnames(colData(Dataset)), collapse = ", "))) }
    stopifnot("DESeqDataSet" %in% class(Dataset))
    
    Metadata <- as.data.table(colData(Dataset))
    switch (class(Metadata[[item]]),
            "factor" = {
                #Hacky implementaton of expand.grid without duplicates or backwards tests
                Grid = data.table()
                FacLevels = as.character(levels(Metadata[[item]]))
                while (length(FacLevels) > 1) {
                    CurLevel = FacLevels[1]
                    FacLevels = FacLevels[2:length(FacLevels)]
                    Grid = rbind(Grid, data.table(CurLevel, FacLevels), fill=T)
                }
                ResultItem = data.table()
                for (i in seq(1, nrow(Grid))){
                    message(paste0("Testing factor '", item, "', class '", class(item), 
                                   "' in contrast ' c( ", item, " ", Grid[i,2], " ", Grid[i,1], " )'...")
                    ) 
                    Contrast = c(item, as.character(Grid[i,2]), as.character(Grid[i,1]))
                    SubSubtable = GrabResult(Dataset, Contrast, paste(Contrast, collapse = "|"), Alpha=Alpha, Filter = filterNonPass)
                    ResultItem = rbind(ResultItem, SubSubtable, fill=TRUE)
                }
            }, #factor
            "numeric" = { 
                ResultItem = GrabResult(Dataset, item, paste(item, "numeric", "numeric", collapse = "|"), Alpha = Alpha, Filter = filterNonPass) 
            }, 
            "integer" = { 
                ResultItem = GrabResult(Dataset, item, paste(item, "integer", "integer", collapse = "|"), Alpha = Alpha, Filter = filterNonPass) 
            }
    )#switch class Metadata
    return(ResultItem)
}

GrabResult <- function(DESeqObject, ContrastOrItem, comparisonStr, Alpha=CONF$MIN_FDR, Filter=TRUE){
    if (!class(DESeqObject)[1]=="DESeqDataSet") { stop("DESeqObject must be of class DESeqDataSet to retrieve results") }
    if (length(ContrastOrItem)==1) {
        if (! ContrastOrItem %in% resultsNames(DESeqObject)) { stop(paste0("ERROR: name must be one of [", paste0(resultsNames(DESeqObject), collapse = ", "), "]")) }
        ResultObj = DESeq2::results(DESeqObject, alpha=Alpha, name = ContrastOrItem)
        ResultObj = DESeqResultToDataTable(ResultObj, Alpha, Filter = Filter)
        ResultObj[,CompVar:=ContrastOrItem]
    } else if (length(ContrastOrItem)==3) {
        ResultObj = DESeq2::results(DESeqObject, alpha=Alpha, contrast = ContrastOrItem)
        ResultObj = DESeqResultToDataTable(ResultObj, Alpha, Filter = Filter)
        ResultObj[,CompVar:=ContrastOrItem[1]]
    } else { stop("ContrastOrItem must be of format 'Parameter' or 'c(Parameter, Denominator, Numerator)'") }
    
    ResultObj[,Comparison:=comparisonStr]
    ResultObj <- merge(ResultObj, Dictionary, by.x="geneID", by.y="Gene stable ID", all.x=T)
    TranscriptCols = which(tolower(colnames(ResultObj)) %like% "transcript stable id")
    if (length(TranscriptCols) > 0) {
        message("Transcript column detected; removing and running unique()...")
        set(ResultObj, j = TranscriptCols, value = NULL)
        ResultObj = unique(ResultObj)
    }
    ResultObj[,`Gene description` := gsub("\\[Source\\:.*\\]$", "", `Gene description`)]
    
    return(ResultObj)
}

GrabAllPossibleResults<-function(Dataset, Filter=TRUE){
    if (!("DESeqDataSet" %in% class(Dataset))) {stop(paste0("Error: '", Dataset, "' is not a DESeqDataSet object."))  }
    Metadata <- as.data.table(colData(Dataset))
    Formula <- all.names(design(Dataset))
    Formula <- Formula[!(Formula %in% c("+", "~"))]
    if (":" %in% Formula) { 
        warning("Formula features an interaction term. GrabAllPossibleResults is not designed to retrieve interaction term results or LRT data. ",
                "Please process these manually.") 
    }
    #full-auto extraction of all the things
    HugeTable <- data.table()
    for (item in Formula){
        if (!(item %in% colnames(Metadata))) {
            message(paste0("'", item, "' is not part of colnames() for the supplied item. It may be a computed variable. Please retrieve results manually."))
            next
        }
        ResData = GetAllResultsForVariable(Dataset, item, filterNonPass = Filter)
        HugeTable = rbind(HugeTable, ResData, fill=TRUE)
    } #for item in Formula
    message("Collection of all Wald contrasts complete. Please take care analysing these results - did you ask the right question? ",
            "(alternatives are e.g. an LRT removing more than one term)")
    return(HugeTable)
} # GrabAllPossibleResults

MakeNormalisedCountTable <- function(CountData, ExpressPct, isGWPM=FALSE){
    #CountData=GWPMGroupModel; ExpressPct=ExpressPct.GWPM; isGWPM = TRUE
    #CountData=SplineModel; ExpressPct=ExpressPct;
    GeneExpression = dcast.data.table(ExpressPct, geneID~Timepoint, value.var = "IsExpressed")
    if (isGWPM == FALSE) {
        GeneExpression[,Pattern:=paste0(`1W`*1, `1M`*1, `3M`*1, AD*1)]
    } else {
        GeneExpression[,Pattern:=paste0(`C`*1, `R`*1, `HF`*1)]
    }
    GeneExpression = GeneExpression[,.(geneID, Pattern)]
    
    #New Method: Quantile normalisation per sample, making them comparable
    #ScaledCounts = aroma.light::normalizeQuantileRank(log(t(counts(CountData))+1), robust = TRUE) 
    ScaledCounts = log(counts(CountData)+1)
    ScaledCounts = ScaledCounts[rowSums(ScaledCounts[,-1])>0,]
    
    NQRankCounts = matrix(nrow = NROW(ScaledCounts), ncol = 0)
    
    if (isGWPM == FALSE) { PossTimes = c("1W", "1M", "3M", "AD") } else { PossTimes = c("C", "R", "HF") }
    
    for (Timepoint in PossTimes) {
        columns = which(colnames(ScaledCounts) %like% Timepoint)
        NQRankCounts = cbind(aroma.light::normalizeQuantileRank(ScaledCounts[,columns], robust=TRUE), NQRankCounts)
    }
    
    a = apply(NQRankCounts, 1, function(x)length(unique(x))>2)
    NQRankCounts.Shadow = NQRankCounts[!a,]
    NQRankCounts = NQRankCounts[a,] #TODO: Keep in different matrix, rbind() back in??
    
    ScaledCounts <- apply(NQRankCounts, 1, function(x) {
        rt <- rntransform(x)
        names(rt) <- names(x)
        rt
    })
    
    ScaledCounts <- t(ScaledCounts)
    ScaledCounts.M = melt(data.table(ScaledCounts, geneID=rownames(ScaledCounts)), id.vars = "geneID")
    ScaledCounts.M = merge(ScaledCounts.M, unique(Dictionary[,.(`Gene stable ID`, `Gene name`)]), 
                           by.x="geneID", by.y="Gene stable ID")
    setnames(ScaledCounts.M, c("variable", "Gene name"), c("Sample", "geneName"))
    Metadata = as.data.table(colData(CountData))
    Metadata[,Sample:=rownames(colData(CountData))]
    
    if (isGWPM == FALSE) { Metadata = Metadata[,.(Sample, ImputedAgeDays, Timepoint)]
    } else { Metadata = Metadata[,.(Sample, Timepoint)] }
    
    ScaledCounts.M = merge(ScaledCounts.M, Metadata, by="Sample")
    ScaledCounts.M[geneName=="", geneName:=geneID] #Ensures there will always be a label available under geneName
    ScaledCounts.M=merge(ScaledCounts.M, GeneExpression, by="geneID", all.x=T) #Transfers information on expression pattern
    
    if (isGWPM == FALSE) { ScaledCounts.M[,Sample:=factor(Sample, levels=unique(ScaledCounts.M[order(ImputedAgeDays), Sample]))]
    } else { ScaledCounts.M[,Sample:=factor(Sample, levels=unique(ScaledCounts.M[order(Timepoint), Sample]))] }

    rm(NQRankCounts, ScaledCounts)
    ScaledCounts.M
}
