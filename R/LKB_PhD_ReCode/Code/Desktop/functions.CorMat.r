library(data.table)

PrepareCorrelationCountMatrix <- function() {
    if (!exists('CountData')) { source("Code/Graphs.Setup.R") }
    if (class(CountData)=="sleuth"){
        Counts.Nt = data.table(kallisto_table(CountData))
        Counts.Nt[,target_id:=substr(target_id, 1, 18)]
        Counts.Nt = merge(Counts.Nt, TxToGene, by.x="target_id", by.y="Transcript stable ID", all.x=T)
        Counts.Nt = Counts.Nt[,.(sum_counts=sum(est_counts)), .(`Gene stable ID`, sample)]
        Counts.Nt = dcast(Counts.Nt, `Gene stable ID`~sample, value.var="sum_counts")
        rownames(Counts.Nt) = Counts.Nt[,`Gene stable ID`]
        Counts.Nt[,`Gene stable ID`:=NULL]
        
        Counts.Nt1 = as.matrix(Counts.Nt)
        rownames(Counts.Nt1)=rownames(Counts.Nt)
        Counts.Nt1 = Counts.Nt1[rowSums(Counts.Nt1)>0,] #important! Not removing 0 rows leads to an error wherein 
    }
    else if (class(CountData)=="DESeqDataSet") { 
        Counts.Nt1 = counts(CountData, normalized=T) 
    } else {
        stop("CountData does not seem to be either a sleuth or DESeq2 object. Cannot continue.")
    }
    Counts.Nt1
}

DoCorrelationMatrix <- function(GeneNameList, isIDList = FALSE, doAdvanced=TRUE) {
    #GeneNameList = c("ENSOARG00000003817"); isIDList = TRUE; doAdvanced = TRUE #works with two but not 1 source gene
    #GeneNameList = c("ENSOARG00000003817", "ENSOARG00000003818"); isIDList = TRUE; doAdvanced = TRUE #works with two but not 1 source gene
    CountMatrix = PrepareCorrelationCountMatrix()
    if (!isIDList) {
        GENE_ID_LIST <- unique(TxToGene[`Gene name` %in% GeneNameList, `Gene stable ID`])
        nFound <- length(GENE_ID_LIST)
        nQuery <- length(GeneNameList) 
        if (nFound < nQuery) { message(paste0("Warning: Could only find ENSEMBL ids for ", nFound, " of ", nQuery, " genes.\n", 
                                              "Submit a list of ENSEMBL IDs and set isIDList=TRUE to bypass.")) }
    } else { GENE_ID_LIST <- GeneNameList }
    
    #TODO: make a loop to allow query-query correlations to be seen!
    if(length(GENE_ID_LIST) == 1) {
        Counts.SourceGenes    <- t(CountMatrix[which(  rownames(CountMatrix) %in% GENE_ID_LIST) , , drop = FALSE]) 
    } else {
        Counts.SourceGenes    <- t(CountMatrix[which(  rownames(CountMatrix) %in% GENE_ID_LIST) , ]) 
    }
    # drop=FALSE stops one-row matrices turning into vectors, which register on dim() but fail corr.test() and cor()
    
    Counts.RemainingGenes <- t(CountMatrix[which(!(rownames(CountMatrix) %in% GENE_ID_LIST)) , ])
    if (doAdvanced == FALSE){
        CorrelationMatrix = cor(Counts.RemainingGenes, Counts.SourceGenes) # p-values?
    } else { CorrelationMatrix = psych::corr.test(Counts.RemainingGenes, Counts.SourceGenes, ci = FALSE, adjust="BH") }
    #make sure to turn ci=FALSE to -massively- speed up processing
    
    CorrelationMatrix
}
#> all(CorMat.Nt == CorrelationMatrix) #[1] TRUE -> Function encapsulation did not change output, function 'works'
#TODO: write test data

PsychCorMatToDT <-function(PsychObj) {
    Corrs = melt(PsychObj$r)
    colnames(Corrs)=c("SourceGene", "CorrelatingGene", "Correlation")
    pVals = melt(PsychObj$p)
    colnames(pVals)=c("SourceGene", "CorrelatingGene", "pValAdj")
    TEMP = as.data.table(merge(Corrs, pVals, by=c("SourceGene", "CorrelatingGene")))
    rm(Corrs, pVals)
    TEMP
}

CheckCorrelationMatrixThreshold <- function(CorrelationMatrix) {
    #CorrelationMatrix=BHB.MeanLength
    if ("psych" %in% class(CorrelationMatrix)) {
        TEMP = PsychCorMatToDT(CorrelationMatrix)
    } else { 
        TEMP = as.data.table(melt(CorrelationMatrix))
        colnames(TEMP)=c("SourceGene", "CorrelatingGene", "Correlation")
    }
    Results = rbindlist(lapply(seq(0.25, 1, 0.05), 
                               function(n) {
                                   if ("psych" %in% class(CorrelationMatrix)) {
                                       return(TEMP[abs(Correlation)>=n,.(nGenes=.N, CorrelationCutoff=n),pValAdj<0.05])
                                   } else { 
                                       TDAT = TEMP[abs(Correlation)>=n,.(nGenes=.N, CorrelationCutoff=n)]
                                       if (nrow(TDAT)==0 ) {return(list(0,0))} else { return(TDAT) }
                                   }
                               }
    ))
    TestPlot = ggplot(Results, aes(x=CorrelationCutoff, y=nGenes))+
        labs(x="Correlation cutoff (%)", y="Number of correlating genes")+scale_x_continuous(labels=scales::percent)
    if ("psych" %in% class(CorrelationMatrix)) { TestPlot = TestPlot+geom_col(aes(fill=pValAdj)) } else { TestPlot = TestPlot+geom_col() }
    return(list("data"=Results, "plot"=TestPlot))
}

CorrelationMatrixToResults <- function(CorrelationMatrix, CorrelationMinimum = CONF$COR_CUTOFF, MakeHumanReadable = TRUE) {
    if ("psych" %in% class(CorrelationMatrix)) { 
        CorMat.Nt2 = PsychCorMatToDT(CorrelationMatrix)
    } else {
        CorMat.Nt2 = as.data.table(CorrelationMatrix)
        CorMat.Nt2[,CorrelatingGene:=rownames(CorrelationMatrix)]
        CorMat.Nt2 = melt(CorMat.Nt2, id.vars="CorrelatingGene")
        setnames(CorMat.Nt2, x("variable", "value"), c("SourceGene", "Correlation"))
    } 
    CorMat.Nt2 = CorMat.Nt2[abs(Correlation)>=CorrelationMinimum] #i'm too lazy to learn how to filter matrix rows, fuck all y'all
    if (MakeHumanReadable == TRUE) {
        CorMat.Nt2 = merge(CorMat.Nt2, Dictionary, by.x="SourceGene", by.y="Gene stable ID", all.x=TRUE) #Dictionary is imported in Setup.Graphs.R; 
        #it's just an ensembl biomart download for Ovis Aris v3.1 containing Gene stable ID, Gene name, Gene description
        setnames(CorMat.Nt2, c("Gene description", "Gene name"), c("SourceDescription", "SourceName"))
        CorMat.Nt2 = merge(CorMat.Nt2, Dictionary, by.x="CorrelatingGene", by.y="Gene stable ID", all.x=TRUE)
        setnames(CorMat.Nt2, c("Gene description", "Gene name"), c("CorrelatingDescription", "CorrelatingName"))
        setcolorder(CorMat.Nt2, c(colnames(CorMat.Nt2)[colnames(CorMat.Nt2) %like% "Source"], "Correlation", 
                                  colnames(CorMat.Nt2)[colnames(CorMat.Nt2) %like% "Correlating"]))
    }
    unique(CorMat.Nt2)
}
