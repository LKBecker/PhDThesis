#WARNING: REQUIRES GGSIGNIF v off GitHub to use annotation_df for faceted plots
library(ggplot2)
library(ggsignif)
library(cowplot)
library(data.table)
library(stringr)
library(psych)

if(!exists("Dictionary")) { Dictionary = fread("Data/biomart/BIOMART-OvAr3.1-GeneIDNameDescrLocation-FINAL.txt") }

ffwrite <- function(OBJECT, FileName=NULL, Folder="./", ...) {
  if (is.null(FileName)) { FileName = "QuickExport" }
  if ( !(substr(Folder, nchar(Folder), nchar(Folder)) == .Platform$file.sep) ) { Folder = paste0(Folder, .Platform$file.sep) }
  if (!FileName=="clipboard") { FileName = paste0(Folder, format(Sys.time(), "%y%m%d-%H%M%S_"), FileName, ".tsv") }
  write.table(OBJECT, file = FileName, quote=F, row.names = F, col.names = (length(colnames(OBJECT)) > 0),  sep="\t", ...)
}

tstamp <- function(DoFileStamp=F) { if(DoFileStamp) {format(Sys.time(), "%y%m%d-%H%M%S_") } else { format(Sys.time(), "[%H:%M:%S] -- ") } }

tmsg <- function(message) { message(paste0(tstamp(), message)) }

se = function(x) sd(x)/sqrt(length(x))

AskYN <- function(question){
	while(TRUE){
		answer <- readline(paste(question, "(y/n): "))
		if (answer %in% c('y', 'Y', 'yes')) { return(TRUE);}
		if (answer %in% c('n', 'N', 'no')) { return(FALSE);}
		message("Please reply with 'y' or 'n'.")
	}
}

SplitDataTableWithMultiRows<-function(DataTable, TargetColumnIndex, Separator=","){
  if(!("data.table" %in% class(DataTable))){ DT<- data.table(DataTable) }
  else { DT <- copy(DataTable) }
  if(class(TargetColumnIndex)=="character") {
    if (!TargetColumnIndex %in% colnames(DT)) { stop(paste0("Error: '", TargetColumnIndex,
                                                            "' is not a valid column in the table submitted to SplitDataTableWithMultiRows")) }
    TargetColumnIndex <- which(colnames(DT) == TargetColumnIndex)
  }
  DTColOrder <- copy(colnames(DT))
  TargetCol <- DTColOrder[TargetColumnIndex]
  DTF<-data.frame(DT)
  SplitRowData<-strsplit(as.character(DTF[,TargetColumnIndex]), Separator, fixed=TRUE)
  SplitRowData<-lapply(SplitRowData, function(x){if(length(x)==0){ return("")} else { return(x)}})
  if(TargetColumnIndex != 1){
    RemainingCols <- DTColOrder
    RemainingCols <- RemainingCols[-TargetColumnIndex]
    setcolorder(DT, c(TargetCol, RemainingCols))
  }
  nreps<-sapply(SplitRowData, function(x){max(length(x), 1)})
  out<-data.table(TEMPCOL= unlist(SplitRowData), DT[ rep(1:nrow(DT), nreps) , -1, with=F] )
  setnames(out, "TEMPCOL", TargetCol)
  if(TargetColumnIndex != 1){ setcolorder(out, DTColOrder) }
  rm(DT, DTF, DTColOrder, TargetCol, SplitRowData, nreps)
  return(out)
}

completeDT <- function(DT, cols, defs = NULL){
    mDT = do.call(CJ, c(DT[, ..cols], list(unique=TRUE))) 
    #Creates a data.table to Join, over cols, with Unique values, and all value combinations (silent repetition?)
    res = DT[mDT, on=names(mDT)] #joins by names
    if (length(defs)) { res[, names(defs) := Map(replace, .SD, lapply(.SD, is.na), defs), .SDcols=names(defs)] } #idk
    #seems to use a call to Map() to call replace() on all .SDs where is.na, replaceing the present value with that in defs?
    #whilst using .SDcols to restrict the columns being edited to those selected
    return(res)
} 

give.n <- function(x){ return(c(ymin = median, label = length(x))) }

GroupNonOverlappingIntervals <- function(DT, start="start", end="end"){
    #DT=GeneMap;start="i.start";end="i.end"
    #DT=SignifList[geneName=="DYSF"];start="xmin_d";end="xmax_d"
    DT2 = copy(DT) #unless this is done early, somehow, even using return(), nothing comes back... unless captured in a variable.
    setnames(DT2, c(start, end), c("TEMPSTART", "TEMPEND"))
    if (nrow(DT2)==0) {          return(data.table()) }
    Coord_Groups = data.table()
    
    if (nrow(DT2)==1){ #just one item, so no need to actually look for overlap(s)
        set(DT2, j = "NoOverlapGroup", value = 1) 
        DT2 = DT2[,.(TEMPSTART, TEMPEND, NoOverlapGroup)] #ensure all outputs are of the same format!
        setnames(DT2, c("TEMPSTART", "TEMPEND"), c(start, end))
        return(DT2) 
    }
    
    DT2 = DT2[,.(TEMPSTART, TEMPEND)]
    currentGroup = 1
    Current_Item = DT2[1]
    DT2 = DT2[-1]
    while(nrow(DT2)>0){
        setkey(Current_Item, TEMPSTART, TEMPEND)
        setkey(DT2, TEMPSTART, TEMPEND)
        Overlap = foverlaps(DT2, Current_Item)
        if (nrow(Overlap[is.na(TEMPSTART)])>0){
            Current_Item = rbind(Current_Item, Overlap[is.na(TEMPSTART), .(TEMPSTART=i.TEMPSTART, TEMPEND=i.TEMPEND)][1]) 
            #grab first non-overlapping, append to current group - gotta rename as foverlaps() prepends "i." to your variable for overlap
            DT2 = DT2[!Current_Item]
        } else {
            Current_Item[,NoOverlapGroup:=currentGroup]
            Coord_Groups = rbind(Coord_Groups, Current_Item)
            Current_Item[,NoOverlapGroup:=NULL]
            DT2 = DT2[!Current_Item]
            currentGroup = currentGroup + 1
            Current_Item = DT2[1]
            DT2 = DT2[-1]
        }
    }
    Current_Item[,NoOverlapGroup:=currentGroup]
    Coord_Groups = rbind(Coord_Groups, Current_Item)
    rm(Current_Item, Overlap, currentGroup)
    setnames(Coord_Groups, c("TEMPSTART", "TEMPEND"), c(start, end))
    return(Coord_Groups)
}

LeafCutterContainerToResultsObj <- function(LeafRData){
    DataObj = load(LeafRData)
    OutputObj = list()
    OutputObj$clusters = as.data.table(DataObj$clusters)
    OutputObj$clusters[,gene:=str_remove_all(gene, "\\<i\\>|\\<\\/i\\>")]
    OutputObj$clusters[,chr:=str_match(coord, "^chr(\\d{1,2}|[X,Y,MT])")[,2]]
}

DoHeatmap<-function(GeneIDList){
  if (!exists('ScaledCounts.M')) { stop("Missing required data: ScaledCounts. Check load() calls and constructor scripts.") }
  GenesData = ScaledCounts.M[geneID %in% GeneIDList] #TODO get genes from peers
  GenesData[geneID %in% Spline.Full[,unique(SheepID)], geneName:=paste0("*", geneName)]
  GenesData[,geneName := factor(geneName, levels = unique(GenesData[order(geneName),geneName]))]
  GenesHeatmapPlot = HEATMAP + geom_raster(data = GenesData)+geom_vline(xintercept = c(7.5, 14.5, 20.5))+
    geom_hline(yintercept = length(GenesData[substr(geneName, 1,1)=="*", unique(geneName)])+0.5, size = 1)
  GenesHeatmapPlot
}

GenerateGeneSummaryPlots <- function(GeneList, isNameList=TRUE, HeatmapThreshold=5, forceUserNames=FALSE, Dataset="LKB.ERCC", Model="Spline",
                                     CrossbarFun=median, UseCommonScale=FALSE, MakeHeatmap=FALSE, UniqueComparisonSymbols=FALSE, PlotTPMs=TRUE){
    #GeneList =c("PRKAA2"="ENSOARG00000007851", "PRKAB1"="ENSOARG00000001411", "PRKAG3"="ENSOARG00000019718","PRKAR2A"="ENSOARG00000016132", "CAMK2B"="ENSOARG00000015725", "CAMK2D"="ENSOARG00000018416");isNameList = TRUE; HeatmapThreshold=5; forceUserNames=FALSE; Dataset = "ERCC+GWPM"; Model="Spline"; UseCommonScale=FALSE;UniqueComparisonSymbols=FALSE;PlotTPMs=TRUE
    #GeneList ="ENSOARG00000016226";isNameList = FALSE; HeatmapThreshold=5; forceUserNames=FALSE; Dataset = "ERCC+GWPM"; Model="Spline"; UseCommonScale=FALSE;UniqueComparisonSymbols=FALSE;PlotTPMs=TRUE
    
    #TODO
    #Edge case on SignifList y_position calculation when there are no TPMs - insert geom_text "No TPM could be calculated"?
    #Add dummy entry at 0,0 to generate a scale, to be abel to display geom_signif.
    #GeneList =c("PPP1R11", "PPP1R1A"); isNameList=TRUE;HeatmapThreshold=5; forceUserNames=FALSE; Dataset = "LKB.ERCC"; Model="Spline"; UseCommonScale=FALSE;UniqueComparisonSymbols=FALSE;PlotTPMs=TRUE
    
    #defensive code ----
    if(!exists('TxCpy')) { 
        TxCpy = fread("Data/biomart/BIOMART-OvAr3.1-GeneIDNameDescrLocation-FINAL.txt") 
        setnames(TxCpy, c("Gene stable ID", "Gene name"), c("geneID", "geneName"))
    }
    
    if( length(GeneList)==0 ) { stop("DoAGeneSummaryThing(): No gene names or IDs were supplied (GeneList is length 0). Please check input and retry.") }
    
    HasNames = !(is.null(names(GeneList)))
    if(!HasNames & forceUserNames) { stop("ForceUserNames of DoAGeneSummaryThing() is set to TRUE but supplies GeneList does not appear to be named.") }
    
    if( length(GeneList)!=length(unique(GeneList)) ) { 
        warning(paste0("There are ", length(GeneList)-length(unique(GeneList)), " redundant items in GeneList. Removing duplicates."))
        GeneList = unique(GeneList)
    }
    
    if (  all(unique(substr(GeneList, 1,7))=="ENSOARG") ) { isNameList = FALSE  } # if there are only geneIDs, it cannot be a list of gene names
    
    if ("ENSOARG" %in% unique(substr(GeneList, 1,7)) & length( unique(substr(GeneList, 1,7)) )>1 ) { 
        stop("Error in DoAGeneSummaryThing(): GeneList contains at least some 'ID' sequences that do not start with 'ENSOARG' whilst others do.\nCheck for gene names/human IDs and retry.")
        #TODO: Write conversion logic, keeping ENSOARG as ENSOARG and converting as many names to geneIDs as possible
    }
    
    if(!(Dataset %in% c("LKB.ERCC", "GWPM", "ERCC+GWPM"))) { stop("Argument 'Dataset' of DoAGeneSummaryThing() must be one of: LKB.ERCC, GWPM, ERCC+GWPM")  }
    switch(Dataset,
           "LKB.ERCC"={ if(!(Model %in% c("Spline", "1WAD", "Both"))) { stop("Argument 'Model' of DoAGeneSummaryThing() must be one of: Spline, 1WAD, Both")  } },
           "ERCC+GWPM"={ if(!(Model %in% c("Spline", "1WAD", "CvHF", "Spline+CvHF", "1WAD+CvHF"))) { 
               stop("Argument 'Model' of DoAGeneSummaryThing() must be one of: Spline, 1WAD, CvHF, Spline+CvHF, 1WAD+CvHF")  
               } }
           )
    #loading data----
    #Best to prepare and use separate "Both" and/or "GWPM" data objects; overwriting the "default" objects like ExpressPct might lead to confusion when it's not explicit
    #e.g. you'll use ExpressPct in another script and have no idea the 'default' expressPct is now actually ExpressPct.GWPM
    if(Dataset %in% c("ERCC+GWPM", "GWPM" )) { load("Data/RNAseq/Output/GWPM.GeneSummaryObjects.rda") }
    if(Dataset %in% c("ERCC+GWPM", "LKB.ERCC")) { load("Data/RNAseq/Output/LKB.ERCC.GeneSummaryObjects.rda") }
    
    GeneSummary = list()
    
    #ID data----
    if(isNameList) { FilterVar = "geneName" } else { FilterVar= "geneID" }
    GeneList.ID = TxCpy[get(FilterVar) %in% GeneList]
    GeneList.ID= unique(GeneList.ID)
    
    stopifnot(nrow(GeneList.ID)>0)
    
    if ( nrow(GeneList.ID) != length(GeneList) ) { 
        warning( sprintf("Only %s of %s items could be matched to an ID number/name. Could not match [ %s ]. Check for errant whitespace and retry?",
                    nrow(GeneList.ID), length(GeneList), paste0(setdiff(GeneList, GeneList.ID[,get(FilterVar)]), collapse = ","))) 
        }
    
    for (ID in GeneList.ID[,geneID]){ 
        if (GeneList.ID[geneID==ID, geneName]=="") { 
            #If there is no gene name found or the user names are forced, use the one supplied by the user. if there is no user-supplied name, use the geneID.
            if (HasNames) {
                if ( forceUserNames ) {
                    GeneList.ID[geneID==ID, geneName:=paste0("[", names(which(GeneList == ID)), "]")]
                } else {
                    GeneList.ID[geneID==ID, geneName:=ifelse(names(which(GeneList == ID))!="", paste0("[", names(which(GeneList == ID)), "]"), ID)]    
                } # if forceUserNames
            } else {
                GeneList.ID[geneID==ID, geneName:=ID]
            } #if HasNames
        } #if geneName=="" 
    } #for ID in GeneList.ID
    
    #n genes expressed----
    switch(Dataset,
           "LKB.ERCC"  = { GeneList.AreExpressed = ExpressPct[geneID %in% GeneList.ID[,geneID]] },
           "GWPM"      = { GeneList.AreExpressed = ExpressPct.GWPM[geneID %in% GeneList.ID[,geneID]] },
           "ERCC+GWPM" = { GeneList.AreExpressed = rbind(ExpressPct.GWPM[geneID %in% GeneList.ID[,geneID]], ExpressPct[geneID %in% GeneList.ID[,geneID]], fill=T)}
    )
    
    if (nrow(GeneList.AreExpressed) > 0){
        GeneList.AreExpressed = dcast(GeneList.AreExpressed, geneID~Timepoint, value.var = "IsExpressed")
        GeneList.AreExpressed = unique(merge(GeneList.AreExpressed, GeneList.ID, by="geneID"))
    }
    
    #Changing----
    switch(Dataset, 
           "LKB.ERCC"  = {  switch(Model, 
                                "Spline" = { GeneList.Change = DEGs.ERCC.SplineDF3[geneID %in% GeneList.ID[,geneID]] },
                                "1WAD" = { GeneList.Change = DEGs.ERCC.1WAD[geneID %in% GeneList.ID[,geneID]] },
                                "Both" = { GeneList.Change = rbind(DEGs.ERCC.SplineDF3[geneID %in% GeneList.ID[,geneID]], DEGs.ERCC.1WAD[geneID %in% GeneList.ID[,geneID]], fill=T) }
                                )
                         },
           "GWPM"      = { GeneList.Change = DEGs.GWPM[geneID %in% GeneList.ID[,geneID]] },
           "ERCC+GWPM" = { switch(Model, 
                                  "CvHF" = { GeneList.Change = DEGs.GWPM[geneID %in% GeneList.ID[,geneID]] },
                                  "1WAD" = { GeneList.Change = DEGs.ERCC.1WAD[geneID %in% GeneList.ID[,geneID]] },
                                  "Spline" = { GeneList.Change = DEGs.ERCC.SplineDF3[geneID %in% GeneList.ID[,geneID]] },
                                  "Spline+CvHF" = { GeneList.Change = rbind(DEGs.GWPM[geneID %in% GeneList.ID[,geneID]], 
                                                                            DEGs.ERCC.SplineDF3[geneID %in% GeneList.ID[,geneID]], fill=T) },
                                  "1WAD+CvHF" = { GeneList.Change = rbind(DEGs.GWPM[geneID %in% GeneList.ID[,geneID]], 
                                                                          GeneList.Change = DEGs.ERCC.1WAD[geneID %in% GeneList.ID[,geneID]], fill=T) }
                                 )
                        }
    )
    
    if ("Transcript stable ID" %in% colnames(GeneList.Change)) { GeneList.Change[,`Transcript stable ID`:=NULL] }
    GeneList.Change = unique(GeneList.Change)
    GeneList.Change = unique(merge(GeneList.Change, GeneList.ID, by="geneID"))
    
    #Graphing All Genes----
    GeneList.GraphChange = data.table()
    
    #counts don't change no matter the model, so we can use the SplineModel DESeqObject to obtain and plot counts for either SplineModel or GroupModel stats
    if (PlotTPMs==FALSE){
        if (!exists("SplineModel")) {
            SplineModel = readRDS("Data/RNAseq/Output/SplineModel.ERCC.ImputeB4Exclude.2excl.maxit5000.batchCutRINOXnsImputedAgeDaysdf3.LRT.RDS")
        }
        if (!exists("GWPMGroupModel")) { GWPMGroupModel = readRDS("Data/RNAseq/Output/GroupGWPM.2excl.maxit5000.BatchGroup.Wald.RDS")  }
        
        MakeCountDT <-function(TEMPCOUNTDATA, GraphFactor) {
            CountData.DT = as.data.table(DESeq2::counts(TEMPCOUNTDATA))
            CountData.DT[,geneID:=rownames(TEMPCOUNTDATA)]
            CountData.DT = melt(CountData.DT, id.vars="geneID")
            setnames(CountData.DT, c("variable", "value"), c("sample", "count"))
            CountData.DT[,group:=str_extract(sample, paste0(GraphFactor, collapse = "|"))]
            CountData.DT
        }
        
        switch(Dataset,
            "LKB.ERCC"  = { GraphFactor = c("1W","1M","3M","AD");                  GeneList.GraphChange = MakeCountDT(SplineModel, GraphFactor);},
            "GWPM"      = { GraphFactor = c("C", "HF", "R");                       GeneList.GraphChange = MakeCountDT(GWPMGroupModel, GraphFactor);},
            "ERCC+GWPM" = { GraphFactor = c("1W","1M","3M","AD", "C", "HF", "R");  GeneList.GraphChange = rbind(MakeCountDT(GWPMGroupModel, c("C", "HF", "R")), 
                                                                                                         MakeCountDT(SplineModel, c("1W","1M","3M","AD")))} 
        )
        rm(MakeCountDT)
    }
    
    if (PlotTPMs==TRUE){
        if (Dataset=="ERCC+GWPM"){
            TPMTable = readRDS("Data/RNAseq/Output/TPMs.GWPMLKBOverlap.RDS") #ERROR
            nMatchesInTPMTable = length(intersect(TPMTable[,geneID], GeneList.ID[,geneID]))
            if(nMatchesInTPMTable != GeneList.ID[,.N]){
                warning(paste0("TPMs must be normalised for LKB and GWPM data; thus only genes appearing in BOTH datasets can be plotted with TPMs.\n Here, only ", 
                               nMatchesInTPMTable, " of ", GeneList.ID[,.N], " requested genes are in the intersection.\n", 
                               "If all genes must be plotted, plot raw counts (not recommended), use PlotTPMs=FALSE")
                        )
            }
            GraphFactor = c("1W","1M","3M","AD", "C", "HF", "R")
            GeneList.GraphChange = TPMTable[geneID %in% GeneList.ID[,geneID]]
        } else {
            switch(Dataset,
                   "LKB.ERCC" = { GeneList.GraphChange = LKB.TPMCountTable[geneID %in% GeneList.ID[,geneID]]; GraphFactor = c("1W","1M","3M","AD"); },
                   "GWPM"     = { GeneList.GraphChange = BHB.TPMCountTable[geneID %in% GeneList.ID[,geneID]]; GraphFactor = c("C", "HF", "R"); }
            )
        }
        GeneList.GraphChange = melt(GeneList.GraphChange, id.vars="geneID")
        setnames(GeneList.GraphChange, c("variable", "value"), c("sample", "TPM"))
        GeneList.GraphChange[,group:=str_extract(sample, paste0(GraphFactor, collapse = "|"))]
        GeneList.GraphChange = GeneList.GraphChange[!is.na(group)]
    }
    
    GeneList.GraphChange[,group:=factor(group, levels = GraphFactor)]
    GeneList.GraphChange = merge(GeneList.GraphChange, GeneList.ID, by="geneID")
    #GeneList.GraphChange[,geneIndex:=NULL]
    
    Comparisons = unique(GeneList.Change[,Comparison])
    if (UniqueComparisonSymbols == TRUE ) { CompSymbols = c("*", "$", "#", "+", "&", "%") } else { CompSymbols = rep("*", length(Comparisons)) }
    SIGNIF_LEVELS = c(0.05, 0.01, 0.001)
    stopifnot(length(Comparisons) <= length(CompSymbols)) #only relevant if UniqueComparisonSymbols == TRUE
    PlotLegend = list()

    for (currentComparison in Comparisons) {
         currentCompSymbol = CompSymbols[which(currentComparison == Comparisons)]
         if (nrow(GeneList.Change[Comparison == currentComparison & pValAdjust < SIGNIF_LEVELS[1]]) > 0) { PlotLegend[[currentComparison]] = currentCompSymbol }
    }
    #GeneList.GraphChange
    
    if(length(PlotLegend)>0) { 
        PlotLegend=data.table(Comparison=factor(names(PlotLegend)), Symbol=unlist(PlotLegend)) 
        #PlotLegend[,SymbolInt := unlist(lapply(PlotLegend[,Symbol], utf8ToInt))]
        #PlotLegend[,SymbolInt := factor(SymbolInt)]
        PlotLegend[,Symbol:=factor(Symbol)]
    
        SignifList = GeneList.Change[,.(geneID,pValAdjust,Comparison,geneName)]
        #These are sort of specific but I don't feel like redoing the models just to make this work...
        SignifList[Comparison!="1WAD" & Comparison!="All", xmax:=str_match(Comparison, "^.*\\|(.*)\\|.*")[,2]]
        SignifList[Comparison!="1WAD" & Comparison!="All", xmin:=str_match(Comparison, "^.*\\|.*\\|(.*)")[,2]]
        #Kludge:
        SignifList[Comparison=="1WAD" | Comparison=="All", xmin:="1W"]
        SignifList[Comparison=="1WAD" | Comparison=="All", xmax:="AD"]
        SignifList = merge(SignifList, PlotLegend, by="Comparison")
        SignifList[pValAdjust < SIGNIF_LEVELS[3], Symbol:= strrep(Symbol, 3)]
        SignifList[pValAdjust >= SIGNIF_LEVELS[3] & pValAdjust < SIGNIF_LEVELS[2], Symbol:= strrep(Symbol, 2)]
        
        if (PlotTPMs==TRUE) {
            SignifList = merge(SignifList, GeneList.GraphChange[,round(max(TPM)*0.9, 0), geneID], by="geneID", all.x=T)
            if (all(is.na(SignifList[,V1]))) { SignifList[, V1 := 0] } #If TPMs cannot be calculated, there will be an empty entry due to no data
            
            #SignifList[is.na(V1), V1:=GeneList.GraphChange[,round(max(TPM)*0.9, 0)]] #If there are no matches, we take a safe value - the max 
            #Commented out because the above two adjustments should avoid any full NA scenarios?
        } else {
            SignifList = merge(SignifList, GeneList.GraphChange[,round(max(count)*0.9, 0), geneID], by="geneID")    
        }
        
        setnames(SignifList, c("V1"), c("y_position"))
        
        #Let's make sure they also don't overlap -
        SignifList[,xmax_d:=as.numeric(factor(xmax, levels=GraphFactor))]
        SignifList[,xmin_d:=as.numeric(factor(xmin, levels=GraphFactor))]
        SignifList[xmax_d<xmin_d, temp := xmax_d]
        SignifList[!is.na(temp), xmax_d:= xmin_d]
        SignifList[!is.na(temp), xmin_d:= temp]
        SignifList[, temp := NULL]
        if (PlotTPMs==TRUE) {
            SignifList = merge(SignifList, SignifList[,GroupNonOverlappingIntervals(.SD, "xmin_d", "xmax_d"), geneName], by=c("geneName", "xmin_d", "xmax_d"), all.x=T)
        } else {
            SignifList = merge(SignifList, SignifList[,GroupNonOverlappingIntervals(.SD, "xmin_d", "xmax_d"), geneName], by=c("geneName", "xmin_d", "xmax_d"))
        }
        
        SignifList[,y_position_fin := y_position + (y_position * (NoOverlapGroup * 0.2))] #Offset 7.5% per non-overlapping comparison
    }
    
    #TODO - add multinamed code to SignifList
    MultiNamed = GeneList.GraphChange[,.N,.(geneID, geneName)]
    MultiNamed = MultiNamed[geneName %in% MultiNamed[,.N, geneName][N>1], geneID]
    GeneList.GraphChange[geneID %in% MultiNamed, geneName:=paste0(geneName, .GRP), .(geneName, geneID)]
    
    if(PlotTPMs==T) {
        BasePlot = ggplot(GeneList.GraphChange, mapping = aes(x=group, y=TPM, color=group)) 
    } else {
        BasePlot = ggplot(GeneList.GraphChange, mapping = aes(x=group, y=count, color=group)) 
    }
    
    GeneList.Plot = BasePlot+geom_boxplot(color="black", outlier.shape = NA)+geom_jitter(width = 0.05)+xlab("")+ylab(ifelse(PlotTPMs, "TPM", "Normalized Counts"))+
        theme(legend.position = "none")+facet_wrap(~geneName, scales = ifelse(UseCommonScale, "fixed", "free_y"), ncol = 3)
    
    if(exists("SignifList")) { 
        GeneList.Plot = GeneList.Plot+geom_signif(data=SignifList, mapping = aes(annotations = Symbol, xmin = xmin, xmax = xmax, y_position = y_position_fin), color="#000000", manual = TRUE)+
            scale_y_continuous(expand = expand_scale(mult = c(0, .075))) 
    }
    
    if(Dataset =="ERCC+GWPM") { GeneList.Plot = GeneList.Plot+geom_vline(xintercept = 4.5) }
    
    #Graphing significant genes only----
    if (nrow(GeneList.GraphChange[geneID %in% GeneList.Change[,geneID]])>0) {
        GeneList.PlotSign = ggplot(GeneList.GraphChange[geneID %in% GeneList.Change[,geneID]], 
                          mapping = aes_string(x="group", y=ifelse(PlotTPMs, "TPM", "count"), color="group"))+
            geom_boxplot(color="black", outlier.shape = NA, show.legend = FALSE)+geom_jitter(width = 0.05, show.legend = FALSE)+
            xlab("")+ylab(ifelse(PlotTPMs, "TPM", "Counts"))+facet_wrap(~geneName, scales = ifelse(UseCommonScale, "fixed", "free_y"))
        
        if(exists("SignifList")) { 
            GeneList.PlotSign = GeneList.PlotSign+ggsignif::geom_signif(data=SignifList, aes(annotations = Symbol, xmin = xmin, xmax = xmax, y_position = y_position_fin), manual = TRUE, color="black")+
                                scale_y_continuous(expand = expand_scale(mult = c(0, .2))) 
        }
        if(Dataset =="ERCC+GWPM") { GeneList.PlotSign = GeneList.PlotSign+geom_vline(xintercept = 4.5) }
    } else {
        GeneList.PlotSign = ggplot() + cowplot::draw_label(label = "No significant results")
    }
    
    #Heatmap----
        if (MakeHeatmap==TRUE){
        switch(Dataset,
               "LKB.ERCC"  = { ScaledForHeatmap = ScaledCounts },
               "GWPM"      = { ScaledForHeatmap = ScaledCounts.GWPM},
               "ERCC+GWPM" = { ScaledForHeatmap = rbind(ScaledCounts, ScaledCounts.GWPM) }
        )
        HeatmapData = ScaledForHeatmap[geneID %in% GeneList.ID[,geneID]]
        if( HeatmapData[substr(geneName, 1, 7) == "ENSOARG", .N] > 0 & !is.null(names(GeneList))){
            HeatmapData[substr(geneName, 1, 7) == "ENSOARG", geneName:=paste0("[", names(GeneList)[match(geneName,GeneList)], "]")]
        }
        for (SignifLevel in c(0.05, 0.01, 0.001)) { HeatmapData[geneID %in% GeneList.Change[pValAdjust<SignifLevel,geneID], geneName:=paste0("*", geneName)] }
        
        HeatmapData[,group:=factor(group, levels = GraphFactor)]
        HeatmapData[,sample:=factor(Sample, levels = HeatmapData[order(group), unique(Sample)])]
        
        Heatmap.Min = HeatmapData[,floor(min(value)*2)/2] #rounded to nearest 0.5
        Heatmap.Max = HeatmapData[,ceiling(max(value)*2)/2] #rounded to nearest 0.5
        
        Heatmap = ggplot(mapping = aes_string(x="sample", y="geneName", fill="value"))+theme(axis.text.x = element_text(angle=90, vjust = 0.3))+
            scale_fill_gradientn(name="Counts scaled\nrelative to\ngene mean", colors=c("#0000FF", "#FFFFFF", "#FF0000"),
                                 # Mon Jul 09 14:29:03 2018 implemented automatic gradient thresholding with squishing
                                 limits=c(max(Heatmap.Min, -HeatmapThreshold), min(Heatmap.Max, HeatmapThreshold)), oob=scales::squish)+
            scale_x_discrete()+ylab("")+theme(axis.text.y = element_text(hjust = 0))
        
        GeneHeatmap = Heatmap + geom_raster(data = HeatmapData)+
            geom_vline(xintercept = switch(Dataset, "LKB.Old"=c(7.5, 14.5, 20.5), "LKB.ERCC"=c(7.5, 14.5, 20.5), "GWPM"=c(7.5,15.5), 
                                           "ERCC+GWPM"=c(7.5, 14.5, 20.5, 28.5, 35.5, 43.5), "Old+GWPM" = c(7.5, 14.5, 20.5, 28.5, 35.5, 43.5)))
        GeneHeatmapSig = Heatmap + geom_raster(data = HeatmapData[geneID %in% GeneList.Change[,geneID]])+
            geom_vline(xintercept = switch(Dataset, "LKB.Old"=c(7.5, 14.5, 20.5), "LKB.ERCC"=c(7.5, 14.5, 20.5), "GWPM"=c(7.5,15.5), 
                                           "ERCC+GWPM"=c(7.5, 14.5, 20.5, 28.5, 35.5, 43.5), "Old+GWPM" = c(7.5, 14.5, 20.5, 28.5, 35.5, 43.5)))
        
        #Heatmap output
        GeneSummary$Heatmap         = GeneHeatmap
        GeneSummary$HeatmapSignif   = GeneHeatmapSig
        GeneSummary$HeatmapData     = HeatmapData
        GeneSummary$HeatmapMeta     = data.table(Min=Heatmap.Min, Max=Heatmap.Max)
        #GeneSummary$HeatmapPatterns = GMHeatMap
    }
    
    #Output----
    GeneSummary$Changing        = GeneList.Change
    GeneSummary$Plot            = GeneList.Plot
    GeneSummary$PlotSignif      = GeneList.PlotSign
    GeneSummary$PlotData        = GeneList.GraphChange
    GeneSummary$Expressed       = GeneList.AreExpressed
    GeneSummary$Genes           = GeneList.ID
    return(GeneSummary)
}
