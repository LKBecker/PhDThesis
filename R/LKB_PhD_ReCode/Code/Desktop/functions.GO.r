library(data.table); library(stringr)
#functions====
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

#TODO - remove SheepID references - replace with GeneIDTerm, supplied by user, 
#which determines which column to collapse multi-mapping genes into/through
RankPANTHERGOTermTable <- function(DT, Term, NLeafNodes=NA, MinLeafRank=-1, FilterObsolete=TRUE, FoldChangeTerm=NULL, geneTerm="SheepID"){
    #DT=copy(a);Term=colIDs[which(names(colIDs)=="BP")];NLeafNodes=NA; MinLeafRank=-1;FilterObsolete=TRUE;FoldChangeTerm="log2FC"
    
    if( !(names(Term) %in% c("BP", "MF", "CC")) ) {
        stop("Term must be one of: 'BP', 'MF', 'CC'")
    }
    switch(names(Term),
           "BP"={ RankTable=GORank$Biological },
           "MF"={ RankTable=GORank$Molecular },
           "CC"={ RankTable=GORank$Cellular },
           { stop(sprintf("Cannot load GO Rank table for term \"%s\". Check DoTheWholePanther() input and retry.", names(Term))) }
    )
    
    if (geneTerm=="") { stop(paste0("Error in RankPANTHERGOTermTable (", deparse(substitute(DT))), "): geneTerm is not defined") }
    if (!(geneTerm %in% colnames(DT))) { stop(paste0("Error in RankPANTHERGOTermTable: ", dQuote(geneTerm), " is not a valid column in ", deparse(substitute(DT)))) }
    if (geneTerm!="") { setnames(DT, geneTerm, "GENETERM") }
    
    tryCatch(expr = {
        DT2   = SplitDataTableWithMultiRows(DT, Term, ";")
        set(DT2, j = DT2[,which(!(colnames(DT2) %in% c(Term, "MappedIDs", "GeneNameSymbolOrtholog", "GENETERM", "Name", FoldChangeTerm)))], value = NULL)
        # new 15/3 -- erase superfluous columns
        if (FilterObsolete) { DT2 = DT2[!(get(Term) %like% "obsolete" )] }
        DT2.E = DT2[get(Term)==""]
        DT2   = DT2[get(Term)!=""]
        
        DT2.E[,Rank:=-1] #So that rbind will have a Rank column
        DT2.E[,GO_ID:=""]
        
        DT2[, GO_ID:=str_match(get(Term), ".*\\((GO:\\d+)\\)$")[,2]]
        DT2[,(Term):=str_replace(get(Term), "\\(GO:\\d+\\)", "")] #removes GO ID from name
        DT2  = merge(DT2,  RankTable,  by.x="GO_ID", by.y="Node", all.x=T)
        
        if(any(is.na(DT2[!is.na(GO_ID),Rank]))) {
            stop(paste0("One or more ", Term, " GO Terms have no associated rank(s). Please manually curate using e.g. https://www.ebi.ac.uk/QuickGO/term/\n",
                        paste(unique(DT2[is.na(Rank), GO_ID]), collapse = ",")))
        }
        
        if( any(DT2[get(Term)=="", Rank] != -1) ) {
            EmptyIDS = DT2[get(Term)=="", GENETERM]
            stop("At least one GO entry of type '", Term, "' has no entry (description is \"\") but a rank that is not -1.\n",
                 "Possible cause: obsolete GO entries from PANTHER, which are given an ID but no name ",
                 "(search for ';(GO:' in your input file to quickly find them).\n",
                 "Please rectify and re-run the function. Affected genes: ", paste(EmptyIDS, collapse = ","))
            #how defensive can my coding get
        }
        
        DT2=DT2[Rank >= MinLeafRank]
        
        if (!is.na(NLeafNodes)){ DT2=DT2[,.SD[Rank %in% .SD[,unique(Rank)][order(Rank, decreasing = T)][1:NLeafNodes]],GENETERM] }
        
        DT2 = rbind(DT2, DT2.E)
        if (geneTerm!="") { setnames(DT2, "GENETERM", geneTerm) }
    }, finally =  { if (geneTerm!="") { setnames(DT, "GENETERM", geneTerm) } })
    return(DT2)
}
#TODO: WrapStr works in QuantifyPANTHERTerms but not in MakePANTHERPlot...
QuantifyPANTHERTerms <- function(DTx, Term, FilterWords=NA, RemoveWords=NA, WrapStrings=0, RemoveMultiMapping=TRUE,
                                 FoldChangeTerm=NULL, FillEmpty=NULL, MultiMapTerm=NA, geneTerm="SheepID"){
    #TODO Write TEST dataset, run TEST dataset
    #DTx=PANTHERObj$Biological.DT; Term="BiologicalProcess";FilterWords=NA; RemoveWords=NA; WrapStrings=20;RemoveMultiMapping=TRUE; FoldChangeTerm="logFC"
    #DTx=TESTDATA; Term="MolecularFunction";FilterWords=NA;RemoveWords=NA;WrapStrings=20;RemoveMultiMapping=TRUE;FoldChangeTerm=NULL
    DT = copy(DTx)
    setnames(DT, Term, "TEMP")
    FilterWords  = ifelse(is.na(FilterWords), NA, paste(FilterWords, collapse = "|")[1])
    if (is.character(FilterWords)) { suppressWarnings( { DT = DT[grepl(FilterWords, TEMP)] } ) }
    if (is.character(RemoveWords)) { DT = DT[!(Term %in% RemoveWords)] }
    
    if (geneTerm=="") { stop(paste0("Error in QuantifyPantherTerms (", deparse(substitute(DT))), "): geneTerm is not defined") }
    if (geneTerm!="") { setnames(DT, geneTerm, "GENETERM") }
    
    if(!is.null(FillEmpty)){ set(x = DT, i = which(DT[,TEMP]==""), j = "TEMP", value = FillEmpty) }
    
    if (!is.null(FoldChangeTerm)) { setnames(DT, FoldChangeTerm, "FOLDCHANGE") }
    
    if(RemoveMultiMapping==TRUE) {
        #if one sheep gene maps to multiple human genes, you'll see (over) enrichment of some Terms (e.g. protein_coding for gene1 and gene2 -> 2x)
        #due to tags of ortholog1 and ortholog2 mixing. Whether this is "valid" is left to the researcher, and can be toggled
        #TODO: Insert simple heuristic for picking among multimappers - e.g. the one with a name (that is not 'unassigned')?
        if (!is.null(FoldChangeTerm)) {
            DT  = unique(DT[,.( TEMP, Rank, FOLDCHANGE), GENETERM])
            DT  = DT[,.N,.(TEMP, Rank, FOLDCHANGE)]
        } else {
            DT  = unique(DT[,.( TEMP, Rank), GENETERM])
            DT  = DT[,.N,.(TEMP, Rank)]
        } #if !is.null(FoldChangeTerm)
    }
    
    if (!is.null(FoldChangeTerm)) {
        tempDT = DT[,.N, TEMP][order(N, decreasing = TRUE)][,TEMP] #temp DT, sums occurences regardless of direction for later ranking
        #180419 changed final [,get(Term)] to [,get]
        DT[,ChangeDir:=FOLDCHANGE>0] #proxy variable for coloring and sorting
        DT  = DT[,.( .N, xFoldChangeSum=sum(FOLDCHANGE) ),.(TEMP, ChangeDir, Rank)]
        if(is.numeric(WrapStrings) & WrapStrings > 0) {
            tempDT = str_wrap(tempDT, WrapStrings)
            DT[,TEMP:=str_wrap(TEMP, WrapStrings)]
        }
        DT[,TEMP:=factor(TEMP, levels=unique(tempDT), ordered = T)]
    } else { #FoldChangeTerm is null
        if(is.numeric(WrapStrings) & WrapStrings > 0) {
            DT[,TEMP:=str_wrap(TEMP, WrapStrings)]
        }
        DT[,TEMP:=factor(TEMP, levels=unique(DT[order(N, decreasing = T)][, TEMP]), ordered = T)]
    }
    
    setnames(DT, c("TEMP"), c(Term))
    return(DT) #View(DT)
}

MakePANTHERGridPlot<-function(PANTHERObj, NTop=100, Terms = c("MF"="MolecularFunction", "BP"="BiologicalProcess", "CC"="CellularComponent"), 
                              SplitByNodeLevel=FALSE, ScaleToFoldChanges=FALSE, WrapStrings=20, TruncateTo=NULL, TruncateStr="[...]", 
                              UnassignedStr="", StackBars=TRUE){
    #PANTHERObj=PANTHERResults; NTop=20; SplitByNodeLevel=FALSE;ScaleToFoldChanges=FALSE;WrapStrings=20;TruncateTo=NULL
    #TruncateStr="[...]";UnassignedStr="";StackBars=TRUE
    stopifnot(!is.null(names(Terms)))
    stopifnot(length(unique(Terms))==3)
    stopifnot(length(unique(names(Terms)))==3)
    stopifnot(sort(names(Terms)) == c("BP", "CC", "MF"))
    
    grid = list(MolePlot=list(), BioPlot=list(), CellPlot=list())
    if(SplitByNodeLevel==TRUE) {
        plotObj = list(list(PANTHERObj$Molecular.Q,  Terms[["MF"]], "Molecular Function", 1),
                       list(PANTHERObj$Biological.Q, Terms[["BP"]], "Biological Process", 2), 
                       list(PANTHERObj$Cellular.Q,   Terms[["CC"]], "Cellular Component", 3) 
        ) #plotObj
        maxRanks = c(PANTHERObj$Molecular.Q[,max(Rank)], PANTHERObj$Biological.Q[,max(Rank)], PANTHERObj$Cellular.Q[,max(Rank)]) 
        #In order of the integers in plotObj[[n]][[4]]
        for (aRank in seq(-1, max(maxRanks))) {
            for (plotObjItem in plotObj) {
                if (aRank <= maxRanks[plotObjItem[[4]]]) {
                    DT = plotObjItem[[1]][Rank==aRank]
                    xPlot = MakePANTHERPlot(DT, plotObjItem[[2]], paste0(plotObjItem[[3]], " Rank ", aRank), NTop, ScaleToFoldChanges,
                                            WrapStrings, TruncateTo, TruncateStr, UnassignedStr, StackBars)
                    grid[[ plotObjItem[[4]] ]][[paste0("Rank", aRank)]] = xPlot
                }
            } #for plotObjItem in plotObj
        } #for Rank in seq(-1, max Rank)
    }#if SplitByNodeLevel==TRUE
    else{
        grid[[1]] <- MakePANTHERPlot(PANTHERObj$Molecular.Q,  Terms[["MF"]], "Molecular Function", NTop, ScaleToFoldChanges,
                                     WrapStrings, TruncateTo, TruncateStr, UnassignedStr, StackBars)
        grid[[2]] <- MakePANTHERPlot(PANTHERObj$Biological.Q, Terms[["BP"]], "Biological Process", NTop, ScaleToFoldChanges,
                                     WrapStrings, TruncateTo, TruncateStr, UnassignedStr, StackBars)
        grid[[3]] <- MakePANTHERPlot(PANTHERObj$Cellular.Q,   Terms[["CC"]], "Cellular Component", NTop, ScaleToFoldChanges,
                                     WrapStrings, TruncateTo, TruncateStr, UnassignedStr, StackBars)
        grid[["ThreeInOne"]] <- plot_grid(grid[[1]], grid[[2]], grid[[3]], ncol = 1, nrow = 3)
    }#SplitByNodeLevel==FALSE
    return(grid)
}

NthMostDTTerm <- function(DT, N, RankColName, IDColName, Descending=TRUE) {
    #DT=DoTheWholePANTHER(Spline.Full[BiologicalProcess %like% "calcium"], foldChangeTerm = "log2FC")$Biological.DT; N=1; RankColName="Rank"
    #Descending=TRUE; IDColName="SheepID";
    setnames(DT, c(RankColName, IDColName), c("RANKCOLUMN.FUNC", "IDCOLUMN.FUNC"))
    setorderv(DT, "RANKCOLUMN.FUNC", ifelse(Descending, -1, 1))
    DT2 <- DT[, .SD[RANKCOLUMN.FUNC==unique(RANKCOLUMN.FUNC)[N]], IDCOLUMN.FUNC]
    setnames(DT, c("RANKCOLUMN.FUNC", "IDCOLUMN.FUNC"), c(RankColName, IDColName))
    
}

NFromBottomPANTHERPlot <- function(PANTHERObj, N=1, RankColName, IDColName, Descending=TRUE) {
    grid = list(MolePlot=list(), BioPlot=list(), CellPlot=list())
    PANTHERObj$Biological.DT
    PANTHERObj$Molecular.DT
    PANTHERObj$Cellular.DT
    
}

#TODO = try 'custom' factor, rip the numerical values (not unique, the whole ones), 
#the 'paste' them onto the truncated, wrapped strings. / order by numeric, get labels.
#What if there's non-unique labels? Compensator / Append a counter?
MakePANTHERPlot<-function(DT, Term, Title, NTop, ScaleToFoldChanges=FALSE, WrapStrings=9999, TruncateTo=9999, TruncateStr="[...]", UnassignedStr = "", 
                          StackBars = TRUE) {
    #DT=copy(PANTHERObj$Biological.Q); Term="BiologicalProcess"; TruncateTo=40; TruncateStr="[...]"; WrapStrings=20; 
    #HasFoldChangeTerm=TRUE; NTop=25;ScaleToFoldChanges=F;Title="DEBUGTITLE"
    setnames(DT, Term, "TEMP")

    #if (HasFoldChangeTerm==TRUE) {
    if ("xFoldChangeSum" %in% colnames(DT)) {
        topTerms = DT[, sum(N), TEMP][order(V1, decreasing = T), TEMP]
        topTerms = topTerms[1:min(nrow(topTerms), NTop)]
        plotData =  DT[TEMP %in% topTerms]
        #plotLabels = plotData[order(TEMP, decreasing = T), unique(TruncLabel)]
        if (ScaleToFoldChanges==TRUE) {
            if ( StackBars == TRUE ) { plotData[,  xFoldChangeSum := abs(xFoldChangeSum)]  }
            # Wed Jun 20 19:48:20 2018 oops if you don't tell them FALSE means Down and only submit TRUE, it'll still count the first value as Down and
            #colors it red. Fixed through using a named vector for explicit association.
            xPlot <- ggplot() + geom_col(data = plotData, mapping = aes_string(x="TEMP", y="xFoldChangeSum", fill="ChangeDir")) +
                ggtitle(Title)+scale_fill_manual(name="Regulated", labels=c(`FALSE`="Down", `TRUE`="Up"), values=c(`FALSE`="#e41a1c", `TRUE`="#4daf4a"))+
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none")+ylab("Sum of log2 fold changes")
        } else {
            if ( StackBars == FALSE ) { plotData[ xFoldChangeSum < 0,  N := -N]  }
            xPlot <- ggplot() + geom_col(data = plotData, mapping = aes_string(x="TEMP", y="N", fill="ChangeDir"))+
                ggtitle(Title)+scale_fill_manual(name="Regulated", labels=c(`FALSE`="Down", `TRUE`="Up"), values=c(`FALSE`="#e41a1c", `TRUE`="#4daf4a"))+
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none")+ylab("Number of genes")
        }# if ScaleToFoldChanges
    } else { 
        topTerms = DT[, sum(N), TEMP][order(V1, decreasing = T), TEMP]
        topTerms = topTerms[1:min(nrow(topTerms), NTop)]
        plotData =  DT[TEMP %in% topTerms] 
        xPlot <- ggplot() + geom_col(data = plotData, mapping = aes_string(x="TEMP", y="N"))+
            ggtitle(Title)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+ylab("Number of genes")
        #TODO fix label
    }
    #TODO: replace with space-seeking function to get nearest whole-word abbreviation?
    #e.g. scale_x_discrete(labels = function(x) ifelse(str_length(x)>40, str_wrap(str_trunc(x, 40, "center"), 30), str_wrap(x, 30)))+coord_flip()
    LabelFunction = function(x) { x = str_trunc(as.character(x), 5, "center", TruncateStr); return(x) }
    xPlot = xPlot+xlab("")
    setnames(DT, "TEMP", Term)
    return(xPlot)
}

DoTheWholePANTHER <- function(DT, NLeafNodes=NA, MinLeafRank=-1, FilterWords=NA, RemoveWords=NA, wrapStrings=20,
                              foldChangeTerm=NULL, FillEmpty=NULL, RemoveMultiMapping=TRUE, geneTerm="SheepID", 
                              colIDs = c("BP"= "BiologicalProcess", "CC"="CellularComponent", "MF"="MolecularFunction")){
    #a = fread("Data/PANTHER/181218_ERCC.ImputeB4Exclude.batchCutRINOX3nsImputedAgeDaysdf3.maxit5000.LRT.GOFullAndSlim.One2OneOrthosOnly.txt")
    #DT=copy(a); geneTerm="MappedIDs";colIDs = c("BP"= "GO.BP.Complete", "CC"="GO.CC.Complete", "MF"="GO.MF.Complete")
    #FilterWords=NA; RemoveWords=NA; wrapStrings=-1; foldChangeTerm=NULL; MinLeafRank=-1; NLeafNodes=NA; RemoveMultiMapping = TRUE; 
    #FillEmpty="[Unannotated Genes]";geneTerm="MappedIDs"
    stopifnot(!is.null(names(colIDs)))
    
    stopifnot(length(unique(colIDs))==3)
    stopifnot(all(colIDs %in% colnames(DT)))
    
    stopifnot(length(unique(names(colIDs)))==3)
    stopifnot(sort(names(colIDs)) == c("BP", "CC", "MF"))
    #why am i trying to make this code idiot-proof when i'm probably the only idiot who'll ever run it?
    
    PANTHERObj = list()
    
    PANTHERObj$Biological.DT = RankPANTHERGOTermTable(DT, colIDs[which(names(colIDs)=="BP")], NLeafNodes, MinLeafRank, FoldChangeTerm = foldChangeTerm, 
                                                      geneTerm = geneTerm)
    PANTHERObj$Biological.Q  = QuantifyPANTHERTerms(PANTHERObj$Biological.DT, colIDs[which(names(colIDs)=="BP")], FilterWords, RemoveWords,
                                                    FoldChangeTerm=foldChangeTerm, WrapStrings = wrapStrings, FillEmpty = FillEmpty,
                                                    RemoveMultiMapping = RemoveMultiMapping, geneTerm = geneTerm)
    
    PANTHERObj$Cellular.DT   = RankPANTHERGOTermTable(DT, colIDs[which(names(colIDs)=="CC")], NLeafNodes, MinLeafRank, FoldChangeTerm = foldChangeTerm, 
                                                      geneTerm = geneTerm)
    PANTHERObj$Cellular.Q    = QuantifyPANTHERTerms(PANTHERObj$Cellular.DT, colIDs[which(names(colIDs)=="CC")], FilterWords, RemoveWords,
                                                    FoldChangeTerm=foldChangeTerm, WrapStrings = wrapStrings, FillEmpty = FillEmpty,
                                                    RemoveMultiMapping = RemoveMultiMapping, geneTerm = geneTerm)
    
    PANTHERObj$Molecular.DT  = RankPANTHERGOTermTable(DT, colIDs[which(names(colIDs)=="MF")], NLeafNodes, MinLeafRank, FoldChangeTerm = foldChangeTerm, 
                                                      geneTerm = geneTerm)
    PANTHERObj$Molecular.Q   = QuantifyPANTHERTerms(PANTHERObj$Molecular.DT, colIDs[which(names(colIDs)=="MF")], FilterWords, RemoveWords,
                                                    FoldChangeTerm=foldChangeTerm, WrapStrings = wrapStrings, FillEmpty = FillEmpty,
                                                    RemoveMultiMapping = RemoveMultiMapping, geneTerm = geneTerm)
    PANTHERObj
}

#Setup: GO rank lists====
GORank = list()
GORank$Biological = readRDS("Data/GORank/180226-1134_GOTravellingSalesmanRank.BP.RDS")
GORank$Biological = rbind(GORank$Biological, data.table(Node=c("GO:0016337", "GO:0044707", "GO:0007067", "GO:0006917", "GO:0106064", 
                                                               "GO:0062009", "GO:0120117", "GO:0150012", "GO:0150020", "GO:0099178",
                                                               "GO:0120162", "GO:0120163", "GO:0062043", "GO:0106074", "GO:0150018",
                                                               "GO:0099173", "GO:0099175", "GO:0140199", "GO:0062026", "GO:0099179",
                                                               "GO:0110091", "GO:0099186", "GO:0106072", "GO:0110061", "GO:0110104",
                                                               "GO:0150024", "GO:0062028", "GO:0099183", "GO:0106071", "GO:0106077",
                                                               "GO:0106101", "GO:0150019", "GO:0150033", "GO:0099180", "GO:0099185",
                                                               "GO:0106089", "GO:0120158", "GO:0061951", "GO:0061952", "GO:0140253",
                                                               "GO:0007126", "GO:0008105", "GO:0061982", "GO:0110095", "GO:0140289",
                                                               "GO:0140290", "GO:0150078", "GO:0150093", "GO:0150094"
), 
Rank=c(3, 1, 3, 4, 5, 
       4, 5, 7, 7, 5,
       4, 3, 5, 4, 4,
       3, 4, 4, 9, 5,
       6, 4, 5, 5, 9,
       3, 6, 7, 6, 9,
       7, 6, 4, 10,5,
       5, 5, 5, 4, 2,
       2, 3, 3, 3, 7,
       7, 5, 3, 3))) 
#manual completion of missed nodes... why missed? Missing from package according to e.g. GOID("GO:0007067")
GORank$Molecular = readRDS("Data/GORank/180226-1119_GOTravellingSalesmanRank.MF.RDS")
GORank$Molecular = rbind(GORank$Molecular, data.table(Node=c("GO:0030695", "GO:0005083", "GO:0061953", "GO:0102092", "GO:0102758",
                                                             "GO:0103068", "GO:0120147", "GO:0120170", "GO:0099186", "GO:0106105",
                                                             "GO:0120153", "GO:0150025", "GO:0106073", "GO:0106078", "GO:0120160",
                                                             "GO:0099184", "GO:0062069", "GO:0062072", "GO:0062076", "GO:0102121",
                                                             "GO:0019204", "GO:0016251", "GO:0106140", "GO:0120146", "GO:0140292"), 
                                                      Rank=c(4, 4, 5, 5, 5,
                                                             4, 5, 3, 4, 5,
                                                             3, 3, 6, 9, 3,
                                                             4, 3, 5, 6, 5,
                                                             -1,3, 3, 4, 5))) 

GORank$Cellular  = readRDS("Data/GORank/180226-1116_GOTravellingSalesmanRank.CC.RDS")
GORank$Cellular = rbind(GORank$Cellular, data.table(Node=c("GO:0016023", "GO:0062023", "GO:0120134", "GO:0120135", "GO:0150034",
                                                           "GO:0150056", "GO:0110070", "GO:0120115", "GO:0140261", "GO:0140268"), 
                                                    Rank=c(4, 3, 6, 6, 4,
                                                           8, 4, 3, 5, 4))) 