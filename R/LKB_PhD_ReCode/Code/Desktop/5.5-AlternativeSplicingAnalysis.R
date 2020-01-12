#TODO: multi-transcript mode. You get erroneous "skipped exons" detected when an alternate transcript has a different first exon...
#TODO: Actually protein coordinates are relative to the splice variant you're looking at.........
#Sequence alignment may be easier. You'll have to do that anyway.
#libraries====
library(data.table); library(stringr); library(ggplot2); library(cowplot); library(dplyr); library(reshape2)
source("Code/Desktop/functions.r")
#graphs theme----
TextSize = 10 #14 for thesis, 20+ for presentation
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
#functions====
#Used to distance exon number labels, should alternate transcripts have different exon numbers of the same region. 
ComputeYOffset = function(x, offset=0.025){
    nItems=length(x)
    if(nItems==1){ return(0.5) }
    return(seq(0.5-(offset * (nItems/2)), 0.5+(offset * (nItems/2)), length.out = nItems))
}

MakeASPlots <- function(GeneMap, GeneGTFData, GENE_TO_PROCESS, transcript, DRAW_MODE, LABEL_INTRONS, LABEL_OVERLAPS){
    stopifnot(all(c("first.exon_number", "last.exon_number", "first.start", "i.start", "last.end", "i.end", "exonsSkippedN") %in% colnames(GeneMap)))
    exon_coords = unique(GeneGTFData[transcript_id == transcript & feature=="exon",.(exon_number, start, end)]) #this includes UTR
    #Toggle full draw / only area of interest draw
    if (is.numeric(DRAW_MODE)) {
        DEArea = exon_coords[exon_number >= (GeneMap[,min(first.exon_number)] - DRAW_MODE) & 
                                 exon_number <= (GeneMap[,min( last.exon_number)] + DRAW_MODE),
                             .(start=min(start), end=max(end))]
    } else {
        switch(DRAW_MODE, 
               "FULL_GENE"   = { DEArea = exon_coords[,.(start=min(start), end=max(end))] },
               "REGION_ONLY" = { DEArea = GeneMap[,.(start=min(first.start, i.start), end=max(last.end, i.end))] }
        )
    }
    stopifnot(exists("DEArea"))
    
    #intron plot====
    tmsg("Preparing to draw AS plot...")
    
    exonsnearDEIntrons = exon_coords[start <= DEArea[,end+500] & end >= DEArea[,start-500] ]
    exonsnearDEIntrons[,mid:=start+((end-start)/2)]
    exonsnearDEIntrons[,exon_instance:=.GRP,mid] #TODO: Insert test for proximity between this and previous group (sort by mid?)
    exonsnearDEIntrons[,y_offset:=ComputeYOffset(seq_len(.N)), exon_instance]
    setorder(exonsnearDEIntrons, exon_instance)
    
    exon_plot = ggplot(exonsnearDEIntrons)+geom_rect(aes(xmin=start, xmax=end, ymin=0.375, ymax=0.625), alpha=1)+ #metafiles can't handle alpha
        geom_text(aes(label=exon_number, x=mid, y=y_offset))+
        geom_hline(yintercept = 0.5, alpha=1)+theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank())
    
    # grab intron coordinates
    Coord_Groups = GroupNonOverlappingIntervals(GeneMap, "i.start", "i.end")
    Coord_Groups = merge(Coord_Groups, GeneMap[,.(i.start, i.end, ID)], by=c("i.start", "i.end"))
    Coord_Groups[, ystart:=0.6+(0.07*NoOverlapGroup)]
    Coord_Groups[, yend  :=0.6+(0.07*NoOverlapGroup)+0.05]
    Coord_Groups[, meanX := i.start+(i.end-i.start)/2]
    Coord_Groups[, meanY := ystart+(yend-ystart)/2]
    
    exon_skipping_introns = GroupNonOverlappingIntervals(GeneMap[exonsSkippedN != 0], "i.start", "i.end")
    if (nrow(exon_skipping_introns)>0) {
        exon_skipping_introns[, ystart:=0.3-(0.07*NoOverlapGroup)]
        exon_skipping_introns[, yend  :=0.3-(0.07*NoOverlapGroup)-0.05]
        exon_skipping_introns[, meanX := i.start+(i.end-i.start)/2]
        exon_skipping_introns[, meanY := ystart+(yend-ystart)/2]
        exon_skipping_introns = merge(exon_skipping_introns, GeneMap[,.(i.start, i.end, ID)], by=c("i.start", "i.end"))
        IntronPlot = exon_plot + geom_rect(data = exon_skipping_introns, aes(xmin=i.start, xmax=i.end, ymin=ystart, ymax=yend), fill="red", alpha=1)
        if(LABEL_OVERLAPS==TRUE) {
            IntronPlot = IntronPlot+geom_text(data = exon_skipping_introns, aes(x=meanX, y=meanY, label=ID))
        }
    } else { IntronPlot = exon_plot }
    
    IntronPlot = IntronPlot + geom_rect(data = Coord_Groups, aes(xmin=i.start, xmax=i.end, ymin=ystart, ymax=yend), fill="blue", alpha=1)+ylab("")+
        xlab("Genomic position")
    if(LABEL_INTRONS==TRUE) {
        IntronPlot = IntronPlot+geom_text(data = Coord_Groups, aes(x=meanX, y=meanY, label=ID))
    }
    #Finalise, add yaxis dual labels
    IntronPlot = ggdraw(IntronPlot+theme(plot.margin = margin(0, 0, 0, 8, unit = "pt"))+
                            geom_label(aes(label = GENE_TO_PROCESS, x=exonsnearDEIntrons[,mean(start)], y=0.5)))+
        cowplot::draw_label(x = 0.015, y = 0.75, angle = 90, label = "DE introns", size = geom.text.size)+
        cowplot::draw_label(x = 0.015, y = 0.35, angle = 90, label = "Exon overlap", size = geom.text.size)

    Coord_Groups.DensityPlot = copy(Coord_Groups)
    Coord_Groups.DensityPlot[,ShiftStart:=data.table::shift(i.start)]
    Coord_Groups.DensityPlot[,Subgroup:=ifelse(i.start-ShiftStart < 50, 1, 2)]
    Coord_Groups.DensityPlot[,ystart:=0.625]
    Coord_Groups.DensityPlot[,yend:=0.625]
    
    ASPLOT = ggplot()
    if (nrow(exon_skipping_introns) > 0) {
        exon_skipping_introns[,ShiftStart:=data.table::shift(i.start)]
        exon_skipping_introns[,Subgroup:=ifelse(i.start-ShiftStart < 50, 1, 2)]
        exon_skipping_introns[,ystart:=0.375]
        exon_skipping_introns[,yend:=0.375]
        
        exonsnearDEIntrons[seq(1,.N, 2), exon_number:=NA]
        ASPLOT = ggplot(exonsnearDEIntrons)+geom_rect(aes(xmin=start, xmax=end, ymin=0.375, ymax=0.625), alpha=1)+
            geom_text(aes(label=exon_number, x=mid, y=y_offset))+
            geom_hline(yintercept = 0.5, alpha=1)+theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank()) +
            geom_curve(data = Coord_Groups.DensityPlot[Subgroup%%2==0], aes(x=i.start, xend=i.end, y=ystart, yend=yend), curvature = -0.5, ncp=200, size=0.7)+
            geom_curve(data = exon_skipping_introns, aes(x=i.start, xend=i.end, y=ystart, yend=yend), curvature = 0.5, ncp=200, size=0.7, color="red")+
            coord_cartesian(ylim=c(-0.25,1.25))+xlab("")+ylab("")
        
    } else {
        message("Cannot draw AS plot - no AS introns.")    
    }
    
    return(list(IntronPlot, ASPLOT))
}

MakePrettyGeneMap <- function(GeneMap, DATASET, MODE){
    switch(DATASET,
           "LKB"= { 
               if (MODE=="LEAFCUTTER") { setnames( GeneMap, c("1W.y", "AD.y"), c("1W", "AD") ) }
               PaperCols = c("ID", "intron", "1W", "AD", "exonsSkipped", "UnorthodoxSplice", "Intron_AA_Start", "Intron_AA_End") 
           },
           "BHB"= { 
               if (MODE=="LEAFCUTTER") { setnames( GeneMap, c("C.y", "HF.y"), c("C", "HF") ) }
               PaperCols = c("ID", "intron", "C", "HF", "exonsSkipped", "UnorthodoxSplice", "Intron_AA_Start", "Intron_AA_End") 
           }
    )
    if (MODE=="LEAFCUTTER") { PaperCols = c(PaperCols, "deltapsi", "FDR") }
    if (MODE=="DESEQ2")     { PaperCols = c(PaperCols, "pValAdjust") }
    
    GeneMap.ForPaper = copy(GeneMap[, PaperCols, with=FALSE ])
    
    switch(DATASET,
           "LKB"= { setnames(GeneMap.ForPaper,  c("ID", "intron", "exonsSkipped",  "UnorthodoxSplice", "Intron_AA_Start", "Intron_AA_End", "1W", "AD"),
                             c("ID", "Intron", "Exons Skipped", "Splices outside exon boundary", "Protein Start (Sheep)", "Protein End (Sheep)", 
                               "Normalised counts (1W)", "Normalised counts (AD)")) },
           "BHB"= { setnames(GeneMap.ForPaper,  c("ID", "intron", "exonsSkipped",  "UnorthodoxSplice", "Intron_AA_Start", "Intron_AA_End", "C", "HF"),
                             c("ID", "Intron", "Exons Skipped", "Splices outside exon boundary", "Protein Start (Sheep)", "Protein End (Sheep)", 
                               "Normalised counts (Control)", "Normalised counts (Heart Failure)")) }
    )
    
    if (MODE=="LEAFCUTTER") { setnames(GeneMap.ForPaper, c("deltapsi", "FDR"), c("deltaPSI", "p-value (Cluster)")); }
    if (MODE=="DESEQ2")     { setnames(GeneMap.ForPaper, c("pValAdjust"), c("B&H adjusted p-value")) }
    
    switch(DATASET,
           "LKB"= { setcolorder(GeneMap.ForPaper, c("ID", "Intron", "Normalised counts (1W)", "Normalised counts (AD)", "Protein Start (Sheep)", 
                                                    "Protein End (Sheep)", "Exons Skipped", "Splices outside exon boundary",
                                                    if (MODE=="DESEQ2") c("B&H adjusted p-value") else c("deltaPSI", "p-value (Cluster)"))) },
           "BHB"= {setcolorder(GeneMap.ForPaper, c("ID", "Intron", "Normalised counts (Control)", "Normalised counts (Heart Failure)", "Protein Start (Sheep)", 
                                                   "Protein End (Sheep)", "Exons Skipped", "Splices outside exon boundary",
                                                   if (MODE=="DESEQ2") c("B&H adjusted p-value") else c("deltaPSI", "p-value (Cluster)"))) }
    )
    #ffwrite(GeneMap.ForPaper, "clipboard") #Export for thesis chapter / paper / spreadsheet
    return(GeneMap.ForPaper)
}

FinaliseASAnalysisAndMakeDataObject <- function(GENE_TO_PROCESS, transcript, GeneGTFData, AltSplIntrons, DATASET, MODE){
    #transcript="ENSOART00000017670"
    GTFData = GeneGTFData[transcript_id==transcript] #resolves multi-transcript problem, sort of.
    
    GTFData.CDS = unique(GTFData[feature=="CDS", .(chr, start, end, gene_id, exon_number)])
    setorder(GTFData.CDS, exon_number)

        GTFData.CDS[,Exon_Len := (end-start)+1]
    stopifnot(GTFData.CDS[,sum(Exon_Len)]%%3 == 0) #check if total length is cleanly divisible by 3 (3 nt to an AA)
    
    GTFData.CDS[,AA_Start:=data.table::shift(cumsum(Exon_Len/3), 1)]    #starts at previous exon's end
    GTFData.CDS[is.na(AA_Start), AA_Start:=1]                           #shift doesn't loop, inserts NA
    GTFData.CDS[,AA_End:=cumsum(Exon_Len/3)]                            #repeat. neater than reordering columns, if maybe slower.
    setkey(GTFData.CDS, chr, start, end)
    
    #Exons / introns may also map to the UTR
    GTFData.UTR = GTFData[start < GTFData[feature=="start_codon", start] | end > GTFData[feature=="stop_codon", end]]
    GTFData.UTR = GTFData.UTR[feature != "gene" & feature != "transcript"] #those would overlap -everything-
    setkey(GTFData.UTR, chr, start, end)

    #Generate Gene Map====
    tmsg("Merging count data and AS data...")
    ColsToMerge = c("intronID", "intron", "chr", "gene_id", "i.start", "i.end")
    if (MODE=="LEAFCUTTER") { 
        ColsToMerge = c(ColsToMerge, "deltapsi", "FDR")
        switch (DATASET,
                "LKB" = { ColsToMerge = c(ColsToMerge, "1W.y", "AD.y", "1W.x", "AD.x") },
                "BHB" = { ColsToMerge = c(ColsToMerge, "C.y", "HF.y", "C.x", "HF.x") }
        )
    }
    if (MODE=="DESEQ2") { 
        ColsToMerge = c(ColsToMerge, "pValAdjust")
        switch (DATASET,
                "LKB" = { ColsToMerge = c(ColsToMerge, "1W", "AD") },
                "BHB" = { ColsToMerge = c(ColsToMerge, "C", "HF") }
        )
    }
    
    ColsToKeep = c(ColsToMerge, "first.start", "first.end", "first.AA_Start", "first.AA_End", "first.exon_number")
    
    GeneMap.First = foverlaps(AltSplIntrons, GTFData.CDS, mult = "first", nomatch = NULL)
    setnames(GeneMap.First, c("start", "end", "AA_Start", "AA_End", "exon_number"), 
             c("first.start", "first.end", "first.AA_Start", "first.AA_End", "first.exon_number"))
    GeneMap.First = GeneMap.First[,..ColsToKeep]
    
    #TODO: if the overlap causes the end of the intron to exactly match the exon's start, check the other end instead?
    #there's a pattern of 'mistakes' not being caught...
    
    ColsToKeep = gsub("first", "last", ColsToKeep)
    
    GeneMap.Last = foverlaps(AltSplIntrons, GTFData.CDS, mult = "last", nomatch = NULL) 
    setnames(GeneMap.Last, c("start", "end", "AA_Start", "AA_End", "exon_number"), 
             c("last.start", "last.end", "last.AA_Start", "last.AA_End", "last.exon_number"))
    GeneMap.Last = GeneMap.Last[,..ColsToKeep]
    
    GeneMap = merge(GeneMap.First, GeneMap.Last, by=ColsToMerge)
    setcolorder(GeneMap, c("gene_id", "intronID", "intron", "chr", "i.start", "i.end", "first.start", "first.end", "last.start", "last.end", 
                           "first.AA_Start", "first.AA_End", "last.AA_Start", "last.AA_End"))
    
    #If we hit the splice site exactly, we can just use the AA coordinate of that
    GeneMap[i.start==first.end, Intron_AA_Start:=first.AA_End]
    GeneMap[i.end==last.start, Intron_AA_End:=last.AA_Start]
    
    #Otherwise, we have to recalculate AA coordinates. Remember to add +1 to get the real number of bp between coords
    GeneMap[is.na(Intron_AA_Start) | is.na(Intron_AA_End), UnorthodoxSplice:=TRUE]
    GeneMap[is.na(UnorthodoxSplice), UnorthodoxSplice := FALSE]
    
    #Intron_AA_Start must depend on the position of the next... exon? 
    GeneMap[UnorthodoxSplice == FALSE & is.na(Intron_AA_Start) & i.start<first.end, Intron_AA_Start:= first.AA_Start + (((i.start - first.start)+1)/3)] 
    GeneMap[UnorthodoxSplice == FALSE & is.na(Intron_AA_End)   & i.end>last.start,  Intron_AA_End:= last.AA_Start + (((i.end - last.start)+1)/3)]
    
    #GeneMap.UTR = foverlaps(AltSplIntrons, GTFData.UTR, nomatch = NULL)
    #if (nrow(GeneMap.UTR) > 0) { GeneMap = rbind(GeneMap, GeneMap.UTR[,.(intron, chr, i.start, i.end)], fill=T) } #insert UTR introns
    rm(GTFData.CDS, GTFData.UTR, GeneMap.First, GeneMap.Last)
    
    UniProtIDFile   = fread("Data/biomart/BIOMART-OvAr3.1-GeneIDTxIDUniProtID-FINAL.txt")
    UniProtIDFile   = unique(UniProtIDFile[,.(`Gene stable ID`, `Transcript stable ID`, `UniProtKB Gene Name ID`)])
    UniProtID = UniProtIDFile[`Transcript stable ID`==transcript][,unique(`UniProtKB Gene Name ID`)]
    rm(UniProtIDFile)
    if(length(UniProtID) > 1) { 
        warning(sprintf("More than one UniProtID found for transcript '%s'. If proceeding to 3D modelling, manually check data the protein is right !", transcript)) 
        UniProtID = UniProtID[1]
    }
    #TODO: this is reductive, is there a better way...?
    GeneMap[,`UniProtKB Gene Name ID`:=UniProtID]
    GeneMap[,exonsSkippedN:= (last.exon_number - first.exon_number)-1]
    GeneMap[,exonsSkipped:=""]
    GeneMap[exonsSkippedN>0, exonsSkipped:=paste(seq(from = first.exon_number+1, to = last.exon_number-1), collapse=", "), intron]
    
    GeneMap[,ID:=sprintf("%s-%d", if (MODE=="LEAFCUTTER") "L" else "D", .I)]

    #GRAPHING====
    GeneMap.AS = data.table()
    if (GeneMap[exonsSkippedN>0,.N]==0) {
        if(GeneMap[UnorthodoxSplice == TRUE,.N]>0) {
            sprintf("No AS introns found for %s, but 'odd' splices found.", GENE_TO_PROCESS) 
        } else {
            sprintf("No AS introns found for %s; DE Introns likely detected due to  DEG status, not AS.", GENE_TO_PROCESS) 
        }
    } else {
        GeneMap.AS = unique(GeneMap[exonsSkippedN>0])
        sprintf("%d AS introns found for %s.", nrow(GeneMap.AS), GENE_TO_PROCESS) 
    }
    return(list(GeneMap, GeneMap.AS))
}

AlternativeSplicingAnalysis <- function(GENE_TO_PROCESS, MODE="LEAFCUTTER", DATASET="LKB", COMPARISON=NULL, 
                                        DRAW_MODE = "REGION_ONLY", USE_FILTER = TRUE, IS_GENE_NAME=TRUE,
                                        LABEL_INTRONS=FALSE, LABEL_OVERLAPS=TRUE){
    #GENE_TO_PROCESS="RYR2";MODE="LEAFCUTTER";DATASET="LKB";COMPARISON=NULL;DRAW_MODE=2;USE_FILTER=TRUE;IS_GENE_NAME=TRUE
    #REGION_ONLY for exons surrounding DE introns only, FULL_GENE for full gene, a number for +/- n exons
    if ( !(DATASET %in% c("LKB", "BHB") )){
        stop("AlternativeSplicingAnalysis(): DATASET must be one of \"LKB\" or \"BHB\"")
    }
    if (!is.logical(USE_FILTER)) { stop("AlternativeSplicingAnalysis(): USE_FILTER must be TRUE or FALSE") }
    if( is.null(COMPARISON)) { COMPARISON = ifelse(DATASET=="LKB", "Timepoint|AD|1W", no = "CvHF") }
    if ( !(COMPARISON %in% c("Timepoint|AD|1W", "CvHF")) ) { 
        stop("To run AlternativeSplicingAnalysis, COMPARISON must be one of \"Timepoint|AD|1W\" or \"CvHF\"") 
    }
    if ( !(DRAW_MODE %in% c("FULL_GENE", "REGION_ONLY") | is.numeric(DRAW_MODE)) ) {
        stop("To run AlternativeSplicingAnalysis, DRAW_MODE must be one of \"FULL_GENE\", 
             \"REGION_ONLY\" or an integer representing the number of exons to draw") 
    }
    if ( !(MODE %in% c("LEAFCUTTER", "DESEQ2") )){ stop("AlternativeSplicingAnalysis(): MODE must be one of \"DESEQ2\" or \"LEAFCUTTER\"") }
    
    if(substr(GENE_TO_PROCESS, 1, str_length("ENSOARG"))=="ENSOARG" & IS_GENE_NAME==TRUE) { 
        tmsg("WARNING: IS_GENE_NAME is set to TRUE, but GENE_TO_PROCESS starts with 'ENSOARG'. Assuming user error and setting IS_GENE_NAME to FALSE.") 
        IS_GENE_NAME = FALSE
    }

    tmsg(sprintf("Beginning analysis.\nTarget gene\t[%s]\nDataset\t[%s]\nAnalysed with\t[%s]\nComparison\t[%s]\n", 
                 GENE_TO_PROCESS, DATASET, MODE, COMPARISON))
    
    tmsg("Loading auxiliary data...")
    #The OvArv3.1 parsed GTF was made using Auxiliary/6-ProcessGTFKeyValuePairs 
    #(which uses sapply()+str_match() onto the feature keyvaluepair string)
    exon_file       = fread("Data/biomart/parsedOvAr3.1GTF.tsv")
    Dictionary      = fread("Data/biomart/BIOMART-OvAr3.1-GeneIDNameDescrLocation-FINAL.txt")
    Dictionary[,`Gene start (bp)`:=as.numeric(`Gene start (bp)`)]
    Dictionary[,`Gene end (bp)`:=as.numeric(`Gene end (bp)`)]
    Dictionary[,`Gene description`:= gsub("\\[Source:.*?\\]", "", `Gene description`)]
    colnames(Dictionary) = c("GeneID", "chr", "start", "end", "geneSymbol", "Description")
    #Dictionary[,chr:=factor(chr, levels = Dictionary[,unique(chr)])]
    setkey(Dictionary, chr, start, end)
    
    if (IS_GENE_NAME == TRUE) {
        GENE_ID = Dictionary[geneSymbol==GENE_TO_PROCESS, GeneID]
        if(length(GENE_ID)<1) { stop(sprintf("Could not determine a sheep gene ID for '%s'. Please determine manually, set IS_GENE_NAME to FALSE, and retry.", GENE_TO_PROCESS)) }
        tmsg(sprintf("Translated gene name \"%s\" to ENSEMBL ID '%s'...", GENE_TO_PROCESS, GENE_ID))
    } else { GENE_ID = GENE_TO_PROCESS }
    
    if (MODE=="DESEQ2"){
        tmsg("Loading DESEQ2 data...")
        #DESeq2--
        switch (DATASET,
            "LKB" = { load("Data/leafcutter/LKB-DESeq2wIntrons-v1-Unfiltered.rda"); 
                      if(USE_FILTER == TRUE) { IntronsResults.Group = IntronsResults.Group[pValAdjust < 0.05] } 
                    },
            "BHB" = { load("Data/leafcutter/BHB-DESeq2wIntrons-v1.rda") }
        )
        
        IntronsResults.Group[,sourceGene:= gsub("-\\d+", "", geneID)]
        AltSplIntrons = IntronsResults.Group[Comparison==COMPARISON & grepl(pattern = paste0("^", GENE_ID), geneID), .(geneID, pValAdjust)]
        AltSplIntrons = IntronMap[AltSplIntrons, on="geneID"]
        setnames(AltSplIntrons, c("Intron.Chr", "Intron.Start", "Intron.Stop", "geneID"), c("chr", "i.start", "i.end", "intronID"))
        AltSplIntrons[,intron:=sprintf("%s:%s-%s", chr, i.start, i.end)]
    } 
    if (MODE=="LEAFCUTTER"){
        tmsg("Parsing LEAFCUTTER results...")
        #Read in data====
        INFILE_CS = if_else(DATASET=="LKB", true = "Data/leafcutter/1WAD-ONLY-m10_cluster_significance.txt", 
                                         false = "Data/leafcutter/BHB-m10-CvHF_cluster_significance.txt")
        cluster_significance = fread(INFILE_CS)
        setnames(cluster_significance, "p.adjust", "FDR")
        
        INFILE_ES = if_else(DATASET=="LKB", true = "Data/leafcutter/1WAD-ONLY-m10_effect_sizes.txt", 
                         false = "Data/leafcutter/BHB-m10-CvHF_effect_sizes.txt")
        effect_sizes = fread(INFILE_ES)
        rm(INFILE_ES, INFILE_CS)
        effectSizesSplit <-  as.data.frame(str_split_fixed(effect_sizes$intron, ":", 4), stringsAsFactors = FALSE )
        names(effectSizesSplit) <- c("chr","start","end","clusterID")
        effectSizes <- cbind( effect_sizes, effectSizesSplit)
        effectSizes$cluster <- paste(effectSizesSplit$chr, effectSizesSplit$clusterID, sep = ":")
        
        all.introns <- merge(x = cluster_significance, y = effectSizes, by = "cluster")
        rm(cluster_significance, effectSizes)
        
        all.introns <- all.introns[ order(all.introns$FDR),]
        #all.introns <- subset( all.introns, FDR <= 0.05 ) #reduces 23747 to 5429
        all.introns$start <- as.numeric(all.introns$start)
        all.introns$end <- as.numeric(all.introns$end)
        all.introns$chr <- gsub("chr", "", all.introns$chr)
        setkey(all.introns, chr, start, end)
        
        stats_overlap = foverlaps(all.introns, Dictionary, by.x=c("chr", "start", "end"), by.y=c("chr", "start", "end"))
        stats_overlap[,length(unique(GeneID))]
        
        stats_overlap.sig = stats_overlap[FDR<0.05]
        switch(DATASET,
            "LKB" = {
                stats_overlap.sig=stats_overlap.sig[,.(intron, clusterID, FDR, deltapsi, `1W`, AD, GeneID, 
                                                       geneSymbol, Description, chr, i.start, i.end)]
                stats_overlap.sig[,AbsDeltaPsi:=abs(deltapsi)]
                rm(effect_sizes, effectSizesSplit, stats_overlap, all.introns)
                AltSplIntrons = stats_overlap.sig[GeneID==GENE_ID, .(intron, chr, i.start, i.end, geneSymbol, `1W`, AD, deltapsi, FDR)]
            },
            "BHB" = {
                stats_overlap.sig=stats_overlap.sig[,.(intron, clusterID, FDR, deltapsi, C, HF, GeneID, 
                                                       geneSymbol, Description, chr, i.start, i.end)]
                stats_overlap.sig[,AbsDeltaPsi:=abs(deltapsi)]
                rm(effect_sizes, effectSizesSplit, stats_overlap, all.introns)
                AltSplIntrons = stats_overlap.sig[GeneID==GENE_ID, .(intron, chr, i.start, i.end, geneSymbol, C, HF, deltapsi, FDR)]
            }
        )
        setnames(AltSplIntrons, "intron", "intronID")
        AltSplIntrons[,intron:=sprintf("%s:%s:%s", chr, i.start, i.end)]
    }
    
    #processing----
    if(nrow(AltSplIntrons)<1) { tmsg("Analysis complete. No alternatively spliced introns detected."); return(NULL) }
    
    tmsg("Processing AS data...")
    
    AltSplIntrons[,overlap_start:=i.start-1]
    AltSplIntrons[,overlap_end:=i.end+1]      #expanding start and stop by 1 each forces them into the nearby exon(s)
    
    if (MODE=="DESEQ2"){
        #Count data is already loaded and merely has to be added to the results table
        Counts = as.data.table(Introns[which(rownames(Introns) %in% AltSplIntrons[,intronID]), ], keep.rownames = "intronID")
        if (nrow(AltSplIntrons)==1) { 
            setnames(Counts, c("V1", "V2"), c("variable", "value")) 
            Counts[,intronID:=AltSplIntrons[,intronID]]
        }
        else { Counts = melt.data.table(Counts, id.vars="intronID") }
    }
    if (MODE=="LEAFCUTTER"){
        #Current 'count data' is actualy % of the cluster, and we have to merge the real counts from the initial matrix.
        #They will be <timepoint>.y whereas the fractional 'leafcutter' counts will be the <timepoint>.x counts
        COUNTFILE = if_else(DATASET=="LKB", "Data/leafcutter/group1WvAD-m10(2)_perind_numers.counts/group1WvAD-m10_perind_numers.counts",
                            "Data/leafcutter/BHB/BHB-m10_perind_numers.counts")
        CountData = fread(COUNTFILE)
        CountData[,intronID:=paste0("chr", intronID)]
        #CountData[,intronID:=gsub(":clu_\\d+_NA$", "", intronID)]    
        Counts = CountData[intronID %in% AltSplIntrons[,intronID]]
        Counts = melt.data.table(Counts, id.vars="intronID")
    }
    EXTRACTSTR = if_else(DATASET=="LKB", true = "1W|AD", "C|HF")
    Counts[,Timepoint:=str_extract(variable, EXTRACTSTR)]
    Counts = Counts[!is.na(Timepoint)]
    Counts = dcast.data.table(Counts[,.(meanCounts=mean(value), n=.N),.(intronID, Timepoint)], intronID~Timepoint, value.var="meanCounts")
    AltSplIntrons = merge(AltSplIntrons, Counts, by="intronID", all.x=T)
    rm(Counts)
    setkey(AltSplIntrons, chr, overlap_start, overlap_end)
    
    GeneGTFData = exon_file[gene_id==GENE_ID]
    #stopifnot(GeneGTFData[transcript_id!="",length(unique(transcript_id))<2]) #TODO multi-transcript processing currently NOT supported
    FinalResults = list()
    #transcript=GeneGTFData[transcript_id!="",unique(transcript_id)]
    for (transcript in GeneGTFData[transcript_id!="",unique(transcript_id)]) {
        tmp = FinaliseASAnalysisAndMakeDataObject(GENE_TO_PROCESS, transcript, GeneGTFData, AltSplIntrons, DATASET, MODE)
        GeneMap          = tmp[[1]]
        GeneMap.ForPaper = MakePrettyGeneMap(GeneMap, DATASET, MODE)
        GeneMap.AS       = tmp[[2]]
        tmp2 = MakeASPlots(GeneMap, GeneGTFData, GENE_TO_PROCESS, transcript, DRAW_MODE, LABEL_INTRONS, LABEL_OVERLAPS)
        FinalResults = append(FinalResults, list(list("Gene"=GENE_TO_PROCESS, "Transcript"=transcript, "GeneMap"=GeneMap, 
                                                      "PublicationMap"=GeneMap.ForPaper, "ASPlot"=tmp2[[1]], "IntronPlot"=tmp2[[2]])) )
    }
    
    #Fin====
    tmsg(sprintf("Complete. Processed data for %d transcripts and output to list.", length(FinalResults)))
    return(FinalResults)
}
#TODO: long genes / genes which have a skip in exon 1 and exon 50 get bad x-scales even in REGION_ONLY.
#redefine 'region' as around each AS cluster? maybe sometihng like facet_wrap() with nrow=1 and a subset of surrounding n exons....

#a=AlternativeSplicingAnalysis("CACNA1C")