#Libraries====
library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)

Dictionary = fread("Data/biomart/BIOMART-OvAr3.1-GeneIDNameDescrLocation-FINAL.txt")
Dictionary[,`Gene start (bp)`:=as.numeric(`Gene start (bp)`)]
Dictionary[,`Gene end (bp)`:=as.numeric(`Gene end (bp)`)]
Dictionary[,`Gene description`:= gsub(" \\[Source:.*?\\]", "", `Gene description`)]
colnames(Dictionary) = c("GeneID", "chr", "start", "end", "geneSymbol", "Description")
setkey(Dictionary, chr, start, end)

AnalyzeASData = function(combination){
    #Read in data====
    stopifnot(combination %in% c("CvHF", "RvC", "HFvR"))
    cluster_significance = fread(sprintf("Data/leafcutter/BHB-m10-%s_cluster_significance.txt", combination))
    setnames(cluster_significance, "p.adjust", "FDR")
    
    effect_sizes = fread(sprintf("Data/leafcutter/BHB-m10-%s_effect_sizes.txt", combination))
    effectSizesSplit <-  as.data.frame(str_split_fixed(effect_sizes$intron, ":", 4), stringsAsFactors = FALSE )
    names(effectSizesSplit) <- c("chr","start","end","clusterID")
    effectSizes <- cbind( effect_sizes, effectSizesSplit)
    effectSizes$cluster <- paste(effectSizesSplit$chr, effectSizesSplit$clusterID, sep = ":")
    
    all.introns <- merge(x = cluster_significance, y = effectSizes, by = "cluster")
    all.introns <- all.introns[ order(all.introns$FDR),]
    all.introns$start <- as.numeric(all.introns$start)
    all.introns$end <- as.numeric(all.introns$end)
    all.introns$chr <- gsub("chr", "", all.introns$chr)
    setkey(all.introns, chr, start, end)
    
    stats_overlap = foverlaps(all.introns, Dictionary, by.x=c("chr", "start", "end"), by.y=c("chr", "start", "end"), verbose = T)
    stats_overlap[,length(unique(GeneID))]
    
    stats_overlap.sig = stats_overlap[FDR<0.05]
    setnames(stats_overlap.sig, unlist(str_split(combination, "v", 2)), c("TEMP1", "TEMP2"))
    stats_overlap.sig=stats_overlap.sig[,.(GeneID, geneSymbol, Description, clusterID, intron, chr, i.start, i.end, FDR, deltapsi, TEMP1, TEMP2)]
    setnames(stats_overlap.sig, c("TEMP1", "TEMP2"), unlist(str_split(combination, "v", 2)))
    stats_overlap.sig[,AbsDeltaPsi:=abs(deltapsi)]
    setorder(stats_overlap.sig, -AbsDeltaPsi)
    stats_overlap.sig
}

CvHF = AnalyzeASData("CvHF")
HFvR = AnalyzeASData("HFvR")
RvC = AnalyzeASData("RvC")
rm(Dictionary)
