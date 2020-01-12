#Libraries====
library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)

#Read in data====
cluster_significance = fread("Data/leafcutter/1WAD-ONLY-m10_cluster_significance.txt")
setnames(cluster_significance, "p.adjust", "FDR")

effect_sizes = fread("Data/leafcutter/1WAD-ONLY-m10_effect_sizes.txt")
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

Dictionary = fread("Data/biomart/BIOMART-OvAr3.1-GeneIDNameDescrLocation-FINAL.txt")
Dictionary[,`Gene start (bp)`:=as.numeric(`Gene start (bp)`)]
Dictionary[,`Gene end (bp)`:=as.numeric(`Gene end (bp)`)]
Dictionary[,`Gene description`:= gsub("\\[Source:.*?\\]", "", `Gene description`)]
colnames(Dictionary) = c("GeneID", "chr", "start", "end", "geneSymbol", "Description")
#Dictionary[,chr:=factor(chr, levels = Dictionary[,unique(chr)])]
setkey(Dictionary, chr, start, end)

stats_overlap = foverlaps(all.introns, Dictionary, by.x=c("chr", "start", "end"), by.y=c("chr", "start", "end"), verbose = T)
stats_overlap[,length(unique(GeneID))]

stats_overlap.sig = stats_overlap[FDR<0.05]
stats_overlap.sig=stats_overlap.sig[,.(intron, clusterID, FDR, deltapsi, `1W`, AD, GeneID, geneSymbol, Description, chr, i.start, i.end)]
stats_overlap.sig[,AbsDeltaPsi:=abs(deltapsi)]
stats_overlap.sig[,length(unique(GeneID))] #1510

stats_overlap.sig[AbsDeltaPsi>=0.1,length(unique(GeneID))] #528
stats_overlap.sig[AbsDeltaPsi>=0.5,length(unique(GeneID))] #12

#TopClusters = stats_overlap.sig[clusterID %in% stats_overlap.sig[AbsDeltaPsi>=0.5,unique(clusterID)]]
#TopClusters = stats_overlap.sig[AbsDeltaPsi>=0.5]

#ggplot(stats_overlap.sig, aes(x=deltapsi))+geom_histogram()
#ggplot(stats_overlap.sig[AbsDeltaPsi>=0.1], aes(x=deltapsi))+geom_histogram()

#Results====
#CACNA1C = unique(stats_overlap.sig[geneSymbol=="CACNA1C"])
#MACF1
#vinculin
#tumor protein p53 inducible protein 11
#troponin I3, cardiac type
#titin
#thyroid hormone receptor associated protein 3
#tensin 1 2 3
#talin2
#syntaxin
#supervillin