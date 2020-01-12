library(data.table)
library(DESeq2)
library(dplyr)
library(foreach)
library(ggplot2)
library(splines)
library(Hmisc)
source("Code/Desktop/functions.R")
source("Code/Desktop/functions.DESeq2.r")
tmsg("Preparing to run 5.4-LeafcutterEsqueViaIntronCounts.R, version 0...")
#Loading in data====
load("Data/leafcutter/DESeq2wIntronsv1.rda")

#ASStats = rbind(IntronsResults.Group[pValAdjust<0.05,.N,Comparison])

MakeSplicePlot <- function(gene, IntronMapData, Timepoints = c("1W", "AD"), evalLen=500, legend_title=gene, main_title=gene, exons_table=NULL){
    #gene = "ENSOARG00000001811"; Timepoints = c("1W", "AD"); IntronMapData = IntronMap; evalLen=500; legend_title=gene; main_title=gene
    new_theme_empty <- theme_bw(base_size = 14)
    new_theme_empty$line <- element_blank()
    new_theme_empty$rect <- element_blank()
    new_theme_empty$strip.text <- element_blank()
    new_theme_empty$axis.text <- element_blank()
    
    SpecificCounts = Introns[,which(colnames(Introns) %in% Metadata[Timepoint %in% Timepoints, SampleName])]
    GeneCounts = t(SpecificCounts[which(rownames(SpecificCounts) %like% gene),]) #y of original function
    #GeneCounts2 = GeneCounts
    #GeneCounts2[,w]
    Types = substr(rownames(GeneCounts), str_length(rownames(GeneCounts))-1, str_length(rownames(GeneCounts)))
    
    GeneCounts.M = as.data.table(melt(GeneCounts))
    
    IntronMapData = IntronMapData[`Gene stable ID`==gene]
    IntronMapData = merge(IntronMapData, GeneCounts.M[,.(counts = mean(value)), Var2], by.x="geneID", by.y="Var2")
    
    max_log = 0.5 * ceiling(2 * log10(1 + max(IntronMapData[,counts]) ) )
    breaks = ifelse(max_log <= 2.5, seq(0, max_log, by = 0.5), seq(0, ceiling(max_log), by = 1)) 
    limits = c(0, max_log)
    IntronMapData[,id := as.factor(.I)]
    IntronCoords = sort( unique(c(IntronMapData[, Intron.Start], IntronMapData[,Intron.Stop])) )
    IntronLengths = IntronCoords[2:length(IntronCoords)] - IntronCoords[1:length(IntronCoords) - 1]
    trans_InLens = log(IntronLengths+1)
    coords = c(0, cumsum(trans_InLens)) #gives intervals for each measured
    names(coords) = IntronCoords
    total_length = sum(trans_InLens)
    my_xlim = c(-0.05 * total_length, 1.05 * total_length)
    first_plot = T
    min_height = 0
    max_height = 0
    
    #plots =
    plots = foreach(tis = Timepoints) %do% { #for each tis in Timepoints
        print(tis)
        group_sample_size = sum(tis == Types) #sums true values, Types is a vector...
        currentGroupCounts = GeneCounts
        allEdges = do.call(rbind, foreach(i = 1:nrow(IntronMapData)) %do% {
            if (IntronMapData$counts[i] == 0) { return(NULL) }
            start = coords[as.character(IntronMapData$Intron.Start[i])]
            end = coords[as.character(IntronMapData$Intron.Stop[i])]
            intron_length = end - start
            intron_height = (1 + sqrt(intron_length)) * ((i%%2) * 2 - 1)
            min_height = min(min_height, intron_height)
            max_height = max(max_height, intron_height)
            #edge describes a smooth bezier curve from start to end, constrained to min and max height
            edge = data.frame(Hmisc::bezier(
                x = c(start, start - 0.1 * total_length, (start + end)/2, end + 0.1 * total_length, end), 
                y = c(0, intron_height * 2/3, intron_height, intron_height * 2/3, 0), 
                evaluation = evalLen
            ))
            edge$Sequence <- log10(1 + IntronMapData$counts[i]) * sin(seq(0, pi, length.out = evalLen))
            edge$log10counts = log10(IntronMapData$counts[i])
            edge$Group <- i
            edge
        })
        
        if (is.na(gene) | !first_plot) { new_theme_empty$plot.title <- element_blank() }
        first_plot = F
        g = ggplot(allEdges) + geom_path(aes(x = x, y = y, group = Group, colour = log10counts, size = 0.1, alpha = 1)) + 
            scale_size(breaks = breaks, labels = format(10^breaks, digits = 0), limits = limits, 
                       range = c(0.3, 10), guide = guide_legend(title = legend_title)) + 
            scale_alpha(guide = "none", range = c(0.1, 1)) + new_theme_empty + 
            scale_color_gradient(breaks = breaks, limits = limits, labels = format(10^breaks, digits = 0), 
                                 low = "blue", high = "red", guide = guide_legend(title = legend_title)) + 
            ylab(paste0(tis, " (n=", group_sample_size, ")")) + xlab("") + xlim(my_xlim) + geom_hline(yintercept = 0, alpha = 0.3)
        g
    }
    
    df = data.frame(x = coords, xend = total_length * (IntronCoords - min(IntronCoords))/(max(IntronCoords) - min(IntronCoords)), y = 0, yend = min_height)
    plots[[length(plots)]] = plots[[length(plots)]] + geom_segment(data = df, aes(x = x, y = y, xend = xend, yend = yend), alpha = 0.1)
    if (!is.null(exons_table)) {
        exons_chr = exons_table[exons_table == intron_meta$chr[1], ]
        exons_here = exons_chr[(min(s) <= exons_chr$start & exons_chr$start <= max(s)) | (min(s) <= exons_chr$end & exons_chr$end <= max(s)), ]
        exons_here$gene_name = factor(exons_here$gene_name)
        gene_heights = min_height - ((1:length(levels(exons_here$gene_name))) - 1) * abs(min_height) * 0.15
        heights = gene_heights[as.numeric(exons_here$gene_name)]
        df = data.frame(x = total_length * (exons_here$start - min(s))/(max(s) - min(s)), xend = total_length * 
                            (exons_here$end - min(s))/(max(s) - min(s)), y = heights, yend = heights)
        if (nrow(exons_here) > 0) 
            plots[[length(plots)]] = plots[[length(plots)]] + 
                geom_segment(data = df, aes(x = x, y = y, xend = xend, yend = yend), alpha = 0.3, size = 5) + 
                geom_hline(yintercept = min_height, alpha = 0.3) + 
                geom_text(data = data.frame(x = my_xlim[1], y = gene_heights, label = levels(exons_here$gene_name)), aes(x, y, label = label))
        invert_mapping = function(pos) if (pos %in% s) 
            coords[as.character(pos)]
        else if (pos < min(s)) { my_xlim[1] }
        else if (pos > max(s)) { my_xlim[2] }
        else {
            w = which(pos < s[2:length(s)] & pos > s[1:(length(s) - 1)])
            stopifnot(length(w) == 1)
            coords[w] + (coords[w + 1] - coords[w]) * (pos - s[w])/(s[w + 1] - s[w])
        }
        if (nrow(exons_here) > 0) {
            df = data.frame(x = sapply(exons_here$start, invert_mapping), xend = sapply(exons_here$end, invert_mapping), y = 0, yend = 0)
            for (i in 1:length(plots)) plots[[i]] = plots[[i]] + geom_segment(data = df, aes(x = x, y = y, xend = xend, yend = yend), alpha = 0.3, size = 5)
        }
    } #Exons
    if (!is.na(main_title)) 
        plots[[1]] = plots[[1]] + ggtitle(main_title)
    do.call(gridExtra::grid.arrange, c(plots, list(ncol = 1)))
}

Top10 = IntronsResults.Group[Comparison=="Timepoint|AD|1W" & pValAdjust < 0.05][order(baseMean, decreasing = T)][1:10]
MakeSplicePlot(gsub("-\\d{1,2}", "", Top10[1, geneID]), IntronMap)
