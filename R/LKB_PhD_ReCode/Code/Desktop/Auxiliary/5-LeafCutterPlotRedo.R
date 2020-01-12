library(data.table)
library(dplyr)
library(foreach)
library(ggplot2)
library(intervals)
library(leafcutter)
library(reshape2)

#load("Data/leafcutter/1WAD.shinyapp.RData") #Example data; this archive is, of course, only exported by leafcutter itself because why should things be easy for you.
load("Data/leafcutter/1WAD_redo.RData")

#in introns table, you can get the delta PSI
clusters=as.data.table(clusters)
intronsDT = as.data.table(introns)
intronsDT[,absDeltaPsi := abs(deltapsi)]
a = merge(intronsDT, clusters[,.(clusterID, FDR)], by="clusterID")
#Clusters is pre-filtered according to FDR; all clusters reported are significant clusters

#below is the code of leafcutter::make_gene_plot, splayed out like a butterfly on cork.

#BIN1 - 0 junctions on chr2. Of 2409. What. 
#JPH2 - 0 junctions on chr13. Of 1153. What.
#CACNA1C - nothing, but ENSOARG00000013089 (CACNA1C) = 6!
#

BIN1Test = introns_to_plot[introns_to_plot$start >= 117267199 & introns_to_plot$end <= 117314475 & introns_to_plot$chr == 2, ] # nothing.

#config====
gene_name="BIN1"  
#Used to filter exons_table; trips error if not all on same chr; used to calculated gene start and end from first/last exon
clusterID = NULL        #optional; cluster_ids (derived from gene_name) is filtered by this assignment if not NULL. result most be >0 row
cluster_list = NULL     #used for FDR and labels if not NULL - optional
#introns = NULL          #Comes with RDA object; used for allEdges$deltaPSI and allEdges$deltaPSI to get... deltaPSI I guess
#introns_to_plot = NULL  #Comes with RDA object; data.frame of 5 columns: chr, start, end, clu(ster), middle (= start + ((end - start)/2) )
#counts                 #Comes with RDA object; data.frame of chr/scaffold:start:end_clusterID versus sample, value being counts (duh)
# exons_table = NULL    #Comes with RDA object; derived from GTF? 5 columns: chr, start, end, strand, gene_name
len = 500
min_exon_length = 0.5 
main_title = NA 
#snp_pos = NA            #Lol if you do SNPs; used for plotting snp specific data etc. for shortness, i removed all that since I don't use any SNP data
#snp = NA                #see snp_pos
summary_func = colSums  #I assume the function used to summarise... (scrolls code) nothing. Is never used again in this function. Ok.
legend_title = "Mean counts" # UNUSED
debug = TRUE            #Turns on verbose mode.
exonMax <- 1000

#grabbing cluster and exon data====
stopifnot(!is.null(exons_table)) ; stopifnot(!is.null(introns_to_plot)) #check for critical data existing

exons <- exons_table[exons_table$gene_name == gene_name, ] #grab exons only related to the gene_name desired
stopifnot(length(unique(exons$chr)) == 1) #all exons must be on the same chromosome or your annotation is wrong

if (debug) { cat(nrow(exons), "exons in gene", gene_name, "\nTesting cluster", clusterID, "\n") }
gene_start <- min(exons$start); gene_end <- max(exons$end) #coordinates of first and last exon, used to filter intron/cluster data
if (debug) { cat("Gene start:", gene_start, " end:", gene_end, "\n") }

introns_to_plot$chr = gsub("^chr", "", introns_to_plot$chr) # simple gsub to remove chr
myChr <- gsub("^chr", "", unique(exons$chr))                # same

#Filters clusters: their introns must be greater than gene_start (i.e. min(exon_start) ) and smaller than gene_end, and on the same chromosome
myclusters <- introns_to_plot[introns_to_plot$start >= gene_start & introns_to_plot$end <= gene_end & introns_to_plot$chr == myChr, ]
cluster_ids <- myclusters$clu #grabs cluster IDs from filtered clusters
if (debug) { cat(nrow(myclusters), "junctions found. Num junctions on chromosome ",  myChr, ":", sum(introns_to_plot$chr == myChr), "\n") }

stopifnot(nrow(myclusters) > 0) #gotta have at least one cluster!

if (!is.null(clusterID)) { stopifnot(clusterID %in% cluster_ids) } #if clusterID is given, must exist

#remove duplicates and create unified coordinate system====
exons <- exons[!duplicated(paste(exons$start, exons$end)), ]    # removes duplicates? duplicated() returns T/F vector of duplicate rows
exons <- exons[exons$end - exons$start <= exonMax, ]            # removes exons that are too long (default: 1000bp)
exons <- exons[order(exons$end), ]                              # orders by final position
myclusters$id <- "cluster"                                      # creates ID column, fills w/ string
all_junctions <- myclusters[, c("start", "end", "clu")]
exons$clu <- row.names(exons)                                   # converts row names to column
all_exons <- exons[, c("start", "end", "clu")]
length_transform <- function(g) log(g + 1)                      # safe log transformation, g+1 means you'll never try to log 0 (infinite result)
m <- reshape2::melt(all_junctions, id.vars = "clu")
s <- unique(m$value)                                            # should be... unique start and stops for all myclusters (introns_to_plot within the gene)
# if (!is.na(snp_pos)) {
#     SNP_pos <- as.numeric(str_split_fixed(snp_pos, ":", 2)[, 2])
#     s <- c(s, SNP_pos)
# }
s <- sort(s)                                                    
d <- s[2:length(s)] - s[1:length(s) - 1]                        # what? this... oh. it's not an exclude, it's an end minus start, each offset by 1...
trans_d <- length_transform(d)                                  # the length of each distance is then log transformed
coords <- c(0, cumsum(trans_d))                                 # and our coordinate system is established as 0 to cumsum of all transformed... lengths?
names(coords) <- s
total_length = sum(trans_d)                                     # hang on how are cumsum and sum different here
#they evaluate to same end result but diff steps! cumulative, after all. 
my_xlim = c(-0.05 * total_length, 1.05 * total_length)          # 5% 'margin' on coords
# if (!is.na(snp_pos)) { snp_coord <- coords[as.character(SNP_pos)] }
exons_here <- exons                                             # this isn't a copy command so isn't this just gonna do this... by reference??
exons_here$gene_name = factor(exons_here$gene_name)
invert_mapping = function(pos, s, coords, xlim) {
    if (pos %in% s) { coords[as.character(pos)] } else if (pos < min(s))  { xlim[1] } else if (pos > max(s)) { xlim[2] } else {
        w = which(pos < s[2:length(s)] & pos > s[1:(length(s) - 1)])
        stopifnot(length(w) == 1)
        coords[w] + (coords[w + 1] - coords[w]) * (pos - s[w])/(s[w + 1] - s[w])
    }
}
exon_df <- data.frame(x = sapply(exons_here$start, FUN = function(X) invert_mapping(pos = X, s = s, coords = coords, xlim = my_xlim)), 
                      xend = sapply(exons_here$end, FUN = function(X) invert_mapping(pos = X, s = s, coords = coords, xlim = my_xlim)), 
                      y = numeric(nrow(exons_here)), yend = numeric(nrow(exons_here)), label = exons_here$gene_name)

upstream_exons <- exons_here[exons_here$start < min(s), ]

if (nrow(upstream_exons) > 0) {
    upstream_s <- c(min(upstream_exons$start), max(upstream_exons$end))
    upstream_coords <- c(-length_transform(upstream_s[2] - upstream_s[1]), 0)
    names(upstream_coords) <- upstream_s
    upstream_df <- data.frame(x = sapply(upstream_exons$start, FUN = function(X) invert_mapping(pos = X, s = upstream_s, coords = upstream_coords, 
                                                                                                xlim = upstream_coords)), 
                              xend = sapply(upstream_exons$end, FUN = function(X) invert_mapping(pos = X, s = upstream_s, coords = upstream_coords, 
                                                                                                 xlim = upstream_coords)), 
                              y = 0, yend = 0, label = upstream_exons$gene_name)
    exon_df[1:nrow(upstream_df), ] <- upstream_df
}
downstream_exons <- exons_here[exons_here$end > max(s), ]
if (nrow(downstream_exons) > 0) {
    downstream_s <- c(min(downstream_exons$start), max(downstream_exons$end))
    downstream_coords <- c(max(coords), max(coords) + length_transform(downstream_s[2] - downstream_s[1]))
    names(downstream_coords) <- downstream_s
    downstream_df <- data.frame(x = sapply(downstream_exons$start, 
        FUN = function(X) invert_mapping(pos = X, s = downstream_s, coords = downstream_coords, xlim = downstream_coords)), 
        xend = sapply(downstream_exons$end, FUN = function(X) invert_mapping(pos = X, s = downstream_s, coords = downstream_coords, 
        xlim = downstream_coords)), y = 0, yend = 0, label = downstream_exons$gene_name)
    exon_df[((nrow(exon_df) + 1) - nrow(downstream_df)):nrow(exon_df), ] <- downstream_df
}

exon_df[(exon_df$xend - exon_df$x) < min_exon_length, ]$xend <- exon_df[(exon_df$xend - exon_df$x) < min_exon_length, ]$x + min_exon_length
my_xlim <- c(min(exon_df$x), max(exon_df$xend))
gene_length <- max(exon_df$xend) - min(exon_df$x)
new_theme_empty <- theme_bw(base_size = 15)
new_theme_empty$line <- element_blank()
new_theme_empty$rect <- element_blank()
new_theme_empty$strip.text <- element_blank()
new_theme_empty$axis.text <- element_blank()
YLIMN = 1000; YLIMP = -1000; curv = 0.35; curveMax = 0.1
curveExponent = 2; yOffset = 0; centreLineWidth = 3; junction_colour <- "#d66464"
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
min_height = 0; max_height = -1; yFactor = 0.8
allEdges = do.call(rbind, foreach(i = 1:nrow(all_junctions)) %do% {
    if (i%%2 == 1) { return(NULL) }
    start = coords[as.character(all_junctions$start[i])]
    end = coords[as.character(all_junctions$end[i])]
    l = end - start
    h = (1 + sqrt(l)) * ((i%%2) * 2 - 1)
    min_height = min(min_height, h)
    max_height = max(max_height, h)
    edge = data.frame(start, end)
    edge$startv <- all_junctions$start[i]
    edge$endv <- all_junctions$end[i]
    edge$start <- start
    edge$end <- end
    edge$clu <- all_junctions$clu[i]
    edge$Group <- i
    edge$xtext <- start + l/2
    edge$ytext <- -l^(yFactor)/2 + 0.5
    edge
    })

allEdgesP = do.call(rbind, foreach(i = 1:nrow(all_junctions)) %do% {
    if (i%%2 == 0) { return(NULL) }
    start = coords[as.character(all_junctions$start[i])]
    end = coords[as.character(all_junctions$end[i])]
    l = end - start
    h = (1 + sqrt(l)) * ((i%%2) * 2 - 1)
    min_height = min(min_height, h)
    max_height = max(max_height, h)
    edge = data.frame(start, end)
    edge$startv <- all_junctions$start[i]
    edge$endv <- all_junctions$end[i]
    edge$start <- start
    edge$end <- end
    edge$clu <- all_junctions$clu[i]
    edge$Group <- i
    edge$xtext <- start + l/2
    edge$ytext <- l^(yFactor)/2 - 0.5
    edge
    })

#deltaPSI here!
if (!is.null(introns)) {
    allEdges$deltaPSI <- introns$deltapsi[match(paste(allEdges$clu, allEdges$startv, allEdges$endv), 
                                                paste(introns$clusterID, introns$start, introns$end))]
    allEdgesP$deltaPSI <- introns$deltapsi[match(paste(allEdgesP$clu, allEdgesP$startv, allEdgesP$endv), 
                                                 paste(introns$clusterID, introns$start, introns$end))]
}

junctions <- rbind(allEdges, allEdgesP)

label_df <- as.data.frame(group_by(junctions, clu) %>% summarise(start = min(start), 
                                                                 end = max(end), middle = start + ((end - start)/2), ytext = max(ytext)))
if (!is.null(cluster_list)) {
    label_df$FDR <- cluster_list$FDR[match(label_df$clu, cluster_list$clusterID)]
    label_df$FDR[is.na(label_df$FDR)] <- "."
    label_df$label <- paste0(label_df$clu, "\n", label_df$FDR)
    label_df$label <- gsub("\n.$", "", label_df$label)
}
n_genes <- seq(1, length(levels(exons_here$gene_name)))
gene_name_df <- data.frame(x = total_length * 0.01, y = YLIMP - 0.1 * YLIMP, label = levels(exons_here$gene_name))
gene_strand <- unique(exons_here$strand)
if (length(gene_strand) > 1) { gene_strand <- NA }
if (!is.na(gene_strand)) {
    strand_df <- exon_df[exon_df$xend - exon_df$x >= 1, ]
    strand_df <- strand_df[!duplicated(paste(strand_df$x, strand_df$xend)), ]
    strand_df$id <- rownames(strand_df)
    arrowFactor = 1
    strand_df$x = strand_df$x + arrowFactor
    strand_df$xend <- strand_df$xend - arrowFactor
    strand_df <- strand_df[strand_df$xend - strand_df$x >= 1, ]
    if (gene_strand == "+") { strand_pos <- "last"; strand_df <- strand_df[1, ]; strand_df <- melt(strand_df[, c("x", "xend", "id")], "id") }
    if (gene_strand == "-") { strand_pos <- "first";strand_df <- strand_df[nrow(strand_df), ]; strand_df <- melt(strand_df[, c("x", "xend", "id")], "id") }
} else { strand_pos <- NULL }
exon_intervals <- Intervals(matrix(data = c(exon_df$x, exon_df$xend), ncol = 2))
intron_intervals <- intervals::interval_complement(exon_intervals)
intron_intervals <- intron_intervals[2:(nrow(intron_intervals) - 1), ]
strand_df <- as.data.frame(intron_intervals)
strand_df <- strand_df[(strand_df$V2 - strand_df$V1) > 0.025 * total_length, ]
strand_df$midpoint <- strand_df$V1 + (strand_df$V2 - strand_df$V1)/2
group <- c()
for (i in 1:nrow(strand_df)) { group <- c(group, rep(i, 2)) }
if (gene_strand == "+") { 
    strand_df <- data.frame(x = c(rbind(strand_df$V1, strand_df$midpoint)),  group = group, y = 0) 
}
if (gene_strand == "-") {
    strand_df <- data.frame(x = c(rbind(strand_df$midpoint, strand_df$V2)), group = group, y = 0)
}
# if (!is.na(snp_pos)) {
#     SNP_df <- data.frame(x = snp_coord - (min_exon_length/2), 
#                          xend = snp_coord + (min_exon_length/2), y = 0, yend = 0, label = snp)
# }
plots <- ggplot()
if (!is.null(cluster_list)) {
    allEdgesP = allEdgesP %>% left_join(label_df %>% select(clu, FDR), by = "clu")
    allEdgesP_significant = allEdgesP %>% dplyr::filter(FDR != ".")
    if (nrow(allEdgesP_significant) > 0) 
        plots <- plots + geom_curve(data = allEdgesP_significant, aes(x = start, xend = end, y = 0, yend = 0, group = Group, 
                                        color = clu, size = curveMax), curvature = curv, lineend = "round", colour = junction_colour)
    allEdges = allEdges %>% left_join(label_df %>% select(clu, FDR), by = "clu")
    allEdges_significant = allEdgesP %>% dplyr::filter(FDR != ".")
    if (nrow(allEdges_significant) > 0) 
        plots <- plots + geom_curve(data = allEdges_significant, aes(x = start, xend = end, y = 0, yend = 0, group = Group, 
                                        color = clu, size = curveMax), curvature = -curv, lineend = "round", colour = junction_colour)
    allEdgesP_nonsignificant = allEdgesP %>% dplyr::filter(FDR == ".")
    if (nrow(allEdgesP_nonsignificant) > 0) {
        plots <- plots + geom_curve(data = allEdgesP_nonsignificant, aes(x = start, xend = end, y = 0, yend = 0, group = Group, 
                                        color = clu, size = curveMax), curvature = curv, lineend = "round", colour = "gray")
    }
    allEdges_nonsignificant = allEdges %>% dplyr::filter(FDR == ".")
    if (nrow(allEdges_nonsignificant) > 0) {
        plots <- plots + geom_curve(data = allEdges_nonsignificant, aes(x = start, xend = end, y = 0, yend = 0, group = Group, 
                                        color = clu, size = curveMax), curvature = -curv, lineend = "round", colour = "gray")
    }
} else {
    if (is.null(clusterID)) {
        plots <- plots + geom_curve(data = allEdgesP, aes(x = start, 
            xend = end, y = 0, yend = 0, group = Group, color = clu, size = curveMax), curvature = curv, lineend = "round", colour = junction_colour) + 
            geom_curve(data = allEdges, aes(x = start, xend = end, y = 0, yend = 0, group = Group, color = clu, size = curveMax), curvature = -curv, 
                        lineend = "round", colour = junction_colour)
    }
    else {
        plots <- plots + geom_curve(data = allEdgesP, aes(x = start, xend = end, y = 0, yend = 0, group = Group, color = factor(clu == clusterID), 
            size = curveMax), curvature = curv, lineend = "round") + geom_curve(data = allEdges, aes(x = start, xend = end, y = 0, yend = 0, group = Group, 
            color = factor(clu == clusterID), size = curveMax), curvature = -curv, lineend = "round") + scale_colour_manual("", breaks = c(TRUE, FALSE), 
            limits = c(TRUE, FALSE), values = c("firebrick2", "gray")) + guides(colour = FALSE)
    }
}
if (!is.null(introns) & !is.null(cluster_list)) {
    plots <- plots + geom_curve(data = allEdgesP_significant, aes(x = start, xend = end, y = 0, yend = 0, group = Group, color = as.factor(deltaPSI > 0), 
            size = curveMax), curvature = curv, lineend = "round") + geom_curve(data = allEdges_significant, aes(x = start, xend = end, y = 0, yend = 0, 
            group = Group, color = as.factor(deltaPSI > 0), size = curveMax), curvature = -curv, lineend = "round") + 
            scale_colour_manual("dPSI", labels = c("down", "up"), values = c("darkturquoise", "firebrick2"))
}
plots <- plots + new_theme_empty + xlab("") + ylab("") + 
    geom_hline(yintercept = 0, size = centreLineWidth, colour = "white") + 
    geom_hline(yintercept = 0, alpha = 0.9, size = 1) + ylim(YLIMN, YLIMP) + scale_size_continuous(limits = c(0, 10), guide = "none") + 
    ggtitle(paste(gene_name_df$label, collapse = "+")) + theme(plot.title = element_text(face = "bold.italic", colour = "black", size = 20)) + 
    geom_segment(data = exon_df, aes(x = x, y = y, xend = xend, yend = yend), alpha = 1, size = 6, colour = "black")
if (!is.null(strand_pos)) {
    plots <- plots + geom_line(data = strand_df, aes(x = x, y = y, group = group), colour = "black", size = 1, 
            arrow = arrow(ends = strand_pos, type = "open", angle = 30, length = unit(0.1, units = "inches")))
}
plots <- plots + geom_segment(data = exon_df, aes(x = xend, xend = xend + 0.05, y = 0, yend = 0), colour = "white", 
                              size = 6) + geom_segment(data = exon_df, aes(x = x, xend = x - 0.05, y = 0, yend = 0), colour = "white", size = 6)
if (!is.null(cluster_list)) {
    label_df$label <- gsub("_[+-]", "", label_df$label)
    label_df$label <- gsub("_", "\n", label_df$label)
    label_df <- label_df[order(label_df$middle), ]
    label_df$labelY <- 0
    label_df$labelY[seq(1, nrow(label_df), 2)] <- YLIMN - 0.2 * YLIMN
    if (nrow(label_df) > 1) { label_df$labelY[seq(2, nrow(label_df), 2)] <- YLIMP - 0.2 * YLIMP }
    
    plots <- plots + geom_segment(data = label_df, aes(x = start, xend = middle, y = 0, yend = labelY), colour = "gray", 
                linetype = 3) + geom_segment(data = label_df, aes(x = end, xend = middle, y = 0, yend = labelY), colour = "gray", 
                linetype = 3) + geom_point(data = label_df, aes(x = middle, y = labelY), colour = "white", size = 22)
    if (!is.null(clusterID)) {
        plots <- plots + geom_text_repel(data = label_df[label_df$clu != 
            clusterID, ], aes(x = middle, y = labelY, label = label), point.padding = NA, direction = "y", segment.alpha = 0) + 
            geom_label(data = label_df[label_df$clu == clusterID, ], aes(x = middle, y = labelY, label = label), 
                       fontface = "bold", label.size = 0.5, label.r = unit(0.3, "lines"), label.padding = unit(0.3, "lines"))
    }
    else { 
        plots <- plots + geom_text_repel(data = label_df, aes(x = middle, y = labelY, label = label), point.padding = NA, direction = "y", segment.alpha = 0)
    }
}
# if (!is.na(snp_pos)) {
#     plots <- plots + geom_segment(data = SNP_df, aes(x = x, y = y, xend = xend, yend = yend), colour = "goldenrod", size = 6.5) + 
#                     geom_text(data = SNP_df, aes(x = x, y = 0.5 * YLIMP, label = label)) + 
#                     geom_segment(data = SNP_df, aes(x = x, xend = x, y = 0, yend = 0.45 * YLIMP), colour = "gray", linetype = 3)
#}
plots

#How many genes do have >0 introns??
exons_table_dt = as.data.table(exons_table)
MaxBoundaries = unique(exons_table_dt[,.(exons_start=min(start), exons_end=max(end), exon_chr=chr), gene_name])
introns_dt = as.data.table(introns_to_plot)
NGenesDT = data.table()
for (i in 1:nrow(MaxBoundaries)) {
    coords=MaxBoundaries[i]
    exon_start = coords[,exons_start]
    exon_end = coords[,exons_end]
    exon_chr = coords[,exon_chr]
    gene_name = coords[,gene_name]
    #message(sprintf("Testing for introns on chr %s, between %d and %d", exon_chr, exon_start, exon_end))
    PossIntrons = introns_dt[start >= exon_start & end <= exon_end & chr == exon_chr]
    #message(sprintf("Found %d introns", PossIntrons[,.N]))
    NGenesDT = rbind(NGenesDT, data.table(geneName = gene_name, nIntrons = PossIntrons[,.N]), fill=T)
}
#6306. That's... weirdly low??
#Most of the high-intron genes are u6 or SNORNA...