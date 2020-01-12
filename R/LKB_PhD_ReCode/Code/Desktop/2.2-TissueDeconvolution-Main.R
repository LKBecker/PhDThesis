#libraries and functions====
library(Biobase); library(bseqsc); library(data.table); library(DESeq2); library(ggplot2)
library(Matrix); library(plyr); library(reshape2); library(xbioc); require(MuSiC)
library(ape)
source("Code/Desktop/functions.r")

tmsg("Loading ExpressionSet objects prepared in file 2.1 ...")
AtrialTissueExpressionData = readRDS("Data/MuSiC/AtrialTissue-whole-RNA-seq-Data.rds")
Skelly.DeLaughterData = readRDS("Data/MuSiC/Skelly+DeLaughter-scRNAData-DESeqNorm.rds")

tmsg("Calculating cell type proportions...")
props <- music_prop(AtrialTissueExpressionData, Skelly.DeLaughterData, clusters="cluster", samples="cell")
tmsg("Exporting extimated cell type proportions...")
saveRDS(props, "Data/MuSiC/AtrialTissue-CellProportionsMkII.rds")

tmsg("Running music_basis() to estimate cell type 'distances'")
Skelly.DeLaughter.MusicBasis <- music_basis(Skelly.DeLaughterData, clusters="cluster", samples="cell", verbose = T)
tmsg("Exporting MuSiC basis data...")
saveRDS(Skelly.DeLaughter.MusicBasis, "Data/MuSiC/Skelly+DeLaughter.MusicBasis.rds")
