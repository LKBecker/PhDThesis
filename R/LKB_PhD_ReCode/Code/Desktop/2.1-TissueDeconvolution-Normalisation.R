#libraries and functions====
library(data.table); library(plyr); library(reshape2); library(xbioc); require(MuSiC)

#Therefore the number produced by music_S() are average library size for all subjects and cell types. 
#Let S_e^1, ..., S_e^N be the average library size of endothelial cells for dataset 1 and S'_e^1, ..., S'_e^M for dataset2. 
#For any subject j in dataset1, divide the average library size of all cell types for subject j by S_e^j. Do the same thing to dataset 2. 
#Then the two new matrices are comparable. We can take the average across subject for all cell types and use that as the library size factor for constructing design matrix.

#load("Data/MuSiC/Incomplete.MusicBasis+Meta.rda")
LibSizes = melt(IncompleteCardiacCellSet.Basis$S)
setDT(LibSizes)
LibSizes = LibSizes[value!="NaN"]
LibSizes = merge(LibSizes, IncompleteCombinedCardiacCells.Meta, by.x="Var1", by.y="cell")
LibSizes[,Var2:=NULL]

#For any subject j in dataset1, divide the average library size of all cell types for subject j by S_e^j. Do the same thing to dataset 2. 
MeanSharedCellLibSize = LibSizes[cluster == "Fibroblasts" | cluster=="Endothelial cells"|cluster=="fibroblast"|cluster=="endothelial", .(meanSize=mean(value), NLibs=.N), .(Individual, cluster)]
MeanSharedCellLibSize.WeightedAvg = MeanSharedCellLibSize[, .(weightedMeanLibSize = sum(meanSize * (NLibs/sum(NLibs)))), Individual]
MeanLibSizePerIndiv = LibSizes[,.(meanLibSize=mean(value)), Individual]
rm(LibSizes)
SizeFactors = merge(MeanSharedCellLibSize.WeightedAvg, MeanLibSizePerIndiv, by="Individual")
SizeFactors[,SizeFactor := weightedMeanLibSize / meanLibSize]
SizeFactors
