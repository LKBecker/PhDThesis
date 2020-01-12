require(data.table)
require(aroma.light)
require(stringr)
library(ggplot2)

load("Data/RNAseq/Output/GWPM.GeneSummaryObjects.rda")
load("Data/RNAseq/Output/LKB.ERCC.GeneSummaryObjects.rda")

rm(Counts.ERCC.M, Counts.GWPM.M, DEGs.ERCC.1WAD, DEGs.ERCC.SplineDF3, DEGs.GWPM, ExpressPct, ExpressPct.GWPM, GeneExpression.S, SamplesPerGroup, ScaledCounts, ScaledCounts.GWPM)

#Combine two datasets into one

AllowedGenes=intersect(BHB.TPMCountTable[,geneID], LKB.TPMCountTable[,geneID])
BHB.TPMCountTable = BHB.TPMCountTable[geneID %in% AllowedGenes]
LKB.TPMCountTable = LKB.TPMCountTable[geneID %in% AllowedGenes]

CountTable = merge(BHB.TPMCountTable, LKB.TPMCountTable, by="geneID")
CountMatrix = as.matrix(CountTable[,-1])
rownames(CountMatrix)=CountTable[,geneID]

#diagnostic plot
CountTable.Melt = melt(CountTable, id.vars="geneID", variable.name = "sample", value.name = "TPM")
RankTable = CountTable.Melt[,.(sample, rank(TPM)), geneID]
RankTable[,group:=str_extract(sample, "C$|R$|HF|1W$|1M$|3M$|AD$")]
RankTable[,group:= factor(group, levels = c("1W", "1M", "3M", "AD", "C", "R", "HF"), ordered = T)]
RankTable[,sample:= factor(sample, levels = RankTable[order(group), unique(sample)], ordered = T)]
ggplot(RankTable, aes(x=sample, y=V2, color=group))+geom_boxplot()+theme(axis.text.x = element_text(angle=90), legend.position = "none")+ylab("Rank of TPMs")+xlab("")

#Utilise robust quartile normalisation to equalize TPMs and (forcibly) make datasets comparable, as a plot of ranks of TPMs per gene and sample has shown GWPM data is almost
#always lower in TPMs
CountMatrix.RQNorm = aroma.light::normalizeQuantileRank(CountMatrix, robust=TRUE)
#rownames(CountMatrix.RQNorm)
CountMatrix.RQNorm = as.data.table(CountMatrix.RQNorm, keep.rownames="geneID")

#I'm fussy and probably have OCD so lets sort the table. Then everything can be NEAT UND TIDY.
ColNames=data.table(col=colnames(CountMatrix))
ColNames[,grp:=str_extract(col, "1W|1M|3M|AD|C|R|HF")]
ColNames[,grp:=factor(grp, levels=c("1W", "1M", "3M", "AD", "C", "R", "HF"))]
ColNames[grp %in% c("C", "HF", "R"), spl:= as.numeric(str_extract(col, "^\\d{1,2}"))]
ColNames[is.na(spl), spl:=as.numeric(substr(col, 1,6))]
setorder(ColNames, grp, spl)
setcolorder(CountMatrix.RQNorm, c("geneID", ColNames[,col]))

saveRDS(CountMatrix.RQNorm, "Data/RNAseq/Output/TPMs.GWPMLKBOverlap.RDS")
