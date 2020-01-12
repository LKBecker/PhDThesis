#DeSeq2 RIP----
#see also: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3398141/
##libraries ----
require(DESeq2)
require(data.table)
require(stringr)
require(ggplot2)
#source("../../../../Code/R/DGE/20170911_DataImputationForDESeqModel.R")
##BEGIN ----
FileNames <- data.table(fileName = dir(path = "Data/RNAseq/Input/CountFiles.ERCC/", pattern = "*.counts", full.names = T))
FileNames[,Timepoint:=str_match(fileName, "Sample\\d{6}-\\d{2}_(\\w{1,2})_")[,2]]
FileNames[,Timepoint:=factor(Timepoint, levels = c("1W", "1M", "3M", "AD"))]
FileNames[,nTimepoint:=c("1W"=1,"1M"=4,"3M"=12,"AD"=78)[Timepoint]]
FileNames[,nTimepoint:=factor(nTimepoint)]
FileNames[,sample:=str_match(fileName, "Sample(\\d{6}-\\d{2})")[,2]]
FileNames[,shortName:= paste(sample, Timepoint, sep="_")]
FileNames[,batch:=str_match(sample, "\\d{4}(\\d{2})-\\d{2}")[,2]]
setcolorder(FileNames, unique(c("shortName", "fileName", colnames(FileNames))))
##DESeq2 data input----
DESeqData.LKB 				<- DESeqDataSetFromHTSeqCount(sampleTable = FileNames,  design = ~1) #Design not needed for count normalisation 
message("--Pre-filtering: removing genes with 0 or 1 counts...")
DESeqData.LKB <- DESeqData.LKB[ rowSums(counts(DESeqData.LKB)) > 1, ]
message("--Estimatine size factors to equalise libraries...")
DESeqData.LKB				<- estimateSizeFactors(DESeqData.LKB) #use ERCC?
DESeqData.LKB.NormalisedCounts <- counts(DESeqData.LKB, normalized=TRUE)
#Begin PCA----
SeqData <- data.table(DESeqData.LKB.NormalisedCounts)
#SeqData[,countSum := rowSums(DESeqData.LKB.NormalisedCounts)]
SeqData[,geneName := rownames(DESeqData.LKB.NormalisedCounts)]
setcolorder(SeqData, c("geneName", colnames(SeqData)[!(colnames(SeqData) %in% c("geneName"))]))

SeqData.M <- melt.data.table(SeqData, id.vars = "geneName") #[,!"countSum"]
SeqData.M[,group:=factor(substring(variable, 11, 12), levels = c("1W", "1M", "3M", "AD"))]

SeqData.Stat.All <- SeqData.M[,.(sumCount=sum(value)), .(geneName)]
a <- quantile(SeqData.Stat.All[,sumCount])
LowerLimit <- a['25%']
UpperLimit <- a['75%']

SeqData.Middle50<- SeqData[geneName %in% SeqData.Stat.All[(sumCount <= UpperLimit) & (sumCount >= LowerLimit), geneName]]

library(FactoMineR)
res.pca = PCA(t(SeqData.Middle50[,2:31]), scale.unit=TRUE, ncp=5, graph=F)

##START GGPLOT CODE TO BE ADJUSTED
PC1 <- res.pca$ind$coord[,1]
PC2 <- res.pca$ind$coord[,2]
labs <- rownames(res.pca$ind$coord)
PCs <- data.frame(cbind(PC1,PC2))
rownames(PCs) <- labs
ggplot(PCs, aes(PC1,PC2, label=rownames(PCs))) + geom_text() # Just showing the individual samples...
#This agrees with DEseq2's PCA

## Put samples and categorical variables (ie. grouping
## of samples) all together
vPC1 <- res.pca$var$coord[,1]
vPC2 <- res.pca$var$coord[,2]
vlabs <- rownames(res.pca$var$coord)
vPCs <- data.frame(cbind(vPC1,vPC2))
rownames(vPCs) <- vlabs
colnames(vPCs) <- colnames(PCs)
# and plot them
pv <- ggplot() +  theme_bw(base_size = 20) #opts(aspect.ratio=1) +
# no data so there's nothing to plot
# put a faint circle there, as is customary
angle <- seq(-pi, pi, length = 50)
df <- data.frame(x = sin(angle), y = cos(angle))
pv <- pv + geom_path(aes(x, y), data = df, colour="grey70")
#
# add on arrows and variable labels
pv <- pv + geom_text(data=vPCs, aes(x=vPC1,y=vPC2,label=rownames(vPCs)), size=4) + xlab("PC1") + ylab("PC2")
pv <- pv + geom_segment(data=vPCs, aes(x = 0, y = 0, xend = vPC1*0.9, yend = vPC2*0.9), arrow = arrow(length = unit(1/2, 'picas')), color = "grey30")
pv # show plot
###END
