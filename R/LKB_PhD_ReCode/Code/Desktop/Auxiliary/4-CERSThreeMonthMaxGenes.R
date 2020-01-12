#Development paper figures
#libraries and functions====
source("Code/Desktop/functions.r")
source("Code/Desktop/functions.DESeq2.r")
#data: input====
tmsg("Loading RNASeq data for graphing, generated using file 1.3")
load("Data/RNAseq/Output/LKB.ERCC.GeneSummaryObjects.rda")

MeanCounts.ERCC = Counts.ERCC.M[,.(meanCounts=mean(value)), .(Timepoint, gene)]
MaxTime = MeanCounts.ERCC[,.SD[meanCounts==max(meanCounts)], gene]
MaxTime = MaxTime[Timepoint=="3M", gene]

ThreeMonthMax.Spline = DEGs.ERCC.SplineDF3[geneID %in% MaxTime]

GenerateGeneSummaryPlots(ThreeMonthMax.Spline[1:20, geneID], isNameList = F)$Plot

GroupModel = readRDS("Data/RNAseq/Output/GroupModel.ERCC.ImputeB4Exclude.2excl.maxit5000.batchCutRINOXTimepoint.Wald.RDS")
OneMThreeM = GrabResult(GroupModel, c("Timepoint", "3M", "1W"), "Group|3M|1W", 0.05, TRUE)
OneMThreeM = OneMThreeM[geneID %in% MaxTime]

length(intersect(OneMThreeM[,geneID], ThreeMonthMax.Spline[,geneID]))   # 623
length(setdiff(ThreeMonthMax.Spline[,geneID], OneMThreeM[,geneID]))     #  98
length(setdiff(OneMThreeM[,geneID], ThreeMonthMax.Spline[,geneID]))     # 166

SplineModel = readRDS("Data/RNAseq/Output/SplineModel.ERCC.ImputeB4Exclude.2excl.maxit5000.batchCutRINOXnsImputedAgeDaysdf3.LRT.RDS")
# for (gene in unique(ThreeMonthMax.Spline[,geneID])) {
#     plot = DESeq2::plotCounts(SplineModel, gene, intgroup = "Timepoint", returnData = TRUE)
#     plot = ggplot(plot, aes(x=Timepoint, y=count, color=Timepoint))+geom_jitter(width = 0.05)+
#         stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.5)+xlab("")+ylab("Counts")+theme(legend.position = "none")
#     ggsave(filename = sprintf("Graphs/CERS.3MMax/%s.png", gene), units = "cm", width = 7, height = 7)
# }
