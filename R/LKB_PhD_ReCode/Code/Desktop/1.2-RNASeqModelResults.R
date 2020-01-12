#Setup functions, libraries and data====
source("Code/Desktop/functions.r")          #Export, pretty-rprint functions
source("Code/Desktop/functions.DESeq2.r")   #Results grabbers
source("Code/Desktop/functions.GO.r")       #PANTHER processing
DictHumHomol = fread("Data/biomart/181218_biomart_SheepHumanIdNameOrthologyType.txt")

#Get Spline.LRT.ERCC results====
SplineModel =readRDS("Data/RNAseq/Output/SplineModel.ERCC.ImputeB4Exclude.2excl.maxit5000.batchCutRINOXnsImputedAgeDaysdf3.LRT.RDS")
Spline.LRT.ERCC = GrabResult(SplineModel, "ns.ImputedAgeDays..df...3.3", "LRT", 0.05, TRUE)
ffwrite(Spline.LRT.ERCC, "DESeq2.Spline.LRT.ERCC", "Data")

#Process Spline  results for GO and output ID list====
DEGs.HumOrthol = merge(Spline.LRT.ERCC, DictHumHomol, by.x="geneID", by.y="Gene stable ID", all.x=T)
message(paste0("Found Orthologs for ", DEGs.HumOrthol[`Human gene stable ID`!="", .N], " of ", 
               Spline.LRT.ERCC[,length(unique(geneID))], " genes."))
DEGs.HumOrthol[`Human gene stable ID`!="",length(unique(`Human gene stable ID`))] #4742
ffwrite(DEGs.HumOrthol[`Human gene name`!="",unique(HumanID)], "DESeq2.Spline.LRT.ERCC.HumanOrthologs.All", "Data/RNAseq/Output/")
DEGs.HumOrthol.O2O = DEGs.HumOrthol[`Human homology type`=="ortholog_one2one"]
#ffwrite(DEGs.HumOrthol.O2O, "DESeq2.Spline.LRT.ERCC.HumanOrthologs.One2One.FullTable", "Data/RNAseq/Output/")
DEGs.HumOrthol.O2O[`Human gene stable ID`!="",length(unique(`Human gene stable ID`))] #4206
ffwrite(DEGs.HumOrthol[`Human gene name`!="",unique(HumanID)], "DESeq2.Spline.LRT.ERCC.HumanOrthologs.One2One", "Data/RNAseq/Output/")

#Get Group.Wald.ERCC results====
GroupModel = readRDS("Data/RNAseq/Output/GroupModel.ERCC.ImputeB4Exclude.2excl.maxit5000.batchCutRINOXTimepoint.Wald.RDS")
Group.Wald.ERCC = GetAllResultsForVariable(GroupModel, "Timepoint", 0.05, TRUE)
DEGs.Group.HumOrthol = merge(Group.Wald.ERCC[,.(geneID, Comparison)], DictHumHomol, by.x="geneID", by.y="Gene stable ID", all.x=T)
#Process Group model to PANTHER====
#For each comparison, how many hits have a human ortholog?
DEGs.Group.HumOrthol[, .SD[`Human gene name`!="",.N] / .N, Comparison] #92-95%
#For each comparison, how many hits have a one2one human ortholog?
DEGs.Group.HumOrthol[, .SD[`Human gene name`!="" & `Human homology type`=="ortholog_one2one",.N] / .N, Comparison] #80-84%

ffwrite(DEGs.Group.HumOrthol[`Human gene name`!="",unique(`Human gene stable ID`)], 
        "DESeq2.ERCC.Group.Wald.AllComparisons.AllHumanOrthologs")
ffwrite(DEGs.Group.HumOrthol[`Human gene name`!="" & `Human homology type`=="ortholog_one2one",unique(`Human gene stable ID`)], 
        "DESeq2.ERCC.Group.Wald.AllComparisons.HumanOrthologsOne2One")
rm(GroupModel)

# #ERCC checkup====
# GroupCounts.ERCC    = DESeq2::counts(SplineModel, normalized=TRUE) #grab normalised counts, regardless of model 
# ERCCData = as.data.table(GroupCounts.ERCC[rownames(GroupCounts.ERCC) %like% "ERCC",])
# ERCCData[,gene:=rownames(GroupCounts.ERCC)[rownames(GroupCounts.ERCC) %like% "ERCC"]]
# rm(GroupCounts.ERCC)
# ERCCData = melt(ERCCData, id.vars="gene")
# setnames(ERCCData, "variable", "sample")
# ERCCData[,TimePoint:=factor(str_extract(sample, "..$"), levels = c("1W", "1M", "3M", "AD"))]
# ggplot(ERCCData, aes(x=TimePoint, y=value, group=gene))+geom_point()+geom_smooth(method = "glm")+facet_wrap(~gene, scales="free_y")
# 
# #Is there a pattern of 7 samples being (aritifically?) too far 'up'?
# #Normalise to ERCC spikein or use DESeq2 standard procedure?
# 
# #Simon Anders argues that normalisation is supposed to iron out differences in sequnecing depth / amount of starting RNA, 
# #whereas spike-ins would _already_ be constant (plus minus differences in pipetting)
# 
# #Seven points above the general 'level' - result of normalisation? e.g. the spikeins were constant indeed 
# #but the samples into which they were spiked were not;
# #when DESeq2 fixed the samples it had to consequently un-even the spikeins e.g. it's a consequnece of WHEN they were added 
# 
# ggplot(ERCCData.M[gene %in% c("ERCC-00113", "ERCC-00084", "ERCC-00025", "ERCC-00081", "ERCC-00002")], 
#        aes(x=TimePoint, y=value, group=gene))+geom_label(aes(label=sample))+geom_smooth(method = "glm")+facet_wrap(~gene, scales="free_y")
