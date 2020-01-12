source("Code/config.r"); source("Code/libraries.r"); source("Code/functions.r"); 
source("Code/Graphs.Setup.R") #We use DESeq2 data here

require(GenABEL)
require(aroma.light)

CONF$N_GENES = 50 #Number of genes per pattern to show 

#TODO: needs function for proper color scaling w/ squash and limits
MakeHeatmap <- function(DT){
  HEATMAP = ggplot(mapping = aes_string(x="Sample", y="geneName", fill="value"))+theme(axis.text.x = element_text(angle=90, vjust = 0.3))+
    scale_fill_gradientn(name="Counts scaled\nrelative to\ngene mean", colors=c("#0000FF", "#FFFFFF", "#FF0000"), limits=c(-5,5))+ #TODO replae -5, 5 with DT-based scaling factor
    scale_x_discrete()+ylab("")+ geom_raster(data = DT, mapping = aes(x=Sample, y=geneName, fill=value))
  HEATMAP
}


#Plotting setup ====
GeneExpression = dcast.data.table(ExpressPct, geneID~Timepoint, value.var = "IsExpressed")
GeneExpression[,Pattern:=paste0(`1W`*1, `1M`*1, `3M`*1, AD*1)]
GeneExpression = GeneExpression[,.(geneID, Pattern)]

#New Method: Quantile normalisation per sample, making them comparable
#ScaledCounts = aroma.light::normalizeQuantileRank(log(t(counts(CountData))+1), robust = TRUE) 
ScaledCounts = log(counts(CountData)+1)
ScaledCounts = ScaledCounts[rowSums(ScaledCounts[,-1])>0,]

NQRankCounts = matrix(nrow = NROW(ScaledCounts), ncol = 0)
for (Timepoint in c("1W", "1M", "3M", "AD")) {
  columns = which(substr(colnames(ScaledCounts), 11, 12) == Timepoint)
  NQRankCounts = cbind(aroma.light::normalizeQuantileRank(ScaledCounts[,columns], robust=TRUE), NQRankCounts)
}

a = apply(NQRankCounts, 1, function(x)length(unique(x))>2)
NQRankCounts.Shadow = NQRankCounts[!a,]
NQRankCounts = NQRankCounts[a,] #TODO: Keep in different matrix, rbind() back in??

ScaledCounts <- apply(NQRankCounts, 1, function(x) {
    rt <- rntransform(x)
    names(rt) <- names(x)
    rt
})

ScaledCounts <- t(ScaledCounts)
ScaledCounts.M = melt(data.table(ScaledCounts, geneID=rownames(ScaledCounts)), id.vars = "geneID")
ScaledCounts.M = merge(ScaledCounts.M, unique(Dictionary[,.(`Gene stable ID`, `Gene name`)]), 
                       by.x="geneID", by.y="Gene stable ID")
setnames(ScaledCounts.M, c("variable", "Gene name"), c("Sample", "geneName"))
Metadata = as.data.table(colData(CountData))
Metadata = Metadata[,.(sample, ImputedAgeDays, Group)]
Metadata[,Sample:=paste0(sample, "_", Group)]
ScaledCounts.M = merge(ScaledCounts.M, Metadata, by="Sample")

ScaledCounts.M[geneName=="", geneName:=geneID] #Ensures there will always be a label available under geneName

ScaledCounts.M=merge(ScaledCounts.M, GeneExpression, by="geneID", all.x=T) #Transfers information on expression pattern
ScaledCounts.M[,Timepoint:=factor(str_match(Sample, "1W|1M|3M|AD"), levels = c("1W", "1M","3M","AD"))]
#ScaledCounts.M[,Sample:=factor(Sample, levels=unique(ScaledCounts.M[order(Timepoint), Sample]))]
ScaledCounts.M[,Sample:=factor(Sample, levels=unique(ScaledCounts.M[order(ImputedAgeDays), Sample]))]

rm(NQRankCounts, ScaledCounts)

if(F){
    PatternHeatmapPlots = list() #Creates container
    for (PlotPattern in GeneExpression[,unique(Pattern)]) {
        plot = #HEATMAP + geom_raster(data = ScaledCounts.M[geneID %in% GeneExpression[Pattern==PlotPattern, geneID][1:CONF$N_GENES]],
               #                     mapping = aes(x=Sample, y=geneName, fill=value))
          MakeHeatmap(ScaledCounts.M[geneID %in% GeneExpression[Pattern==PlotPattern, geneID][1:CONF$N_GENES]]) #TODO select biggest LFC / pVal?
        PatternHeatmapPlots[[PlotPattern]] = plot
    }
    PatternHeatmapPlots$`0001`
}
