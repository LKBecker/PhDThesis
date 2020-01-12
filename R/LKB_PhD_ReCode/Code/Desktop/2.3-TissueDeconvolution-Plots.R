library(Biobase); library(bseqsc); library(data.table); library(ggplot2)
library(Matrix); library(plyr); library(reshape2); library(xbioc); require(MuSiC)
library(ape); library(ggtree); library(RColorBrewer)


Skelly.DeLaughter.MusicBasis = readRDS("Data/MuSiC/Skelly+DeLaughter.MusicBasis.rds")
Skelly.DeLaughter.MusicProps = readRDS("Data/MuSiC/AtrialTissue-CellProportionsMkII.rds")

#Which genes are most abundant in myocardial cells?
Myocardial = as.data.table(Skelly.DeLaughter.MusicBasis$M.theta, keep.rownames = T)
Myocardial.M=melt(Myocardial, id.vars = "rn")
Myocardial=Myocardial.M[,.SD[value==max(value), variable],rn]
Myocardial[V1=="myocardial", .N] #909
#Appears to need more normalisation

Skelly.DeLaughter.Distances1 <- dist( t( log(Skelly.DeLaughter.MusicBasis$Disgn.mtx + 1e-6) ), method = "euclidean") #Extract design matrix
Skelly.DeLaughter.Clust1 <- hclust(Skelly.DeLaughter.Distances1, method = "complete" ) # Hierarchical clustering using Complete Linkage
plot(Skelly.DeLaughter.Clust1, cex = 0.6, hang = -1, main = 'Cluster log(Design Matrix)') # Plot the obtained dendrogram
rm(Skelly.DeLaughter.Distances1)

Skelly.DeLaughter.Phylo <- as.phylo(Skelly.DeLaughter.Clust1)
ggtree(Skelly.DeLaughter.Phylo,size=1) + xlim(c(NA,550)) + geom_tiplab(size=4,offset = 5, lineheight = 0.8) + geom_treescale()
#ggsave("tree-basic.pdf",plot=ggtobj,width=4,height=9)

Skelly.DeLaughter.Distances2 <- dist( t( log( Skelly.DeLaughter.MusicBasis$M.theta + 1e-8 ) ), method = "euclidean") 
#Extract average relative abundance
Skelly.DeLaughter.Clust2 <- hclust(Skelly.DeLaughter.Distances2, method = "complete") # Hierarchical clustering using Complete Linkage
plot(Skelly.DeLaughter.Clust2, cex = 0.6, hang = -1, main = 'Cluster log(Mean of RA)') # Plot the obtained dendrogram

MuSiC.CellProps = as.data.table(Skelly.DeLaughter.MusicProps$Est.prop.weighted, keep.rownames = T)
setnames(MuSiC.CellProps, "rn", "sample")
MuSiC.CellProps[,Timepoint:=str_extract(sample, "..$")]
MuSiC.CellProps.M = melt(MuSiC.CellProps, id.vars = c("sample", "Timepoint"))
MuSiC.CellProps.M[,Timepoint:=factor(Timepoint, levels=c("1W","1M","3M","AD"))]
se = function(x) sd(x)/sqrt(length(x))
CellPropStats = MuSiC.CellProps.M[,.(MeanProp=mean(value), StdErr=se(value)), .(variable, Timepoint)]
CellPropStats[,Timepoint:=factor(Timepoint, levels=c("1W","1M","3M","AD"))]
colourCount = MuSiC.CellProps.M[,length(unique(variable))]
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

ggplot(MuSiC.CellProps.M, aes(x=Timepoint, y=value, fill=factor(variable)))+geom_boxplot()+geom_jitter(alpha=0.5, width = 0.2)+
    xlab("")+ylab("Inferred % of tissue")+scale_fill_manual(values = getPalette(colourCount))+facet_wrap(~variable)


stop("The below doesn't yet work; can't get clusters.")

CellProps.GeneWeights <- data.table(melt(Skelly.DeLaughter.MusicProps[["Weight.gene"]],varnames = c("gene_name","sample"),value.name = "prop")) 
CellProps.GeneWeights = CellProps.GeneWeights[!is.na(prop)]
#"Weight.gene contains the 'marker' assignments

#mp <- merge(mp,refined_clusters[,.(cluster_id,cluster_group,cluster_colour)],by="cluster_id")
#mp <- mp[,.(prop=sum(prop)),.(cluster_group, sample_id, cluster_colour)]
#mp[,sample_type:=ifelse(sample_id %like% "POZT|N","Biopsy","Nephrectomy")]
#mp[,study:=ifelse(sample_id %like% "-","TCGA","MTG")]
CellProps.GeneWeights[,old_prop:=prop]
CellProps.GeneWeights[,prop:=prop/sum(prop),sample]
CellProps.GeneWeights[,sum(prop),sample]
mat <- acast(CellProps.GeneWeights,cluster_group ~ sample, value.var = "prop")
hc <- hclust(dist(mat))
plot(hc)
ordered_samples <- hc$labels[hc$order]
