#Development paper figures
#libraries and functions====
library(extrafont); library(ggsignif)
source("Code/Desktop/functions.r")
source("Code/Desktop/functions.DESeq2.r")
source("Code/Desktop/functions.GO.r")
#Setup: Graphs theme====
TextSize = 10 #14 for thesis, 20+ for presentation
geom.text.size =  TextSize #/ ggplot2::.pt #geom_text uses "mm", ggplot theme uses point... yikes.
geom.text.size.small = geom.text.size / 2
theme_set(theme_cowplot(line_size = 1))
theme_update( text = element_text(family = windowsFont("Arial"), size = geom.text.size),  line = element_line(size=unit(0.7, "cm")), 
              legend.text.align = 1, plot.margin=margin(10,2,-15,0, "pt")
              )
theme_update(
    axis.ticks.length=unit(-0.2, "cm"), strip.background = element_rect(fill = NA, colour = "#000000", linetype = "solid"),
    axis.text.x = element_text(family = "Arial", size = geom.text.size, margin=margin(8,0,0,0,"pt")), 
    axis.text.y = element_text(family = "Arial", size = geom.text.size, margin=margin(0,10,0,0,"pt")), 
    strip.text = element_text(family = "Arial", size = geom.text.size, margin = margin(2,2,2,2, "pt"))
)
#Make sure to use size=geom.text.size/ggplot2::.pt for geom_label
#Setup: load RNA data====
load("Data/RNAseq/Output/LKB.ERCC.GeneSummaryObjects.rda")
GroupModel          = readRDS("Data/RNAseq/Output/GroupModel.ERCC.ImputeB4Exclude.2excl.maxit5000.batchCutRINOXTimepoint.Wald.RDS")
OneWeekAdult        = GrabResult(GroupModel, c("Timepoint", "AD", "1W"), "1W-AD", 0.05, TRUE)
SplineModel         = readRDS("Data/RNAseq/Output/SplineModel.ERCC.ImputeB4Exclude.2excl.maxit5000.batchCutRINOXnsImputedAgeDaysdf3.LRT.RDS")
SplineResult.Filter = GrabResult(SplineModel, "ns.ImputedAgeDays..df...3.3", "ns.Timepoint.df3", Alpha = 0.05, Filter = TRUE)

#Figure 5.1====
Fig5.1AData = fread("Data/CERS/CERS_ICaL_Development.txt")
Fig5.1AData[,Group:=factor(Group, levels=c("1W", "1M", "3M", "AD"))]
setorder(Fig5.1AData, "Group")
CERS.aov <- aov(Peak1 ~ Group, data = Fig5.1AData) ; summary(CERS.aov) #P < 0.001
TukeyHSD(CERS.aov)

Fig5.1AData.Avgs = Fig5.1AData[,.(meanPeak=abs(mean(Peak1)), sePeak=abs(se(Peak1)), n=.N), Group]
Fig5.1A = ggplot(Fig5.1AData.Avgs, aes(x=Group))+geom_point(aes(y=meanPeak), size=2)+
    geom_errorbar(aes(ymax=meanPeak+sePeak, ymin=meanPeak-sePeak, width=0.5), size=1)+
    scale_y_continuous(limits = c(0, 175), expand = c(0,0), breaks = seq(0, 175, length.out = 6))+
    ylab("Peak current (pA)")+xlab("")+
    geom_text(data = data.table(Group=c(3,4), meanPeak=c(150, 160), label=c("$","$")), size=geom.text.size/ggplot2::.pt, aes(y=meanPeak, label=label))

#ENSOARG00000017677 | CACNG3 is not in this data - not found -at all-. 
#According to expression atlas, no one else has detected it in heart so that's fine.  
Fig5.1B = GenerateGeneSummaryPlots(c("ENSOARG00000013089", "ENSOARG00000019143", "ENSOARG00000015836", "ENSOARG00000000883"), 
                                 isNameList = F, Dataset = "LKB.ERCC", PlotTPMs = F)$Plot+facet_wrap(~geneName, ncol = 2, scales = "free_y")+
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.25)))
Fig5.1B
#fig5.1C
PerUnit = unique(Dictionary[grepl("^CACN", `Gene name`), .(geneID=`Gene stable ID`, geneName=`Gene name`)])
PerUnit[,UnitType:=str_match(geneName, "CACN(A2D|A|B|G)\\d{1}")[,2]]

CACNA   = GenerateGeneSummaryPlots(PerUnit[UnitType=="A", geneID], UseCommonScale = T, PlotTPMs = F)$Plot  +facet_wrap(~geneName, nrow = 1)+ylab("")+
    theme(axis.text.x = element_blank(), plot.margin = margin(r=1, b=5, l=20))+scale_y_continuous(expand = expand_scale(mult = c(0, 0.2)))
CACNB   = GenerateGeneSummaryPlots(PerUnit[UnitType=="B", geneID], UseCommonScale = T, PlotTPMs = F)$Plot  +facet_wrap(~geneName, nrow = 1)+ylab("")+
    theme(axis.text.x = element_blank(), plot.margin = margin(r=1, b=5, l=20))
CACNG   = GenerateGeneSummaryPlots(PerUnit[UnitType=="G", geneID], UseCommonScale = T, PlotTPMs = F)$Plot  +facet_wrap(~geneName, nrow = 1)+ylab("")+
    theme(axis.text.x = element_blank(), plot.margin = margin(r=1, b=5, l=20))+scale_y_continuous(expand = expand_scale(mult = c(0, 0.2)))
CACNA2D = GenerateGeneSummaryPlots(PerUnit[UnitType=="A2D", geneID], UseCommonScale = T, PlotTPMs = F)$Plot+facet_wrap(~geneName, nrow = 1)+ylab("")+
    theme(plot.margin = margin(r=1, b=-5, l=20))

Fig5.1C = plot_grid(CACNA, CACNB, CACNG, CACNA2D, ncol = 1, rel_widths = 0.9)+
    cowplot::draw_label("Normalised Counts", 0.015, 0.5, angle = 90, size = geom.text.size) 
Fig5.1C # 550 * 450

#Figure 5.2: Phosphorylation - does basal phosphorylation change in development? ====
#Ratio of phosphorylation genes to LTCC
Dictionary = fread("Data/biomart/BIOMART-OvAr3.1-GeneIDNameDescrLocation-FINAL.txt")
Dictionary[,`Gene description`:=gsub("\\[Source:.*\\]$", "", `Gene description`)]
load("Data/RNASeq/Output/LKB.ERCC.GeneSummaryObjects.rda")

b = fread("Data/190402-LTCCChapter-PhosphoGraphGenes.txt")
# How we got this:
#   - grab all members from the families: AKAP, ADRB, PRKA, PRKG, PPP, PRKC, PRKG, CAMK2, PPP2A, PPP2B
#   - take only genes with more than 500 counts
#   - multiply meanCounts by abs(log2FC)
#   - order (descending) by that, take top 20

c = Counts.ERCC.M[gene %in% b[,`Gene stable ID`]]
d = Counts.ERCC.M[gene == "ENSOARG00000013089"] # Counts for CACNA1C

c = merge(c,d, by=c("sample", "Timepoint"))
set(c, j=which(grepl("Is_Expr", colnames(c))), value = NULL)
c[,Ratio:=value.x/value.y]

c = merge(c, b, by.x="gene.x", by.y="Gene stable ID", all.x=T)
c[`Gene name`=="", `Gene name`:=gene.x]
c[,Timepoint:=factor(Timepoint, levels=c("1W", "1M", "3M", "AD"))]
c[,Order:=sprintf("%d%s", Group, `Gene name`)]
c[,`Gene name`:=factor(`Gene name`, levels=c[order(Order), unique(`Gene name`)])]

TukeyHSD(aov(formula = Ratio~Timepoint, data=c[`Gene name`=="[CAMK2D]"]))

d = as.data.frame(c)
by(d, d$gene.x, function(x) { TukeyHSD(aov(formula = Ratio~Timepoint, data=x)) }, simplify = F)

e=unique(c[`Gene name` %in% c("PPP1CA", "CAMK2B", "[CAMK2D]"), .(Ratio=mean(Ratio), RatioSE = se(Ratio), `Gene name`), .(gene.x, Timepoint)])
# ggplot(mapping = aes(x=Timepoint, y=Ratio, color=Timepoint))+geom_point(data = e, size=5)+
#     geom_errorbar(data = e, aes(ymin=Ratio-RatioSE, ymax=Ratio+RatioSE), size=2)+facet_wrap(~`Gene name`, scales="free_y", ncol = 3)

BasePlot = ggplot(mapping = aes(x=Timepoint, y=Ratio, color=Timepoint))+theme(panel.spacing.y = unit(0, "pt"), 
            strip.text.x = element_text(margin = margin(1, 0, 1, 0, "pt")), axis.text.x.bottom = element_blank(), legend.position = "none", 
            plot.margin=margin(0,1,-10,10, "pt"))+ylab("")+xlab("")+facet_wrap(~`Gene name`, scales="free_y", ncol = 3)

Plot1Data.KMD = c[Group==1,.(Ratio=mean(Ratio), RatioSE = se(Ratio), `Gene name`), .(gene.x, Timepoint)]
#Plot1= BasePlot+geom_point(data = Plot1Data.KMD)+geom_errorbar(data = Plot1Data.KMD, aes(ymin=Ratio-RatioSE, ymax=Ratio+RatioSE), size=2)#+

Plot1 = BasePlot+geom_boxplot(data = c[Group==1], color="#e87777", outlier.shape = NA)+geom_jitter(data = c[Group==1], alpha=1, width=0.2)+
theme(plot.margin = margin(0,0,0,20))
Plot2 = BasePlot+geom_boxplot(data = c[Group==2], color="#5799db", outlier.shape = NA)+geom_jitter(data = c[Group==2], alpha=1, width=0.2)+
    theme(plot.margin = margin(0,0,0,20))
Plot3 = BasePlot+geom_boxplot(data = c[Group==3], color="#b88bed", outlier.shape = NA)+geom_jitter(data = c[Group==3], alpha=1, width=0.2)+
    theme(plot.margin = margin(0,0,0,20))
Plot4 = BasePlot+geom_boxplot(data = c[Group==4], color="#84c667", outlier.shape = NA)+geom_jitter(data = c[Group==4], alpha=1, width=0.2)+
    theme(plot.margin = margin(0,0,0,20))
Plot5 = BasePlot+geom_boxplot(data = c[Group==5], color="#A0A0A0", outlier.shape = NA)+geom_jitter(data = c[Group==5], alpha=1, width=0.2)+
    theme(axis.text.x.bottom = element_text(margin = margin(8,1,-5,1)), plot.margin = margin(0,0,0,20))
Fig5.2 = plot_grid(Plot1, Plot2, Plot3, Plot4, Plot5, ncol = 1, rel_heights = c(1, 0.6, 0.3, 0.3, 0.3))+
    cowplot::draw_label(x=0.01, y=0.5, angle = 90, "Ratio of gene counts to CACNA1C counts", size = geom.text.size)
Fig5.2
# d = dcast(c, `Gene name`~Timepoint, value.var="Ratio", fun.aggregate = mean)
# d[,RatioDiff:=`AD`-`1W`]

#figure 5.3 - alternative splicing====
source("Code/Desktop/5.5-AlternativeSplicingAnalysis.R")
LTCC_AS = AlternativeSplicingAnalysis("CACNA1C", "LEAFCUTTER", "LKB")
#LTCC_AS[[1]]$PublicationMap
Fig5.3A = LTCC_AS[[1]]$ASPlot
Fig5.3B = LTCC_AS[[1]]$IntronPlot

#Fig4 is colored text from an ENSEMBL protein alignment
#Fig5 is PyMOL
#Fig6 is PyMOL
#Table 3 - Transcription factors====
IPA.TF = fread("Data/IPA/1WAD-TranscriptionFactor-IPA-HHone2one-NewDEG.txt", sep = "\t", dec = ",")
IPA.TF[,V9:=NULL]
IPA.TF = IPA.TF[`Molecule Type`=="transcription regulator"]
#IPA.TF = merge(IPA.TF, DEGs.HumOrthol.All, by.x="Upstream Regulator", by.y="Human gene name", all.x=T)
#IPA.TF.Long = SplitDataTableWithMultiRows(IPA.TF, "Target molecules in dataset", ",")
Table3A=IPA.TF[`Target molecules in dataset` %like% "CACN" & !is.na(`Activation z-score`)] #TODO Check which subunit!
Table3A=SplitDataTableWithMultiRows(Table3A, "Target molecules in dataset")
Table3A=merge(Table3A, SplineResult.Filter[,.(`Gene name`, `log2 Fold Change of Target`=log2FC, `p-value (target)`=pValAdjust)], 
              by.x="Target molecules in dataset", by.y="Gene name", all.x=T)
setnames(Table3A, "Upstream Regulator", "Transcription Factor")
Table3A=Table3A[substr(`Target molecules in dataset`, 1, 4)=="CACN",
    .(`Transcription Factor`, `Activation z-score`, `p-value of overlap`, `Target molecules in dataset`, 
      `log2 Fold Change of Target`, `p-value (target)`)]
ffwrite(Table3A, "clipboard", dec=",")
#Export via ggsave as EMF====
if(F){
    #Warning, inserting these into a table will auto-scale them, but the proportions will be as dictated by width and height
    tmsg("exporting hires PNGs via cairo (consider EMF via ggsave?)...")
    ggsave("Graphs/CAPaper.Fig1A.ggsave.emf", Fig1A, device = "emf", units = "cm", width=6.5, height=5.4) #Fig1A
    ggsave("Graphs/CAPaper.Fig1B.ggsave.emf", Fig1B, device = "emf", units = "cm", width=8.5, height=5.5) #Fig1B
    ggsave("Graphs/CAPaper.Fig1C.ggsave.emf", Fig1C, device = "emf", units = "cm", width=16, height=11) #Fig1C
    ggsave("Graphs/CAPaper.Fig2.ggsave.emf", Fig2, device = "emf", units = "cm", width=12, height=22) #Fig4A
}
#cleanup====
#rm(list = (ls()[!grepl("^Fig|^Table", ls())]))

#gene check====
GenesOfInterest = Dictionary[substr(`Gene name`, 1, 3)=="PDE"]
GenesOfInterest = rbind(GenesOfInterest, Dictionary[substr(`Gene name`, 1, 4)=="AKAP"]) #AKAPs
GenesOfInterest = rbind(GenesOfInterest, Dictionary[substr(`Gene name`, 1, 4)=="ADRB"]) #Adrenoreceptor
GenesOfInterest = rbind(GenesOfInterest, Dictionary[substr(`Gene name`, 1, 4)=="PRKA"]) #Protein kinase A
GenesOfInterest = rbind(GenesOfInterest, Dictionary[substr(`Gene name`, 1, 4)=="PRKG"]) #Protein kinase G
GenesOfInterest = rbind(GenesOfInterest, Dictionary[substr(`Gene name`, 1, 3)=="PPP"])
GenesOfInterest = rbind(GenesOfInterest, Dictionary[substr(`Gene name`, 1, 4)=="PRKC"])
GenesOfInterest = rbind(GenesOfInterest, Dictionary[substr(`Gene name`, 1, 4)=="PRKG"])
GenesOfInterest = rbind(GenesOfInterest, Dictionary[substr(`Gene name`, 1, 5)=="CAMK2"]) #CAMK2
GenesOfInterest = rbind(GenesOfInterest, Dictionary[substr(`Gene name`, 1, 5)=="PPP2A"])
GenesOfInterest = rbind(GenesOfInterest, Dictionary[substr(`Gene name`, 1, 5)=="PPP2B"])
GenesOfInterest = rbind(GenesOfInterest, Dictionary[substr(`Gene name`, 1, 4)=="CACN"])
GenesOfInterest[,`Transcript stable ID`:=NULL]
GenesOfInterest = unique(GenesOfInterest[,.(`Gene stable ID`, `Gene name`, `Gene description` )])

PrimaryDEGs.Spline = SplineResult.Filter[geneID %in% GenesOfInterest[,`Gene stable ID`]]
tmsg(sprintf("Of %d LTCC-related genes, %d are differentially expressed in the Spline model.", nrow(GenesOfInterest), nrow(PrimaryDEGs.Spline)))
#ffwrite(PhosphoDEGs.Spline, "clipboard")

PrimaryDEGs.1WAD = OneWeekAdult[geneID %in% GenesOfInterest[,`Gene stable ID`]]
tmsg(sprintf("Of %d LTCC-related genes, %d are differentially expressed in the Group 1W/AD model.", nrow(GenesOfInterest), nrow(PrimaryDEGs.1WAD)))

PrimaryDEGs.1WADvSpline = intersect(PrimaryDEGs.1WAD[,geneID], PrimaryDEGs.Spline[,geneID])
tmsg(sprintf("Spline and 1WAD models agree on %d of %d genes.", length(PrimaryDEGs.Spline), length(unique(c(PrimaryDEGs.1WAD[,geneID], PrimaryDEGs.Spline[,geneID])))))

PhosphoDEGs.Both = merge(PrimaryDEGs.1WAD, PrimaryDEGs.Spline, 
                         by=c("geneID", "Chromosome/scaffold name", "Gene start (bp)", "Gene end (bp)", "Gene name", "Gene description"), all=TRUE)
set(PhosphoDEGs.Both, j=which(colnames(PhosphoDEGs.Both) %in% c("Chromosome/scaffold name", "Gene start (bp)", "Gene end (bp)")), value=NULL)
Fig4C.Alt1 = GenerateGeneSummaryPlots(GenesOfInterest[,`Gene stable ID`], isNameList = TRUE, Model="Both")
Fig4C.AltData = Fig4C.Alt1$PlotSignif$data
Fig4C.Alt1.NS = Fig4C.Alt1$Plot
Fig4C.Alt1 = Fig4C.Alt1$PlotSignif
#ffwrite(Dictionary[`Gene stable ID` %in% Fig4C.AltData[,unique(geneID)]], "clipboard") #For characterisation in Excel and Google Scholar
