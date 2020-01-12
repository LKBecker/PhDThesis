#RyR chapter figures
#"Does Sr-based calcium handling change in postnatal development?"
#libraries and functions====
library(extrafont)
source("Code/Desktop/functions.DESeq2.r")
source("Code/Desktop/functions.GO.r")
source("Code/Desktop/functions.r")
#Setup: Graphs theme====
TextSize = 10 #10 for thesis, 20+ for presentation
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
#Setup: load RNA data====
load("Data/RNAseq/Output/LKB.ERCC.GeneSummaryObjects.rda")
SplineModel= readRDS("Data/RNAseq/Output/SplineModel.ERCC.ImputeB4Exclude.2excl.maxit5000.batchCutRINOXnsImputedAgeDaysdf3.LRT.RDS")

#Gene lists====
# SR_FUNCTION_GENES = c("RYR2", "JPH2", "ATP2A1", "ATP2A2", "ATP2A3", "CASQ1", "CASQ2", 
#                       "PLN", "SRL", "HRC", "CALU", "CALR", "HSPA5", "CANX", "ASPH", "RCN1", "SDF4", "RCN2", "RCN3")
# SR_FUNCTION_GENES_ID = Dictionary[`Gene name` %in% SR_FUNCTION_GENES, .(`Gene name`, `Gene stable ID`)][order(`Gene name`)]
# SR_FUNCTION_GENES_ID = rbind(SR_FUNCTION_GENES_ID, data.table("Gene name"=c("[RCN3]", "[TRDN]"), 
#                                         "Gene stable ID"=c("ENSOARG00000013207", "ENSOARG00000008072")))

SR_TETHER_GENES = c("RYR2", "JPH2", "PLN")
SR_TETHER_GENES = Dictionary[`Gene name` %in% SR_TETHER_GENES, .(`Gene name`, `Gene stable ID`)][order(`Gene name`)]
SR_TETHER_GENES = rbind(SR_TETHER_GENES, data.table("Gene name"="TRDN", "Gene stable ID"="ENSOARG00000008072"))
SR_TETHER_GENES_ID = SR_TETHER_GENES[,`Gene stable ID`]
names(SR_TETHER_GENES_ID)=SR_TETHER_GENES[,`Gene name`]

SERCAs = c("ATP2A1", "ATP2A2", "ATP2A3")
SERCAs = Dictionary[`Gene name` %in% SERCAs, .(`Gene name`, `Gene stable ID`)][order(`Gene name`)]
SERCA_ID = SERCAs[,`Gene stable ID`]
names(SERCA_ID)=SERCAs[,`Gene name`]

#Figure 6.1 is Dr. Charlotte Smith's Figure from p 204 of their thesis, PROPERLY REREFENCED AS A REPRODUCTION====
#Figure 6.2 - RyR and associated SR DEGs====
Fig6.2 = GenerateGeneSummaryPlots("RYR2", Dataset = "ERCC+GWPM", Model = "1WAD+CvHF")$Plot+
    scale_y_continuous(expand = expand_scale(mult = c(0, .1)))+facet_wrap(~geneName, ncol = 2, scales="free_y") 
Fig6.2
#RYR AS analysis - 1 hit, small
#source("Code/Desktop/5.5-AlternativeSplicingAnalysis.R")
#RYR2.AS = AlternativeSplicingAnalysis("RYR2", MODE = "LEAFCUTTER", DATASET = "LKB", DRAW_MODE = 2, USE_FILTER = TRUE)
#RYR2.AS[[1]]$ASPlot
#Nothing exciting here
#JPH2.AS = AlternativeSplicingAnalysis("JPH2", MODE = "LEAFCUTTER", DATASET = "LKB", DRAW_MODE = "REGION_ONLY", USE_FILTER = TRUE)
#Nothing.

#Figure 6.3 - Control and grouping of the RyR ====
Fig6.3 = GenerateGeneSummaryPlots(c("JPH2"="ENSOARG00000003672", "TRDN"="ENSOARG00000008072", "[Junctin]"="ENSOARG00000016594", 
                                   "CASQ1"="ENSOARG00000008405", "CASQ2"="ENSOARG00000020163"), Dataset = "ERCC+GWPM", Model = "1WAD+CvHF")$Plot+
    scale_y_continuous(expand = expand_scale(mult = c(0, .15)))+facet_wrap(~geneName, ncol = 3, scales="free_y")
Fig6.3

#Figure 6.4 - Phosphorylation of the RyR ====
Fig6.4 = GenerateGeneSummaryPlots(c("PRKAA2"="ENSOARG00000007851", "PRKAB1"="ENSOARG00000001411", "PRKAG3"="ENSOARG00000019718", 
                                  "PRKAR2A"="ENSOARG00000016132", "CAMK2B"="ENSOARG00000015725", "CAMK2D"="ENSOARG00000018416"),
                                  Dataset = "ERCC+GWPM", Model = "1WAD+CvHF")$Plot+
    scale_y_continuous(expand = expand_scale(mult = c(0, .2)))+facet_wrap(~geneName, ncol = 3, scales="free_y")
Fig6.4

Fig6.5 = GenerateGeneSummaryPlots(Dictionary[grepl("^PPP\\dC", `Gene name`), `Gene stable ID`],  Dataset="ERCC+GWPM", Model="Spline+CvHF", PlotTPMs = T)$Plot+
    scale_y_continuous(expand = expand_scale(mult = c(0, .2)))
Fig6.5

#Figure 6.6 - Ca buffering ====
#Cytosolic:
#Troponin C, SERCA, Myosin = ~90% of cytosolic buffering in the ventricle 40. 
#sarcalumenin 

#CYTOSOLIC_BUFFERS
Myosins = GenerateGeneSummaryPlots(Dictionary[grepl("^myosin", `Gene description`) & grepl("^MY", `Gene name`), `Gene stable ID`])$PlotData
Myosins = Myosins[,mean(TPM), .(geneID, group)]
#Fig5A = GenerateGeneSummaryPlots(Myosins[,sum(V1), geneID][V1>1000, geneID], Dataset="LKB.ERCC", Model="Spline")$PlotSignif 
#Fig5Av2 = GenerateGeneSummaryPlots(Myosins[,sum(V1), geneID][V1>1000, geneID], Dataset="GWPM", Model="CvHF")$PlotSignif 
#"Myosin genes with a total average count of at least 200, significant changes only"

Troponins = GenerateGeneSummaryPlots(Dictionary[grepl("^troponin", `Gene description`), `Gene stable ID`])$PlotData
Troponins = Troponins[,mean(TPM), .(geneID, group)]
#Fig5B = GenerateGeneSummaryPlots(Troponins[,sum(V1), geneID][V1>400, geneID])$Plot+scale_y_continuous(expand = expand_scale(mult = c(0, .09)))
#Fig5Bv2 = GenerateGeneSummaryPlots(Troponins[,sum(V1), geneID][V1>400, geneID], Dataset="GWPM", Model="CvHF")$Plot+
#    scale_y_continuous(expand = expand_scale(mult = c(0, .09)))

#"Myosin genes with a total average count of at least 200, significant changes only"

SERCA = GenerateGeneSummaryPlots(Dictionary[grepl("^ATP\\dA\\d", `Gene name`), `Gene stable ID`])$PlotData
SERCA = SERCA[,mean(TPM), .(geneID, group)]
#Fig5C = GenerateGeneSummaryPlots(SERCA[,sum(V1), geneID][V1>400, geneID])$Plot+
#    scale_y_continuous(expand = expand_scale(mult = c(0, .09)))
#Fig5Cv2 = GenerateGeneSummaryPlots(SERCA[,sum(V1), geneID][V1>400, geneID], Dataset="GWPM", Model="CvHF")$Plot+
#    scale_y_continuous(expand = expand_scale(mult = c(0, .09)))
#"Myosin genes with a total average count of at least 200, significant changes only"

#Sarcalumenin = GenerateGeneSummaryPlots("ENSOARG00000003074")$Plot
#Sarcalumeninv2 = GenerateGeneSummaryPlots("ENSOARG00000003074", Dataset="GWPM", Model="CvHF")$Plot+scale_y_continuous(expand = expand_scale(mult = c(0, .1)))

#SR Lumen, List from Dr. Katharine Dibb:
CA_BUFFER_GENES = c("SRL", "HRC", "CALU", "CALR", "HSPA5", "CANX", "ASPH", "RCN1", "SDF4", "RCN2", "RCN3", "CASQ1", "CASQ2")
CA_BUFFER_GENES = Dictionary[`Gene name` %in% CA_BUFFER_GENES, .(`Gene name`, `Gene stable ID`)][order(`Gene name`)]
CA_BUFFER_GENES = rbind(CA_BUFFER_GENES, data.table("Gene name"=c("RCN3"), "Gene stable ID"=c("ENSOARG00000013207")))
CA_BUFFER_GENES_ID = CA_BUFFER_GENES[,`Gene stable ID`]
names(CA_BUFFER_GENES_ID)=CA_BUFFER_GENES[,`Gene name`]

SR_BUFFER_DATA = GenerateGeneSummaryPlots(CA_BUFFER_GENES_ID, isNameList = F, UseCommonScale = T, Dataset="ERCC+GWPM", Model="Spline+CvHF", PlotTPMs = T)

#SR_BUFFER_DATA = SR_BUFFER_DATA[,.(meanTPM=mean(count), seTPM=se(count)), .(geneName, group)]
SR_BUFFER_STATS = SR_BUFFER_DATA$Changing
SR_BUFFER_STATS[,Symbol:=""]
for (SIGNIF_LEVEL in c(0.05, 0.01, 0.001)) { SR_BUFFER_STATS[pValAdjust <= SIGNIF_LEVEL, Symbol:=paste0(Symbol, "*")] }

#SR_BUFFER_DATA_PLOT = SR_BUFFER_DATA$PlotData[group %in% c("1W", "AD")][,.(meanTPM=mean(TPM), seTPM=se(TPM)), .(geneName, group)]
SR_BUFFER_DATA_PLOT = SR_BUFFER_DATA$PlotData[group %in% c("1W", "AD", "C", "HF", "R")][,.(meanTPM=mean(TPM), seTPM=se(TPM)), .(geneName, group)]
SR_BUFFER_DATA_PLOT[,group:=factor(group, levels=c("1W", "AD", "C", "HF", "R"))]

#Fig6.6 - Ca2+ buffers over the ages====
#Calculation
#SR_CALC = SR_BUFFER_DATA$PlotData[group %in% c("1W", "AD"),sum(TPM), .(group, sample)]
SR_CALC = SR_BUFFER_DATA$PlotData[group %in% c("1W", "AD", "C", "HF", "R"), sum(TPM), .(group, sample)]

t.test(SR_CALC[group=="1W", V1], SR_CALC[group=="AD", V1])
t.test(SR_CALC[group=="AD", V1], SR_CALC[group=="C", V1]) #NS
t.test(SR_CALC[group=="C", V1], SR_CALC[group=="HF", V1]) #0.0003

SR_ANOVA = aov(V1~group, SR_CALC)
summary(SR_ANOVA) #***
TukeyHSD(SR_ANOVA)

Fig6.6 = ggplot(SR_BUFFER_DATA_PLOT, mapping = aes(x=group, y=meanTPM, fill=geneName))+geom_col(color="#000000", position = "stack")+
    geom_signif(comparisons = list(c("1W", "AD"), c("AD", "C"), c("C", "HF"), c("HF", "R")), annotations = c("***",  "0.99", "***", "**"), 
                y_position = c(3200, 2750, 4000, 4150))+theme(legend.position = "none")
Fig6.6+coord_flip()+theme(plot.margin = margin(5,0,5,0))+xlab("")+geom_vline(xintercept = 2.5)

#Fig6.7 ====
SignifData = merge(SR_BUFFER_DATA$Changing[,.(geneName, Comparison, Symbol)], 
                   SR_BUFFER_DATA_PLOT[group %in% c("C", "HF"),.(y_pos_HF = 20+round(max(meanTPM*1.1))), geneName], by="geneName", all.y=T)
SignifData = merge(SignifData, SR_BUFFER_DATA_PLOT[group %in% c("1W", "AD"),.(y_pos_Dev = 20+round(max(meanTPM)*1.1)), geneName], by="geneName", all.y=T)
SignifData = SignifData[!grepl("^Group\\|R", Comparison)]
SignifData[,geneName:=factor(geneName, levels=c("[RCN3]", "ASPH", "CALR", "CALU", "CANX", "CASQ1", "CASQ2", "HRC", "HSPA5", "RCN2", "SDF4", "SRL"))]
setorder(SignifData, geneName)

SignifData = merge(SignifData, data.table(geneName=SignifData[,unique(geneName)], xpos=seq(1, length(SignifData[,unique(geneName)]))), by="geneName")

SignifData[Comparison=="All", from := "1W"]
SignifData[Comparison=="All", to := "AD"]
SignifData[Comparison!="All" & !is.na(Comparison),c("Delete", "to", "from") := tstrsplit(Comparison , "|", fixed=TRUE)]
SignifData[,Delete:=NULL]
SignifData[is.na(Comparison), c("from", "to") := ""]

XPOS = c("1W"=-0.35, "AD"=-0.12, "C"=0.12, "HF"=0.35)

SignifData[,xmin := xpos + XPOS[from]]
SignifData[,xmax := xpos + XPOS[to]]
SignifData[, y_pos:= ifelse(from=="C", y_pos_HF, y_pos_Dev)]
SignifData = SignifData[,.(geneName, xmin, xmax, y_pos, Symbol)]
SignifData = SignifData[!is.na(Symbol)]

#the ggsignif github suggests adding a dummy group
SignifData[,dummyGroup:= .I]

Fig6.7 = ggplot(SR_BUFFER_DATA_PLOT, aes(x=geneName, y=meanTPM, group=geneName, fill=group))+geom_col(position = position_dodge2(width = 0.2), width=0.9)+
    geom_errorbar(aes(ymin=meanTPM-seTPM, ymax=meanTPM+seTPM), position = position_dodge2(width = 0.2), width=0.9)+
    theme(axis.ticks.length=unit(0.25, "cm"), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2))+
    scale_fill_manual(values = c("#f8766d", "#c67cff", "#00B6EB", "#A58AFF", "#f0f0f0"), name="Time point")+
    scale_y_continuous(expand = expand_scale(mult = c(0, .075)))+
    #geom_vline(xintercept = seq(1, 12, 1))+
    geom_signif(data=SignifData, mapping = aes(annotations = Symbol, xmin = xmin, xmax = xmax, y_position = y_pos, fill=NA, group=dummyGroup), 
                color="#000000", manual = TRUE, tip_length = 0.01)

#stat_summary(geom="crossbar", fun.y = mean, fun.ymax = mean, fun.ymin = mean)
Fig6.7 # 950*300 px

#ARCHIVE: Figure 6.7 as %====
# SR_BUFFER_DATA_PLOT[,sumAge:=sum(meanTPM), group]
# SR_BUFFER_DATA_PLOT[,percent:= meanTPM / sumAge]
# 
# Fig6 = ggplot(SR_BUFFER_DATA_PLOT, mapping = aes(x=group, y=percent, fill=geneName))+geom_col(color="#000000", position = "stack")+
#     scale_y_continuous(labels = scales::percent, breaks = c(0, 0.5, 1))+theme(legend.position = "none", 
#     plot.margin = margin(0,5,8,5))+ylab(expression(paste("Percentage of these Ca"^"2+", "buffering genes' TPM")))+xlab("")
# Fig6

#figure 6.10 - NCX, SERCA, PLB - Ca Destinations. Max 550*200 px a ====
Fig6.10A = GenerateGeneSummaryPlots(c("PLN"="ENSOARG00000019971", "SLN"="ENSOARG00000016955"), Dataset="ERCC+GWPM", Model="1WAD+CvHF", 
                                    PlotTPMs = F, UseCommonScale = T)$Plot+scale_y_continuous(expand = expand_scale(mult=c(0, 0.15)))

Fig6.10B = GenerateGeneSummaryPlots(c("PLN"="ENSOARG00000019971", "ATP2A1"="ENSOARG00000002132", "ATP2A2"="ENSOARG00000016226", "ATP2A3"="ENSOARG00000017976"),
                                 Dataset="ERCC+GWPM", Model="1WAD+CvHF")$Plot+facet_wrap(~geneName, ncol = 3, scales = "free_y")+
    scale_y_continuous(expand = expand_scale(mult=c(0, 0.2)))

Fig6.10C.NCX = GenerateGeneSummaryPlots(c("SLC8A1"="ENSOARG00000007758", "SLC8A2"="ENSOARG00000011168", "SLC8A3"="ENSOARG00000021191"),
                                        Dataset="ERCC+GWPM", Model="1WAD+CvHF")$Plot+facet_wrap(~geneName, ncol = 3, scales = "free_y")+
    scale_y_continuous(expand = expand_scale(mult=c(0, 0.2)))
Fig6.10C.NCX

#Fig6.11 Calculations====
#Ratio of SLN/PLN to SERCAaaaaaa
RyRModData.Counts = GenerateGeneSummaryPlots(c("SERCA2"="ENSOARG00000016226", "PLN"="ENSOARG00000019971", "SLN"="ENSOARG00000016955"), 
                                             Dataset="ERCC+GWPM", Model="1WAD+CvHF", PlotTPMs = FALSE)$PlotData
#We have to plot counts as SLN, PLN don't have enough insert size to calculate TPMs and thus are eliminated from the count table there.
RyRModRatios = dcast.data.table(RyRModData.Counts, sample~geneName, value.var = "count")
RyRModRatios[,`SERCA : PLN`:= ATP2A2/PLN]
RyRModRatios[,`SERCA : PLN+SLN`:= ATP2A2/(PLN+`[SLN]`)]
RyRModRatios[,`SERCA : SLN`:= ATP2A2/`[SLN]`]
RyRModRatios[,group:= str_extract(sample, "1W$|1M$|3M$|AD$|HF$|C$|R$")]
RyRModRatios[,group:= factor(group, levels = c("1W","1M","3M","AD","C","HF","R"))]

RyRModRatios.M = melt.data.table(RyRModRatios[,.(sample, group, `SERCA : PLN`, `SERCA : PLN+SLN`, `SERCA : SLN`)], id.vars = c("sample", "group"))

SERCAPLN.Anova.LKB = aov(lm(`SERCA : PLN`~group, RyRModRatios[group %in% c("1W","1M","3M","AD")]))
summary(SERCAPLN.Anova.LKB) #NS

SERCAPLN.Anova.BHB = aov(lm(`SERCA : PLN`~group, RyRModRatios[group %in% c("C","HF","R")]))
summary(SERCAPLN.Anova.BHB) #***
#TukeyHSD(SERCAPLN.Anova.BHB)

SERCASLN.Anova.LKB = aov(lm(`SERCA : SLN`~group, RyRModRatios[group %in% c("1W","1M","3M","AD")]))
summary(SERCASLN.Anova.LKB) #***
TukeyHSD(SERCASLN.Anova.LKB)

SERCASLN.Anova.BHB = aov(lm(`SERCA : SLN`~group, RyRModRatios[group %in% c("C","HF","R")]))
summary(SERCASLN.Anova.BHB) #**
#TukeyHSD(SERCASLN.Anova.BHB)

SERCAPLNSLN.Anova.LKB = aov(lm(`SERCA : PLN+SLN`~group, RyRModRatios[group %in% c("1W","1M","3M","AD")]))
summary(SERCAPLNSLN.Anova.LKB) #***
#TukeyHSD(SERCAPLNSLN.Anova.LKB)

SERCAPLNSLN.Anova.BHB = aov(lm(`SERCA : PLN+SLN`~group, RyRModRatios[group %in% c("C","HF","R")]))
summary(SERCAPLNSLN.Anova.BHB) #NS

RyRModRatios.M2 = RyRModRatios.M[, .(meanRatio=mean(value), seRatio=se(value)), .(variable, group)]
RyRModRatios.M2[,group:= factor(group, levels = c("1W","1M","3M","AD","C","HF","R"))]

Fig6.11.Bottom = ggplot(RyRModRatios.M2, aes(x=group, y=meanRatio))+geom_point(aes(color=variable, group=variable), size=2)+
    geom_signif(xmin = c(1),   xmax = c(4),   y_position = 12,   annotations = c("***"),       tip_length = 0.005, color="#00BA38")+
    geom_signif(xmin = c(5,6), xmax = c(6,7), y_position = 15, annotations = c("***", "**"), tip_length = 0.005, color="#F8766D")+
    geom_signif(xmin = c(5),   xmax = c(6),   y_position = 23,   annotations = c("**"),        tip_length = 0.005, color="#00bfc4")+
    geom_errorbar(aes(ymax=meanRatio+seRatio, ymin=meanRatio-seRatio), color="black", width=0.5)+geom_vline(xintercept = 4.5, size=1)+
    scale_y_continuous(limits = c(3, 24), breaks = seq(3,24,7), expand = c(0,1.3))+
    theme(legend.title = element_blank(), legend.position = "none", plot.margin = margin(5, 15, 0, 5))+labs(x="", y="Ratio of counts")

Fig6.11.Top = ggplot(RyRModRatios.M2, aes(x=group, y=meanRatio))+geom_point(aes(color=variable, group=variable), size=2)+
    geom_signif(xmin = c(1), xmax = c(4), y_position = 77, annotations = c("***"), tip_length = 0.03, color="#00bfc4")+
    geom_errorbar(aes(ymax=meanRatio+seRatio, ymin=meanRatio-seRatio), color="black", width=0.5)+geom_vline(xintercept = 4.5, size=1)+
    scale_y_continuous(limits = c(25, 85), breaks = seq(25,85,30), expand = c(0,20))+
    theme(legend.title = element_blank(), legend.position = "none", axis.text.x = element_blank(), axis.ticks.length.x = unit(0, "cm"),
          plot.margin = margin(5, 15, 0, 5))+labs(x="", y="")

Fig6.11.Legend = cowplot::get_legend(ggplot(RyRModRatios.M2, aes(x=group, y=meanRatio))+geom_point(aes(color=variable, group=variable), size=2)+
                                       theme(legend.position = "right", legend.spacing.x = unit(-.5, 'cm')))

Fig6.11 = plot_grid(plot_grid(Fig6.11.Top, Fig6.11.Bottom, ncol=1, rel_heights = c(.5, 1)), plot_grid(NULL, Fig6.11.Legend), ncol = 2, align = "h", axis = "l", 
          rel_widths = c(0.9, 0.2)) # 500 * 375

#Fig6.12 is ~hand drawn~

#misc====
# Kat.SERCAPlot = GenerateGeneSummaryPlots(SERCA_ID, isNameList = F, UseCommonScale = T, Model = "1WAD")$Plot
# #ATP2A (SERCA2) dominates the field, NS change but on average down
# Kat.SERCA.Bottom = Kat.SERCAPlot+coord_cartesian(ylim = c(0, 500))+theme(strip.background = element_blank(), strip.text = element_blank(),
#                                                                          plot.margin = margin(20,5,5,15))
# Kat.SERCA.Top = Kat.SERCAPlot+coord_cartesian(ylim = c(100000, 250000))+theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
#                                                                               axis.line.x = element_blank())
# Fig3 = plot_grid(Kat.SERCA.Top, Kat.SERCA.Bottom, ncol = 1, rel_heights = c(0.5, 0.5)) 
# 

#TRDN.AS = AlternativeSplicingAnalysis("ENSOARG00000008072", "LEAFCUTTER", IS_GENE_NAME = F)
#TRDN.AS[[1]]$ASPlot
#TRDN.AS[[1]]$PublicationMap
#v low counts

#FK506 binding protein - interacts w/ RYR 
#"ENSOARG00000018587" #significant decrease but not in TPM table due to fragment length ~ gene length
#GenerateGeneSummaryPlots("ENSOARG00000017806", PlotTPMs = F)$Plot+facet_wrap(~geneName, ncol = 2, scales = "free_y")+
#    scale_y_continuous(expand = expand_scale(mult=c(0, 0.2)))

#KMD buffer proportion====
#"Can you look statistically at how the complex changes ie triadin, CSQ and junction?"
#"What's the proportion in NB and how does that change in adult?"
BufferData = GenerateGeneSummaryPlots(c("TRDN"="ENSOARG00000008072", "[Junctin]"="ENSOARG00000016594", 
                           "CASQ1"="ENSOARG00000008405", "CASQ2"="ENSOARG00000020163"), Model = "1WAD")$PlotData[group %in% c("1W", "AD")]
BufferData.Avg = BufferData[, .(meanTPM = mean(TPM), sdTPM=sd(TPM)), .(geneName, group)]

#Fisher's Test
BufferData.Counts = GenerateGeneSummaryPlots(c("TRDN"="ENSOARG00000008072", "[Junctin]"="ENSOARG00000016594", 
                                        "CASQ1"="ENSOARG00000008405", "CASQ2"="ENSOARG00000020163"), 
                                      Model = "Spline", PlotTPMs = FALSE)$PlotData[group %in% c("1W", "AD")]
BufferData.Counts = BufferData.Counts[, .(meanCounts = mean(count)), .(geneName, group)]
FisherTable.Buffers = as.matrix(dcast.data.table(BufferData.Counts, geneName~group, value.var = "meanCounts")[,2:3])
rownames(FisherTable.Buffers)=dcast.data.table(BufferData.Counts, geneName~group, value.var = "meanCounts")[,geneName]
fisher.test(FisherTable.Buffers, workspace = 8e+6)

#Plot
#ggplot(BufferData.Avg, aes(x=group, group=geneName, y=meanTPM, fill=geneName))+geom_col(width=0.5, color="#000000")+
#    geom_signif(comparisons = list(c("1W","AD")), annotations = c("***"), map_signif_level=TRUE, y_position = 800)+
#    scale_y_continuous(limits = c(0,850), breaks = seq(0,850,length.out = 3), expand=c(0,0))
#ggplot(BufferData.Avg, aes(x=group, group=geneName, y=meanTPM, color=geneName))+geom_point()+geom_errorbar(aes(ymax=meanTPM+sdTPM, ymin=meanTPM-sdTPM), width = 0.75)

#Phosphorylation verification====
#LTCCChapterPhosphoGenes = fread("Data/190402-LTCCChapter-PhosphoGraphGenes.txt")
# How we got this:
#   - grab all members from the families: AKAP, ADRB, PRKA, PRKG, PPP, PRKC, PRKG, CAMK2, PPP2A, PPP2B
#   - take only genes with more than 500 counts
#   - multiply meanCounts by abs(log2FC)
#   - order (descending) by that, take top 20
#Phospho = GenerateGeneSummaryPlots(LTCCChapterPhosphoGenes[,`Gene stable ID`])
#Phospho$Plot
#ggplot(Phospho$PlotData, aes(x=group, y=TPM))+geom_boxplot(outlier.shape = NA)+geom_jitter(aes(color=group), width = 0.1)+facet_wrap(~geneName, scales="free_y")
#GenerateGeneSummaryPlots(Dictionary[grepl("^PPP1C", `Gene name`), `Gene stable ID`])$Plot
#GenerateGeneSummaryPlots(Dictionary[grepl("^PPP1R", `Gene name`), `Gene stable ID`])$PlotSignif
