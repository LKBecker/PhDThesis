#T-tubule paper figures
#"How do T-tubules change in development?"
#libraries and functions====
library(extrafont)
library(stringr)
source("Code/Desktop/functions.r")
source("Code/Desktop/functions.DESeq2.r")
source("Code/Desktop/functions.GO.r")
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
GroupModel = readRDS("Data/RNAseq/Output/GroupModel.ERCC.ImputeB4Exclude.2excl.maxit5000.batchCutRINOXTimepoint.Wald.RDS")
SplineModel= readRDS("Data/RNAseq/Output/SplineModel.ERCC.ImputeB4Exclude.2excl.maxit5000.batchCutRINOXnsImputedAgeDaysdf3.LRT.RDS")

Spline.ERCC.df3 = GrabResult(SplineModel, "ns.ImputedAgeDays..df...3.3", "Spline.df3", 0.05, TRUE) 
Group.ERCC.All = GetAllResultsForVariable(GroupModel, "Timepoint", 0.05, TRUE) #DO NOT SET IT TO 0.999 JFC

#Figure Fig4.1 - TT rise but BIN1, CAV3, DYSF do not====
CharlotteTTData <- fread("Data/CERS/CERS_HalfDistanceDevelopment.txt")
CharlotteTTData[,Age:=factor(Age, levels = c("1W", "1M", "3M", "AD"))]
CharlotteTTData[,S1:=as.Date.character(S1, "%y%m%d")]
CharlotteTTData[,DOB:=as.Date.character(DOB, "%d/%m/%Y")]
CharlotteTTData2 = unique(CharlotteTTData[,.(meanFA=mean(`F.A.`), sdFA = sd(`F.A.`), Age), Animal])
TukeyHSD(aov(formula = meanFA~Age, data=CharlotteTTData2))

#CERS data shows that t-tubules increase up to 3M
Fig4.1A = ggplot(CharlotteTTData2, aes(x=Age, y=meanFA, color=Age))+geom_boxplot(outlier.shape = NA, color="black", width=0.5)+
    geom_jitter(width = 0.05)+xlab("")+ylab("Fractional Area of t-tubules")+theme(legend.position = "none")+
    scale_y_continuous(limits = c(0,0.18), breaks = seq(0,0.18,0.09), expand = c(0,0))+
    geom_signif(comparisons=list(c("1W", "AD"), c("1W", "3M")), map_signif_level = TRUE, y_position = c(0.14, 0.16), color="black")

#This suggess genes involved in the support, maintenance or assembly of t-tubules should likewise increase.
#We know several t-tubules involved in formation, e.g. BIN1, JPH2, CAV3
Fig4.1B = GenerateGeneSummaryPlots(c("BIN1"="ENSOARG00000016188", "DYSF"="ENSOARG00000011500", "CAV3"="ENSOARG00000007515"), Dataset="ERCC+GWPM", Model = "1WAD+CvHF")$Plot
Fig4.1B+scale_y_continuous(expand = expand_scale(mult = c(0, .1))) 

# Fig4.1BData = GenerateGeneSummaryPlots(c("DYSF"="ENSOARG00000011500"), Model = "1WAD")$PlotData
# DYSF.TT.Correlation = merge(Fig4.1BData, CharlotteTTData2[,mean(meanFA), Age], by.x="group", by.y="Age", all.x=T)
# cor.test(DYSF.TT.Correlation[,V1], DYSF.TT.Correlation[,count])

#Fig4.1C - western blot of AMPII/BIN1
CERSBin1WesternData = fread("Data/CERS.WesternBlotData")
CERSBin1WesternData = melt(CERSBin1WesternData, id.vars="Sample")
CERSBin1WesternData[,variable:=factor(variable, levels=c("1 week", "1 month", "3 months", "Adult"), labels=c("1W", "1M", "3M", "AD"))]

##Anova -> Tukey
#CERSBin1Stats = rbindlist(lapply(CERSBin1WesternData[,unique(Sample)], function(x) { y = TukeyHSD(aov(formula = value~variable, data=CERSBin1WesternData[Sample==x])); data.table(y$variable, Sample=x, Comparison=rownames(y$variable)) }))
#CERSBin1Stats = CERSBin1Stats[Comparison=="AD-1W"]
##Kruskal-Wallis
#CERSBin1Stats = rbindlist(lapply(CERSBin1WesternData[,unique(Sample)], function(x) as.data.table(kruskal.test(formula = value~variable, data=CERSBin1WesternData[Sample==x])$p.value)[,Sample:=x]))
#One-way ANOVA
CERSBin1Stats = rbindlist(lapply(CERSBin1WesternData[,unique(Sample)], function(x) { y = aov(formula = value~variable, data=CERSBin1WesternData[Sample==x]); data.table("p.value"=summary(y)[[1]][["Pr(>F)"]][[1]], Sample=x) }))
setkey(CERSBin1Stats, "Sample")
#setnames(CERSBin1Stats, "V1", "p.value") #kruskal
#setnames(CERSBin1Stats, "p adj", "p.value") #Tukey
CERSBin1Stats = CERSBin1Stats[CERSBin1WesternData[,.(y_position_fin=round(max(value)*1.1, 1)), Sample]]
#CERSBin1Stats[,xmin:=substr(Comparison, 4, 5)] #Tukey
#CERSBin1Stats[,xmax:=substr(Comparison, 1, 2)] #Tukey
CERSBin1Stats[,xmin:="1W"] #ANOVA
CERSBin1Stats[,xmax:="AD"] #ANOVA
CERSBin1Stats[p.value<0.05, Symbol:="*"]
CERSBin1Stats[p.value<0.01, Symbol:="**"]
CERSBin1Stats[p.value<0.001, Symbol:="***"]
CERSBin1Stats[p.value>0.05, Symbol:="NS"]

Fig4.1C = ggplot(CERSBin1WesternData, aes(x=variable, y=value, color=Sample))+geom_point(position = position_dodge(width=0.66), show.legend = FALSE)+
                    stat_summary(geom = "crossbar", fun.ymax = mean, fun.ymin = mean, fun.y = mean, width=0.66, position = position_dodge(), mapping = aes(group=Sample), color="black")+
    scale_color_brewer(palette=6, type = "qual")+facet_wrap(~Sample, scales = "free_y")+
    geom_signif(data = CERSBin1Stats, aes(annotations = Symbol, xmin = xmin, xmax = xmax, y_position = y_position_fin), manual = TRUE, color="black")+
    scale_y_continuous(expand = expand_scale(mult = c(0, .15))) 

#TukeyHSD(aov(formula = value~variable, data=CERSBin1WesternData[Sample=="AMPII Total"]))

#Cavins make caveolae
# GenerateGeneSummaryPlots(Dictionary[grepl("^CAVIN", `Gene name`), `Gene stable ID`], Dataset="ERCC+GWPM", Model = "CvHF")$Plot #CAVIN2 decreases at 1M, but... No idea if that's meaningful?
# ggplot(GenerateGeneSummaryPlots(Dictionary[grepl("^CAV\\d+", `Gene name`), `Gene stable ID`], Dataset="ERCC+GWPM", Model = "CvHF")$PlotData, aes(x=group, y=count, color=group))+
#     geom_boxplot(color="black")+geom_jitter(width=0.05)+facet_wrap(~geneName)
#None of these changes are significant. The counts suggest that BIN1 goes down, CAV3 shows the *opposite* pattern to TT density, JPH2 is all over the place. 
#Either, these genes are not required to increase in expression to increase TT density (the addition is cumulative or a separate mechanism exists)
#Or other genes are involved up to 3M

#Fig4 is alt splicing plot of DYSF
source("Code/Desktop/5.5-AlternativeSplicingAnalysis.R")
Fig4.2 = AlternativeSplicingAnalysis("DYSF", MODE = "DESEQ2")
Fig4.2[[1]]$IntronPlot

Fig4.3 = GenerateGeneSummaryPlots(c("MG29"="ENSOARG00000019213", "Obscurin"="ENSOARG00000005438", "TPM3"="ENSOARG00000001845", "MTM"="ENSOARG00000008921",
                                              "TCAP"="ENSOARG00000011575", "TRDN"="ENSOARG00000008072", "ENSOARG00000013187"), 
                                            Dataset="ERCC+GWPM", Model = "1WAD+CvHF")$PlotSignif+scale_y_continuous(expand = expand_scale(mult = c(0, .2))) 
#"PARVG"="ENSOARG00000019286", 
Fig4.3

#ARCHIVE: Figure2====
#T-tubules peak in density at 3Months (the following decrease is relative to cell size)
MeanCounts = Counts.ERCC.M[,.(meanCount=mean(value)), .(gene, Timepoint)]
MaxCountTimepoints = MeanCounts[,.SD[meanCount==max(meanCount), Timepoint], gene]
#These are the genes that peak in mean counts at 3M:
ThreeMonthMax = MaxCountTimepoints[V1=="3M", unique(gene)]
#The 1WAD model might not capture their significance, so a 1W3M model was generated
DEGs.1W3M = GrabResult(GroupModel, "Timepoint_3M_vs_1W", "1W v 3M", 0.05, T)
sprintf("Of %d genes with max mean counts at 3M, %d are DEGs in 1W3M and %d are DEGs in the Spline. %d of these DEGs are in both models.", 
        length(ThreeMonthMax), DEGs.1W3M[geneID %in% ThreeMonthMax, .N], DEGs.ERCC.SplineDF3[geneID %in% ThreeMonthMax, .N],
        length(intersect(DEGs.1W3M[geneID %in% ThreeMonthMax, geneID], DEGs.ERCC.SplineDF3[geneID %in% ThreeMonthMax, geneID]))
)
#The list of genes peaking at 3M was cross-references with both the 1W3M and spline models, creating a list of DEGs which peak at 3M.
DEGs.3MMax.SplineOr1W3M = rbind(DEGs.1W3M[geneID %in% ThreeMonthMax], DEGs.ERCC.SplineDF3[geneID %in% ThreeMonthMax])
GenerateGeneSummaryPlots(DEGs.3MMax.SplineOr1W3M[baseMean > 500 & log2FC > 0.5,unique(geneID)], Model = "Both")$Plot #166 genes... filter.

#Some highlights: 
#EDNRA follows pattern but appears to be expressed inside t-tubule....
#LAMB2 - laminin, colocalizes with t-tubule (Klietsch 1993), though v dispersed; like EDNRA might be 'side effect' of TT development
#TUBA8 - alpha tubulin, microtubule forming, peaks when TT does - role of MT in tt formation?
#SRL - drops after TT system reaches max density... less buffering needed?
#OBSCN - cytoskeletal calmodulin and titin-interacting RhoGEF 
    # Rhos - Tho/Rac/Cdc42 -> 
        # actin remodelling
        # cell cycle
        # cell junctions
#FNBP1L (TOCA1) - decreases in synch w/ TT!!
#TOCA has F-BAR domain (PIP2 membranes), HR1 (cdc42), SH3 (dynamin, n-wasp). BIN1 does not have the HR1...?
#https://www.ncbi.nlm.nih.gov/pubmed/23872330
#NWASP also binds cdc42
#The above-mentioned BAR domains, such as in the amphiphysin, endophilin and SNX9 proteins, have SH3 domains that bind to dynamin and N-WASP, 
#providing the connection between the membrane shape and the actin polymerization machinery (17–19). 
#Amphiphysin induces N-WASP-dependent actin polymerization that is dependent on liposomes, as demonstrated using a cytosolic fraction from 
#amphiphysin knockout mice (20).
#https://academic.oup.com/jb/article/148/1/1/886794#15340274

#OBSCN could be multifunctional, but could provide a link between progressive f-actin remodelling and t-tubule density increase. 
#Endosomes require cdc42 and ARP2/3 complex...
#"There is agreement that the force that is required for membrane extension is derived from such F-actin assembly, as discussed in the following sections." (Sit, 2011, J Cell Sci)

#ARCHIVE: Figure4 Fatty acid production?====
# PhosphatidylGenes = Dictionary[`Gene description` %like% "phosphatidyl"]
# Fig4.5 = GenerateGeneSummaryPlots(PhosphatidylGenes[,`Gene stable ID`], Dataset="ERCC+GWPM", Model = "CvHF")
# Fig4.5$PlotSignif
# 
# Fig4.5Data = rbind(Spline.ERCC.df3[geneID %in% PhosphatidylGenes[,`Gene stable ID`]], 
#                  Group.ERCC.All[geneID %in% PhosphatidylGenes[,`Gene stable ID`]])
# PIP2 = Fig4.5Data[`Gene description` %like% "phosphatidylinositol-4,5-"]
# GenerateGeneSummaryPlots(PIP2[,geneID], Dataset="ERCC+GWPM", Model = "CvHF")$Plot #overall loss of PIP2 kinases; more PIP2 or less PI3K signalling or......
# #Increase in HF -> Less PIP2, then reversion.

#Processing GO annotations====
DictOrtholog = fread("Data/biomart/181205_SheepMouseHuman_geneIDs_geneNames_OrthologTypes.txt")
setkey(DictOrtholog, "Human gene stable ID")
HumanUniProt = fread("Data/biomart/190526-HumanUniProtToGeneID.txt")
setkey(HumanUniProt, `UniProtKB Gene Name ID`)

ProcessGOEntries = function(Genes){
    cat(sprintf("Processing %d UniProt IDs...\n", length(Genes)))
    
    Prots = HumanUniProt[Genes]
    cat(sprintf("%d UniProt IDs correspond to a stable human Gene ID\n", Prots[!is.na(`Gene stable ID`),.N]))
    
    Prots = Prots[!is.na(`Gene stable ID`), unique(`Gene stable ID`)] #49 unique gene stable IDs
    cat(sprintf("%d unique human Gene IDs.\n", nrow(Prots)))

    Prots.OA = DictOrtholog[Prots][!is.na(`Gene stable ID`), unique(`Gene stable ID`)] #50 due to some one2many mappings and a few losses
    cat(sprintf("%d gene IDs after mapping to sheep.\n", nrow(Prots.OA)))

    GenerateGeneSummaryPlots(Prots.OA, Dataset="ERCC+GWPM", Model = "CvHF")$PlotSignif #14 genes significant
}

GO.Tubule.Human = fread("Data/AmiGO/190526-TTubule.txt")[substr(V1, 1,7)=="UniProt",unique(gsub("^UniProtKB:", "", V1))] #253 IDs
ProcessGOEntries(GO.Tubule.Human) #meh

GO.TTorg = fread("Data/AmiGO/190526-TTubuleOrganisation.txt")[substr(V1, 1,7)=="UniProt",unique(gsub("^UniProtKB:", "", V1))] #27 IDs
ProcessGOEntries(GO.TTorg) #none Significant

GO.Invagination.Human = fread("Data/AmiGO/190526-MembraneInvagination.txt")[substr(V1, 1,7)=="UniProt",unique(gsub("^UniProtKB:", "", V1))] #253 IDs
ProcessGOEntries(GO.Invagination.Human) #TODO <- interesting!!

#Membrane organisation====
#UniProt / PROSITE searches
#https://www.uniprot.org
#database:(type:prosite ps51021) #bar domain
#database:(type:prosite PS51741) #F-bar domain
#taxonomy:"Mammalia [40674]" #mammal
#taxonomy:"Homo sapiens (Human) [9606]" #human
#database:(type:prosite ps50002) #SH3 domain
#database:(type:prosite ps50238) #RHOGAP
#database:(type:prosite PS50010) #DBL homology  
#Biochemical data have established the role of the conserved DH domain in Rho GTPase interaction and activation, and the role of the tandem PH domain in intracellular targeting 
#and/or regulation of DH domain function. 
#BAR, MAMMAL, SH3, RHOGAP

#Figure 6
#https://prosite.expasy.org/cgi-bin/prosite/mydomains/
#BIN1
#29  ,276 ,   1,1, BAR
#520 ,593 ,   3,2, SH3

#TOCA1
#  1, 263, 1, 1, FBAR
# 397, 474, 2, 3, REM1
# 538, 599, 3, 2, SH3
# 
# 245, 535, 1
# 541, 597, 1

BIN1SH3MammalProts = fread("Data/UniProtKB/uniprot-ps51021+BARdomain-taxonomy+Mammalia-SH3domain.tab")
BIN1SH3MammalProts = BIN1SH3MammalProts[`Protein names` != "Uncharacterized protein"]
BIN1SH3MammalProts = BIN1SH3MammalProts[!(`Protein names` %like% "(Fragment)")]
BIN1SH3MammalProts = BIN1SH3MammalProts[!(`Protein names` %like% "isoform")]
BIN1SH3MammalProts = BIN1SH3MammalProts[!(`Protein names` %like% "cDNA")]

BIN1SH3MammalProts.HumanPN = BIN1SH3MammalProts[Organism=="Homo sapiens (Human)"]
#Comparing two tables, the human PN appears to cover p much all there is. 
#This gives us: 
#BIN1/AMPH (BIN1), 
#Endophilin-A(1-3)|B(1-3) (SH3GL[1/2] CNSA2 SH3D2A SH3P4 SH3GLB[1/2]), 
#DNMBP (ARHGEF36 / TUBA), 
#Rho guanine nucleotide exchange factor 37/38 (ARHGEF37|38)
#SH3-containing Grb-2-like 1 protein
#TOCA1

BARSH3Genes = Dictionary[grepl("^(BIN1|AMPH)", `Gene name`)]
BARSH3Genes = rbind(BARSH3Genes, Dictionary[grepl("^SH3(GL|D2|P4)|CNSA", `Gene name`)])
BARSH3Genes = rbind(BARSH3Genes, Dictionary[grepl("^ARHGEF3[6-8]", `Gene name`)])
BARSH3Genes = rbind(BARSH3Genes, Dictionary[grepl("^DNMBP", `Gene name`)])
Fig4.5 = GenerateGeneSummaryPlots(BARSH3Genes[,`Gene stable ID`], Dataset="ERCC+GWPM", Model="1WAD+CvHF")$Plot+scale_y_continuous(expand = expand_scale(mult = c(0, .3)))
Fig4.5

BARMammalProts = fread("Data/UniProtKB/uniprot-prosite+ps51021-taxonomy+Mammalia+40674-BAR+Domain.tab")
BARMammalProts = BARMammalProts[Organism=="Homo sapiens (Human)"] #filter...
BARMammalProts = BARMammalProts[`Ensembl transcript`!=""]
BARMammalProts = SplitDataTableWithMultiRows(BARMammalProts, "Ensembl transcript", ";")
BARMammalProts[,`Ensembl transcript`:=gsub(" \\[.{6}-\\d{1,2}\\]", "", `Ensembl transcript`)]
setkey(BARMammalProts, `Ensembl transcript`)
DictHuman = fread("Data/biomart/190526-HumanUniProtToGeneID.txt")
DictOrtholog[BARMammalProts]

#Fig4.6 Made in JalView
#Fig4.7, Fig4.8 made in PyMOL

#Fig 4.9 F-BAR + SH3 genes====
FBAR.SH3 = c(
    "ENSOARG00000001033", "ENSOARG00000001509", "ENSOARG00000001650", "ENSOARG00000003915", "ENSOARG00000004206", 
    "ENSOARG00000005162", "ENSOARG00000005293", "SRGAP2"="ENSOARG00000005343", "ENSOARG00000006122", "ENSOARG00000006995", 
    "ENSOARG00000007281", "ENSOARG00000008209", "ENSOARG00000010682", "ENSOARG00000016857", "ENSOARG00000019085"
)

Fig4.9 = GenerateGeneSummaryPlots(FBAR.SH3, isNameList = FALSE, Dataset = "ERCC+GWPM", Model = "1WAD+CvHF")$PlotSignif+
    scale_y_continuous(expand = expand_scale(mult = c(0, .2)))+facet_wrap(~geneName, ncol = 3, scales="free_y")
Fig4.9 #550 x 700

#Fig 4.10 and 4.11 Repair complex====
Fig4.10A = GenerateGeneSummaryPlots(Dictionary[grepl("^EHD", `Gene name`), `Gene stable ID`], Dataset = "ERCC+GWPM", Model = "1WAD+CvHF")$Plot+
    facet_wrap(~geneName, ncol = 2, scales="free_y")+scale_y_continuous(expand = expand_scale(mult = c(0, .15)))

    scale_y_continuous(expand = expand_scale(mult = c(0, .15)))

Fig4.11 = GenerateGeneSummaryPlots(Dictionary[grepl("^annexin", `Gene description`), `Gene stable ID`], Dataset = "ERCC+GWPM", Model = "1WAD+CvHF")$Plot+
    scale_y_continuous(expand = expand_scale(mult = c(0, .15)))+facet_wrap(~geneName, ncol = 4, scales="free_y")

#Fig4.12A Dynactin====
Fig4.12A = GenerateGeneSummaryPlots(Dictionary[grepl("^DYN", `Gene name`), `Gene stable ID`], Dataset="ERCC+GWPM", 
                                     Model = "1WAD+CvHF")$PlotSignif+scale_y_continuous(expand = expand_scale(mult = c(0, .15)))

#Fig4.12B CLIP170 and related EB1====
Fig4.12B = GenerateGeneSummaryPlots(c("ENSOARG00000011159", "EB1"="ENSOARG00000015417"),Dataset="ERCC+GWPM", Model = "1WAD+CvHF")$Plot+scale_y_continuous(expand = expand_scale(mult = c(0, .1)))

#FBARSH3Plot2 = GenerateGeneSummaryPlots(FBAR.SH3[8:14])$PlotSignif#+scale_y_continuous(expand = expand_scale(mult = c(0, .1))) 
#ggsave("Graphs/TTPaper.Fig10-1.ggsave.emf", FBARSH3Plot, device = "emf", units = "cm", width=15, height=20, pointsize=10, res=1024) #Fig1A
# ggsave("Graphs/TTPaper.Fig10-2.ggsave.emf", FBARSH3Plot2, device = "emf", units = "cm", width=15, height=20, pointsize=10, res=1024) #Fig1A

#Dynamin====
#is BIN1 a dynamin 'trap'? It accumulates curvature and PIP2, both things dynamin likes, but also has a unique Dynamin-scission-killing SH3 domain
GenerateGeneSummaryPlots(Dictionary[grepl("dynamin", `Gene description`), `Gene stable ID`], Dataset="ERCC+GWPM", Model = "CvHF")$Plot #global fukken decrease
#Are 1W levels sufficient to run the full production cycle, and a gradual decrease is acceptable as long as tubule levels 'compound'?

#more BAR genes
GenerateGeneSummaryPlots("ENSOARG00000018352", Dataset="ERCC+GWPM", Model = "CvHF")$Plot #Epsin2, membrane bender and PIP2 lover; BAR domain haver
#https://www.ncbi.nlm.nih.gov/pubmed/16938488?dopt=Abstract&holding=npg

#Molecular rulers
GenerateGeneSummaryPlots(c("Obscurin"="ENSOARG00000005438", "Titin"="ENSOARG00000017195", "Nebulin"="ENSOARG00000008659", "ENSOARG00000020239"), 
                         forceUserNames = T, Dataset="ERCC+GWPM", Model = "CvHF")$Plot
#Obscuring peaks at 3M, is associated with 'auxiliary' cytoskeleton (extra-sarcomeric cytoskeleton) {https://jcs.biologists.org/content/joces/122/15/2640.full.pdf?with-ds=yes}4
GenerateGeneSummaryPlots(Dictionary[grepl("^SYNE", `Gene name`), `Gene stable ID`], Dataset="ERCC+GWPM", Model = "CvHF")$Plot #nothing

#Microtubule genes
GenerateGeneSummaryPlots(c("IQGAP1"="ENSOARG00000013053", "cdc42"="ENSOARG00000008088", "RAC1"="ENSOARG00000002021", "clip170"="ENSOARG00000011159"), Dataset="ERCC+GWPM", Model = "CvHF")$Plot 
#MT complex generally decreases in pn development - more activation required / less reserves?


#Arp2/3 complex (branch vs straight MT)
GenerateGeneSummaryPlots(c("ENSOARG00000019990", "ENSOARG00000011780", "ENSOARG00000018781", "ENSOARG00000018731", 
                           "ENSOARG00000019495", "ENSOARG00000016761", "ENSOARG00000002017"), Dataset = "GWPM", Model = "CvHF")$Plot

#MG29 AND JPH2 = apoptosis
GenerateGeneSummaryPlots(c("ENSOARG00000019213", "ENSOARG00000003672"), Dataset = "ERCC+GWPM", Model = "CvHF")$Plot 

#dystrophin - CvHF ++, 1WAD NSDa    
GenerateGeneSummaryPlots(c("ENSOARG00000018256"), Dataset = "ERCC+GWPM", Model = "1WAD")$Plot 

GenerateGeneSummaryPlots(c("ENSOARG00000017794", "ENSOARG00000006511", "ENSOARG00000008088"), Dataset = "ERCC+GWPM", Model = "1WAD+CvHF")$Plot

# dystrophin family dynamin, membrane scissiojn protein: decrease in2, all others NS
GenerateGeneSummaryPlots(c("EHD1"="ENSOARG00000010447"), Dataset = "ERCC+GWPM", Model = "1WAD")$Plot

GenerateGeneSummaryPlots(c("PLCB1"="ENSOARG00000009195"), Dataset = "LKB.ERCC", Model = "1WAD")$Plot 

#we are aparently running a GSEA now====
TrainingData = GrabResult(GroupModel, c("Timepoint", "AD", "1W"), comparisonStr = "1W v AD", Alpha = 0.05, Filter = FALSE)
TrainingData = TrainingData[,.(`Gene name`,pValAdjust,t_stat,log2FC)]
TrainingData[,abs_t:=abs(t_stat)]                          #generate column of absolute t_stat 
setorder(TrainingData, -abs_t)                             #order by absolute t_stat column, in decreasing order
TrainingData <- unique(TrainingData,by="Gene name")   #each gene name / t_stat combination may only exist once
TrainingData <- TrainingData[!is.na(t_stat)]          #exclude values without a t_stat

TrainingSet         <- TrainingData[,t_stat]               #new vector of only t_stat values, ordered by absolute value, with no gaps
names(TrainingSet)  <- TrainingData[,`Gene name`]          #assigns (gene)names, which generates FGSEA-ready tstat vector
saveRDS(TrainingSet, "TTEndothelin-TrainingSet-1WAD.rds")
#then process on the cluster using do_fGSEA, reproduced here:
# library(data.table)
# library(reshape2)
# library(fgsea)
# 
# TrainingSet = readRDS("TTEndothelin-TrainingSet.rds")
# to_test = gmtPathways("LKB_Phospho_Cust.gmt") 
# 
# res <- fgsea(pathways=to_test, stats=TrainingSet, nperm=1e7, minSize = 5, nproc = 8) #runs actual FGSEA
# res[,leading_edge:=paste0(unlist(leadingEdge),collapse=", "), pathway]  		#takes data from that awful c() mess FGSEA leaves behind and prettifies
# res[,leadingEdge:=NULL]                                                 		#deletes old leading edge
# setorder(res,padj)                                                      		#orders by pvalue, descending(?)
# fwrite(res,sprintf("%s-fgsea-phospho.csv", format(Sys.time(), "%y%m%d-%H%M")))	#exports results to file. I forgot sprintf is a thing.
#

#dysferlin AS share-of-introns====
IntronData=fread("Code/Desktop/leafcutter_Hack/juncs_perIntronCounts-Raw.out", sep = "\t")
#We observed AS in exon 16...
#according to ENSEMBL: ENSOART00000012517.1 is the sole sheep transcript... it's exon 16 is at...
# Intron 15-16	93359312	93366564
# Exon 16		93366565	93366606
# Intron 16-17	93366607	93368165
#Checkin
DYSFData=IntronData[Intron.Chr==3 & Intron.Start <= 93366610 & Intron.Stop >= 93359310]
#Selects for introns on DYSF's chromosome that start before Intron16-17 starts and end after Intron15-16 has started
#with coords rounded down/up for overlap room
DYSFData[,Intron:=paste(Intron.Chr, Intron.Start, Intron.Stop, sep = ":")]
set(DYSFData, j = which(colnames(DYSFData) %in% c("Intron.Chr", "Intron.Start", "Intron.Stop")), value = NULL)

DYSFData = melt(DYSFData, id.vars = "Intron")
DYSFData[,Timepoint:=str_extract(variable, "1W|1M|3M|AD")]
DYSFData[,Skips16 := substr(Intron, 12, 20)!="93366565"] #If the intron overlaps exon 16 but does not end when exon16 should start, it must skip it

#Introns that don't stop at the canonical exon 16 start site
DYSF1 = DYSFData[Skips16 == TRUE, .(NoExon16Counts=sum(value)), .(Timepoint, variable)]

#All introns
DYSF2 = DYSFData[, .(TotalCounts=sum(value)), .(Timepoint, variable)]
#Put together...
DYSFMaths = merge(DYSF1, DYSF2, by=c("variable", "Timepoint"))
DYSFMaths[,Timepoint:=factor(Timepoint, levels = c("1W", "1M", "3M", "AD"))]

#Calculate % of no-exon-16 intron reads
DYSFMaths[,NoExon16Fraction := NoExon16Counts / TotalCounts]
setorder(DYSFMaths, Timepoint)
DYSFMaths[,mean(NoExon16Fraction), Timepoint]

#BONUS: Thyroid hormone ====
#GO:0038194 - thyroid-stimulating hormone signaling pathway (BP)
#Load data
PANTHERData = fread("Data/PANTHER/181218_ERCC.ImputeB4Exclude.batchCutRINOX3nsImputedAgeDaysdf3.maxit5000.LRT.GOFullAndSlim.One2OneOrthosOnly.txt")
if (!exists("DictOrtholog")) { DictOrtholog = fread("Data/biomart/181205_SheepMouseHuman_geneIDs_geneNames_OrthologTypes.txt") }
setkey(DictOrtholog, "Human gene stable ID")

#Process PANTHERDB.org output table through utility function:
PANTHERData2 = DoTheWholePANTHER(PANTHERData, colIDs = c("BP"="GO.BP.Complete", "MF"="GO.MF.Complete", "CC"="GO.CC.Complete"), geneTerm = "GeneID")
#Extract relevant results:
ThyroidHormoneRelatedGenes = PANTHERData2$Biological.DT[GO.BP.Complete %like% "thyroid", .(MappedIDs, GeneNameSymbolOrtholog, GO.BP.Complete)]
setkey(ThyroidHormoneRelatedGenes, "MappedIDs")
#Grab human ID:
ThyroidHormoneRelatedGenes =  DictOrtholog[ThyroidHormoneRelatedGenes]
ThyroidHormoneRelatedGenes =  ThyroidHormoneRelatedGenes[,.(`Gene stable ID`, `Human gene name`, `Gene name`, GeneNameSymbolOrtholog, GO.BP.Complete)]
#Use second utility function to plot:
Thyroid.Plots = GenerateGeneSummaryPlots(ThyroidHormoneRelatedGenes[,`Gene stable ID`], isNameList = F)
