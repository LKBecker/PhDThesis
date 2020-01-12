#===================================#
#       Lorenz' RNA-seq analysis    #
#           v1.0  August 2019       #
#           for the Dibb Lab        #
#       University Of Manchester    #
#===================================#
#libraries and functions====
library(extrafont)
source("Code/Desktop/functions.DESeq2.r")
source("Code/Desktop/functions.GO.r")
source("Code/Desktop/5.5-AlternativeSplicingAnalysis.R")
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
GroupModel = readRDS("Data/RNAseq/Output/GroupModel.ERCC.ImputeB4Exclude.2excl.maxit5000.batchCutRINOXTimepoint.Wald.RDS")
SplineModel= readRDS("Data/RNAseq/Output/SplineModel.ERCC.ImputeB4Exclude.2excl.maxit5000.batchCutRINOXnsImputedAgeDaysdf3.LRT.RDS")

Spline.ERCC.df3 = GrabResult(SplineModel, "ns.ImputedAgeDays..df...3.3", "Spline.df3", 0.05, TRUE) 
Group.ERCC.All = GetAllResultsForVariable(GroupModel, "Timepoint", 0.05, TRUE) #DO NOT SET IT TO 0.999 JFC

## HOW TO USE ## ====
# Use RStudio and ensure the project LKB_PhD_ReCode is loaded, or the working directory is set to the right target using setwd()
# Use the below code to load the needed libraries and function, as well as the RNA data
#
# The functions you want to use are GenerateGeneSummaryPlots() and AlternativeSplicingAnalysis()

##GenerateGeneSummaryPlots()====
#
NCX.DEG = GenerateGeneSummaryPlots(c("NCX (SLC8A1)"="ENSOARG00000007758"), Dataset = "ERCC+GWPM", Model = "1WAD")$Plot
NCX.DEG


##AlternativeSplicingAnalysis()====
#
NCX.AS = AlternativeSplicingAnalysis("ENSOARG00000007758", MODE = "LEAFCUTTER", DATASET = "LKB", DRAW_MODE = 3, USE_FILTER = TRUE)
NCX.AS[[1]]$ASPlot

#Export for excel / word / other====
#ffwrite(NCX.AS[[1]]$PublicationMap, "clipboard")
#ggsave()