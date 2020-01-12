#Development paper figures
#libraries and functions====
source("Code/Desktop/functions.r")
source("Code/Desktop/functions.DESeq2.r")
library(data.table)
library(extrafont)
library(ggplot2)
library(ggsignif)
library(cowplot)
library(DESeq2)
library(RColorBrewer)

#Theme====
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
#data: input====
tmsg("Loading RNASeq data for graphing, generated using file 1.3")
load("Data/RNAseq/Output/LKB.ERCC.GeneSummaryObjects.rda")
DEGs.ERCC.SplineDF3[,Comparison:="Spline.df3"] #For graph purposes

GroupModel = readRDS("Data/RNAseq/Output/GroupModel.ERCC.ImputeB4Exclude.2excl.maxit5000.batchCutRINOXTimepoint.Wald.RDS")
SplineModel= readRDS("Data/RNAseq/Output/SplineModel.ERCC.ImputeB4Exclude.2excl.maxit5000.batchCutRINOXnsImputedAgeDaysdf3.LRT.RDS")

#Ortholog data / translate to human
tmsg("generating human ortholog tables for functional enrichment...")
DictOrtholog = fread("Data/biomart/181205_SheepMouseHuman_geneIDs_geneNames_OrthologTypes.txt")
DictOrtholog[`Human gene name`!="", .N] #28521
DictOrtholog = DictOrtholog[`Human homology type`=="ortholog_one2one"]
DictOrtholog[`Human gene name`!="", .N] #17064

#Group model results
OneWeekOneMonth  <- GrabResult(GroupModel, c("Timepoint", "1M", "1W"), "1W-1M", 0.05, TRUE)
OneMonth3Months  <- GrabResult(GroupModel, c("Timepoint", "3M", "1M"), "1M-3M", 0.05, TRUE)
ThreeMonthsAdlt  <- GrabResult(GroupModel, c("Timepoint", "AD", "3M"), "3M-AD", 0.05, TRUE)
OneWeekAdult	 <- GrabResult(GroupModel, c("Timepoint", "AD", "1W"), "1W-AD", 0.05, TRUE)

#ffwrite(OneWeekAdult[baseMean>10][order(abs(log2FC), decreasing = T)][1:10], "clipboard")
#ffwrite(OneWeekAdult[baseMean>10][order(baseMean, decreasing = T)][1:10], "clipboard")
OneWAD.baseMeanTop10 = c("ENSOARG00000000032", "ENSOARG00000006257", "ENSOARG00000020185", "ENSOARG00000001811", "ENSOARG00000003713",
                         "ENSOARG00000005748", "ENSOARG00000015908", "ENSOARG00000007666", "ENSOARG00000016476", "ENSOARG00000011575")
Fig3.4 = GenerateGeneSummaryPlots(OneWAD.baseMeanTop10, isNameList = F)$Plot+facet_wrap(~geneName, ncol = 4, scales="free_y")+
    scale_y_continuous(expand = expand_scale(mult = c(0, .2)))
ffwrite(DEGs.ERCC.SplineDF3[baseMean>10][order(abs(log2FC), decreasing = T)][1:10], "clipboard")
ffwrite(DEGs.ERCC.SplineDF3[baseMean>10][order(baseMean, decreasing = T)][1:10], "clipboard")

setkey(DictOrtholog, "Gene stable ID")
SplineHuman = DictOrtholog[DEGs.ERCC.SplineDF3[,.(geneID, foldChange)], .(geneID, foldChange, `Human gene stable ID`, `Human homology type`)]
SplineHuman = unique(SplineHuman[`Human homology type`=="ortholog_one2one"])
ffwrite(SplineHuman[, .(foldChange, `Human gene stable ID`)], "SplineDF3.one2oneHumanOrtho.txt")

#Test1 all.equal(DEGs.ERCC.1WAD[,geneID], OneWeekAdult[,geneID]) #true

Full.Chrono <- rbind(OneMonth3Months, OneWeekOneMonth, ThreeMonthsAdlt) 
#Full contains Transcript stable ID, which it is best to wipe (duplication issues)
Full.Chrono = unique(Full.Chrono)
Full.Chrono = Full.Chrono[!grepl("^ERCC", geneID)]
Full.Chrono[,Regulated:=ifelse(log2FC>0, "Up", "Down")]
Full.Chrono[,Comparison:=gsub("\\-", " vs ", Comparison)]
Full.Chrono[,.N]
###ERCC: 2080, same as DESeq2Results.Group[Comparison %in% c("Group|1M|1W", "Group|3M|1M", "Group|AD|3M") & pValAdjust<0.05, .N]
Full.Chrono[,length(unique(geneID))]                        
#ERCC: 1931, same as DESeq2Results.Group.Chronological[, length(unique(geneID))], 2111 w/ spline
Full.Chrono[!grepl("^ERCC", geneID),length(unique(geneID))] #1926, e.g. 5 genes in Full Chrono are ERCCs


ChronoHuman = DictOrtholog[Full.Chrono[,.(geneID, foldChange, Comparison)], .(geneID, foldChange, `Human gene stable ID`, `Human homology type`, Comparison)]
ChronoHuman = ChronoHuman[`Human homology type`=="ortholog_one2one"]
ffwrite(ChronoHuman[, .(foldChange, `Human gene stable ID`, Comparison)], "Group.Chrono.one2oneHumanOrtho.txt")
    
tmsg("Analysing expression patterns...")
ExpressPct = ExpressPct[!grepl("^ERCC", geneID)]
GeneExpression = dcast.data.table(ExpressPct, geneID~Timepoint, value.var = "IsExpressed")
GeneExpression[,sPattern:=paste0(`1W`*1, `1M`*1, `3M`*1, `AD`*1)]

#Patterns of n genes per timepoint:
DiffExprPatternDigest=GeneExpression[,.(`n Genes`=.N),.(`1M`, `1W`, `3M`, AD)]
#DiffExprPatternDigest[,sum(`n Genes`)] #20204, all genes with counts >1 are accounted for; 19987 ERCC (shift of counts); 19901 (removing ERCC)
setcolorder(DiffExprPatternDigest, c("n Genes", "1W", "1M", "3M", "AD"))
setorder(DiffExprPatternDigest, `n Genes`)
DiffExprPatternDigest[,Pattern:=.I]
DiffExprPatternDigest[,nSwitches:=sum(abs(diff(c(`1W`,`1M`,`3M`,`AD`)))),Pattern]
DiffExprPatternDigest[,sPattern:=paste0(`1W`*1, `1M`*1, `3M`*1, `AD`*1)]
DiffExprPatternDigest.M = melt.data.table(DiffExprPatternDigest, id.vars=c("Pattern", "nSwitches"))
#Mean Expression
GeneMeans = Counts.ERCC.M[,.(MeanExpr=mean(value)),.(gene, Timepoint)]
GeneMeans = dcast(GeneMeans, gene ~ Timepoint, value.var = "MeanExpr")
GeneMeans.Named = merge(GeneMeans, unique(Dictionary[,.(`Gene stable ID`, `Gene name`)]), by.x="gene", by.y="Gene stable ID", all.x=T)

#
OnOffGenes = Counts.ERCC.M[gene %in% GeneExpression[sPattern != "1111" & sPattern !="0000", geneID]]
OnOffGenes = OnOffGenes[,mean(value), .(gene, Timepoint)]
OnOffGenes = dcast(OnOffGenes, gene~Timepoint, value.var = "V1")
setcolorder(OnOffGenes, c("gene", "1W", "1M", "3M", "AD"))
sprintf("Detected %d genes with expression patterns that are not 'always on' or 'always off'", OnOffGenes[,.N])
sprintf("Of these, %d are detected as DEG by the SPLINE model", length(intersect(OnOffGenes[,gene], DEGs.ERCC.SplineDF3[,geneID])))
sprintf("Of these, %d are detected as DEG by the GROUP model", length(intersect(OnOffGenes[,gene], Full.Chrono[,geneID])))
sprintf("A total of %d gene(s) are NOT detected by either model as DEG", length(setdiff(OnOffGenes[,gene], unique(c(Full.Chrono[,geneID], DEGs.ERCC.SplineDF3[,geneID])))))
OnOffGenes[,ExprSum:=`1W`+`1M`+`3M`+`AD`]

Overlap = length(intersect(DEGs.ERCC.SplineDF3[,geneID], OneWeekAdult[,geneID]))
sprintf( "There are %d genes in the Spline model, of which %d are also in the 1W v AD model, leaving %d unique genes; an overlap of %f%%", DEGs.ERCC.SplineDF3[,.N], 
         Overlap, length(setdiff(DEGs.ERCC.SplineDF3[,geneID], OneWeekAdult[,geneID])), (Overlap / DEGs.ERCC.SplineDF3[,.N])*100 )
SplineUnique = DEGs.ERCC.SplineDF3[geneID %in% setdiff(DEGs.ERCC.SplineDF3[,geneID], OneWeekAdult[,geneID])]
#setkey(SplineUnique, "geneID")
#setkey(DictOrtholog, "Gene stable ID")
#ffwrite(DictOrtholog[SplineUnique][!is.na(`Human gene stable ID`) & `Human homology type`=="ortholog_one2one", `Human gene stable ID`], "SplineOnly.HumanOrtho.One2One")

sprintf( "There are %d genes in the 1W v AD model, of which %d are also in the Spline model, leaving %d unique genes; an overlap of %f%%", OneWeekAdult[,.N], 
         Overlap, length(setdiff(OneWeekAdult[,geneID], DEGs.ERCC.SplineDF3[,geneID])), (Overlap / OneWeekAdult[,.N])*100 )
`1WADUnique` = OneWeekAdult[geneID %in% setdiff(OneWeekAdult[,geneID], DEGs.ERCC.SplineDF3[,geneID])]
#setkey(`1WADUnique`, "geneID")
#ffwrite(DictOrtholog[`1WADUnique`][!is.na(`Human gene stable ID`) & `Human homology type`=="ortholog_one2one", `Human gene stable ID`], "1WADOnly.HumanOrtho.One2One")

OnOffGenes.Neither = OnOffGenes[!(gene %in% Full.Chrono[,geneID]) & !(gene %in% DEGs.ERCC.SplineDF3[,geneID])]
OnOffGenes.Neither=melt(OnOffGenes.Neither, id.vars = "gene")
OnOffGenes.Neither[,mean(value), variable]
#Mostly low values - likely noise or genes that do not meet the 50% of samples threshold. 

OnOffGenes.DEG = OnOffGenes[!(gene %in% OnOffGenes.Neither[,gene])]
#rm(OnOffGenes)

Over50Genes = OnOffGenes[,.SD[max(`1W`,`1M`,`3M`,`AD`)>=50], gene][,gene] #50 or more counts at ANY timepoint and a 0 at any timepoint
GeneExpression[geneID %in% Over50Genes, .N, sPattern][order(N)]

#Graphing those genes which turn on or off and have at least 50 counts at any timepoint
Over50DiffExprPattern=GeneExpression[geneID %in% Over50Genes,.(`n Genes`=.N),.(`1M`, `1W`, `3M`, AD)]
setcolorder(Over50DiffExprPattern, c("n Genes", "1W", "1M", "3M", "AD"))
setorder(Over50DiffExprPattern, `n Genes`)
Over50DiffExprPattern[,Pattern:=.I]
Over50DiffExprPattern[,nSwitches:=sum(abs(diff(c(`1W`,`1M`,`3M`,`AD`)))),Pattern]
Over50DiffExprPattern[,sPattern:=paste0(`1W`*1, `1M`*1, `3M`*1, `AD`*1)]
Over50DiffExprPattern.M = melt.data.table(Over50DiffExprPattern, id.vars=c("Pattern", "nSwitches"))

Fig3.1B = ggplot(Over50DiffExprPattern.M[!(variable %in% c("n Genes", "sPattern"))] )+geom_raster(aes(fill=factor(value), x=variable, y=Pattern))+
    geom_text(data=Over50DiffExprPattern.M[variable=="n Genes"], aes(label=value, x=variable, y=Pattern), size=geom.text.size/ggplot2::.pt)+xlab("")+ylab("")+
    theme(legend.position = "none", axis.line = element_blank(), axis.ticks = element_blank(), plot.margin = unit(c(10,30,0,30), "pt"))+
    scale_x_discrete(limits=c("n Genes", "1W", "1M", "3M", "AD"))+
    scale_y_continuous(breaks=seq(1, nrow(Over50DiffExprPattern.M), 1), labels = NULL, expand=c(0,0))+
    scale_fill_manual(values = c("#ffffff", "#a0a0a0"))+geom_hline(yintercept=seq(1.5, nrow(Over50DiffExprPattern)+0.5, 1))

Over50.Explore = merge(GeneExpression[geneID %in% Over50Genes], Dictionary[`Gene stable ID` %in% Over50Genes], by.x="geneID", by.y="Gene stable ID")

EarlyOnGenes = Counts.ERCC.M[gene %in% GeneExpression[AD==F & (`1W`==T | `1M`==T | `3M`==T), geneID]]
EarlyOnGenes = EarlyOnGenes[,mean(value), .(gene, Timepoint)]
EarlyOnGenes = dcast(EarlyOnGenes, gene~Timepoint, value.var = "V1")
EarlyOnGenes[,ExprSum:=`1W`+`1M`+`3M`+`AD`]
EarlyOnGenes = merge(EarlyOnGenes, DictOrtholog, by.x="gene", by.y="Gene stable ID", all.x=T)

LateOnGenes = Counts.ERCC.M[gene %in% GeneExpression[AD==T & (`1W`==F | `1M`==F | `3M`==F), geneID]]
LateOnGenes = LateOnGenes[,mean(value), .(gene, Timepoint)]
LateOnGenes = dcast(LateOnGenes, gene~Timepoint, value.var = "V1")
LateOnGenes[,ExprSum:=`1W`+`1M`+`3M`+`AD`]
LateOnGenes = merge(LateOnGenes, DictOrtholog, by.x="gene", by.y="Gene stable ID", all.x=T)

ExpressedGenes=Counts.ERCC.M[Is_Expr==TRUE,.N,.(sample, Timepoint)] 
#Fetch the number of expressed genes (as defined by CONF$EXPRESSED_MIN_COUNTS) for each sample
ExpressedGenes[,Timepoint:=factor(Timepoint, levels=c("1W","1M","3M","AD"), ordered = T)]

#data: figure 1====
tmsg("Generating Fig 1 data...")
Full.All <- rbind(OneMonth3Months, OneWeekOneMonth, ThreeMonthsAdlt, DEGs.ERCC.SplineDF3[pValAdjust<0.05]) 
#Full contains Transcript stable ID, which it is best to wipe (duplication issues)
Full.All = unique(Full.All)
Full.All = Full.All[!grepl("^ERCC", geneID)]
Full.All[,Regulated:=ifelse(log2FC>0, "Up", "Down")]
Full.All[,Comparison:=gsub("\\-", " vs ", Comparison)]

OneWeekAdult[,Comparison:=gsub("\\-", " vs ", Comparison)]
OneWeekAdult[,Regulated:=ifelse(log2FC>0, "Up", "Down")]

Full2 = merge(Full.All, GeneMeans, by.x="geneID", by.y="gene") #WARNING seems to act as unique()
set(Full2, j = Full2[, which(colnames(Full2) %in% c("ChangeFrom", "1W", "1M", "3M", "AD"))], value=NULL)

tmsg("summarising group- and spline models...")
#FullStats <- Full2[,.(Log2FC=sum(log2FC), NRegs = .N, SumBaseMean=sum(GeneChange)),.(Comparison, Regulated)] 
#This is a different BaseMean than DESeq2's baseMean;
FullStats <- Full2[,.(Log2FC=sum(log2FC), NRegs = .N),.(Comparison, Regulated)] #This is a different BaseMean than DESeq2's baseMean;
FullStats[,SumN:=sum(NRegs), Comparison]
FullStats[,Pct:=abs(NRegs/SumN)]
FullStats[Regulated=="Down", Pct:=-Pct]
FullStats[Regulated=="Down", NRegs:=-NRegs]
#FullStats[,Comparison:=factor(Comparison, levels = c("1W vs 1M", "1M vs 3M", "3M vs AD"))]
FullStats[,Comparison:=factor(Comparison, levels = c("1W vs 1M", "1M vs 3M", "3M vs AD", "Spline.df3"))]

#Patterns of when which genes switch expression how

Full.Patterns=dcast.data.table(Full.Chrono[,.(Comparison, Regulated, geneID)], formula = geneID~Comparison, fill = NA, value.var = "Regulated", 
                    fun.aggregate = c)
Full.Patterns=Full.Patterns[,.(`n Genes`=.N), .(`1W vs 1M`, `1M vs 3M`, `3M vs AD`)]
Full.Patterns[,sum(`n Genes`)] #1931

Full.Patterns[(`1W vs 1M`=="Up"|is.na(`1W vs 1M`)) & (`1M vs 3M`=="Up"|is.na(`1M vs 3M`)) & is.na(`3M vs AD`), Type:="Early Rise"]
Full.Patterns[is.na(Type) & (`1W vs 1M`=="Down"|is.na(`1W vs 1M`)) & (`1M vs 3M`=="Down"|is.na(`1M vs 3M`)) & is.na(`3M vs AD`), Type:="Early Fall"]
Full.Patterns[is.na(Type) & is.na(`1W vs 1M`) & (`1M vs 3M`=="Up"|is.na(`1M vs 3M`)) & `3M vs AD`=="Up", Type:="Late Rise"]
Full.Patterns[is.na(Type) & is.na(`1W vs 1M`) & (`1M vs 3M`=="Down"|is.na(`1M vs 3M`)) & `3M vs AD`=="Down", Type:="Late Fall"]
Full.Patterns[is.na(Type) & `1W vs 1M`=="Down" & `3M vs AD`=="Down", Type:="Global Fall"]
Full.Patterns[is.na(Type) & `1W vs 1M`=="Up" & `3M vs AD`=="Up", Type:="Global Rise"]
Full.Patterns[is.na(Type), Type:="Mixed"]
Full.Patterns[,sum(`n Genes`)]
#still 1931 (ERCC), so the 'movement' of some from early to late must be a result of reassignment through the above

FullPatterns2 = melt(Full.Patterns, id.vars="Type")[variable!="n Genes"]
FullPatterns2 = FullPatterns2[,.SD[,.N,value][order(N, decreasing = T)][1],.(Type, variable)][,.(Type, variable, value)]
FullPatterns2 = dcast.data.table(FullPatterns2, Type~variable, value.var="value")
FullPatterns3 = merge(FullPatterns2, Full.Patterns[,.("n Genes"=sum(`n Genes`)), Type], by="Type")
FullPatterns3[,sum(as.numeric(`n Genes`))] #1931 (ERCC)
setorder(FullPatterns3, `n Genes`)
FullPatterns3[,Pattern:=.I]
FullPatterns3=melt.data.table(FullPatterns3, id.vars="Pattern")
FullPatterns3[is.na(value),value:="No change"]
FullPatterns3[variable=="n Genes", sum(as.numeric(value))] #2104 (no ERCC) / 1931 (ERCC)

#TODO: any above 100 counts ever?

rm(FullPatterns2)

#DEG summary stats
OneWeekADStats = merge(OneWeekAdult, GeneMeans, by.x="geneID", by.y="gene") #WARNING seems to act as unique()
OneWeekADStats = OneWeekADStats[,.(Log2FC=sum(log2FC), NRegs = .N),.(Comparison, Regulated)]
OneWeekADStats[,SumN:=sum(NRegs), Comparison]
OneWeekADStats[,Pct:=abs(NRegs/SumN)]
OneWeekADStats[Regulated=="Down", Pct:=-Pct]
OneWeekADStats[Regulated=="Down", NRegs:=-NRegs]

AnyDEG.HH = unique(merge(Full.All, DictOrtholog, by.x="geneID", by.y="Gene stable ID", all.x=T))
AnyDEG.HH[,.N, .(is.na(`Human gene name`), Comparison)][, .(`% DEGs lost due to no 1:1 ortholog`=.SD[is.na==TRUE, N]/sum(N)), Comparison] #15.5% loss one2one, again
#Rates of loss roughly equal between models.

if(FALSE){
    #Run this manually if the graph needs to be totally redone, not just re-drawn - this makes the input files for toppFun
    # ffwrite(AnyDEG.HH[Comparison=="Group|1M|1W" & !is.na(`Human gene name`), unique(`Human gene name`)], "Group1M1W-HHone2one-NewDEG")
    # ffwrite(AnyDEG.HH[Comparison=="Group|3M|1W" & !is.na(`Human gene name`), unique(`Human gene name`)], "Group3M1W-HHone2one-NewDEG")
    # ffwrite(AnyDEG.HH[Comparison=="Group|AD|1W" & !is.na(`Human gene name`), unique(`Human gene name`)], "GroupAD1W-HHone2one-NewDEG")
    # ffwrite(AnyDEG.HH[Comparison=="Group|3M|1M" & !is.na(`Human gene name`), unique(`Human gene name`)], "Group3M1M-HHone2one-NewDEG")
    # ffwrite(AnyDEG.HH[Comparison=="Group|AD|3M" & !is.na(`Human gene name`), unique(`Human gene name`)], "GroupAD3M-HHone2one-NewDEG")
    # ffwrite(AnyDEG.HH[Comparison=="Spline.df3" & !is.na(`Human gene name`),unique(`Human gene name`)], "SplineAll-HHone2one-NewDEG")
    #Now run those files through https://toppgene.cchmc.org/enrichment.jsp and save output, replacing filenames below
} 

rm(OneMonth3Months, ThreeMonthsAdlt, OneWeekOneMonth)

#data: figure 2====
tmsg("processing toppFun results...")
W1M1.Up     = fread("Data/ToppFun/190429_ToppFun.Group.HumanOne2One.1W1M.Up.txt")[,Direction:="Up"][,Group:="1W v 1M"]
W1M1.Down   = fread("Data/ToppFun/190429_ToppFun.Group.HumanOne2One.1W1M.Down.txt")[,Direction:="Down"][,Group:="1W v 1M"]
W1M1 = rbind(W1M1.Up, W1M1.Down)
M1M3.Up     = fread("Data/ToppFun/190429_ToppFun.Group.HumanOne2One.1M3M.Up.txt")[,Direction:="Up"][,Group:="1M v 3M"]
M1M3.Down   = fread("Data/ToppFun/190429_ToppFun.Group.HumanOne2One.1M3M.Down.txt")[,Direction:="Down"][,Group:="1M v 3M"]
M1M3 = rbind(M1M3.Up, M1M3.Down)
M3AD.Up     = fread("Data/ToppFun/190429_ToppFun.Group.HumanOne2One.3MAD.Up.txt")[,Direction:="Up"][,Group:="3M v AD"]
M3AD.Down   = fread("Data/ToppFun/190429_ToppFun.Group.HumanOne2One.3MAD.Down.txt")[,Direction:="Down"][,Group:="3M v AD"]
M3AD = rbind(M3AD.Up, M3AD.Down)
SPLN.Up     = fread("Data/ToppFun/190429_ToppFun.SplineDF3.HumanOne2One.Up.txt")[,Direction:="Up"][,Group:="Spline"]
SPLN.Down   = fread("Data/ToppFun/190429_ToppFun.SplineDF3.HumanOne2One.Down.txt")[,Direction:="Down"][,Group:="Spline"]
SPLN = rbind(SPLN.Up, SPLN.Down)
rm(W1M1.Up, W1M1.Down, M1M3.Up, M1M3.Down, M3AD.Up, M3AD.Down, SPLN.Up, SPLN.Down)

#Reduce to only GO results (Exclude e.g. transcription factor motifs, whose results name weird names)
W1M1.GO = W1M1[substr(Category,1,2)=="GO"]
M1M3.GO = M1M3[substr(Category,1,2)=="GO"]
M3AD.GO = M3AD[substr(Category,1,2)=="GO"]
#W1AD.GO2= copy(W1AD.GO)
SPLN.GO = SPLN[substr(Category,1,2)=="GO"]
rm(W1M1, M1M3, M3AD, SPLN)

# #Grab the Top 10 (by q value) for each of the 3 Categories
W1M1.GO2 = W1M1.GO[order(`q-value FDR B&H`),.SD[1:10], .(Category, Direction)]
M1M3.GO2 = M1M3.GO[order(`q-value FDR B&H`),.SD[1:10], .(Category, Direction)]
M3AD.GO2 = M3AD.GO[order(`q-value FDR B&H`),.SD[1:10], .(Category, Direction)]
SPLN.GO2 = SPLN.GO[order(`q-value FDR B&H`),.SD[1:10], .(Category, Direction)]

#And bind
All = rbind( W1M1.GO2, M1M3.GO2, M3AD.GO2, SPLN.GO2 )
TopTerms=All[,unique(Name)]
rm( W1M1.GO2, M1M3.GO2, M3AD.GO2, SPLN.GO2 )
All = rbind( W1M1.GO[Name %in% TopTerms], M1M3.GO[Name %in% TopTerms], M3AD.GO[Name %in% TopTerms], SPLN.GO[Name %in% TopTerms] )
rm( W1M1.GO, M1M3.GO, M3AD.GO, SPLN.GO )

All[,GraphVal:= - log10(`q-value FDR B&H`)]
All[,Group:=factor(Group, levels=c("1W v 1M", "1M v 3M", "3M v AD", "Spline"))]

BioData = All[Category=="GO: Biological Process"]
BioData = rbind(BioData, data.table(Name=BioData[,Name], Direction="Up", DirGroup="3M v AD.Up", Group="3M v AD", GraphVal=NA), fill=T) #There's no 'natural' 3M v AD up BPs, so we add a blank row

#Category sort
BioData[Name %in% c("cell activation", "cell cycle", "cell cycle process", "cell division", "chromosome segregation", 
                    "mitotic cell cycle", "mitotic cell cycle process", "nuclear chromosome segregation", "nuclear division", "sister chromatid segregation"),OverCategory:="Cell Cycle"]
BioData[Name %in% c("apoptotic process involved in development", "apoptotic process involved in heart morphogenesis", "apoptotic process involved in morphogenesis", 
                    "blood vessel development", "blood vessel morphogenesis", "cardiovascular system development", "circulatory system development", "cell motility", 
                    "localization of cell", "movement of cell or subcellular component", "organelle fission", "positive regulation of cell differentiation", 
                    "positive regulation of developmental process", "regulation of cell adhesion", "regulation of cellular component movement", 
                    "skin morphogenesis", "vasculature development", "vasculogenesis"), OverCategory:="Tissue/organ development"]
BioData[Name %in% c("collagen catabolic process", "collagen fibril organization", "collagen metabolic process", 
                    "extracellular matrix organization", "extracellular structure organization"),OverCategory:="Extracellular Matrix"]
BioData[Name %in% c("actomyosin structure organization", "contractile actin filament bundle assembly", "positive regulation of actin filament bundle assembly", 
                    "regulation of actin filament bundle assembly", "regulation of stress fiber assembly", "stress fiber assembly"),OverCategory:="Muscle fiber development"]
BioData[Name %in% c("leukocyte activation", "leukocyte aggregation", "leukocyte cell-cell adhesion", "lymphocyte activation", 
                    "lymphocyte aggregation", "regulation of immune system process", "T cell activation", "T cell aggregation"),OverCategory:="Immune system"]
BioData[Name %in% c("cellular response to amino acid stimulus", "enzyme linked receptor protein signaling pathway", 
                    "platelet activation", "protein heterotrimerization", "regeneration", "response to endogenous stimulus", "response to mechanical stimulus", 
                    "response to progesterone", "small GTPase mediated signal transduction"),OverCategory:="Other"]
BioData[,OverCategory:=factor(OverCategory, c("Cell Cycle", "Tissue/organ development", "Muscle fiber development", "Extracellular Matrix", "Immune system", "Other"))]
BioData[,Name := factor(Name, levels=BioData[order(OverCategory), unique(Name)])]

BioData.Up = BioData[Direction == "Up" | is.na(Direction)]
BioData.Up[,DirGroup:=factor(paste0(Group, ".", Direction), levels=c("Spline.Up", "1W v 1M.Up", "1M v 3M.Up", "3M v AD.Up"))]
BioData.Up = completeDT(BioData.Up, c("Name", "DirGroup"), defs=c(GraphVal=NA)) #Filling out the DT to make absent values EXPLICITLY NA,

BioData.Down = BioData[Direction=="Down" | is.na(Direction)]
BioData.Down[,DirGroup:=factor(paste0(Group, ".", Direction), levels=c("Spline.Down", "1W v 1M.Down", "1M v 3M.Down", "3M v AD.Down"))]
BioData.Down = rbind(BioData.Down, data.table(Name=setdiff(BioData[,unique(Name)], BioData.Down[,unique(Name)]), Direction="Down", DirGroup=BioData.Down[,unique(DirGroup)], 
                                              Group=BioData.Down[,unique(Group)], GraphVal=NA), fill=T) #There's no 'natural' 3M v AD up BPs, so we add a blank row
BioData.Down = completeDT(BioData.Down, c("Name", "DirGroup"), defs=c(GraphVal=NA)) #Filling out the DT to make absent values EXPLICITLY NA,

#because ggplot::scale_fill_gradient() has no attribute to explicitly color absent/inferred data
MolData = All[Category=="GO: Molecular Function"]
MolData = completeDT(MolData, c("Name", "Group"), defs=c(GraphVal=NA)) #Filling out the DT to make absent values EXPLICITLY NA,
#TODO: truncate MolData labels
MolData[,Name:=str_trunc(Name, 50, "center")]
rm(All)

#data: figure 3 (Transcription factors)====
#IPA.Input = GroupResults.HH[Comparison=="Timepoint|AD|1W" & !is.na(`Human gene name`), .(`Human gene name`,log2FC, pValAdjust)]
#ffwrite(IPA.Input, "IPA.TransFactorAnalysis.1WAD")
tmsg("processing IPA results...")
# GroupResForIPA = GroupResults[Comparison=="Timepoint|AD|1W",.(geneID, log2FC, pValAdjust)]
# GroupResForIPA = merge(GroupResForIPA, DictOrtholog[,.(`Gene stable ID`, `Human gene name`)], by.x="geneID", by.y="Gene stable ID", all.x=T)
# GroupResForIPA = GroupResForIPA[!is.na(`Human gene name`)]

IPA.TF.1WAD = fread("Data/IPA/1WAD-TranscriptionFactor-IPA-HHone2one-NewDEG.txt", sep = "\t", dec = ",")
IPA.TF.1WAD[,V9:=NULL]
IPA.TF.1WAD = IPA.TF.1WAD[`p-value of overlap`<0.05]
IPA.TF.1WAD = IPA.TF.1WAD[`Molecule Type`=="transcription regulator"]
IPA.TF.1WAD[!is.na(`Activation z-score`),AbsZScore:=abs(`Activation z-score`)]
#IPA.TF.1WAD = merge(IPA.TF.1WAD, GroupResForIPA, by.x="Upstream Regulator", by.y="Human gene name", all.x=T)
#View(IPA.TF.1WAD[,.(`Expr Log Ratio`,log2FC)]) #values are sort of... distorted? Difficulty in IPA matching IDs?

#IPA/ToppFun overlap
W1AD.GO2 = fread("Data/ToppFun/181214_ToppFun_AD1W_HSone2one_newDEG.txt")
W1AD.GO2 = W1AD.GO2[substr(Category,1,2)=="GO"]
W1AD.GO2 = sort(unique(SplitDataTableWithMultiRows(W1AD.GO2, "Hit in Query List", ",")[,`Hit in Query List`]))

IPA.TF.1WAD.LOverlap = IPA.TF.1WAD[,.(totalLen=str_count(`Target molecules in dataset`, ",")+1, 
                                      nOverlap=length(intersect(unlist(strsplit(`Target molecules in dataset`, ",")), W1AD.GO2)) ), `Upstream Regulator`]
IPA.TF.1WAD.LOverlap[,Percent:=(nOverlap/totalLen)*100]

TP53Targets = IPA.TF.1WAD[`Upstream Regulator`=="TP53",unlist(strsplit(`Target molecules in dataset`, ","))]

IPA.TF.Spline = fread("Data/IPA/Spline.DF3-TranscriptionFactor-IPA-HHone2one-NewDEG.txt", sep = "\t")
IPA.TF.Spline[,V9:=NULL]
IPA.TF.Spline = IPA.TF.Spline[`p-value of overlap`<0.05]
IPA.TF.Spline = IPA.TF.Spline[`Molecule Type`=="transcription regulator"]
IPA.TF.Spline[!is.na(`Activation z-score`),AbsZScore:=abs(`Activation z-score`)]
IPA.TF.Spline = merge(IPA.TF.Spline, DEGs.ERCC.SplineDF3, by.x="Upstream Regulator", by.y="Gene name", all.x=T)

IPA.TF.1WAD.Expressed = GenerateGeneSummaryPlots(c("ENSOARG00000011054","E2F1"="ENSOARG00000008548","ENSOARG00000017221","ENSOARG00000009655","ENSOARG00000009493",
                                                   "ENSOARG00000017426","ENSOARG00000008027","ENSOARG00000006421","ENSOARG00000003184","ENSOARG00000003450",
                                                   "TBX2"="ENSOARG00000016407", "ENSOARG00000015228", "ENSOARG00000006614", "ENSOARG00000004153"))
#IPA.TF.1WAD.Expressed$Plot

#NKX2-5, TBX5, GATA4 all have a strange bump at 1 months (a decrease)
#https://www.sciencedirect.com/science/article/pii/S0022282803000026 ?
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC309615/ MEF2C ANKRD1 HAND1 but not ANPA ANPB SCN5A TAGLN(more or less)
WTF = GenerateGeneSummaryPlots(c("ENSOARG00000001511", "ENSOARG00000015228", "ENSOARG00000006614", "ENSOARG00000004153", "ENSOARG00000019963", "ENSOARG00000001512",
                                 "ENSOARG00000015908", "ENSOARG00000009843", "ENSOARG00000016099", "ENSOARG00000017779", "ENSOARG00000004174", "ENSOARG00000004347", 
                                 "ENSOARG00000016462", "ENSOARG00000015496", "ENSOARG00000012149", "ENSOARG00000007915"))$Plot

#IPA.TF.Long = SplitDataTableWithMultiRows(IPA.TF, "Target molecules in dataset", ",")

#data: figure 4 (leafcutter)====
combinations = c("1W1M", "1M3M", "3MAD")
LoadResult <- function(token, sourcepath="../leafcutter/"){
    #token="1W1M"; sourcepath="../leafcutter/"
    LeafCutterObj = list()
    ClusterFile = fread(paste0(sourcepath, token, "_cluster_significance.txt"))
    IntronFile  = fread(paste0(sourcepath, token, "_effect_sizes.txt"))
    IntronFile[,location:=str_extract(intron, "^chr[\\d{1,2}MTXY]\\:\\d+\\:\\d+")]
    IntronFile[,cluster:=str_extract(intron, "clu_\\d+_.*$")]
    IntronFile[,cluster:=paste0(str_extract(location, "^chr[\\d{1,2}MTXY]"), ":", cluster)]
    #IntronFile  
    ClusterFile
}

W1M1= LoadResult(combinations[[1]])[,Comparison:="1W vs 1M"]
M1M3= LoadResult(combinations[[2]])[,Comparison:="1M vs 3M"]
M3AD= LoadResult(combinations[[3]])[,Comparison:="3M vs AD"]
W1AD = LoadResult("1WAD")[,Comparison:="1W vs AD"]
Fin = rbind(W1M1, M1M3, M3AD, W1AD, fill=T)
Fig4AData = Fin[p.adjust<0.05,.N, Comparison]
Fig4AData[,Comparison:=factor(Comparison, levels=c("1W vs 1M", "1M vs 3M", "3M vs AD", "1W vs AD"))]

tmsg("data I/O complete. Proceeding to construction of graphs...")
#figure 1: overview over DEGs====
#question: how does gene expression change in the developing sheep atrium?

tbl = table(ExpressPct$IsExpressed, ExpressPct$Timepoint)
chisq.test(tbl) 
#p=0.03092 w  ERCC and removing all where rowSum(counts) == 0 (also 0.13 when not removing)
#Gene expression changes over time when applying the 50% filter

#Fisher's test - are the number of DEGs per comparison changing?
ExampleCount = fread("Data/RNAseq/Input/CountFiles.ERCC/Sample161206-01_1W_OAv3.1+ERCC_aligned_ENSEMBL_Union_RevStranded.counts")
DetectableGenes = ExampleCount[substr(V1,1,1)!="_" & substr(V1,1,4)!="ERCC",.N] #number of all genes in the reference genome.
rm(ExampleCount)
FisherStats     <- Full2[,.N,Comparison]  # the number of genes changed per timepoint

#1W/1M vs 1M/3M:
FisherTable     <- matrix(ncol=2, nrow=2)
FisherTable[1,] <- FisherStats[,N][2:3] #Number of Changing Genes 1W v 1M and 1M v 3M
FisherTable[2,] <- DetectableGenes-FisherTable[1,]
fisher.test(FisherTable) #0.0007, so yes

#1M/3M vs 3M/AD:
FisherTable     <- matrix(ncol=2, nrow=2)
FisherTable[1,] <- FisherStats[,N][3:4]
FisherTable[2,] <- DetectableGenes-FisherTable[1,]
fisher.test(FisherTable) # < 2.2e-16, extremely yes

#Fig3.1A - which genes turn on/off?
Fig3.1A = ggplot(DiffExprPatternDigest.M[!(variable %in% c("n Genes", "sPattern"))] )+geom_raster(aes(fill=factor(value), x=variable, y=Pattern))+
    geom_text(data=DiffExprPatternDigest.M[variable=="n Genes"], aes(label=value, x=variable, y=Pattern), size=geom.text.size/ggplot2::.pt)+xlab("")+ylab("")+
    theme(legend.position = "none", axis.line = element_blank(), axis.ticks = element_blank(), plot.margin = unit(c(10,30,0,30), "pt"))+
    scale_x_discrete(limits=c("n Genes", "1W", "1M", "3M", "AD"))+
    scale_y_continuous(breaks=seq(1, nrow(DiffExprPatternDigest.M), 1), labels = NULL, expand=c(0,0))+
    scale_fill_manual(values = c("#ffffff", "#a0a0a0"))+geom_hline(yintercept=seq(1.5, nrow(DiffExprPatternDigest)+0.5, 1))

#Column graph or Table giving the numbers of DEGs per model
Table1B=rbind(FullStats, OneWeekADStats)

Fig3.2A.All.Btm = ggplot(Table1B, aes(x=Comparison, y=NRegs, fill=Regulated))+geom_bar(stat="identity", width=0.5)+geom_hline(yintercept = 0)+
    scale_y_continuous(labels = abs, expand = c(0,0), limits = c(-3000,3000), breaks = seq(-2750, 2250, 500))+ylab("")+
    scale_fill_manual(values=c("#e41a1c", "#4daf4a"))+theme(legend.position = "none", plot.margin = unit(c(0,0,0,4), "pt"),
        axis.title.y = element_text(margin = unit(c(0,5,0,0), "pt")), axis.text.x = element_text(angle=35, hjust = 1))+
    coord_cartesian(ylim = c(-2750, -2250))

Fig3.2A.All.Med = ggplot(Table1B, aes(x=Comparison, y=NRegs, fill=Regulated))+geom_bar(stat="identity", width=0.5)+geom_hline(yintercept = 0)+
    scale_y_continuous(labels = abs, expand = c(0,0), limits = c(-3000,3000), breaks = seq(-750, 750, 750))+ylab("Number of changes")+
    scale_fill_manual(values=c("#e41a1c", "#4daf4a"))+theme(legend.position = "none", plot.margin = unit(c(0,0,0,7), "pt"),
        axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank())+xlab("")+
    geom_signif(y_position = c(660), comparisons = list(c("1W vs 1M", "1M vs 3M"), c("1M vs 3M", "3M vs AD")), annotation=c("***", "***"), 
                tip_length=0.01, color="#000000", size=1, textsize = geom.text.size)+coord_cartesian(ylim = c(-750, 750))

Fig3.2A.All.Top = ggplot(Table1B, aes(x=Comparison, y=NRegs, fill=Regulated))+geom_bar(stat="identity", width=0.5)+geom_hline(yintercept = 0)+
    scale_y_continuous(labels = abs, expand = c(0,0), limits = c(-3000,3000), breaks = seq(2250, 2750, 500))+ylab("")+
    scale_fill_manual(values=c("#e41a1c", "#4daf4a"))+theme(legend.position = "none", plot.margin = unit(c(10,0,0,0), "pt"),
        axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank())+xlab("")+
    coord_cartesian(ylim = c(2250, 2750))

Fig3.2A = plot_grid(Fig3.2A.All.Top, Fig3.2A.All.Med, Fig3.2A.All.Btm, rel_heights = c(0.1, 0.2, 0.2), ncol = 1)

Fig3.2AAlt.All.Btm = ggplot(Table1B, aes(x=Comparison, y=Log2FC, fill=Regulated))+geom_bar(stat="identity", width=0.5)+geom_hline(yintercept = 0)+
    scale_y_continuous(labels = abs, expand = c(0,0), limits = c(-4000, 3500), breaks = seq(-4000, -2250, 1750))+ylab("")+
    scale_fill_manual(values=c("#e41a1c", "#4daf4a"))+theme(legend.position = "none", plot.margin = unit(c(0,0,0,4), "pt"),
        axis.title.y = element_text(margin = unit(c(0,5,0,0), "pt")), axis.text.x = element_text(angle=35, hjust = 1))+
    coord_cartesian(ylim = c(-2250, -4000))

Fig3.2AAlt.All.Med = ggplot(Table1B, aes(x=Comparison, y=Log2FC, fill=Regulated))+geom_bar(stat="identity", width=0.5)+geom_hline(yintercept = 0)+
    scale_y_continuous(labels = abs, expand = c(0,0), limits = c(-3500,3500), breaks = seq(-750, 750, 750))+ylab("Number of changes")+
    scale_fill_manual(values=c("#e41a1c", "#4daf4a"))+theme(legend.position = "none", plot.margin = unit(c(0,0,0,7), "pt"),
                                                            axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank())+xlab("")+
    geom_signif(y_position = c(660), comparisons = list(c("1W vs 1M", "1M vs 3M"), c("1M vs 3M", "3M vs AD")), annotation=c("***", "***"), 
                tip_length=0.01, color="#000000", size=1, textsize = geom.text.size)+coord_cartesian(ylim = c(-750, 750))

Fig3.2AAlt.All.Top = ggplot(Table1B, aes(x=Comparison, y=Log2FC, fill=Regulated))+geom_bar(stat="identity", width=0.5)+geom_hline(yintercept = 0)+
    scale_y_continuous(labels = abs, expand = c(0,0), limits = c(-3000,4000), breaks = seq(1500, 4000, 2500))+ylab("")+
    scale_fill_manual(values=c("#e41a1c", "#4daf4a"))+theme(legend.position = "none", plot.margin = unit(c(10,0,0,0), "pt"),
                                                            axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank())+xlab("")+
    coord_cartesian(ylim = c(1500, 4000))

Fig3.2AAlt = plot_grid(Fig3.2AAlt.All.Top, Fig3.2AAlt.All.Med, Fig3.2AAlt.All.Btm, rel_heights = c(0.1, 0.2, 0.2), ncol = 1)

#Venn diagram showing the overlap/difference in DEGs of the different models
Fig3.2B.Data = rbind(Full.All[pValAdjust<0.05,.(geneID, Comparison, Log2FCAbs)],
                  OneWeekAdult[pValAdjust<0.05,.(geneID, Comparison, Log2FCAbs)])
Fig3.2B.Data = Fig3.2B.Data[order(Comparison)]
#ffwrite(Fig3.2B.Data[substr(geneID, 1,4)!="ERCC"], "VennDiagramData", dec=",")

#ffwrite(rbind(OneWeekAdult, Full.All[pValAdjust<0.05]), paste0("VennV2.", CONF$SPLINE_DF), dec=",")
#Using http://bioinformatics.psb.ugent.be/webtools/Venn/

#Fig1E is a venn diagram again, this time 1W1M, 1M3M, 3MAD
#ffwrite(Full.All[pValAdjust<0.05], "Venn.df3.Groups", dec=",")

#'Heatmap' of regulation patterns (regulation, not on/off)
Fig3.3B = ggplot(FullPatterns3[variable!="n Genes" & variable!="Type"])+geom_raster(aes(fill=factor(value), x=variable, y=Pattern))+
    scale_fill_manual(name="Gene regulation", values = c("#e41a1c", "#494949", "#4daf4a"))+
    theme(legend.position = "none", axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(angle=90, hjust=1))+
    ylab("")+xlab("")+scale_x_discrete(limits=c("Type", "n Genes", "1W vs 1M", "1M vs 3M", "3M vs AD"))+
    scale_y_continuous(breaks=seq(1, nrow(Full.Patterns), 1), labels = NULL, expand = c(0,0))+
    geom_text(data=FullPatterns3[variable=="Type"], aes(label=value, x="Type", y=Pattern), size=geom.text.size/ggplot2::.pt)+
    geom_text(data=FullPatterns3[variable=="n Genes"], aes(label=value, x="n Genes", y=Pattern), size=geom.text.size/ggplot2::.pt)+
    geom_hline(yintercept=seq(0.5, FullPatterns3[,max(Pattern)]+0.5, 1))

#figure 2: Prosser style GO analysis====
Fig3.6.DummyPlot = ggplot(BioData.Up, aes(x=DirGroup, y=Name, fill=GraphVal))+geom_raster()+scale_fill_gradientn(limits=c(0, 80), breaks=c(0,80), name="-log10 qval", 
    colours = c("#ff0000", "#00ff00"), na.value = "#555555")+geom_point(data = BioData.Up, aes(size="NA"), shape = NA)+scale_colour_manual(values=NA) + 
    guides(size=guide_legend("No data", override.aes=list(shape=15, size = 10, colour = "#555555"), keywidth = 3, keyheight = 1))+theme(legend.position = "left")
Fig3.6.Legend = get_legend(Fig3.6.DummyPlot)
rm(Fig3.6.DummyPlot)

Fig3.6.PanelUp = ggplot(BioData.Up, aes(x=DirGroup, y=Name, fill=GraphVal))+geom_raster()+scale_fill_gradientn(limits=c(0, 80), breaks=c(0,80), name="-log10 qval", 
    colours = c("#ff0000", "#00ff00"), na.value = "#555555")+xlab("")+ylab("")+theme(axis.text.x = element_text(angle=90, hjust = 1), legend.position = "none",
    axis.ticks.length = unit(3, "pt"))+geom_hline(yintercept = c(10.5, 28.5, 34.5, 39.5, 47.5), size=1, color="#999999")+
    geom_vline(xintercept = seq(1.5, 7.5), size=1, color="#a0a0a0")

Fig3.6.PanelDown = ggplot(BioData.Down, aes(x=DirGroup, y=Name, fill=GraphVal))+geom_raster()+scale_fill_gradientn(limits=c(0, 80), breaks=c(0,80),
    name="-log10 qval", colours = c("#ff0000", "#00ff00"), na.value = "#555555")+xlab("")+ylab("")+
    theme(axis.text.x = element_text(angle=90, hjust = 1), legend.position = "none", axis.text.y = element_blank(), axis.ticks.length = unit(3, "pt"),
          plot.margin = margin(0,0,0,0))+
    geom_hline(yintercept = c(10.5, 28.5, 34.5, 39.5, 47.5), size=1, color="#999999")+ geom_vline(xintercept = seq(1.5, 7.5), size=1, color="#a0a0a0")
    
Fig3.6 = plot_grid(Fig3.6.PanelUp, Fig3.6.PanelDown, nrow=1, rel_widths = c(0.8, 0.2), axis = "tb", align = "h")
#Fig3.6 = plot_grid(Fig3.6.Legend, Fig3.6, rel_widths = c(0.08, 0.92))

Fig3.6.Alt = ggplot(MolData, aes(x=Group, y=Name, fill=GraphVal))+geom_raster()+	scale_fill_gradientn(limits=c(0, 15), breaks=c(0,15),
        name="-log10 qval", colours = c("#ff0000", "#00ff00"), na.value = "#555555")+
    theme(axis.text.x = element_text(angle=90, hjust = 1), legend.position = "left")+xlab("")+ylab("")

#TODO: check out tubulin binding, microtubule motor actibity in 1W v 1M dataset...

#figure 3: transcription factor analysis====
IPA.TF.1WAD=IPA.TF.1WAD[!is.na(`Activation z-score`)]
Table3.Suppl =IPA.TF.1WAD[order(abs(`Activation z-score`), decreasing = T)][,.(`Upstream Regulator`, `Expr Log Ratio`, `p-value of overlap`, `Activation z-score`, `Target molecules in dataset`)]
Table3=Table3.Suppl[1:20]

#Freq histogram of top terms
TFactors=Table3[,`Upstream Regulator`]

#Top terms for each TF:
if(F){
    #specific
    TFCountTable=data.table()
    for (Factor in TFactors) {
        File = sprintf("Data/ToppFun/TFanalysis/181217_TF-GO-%s.txt", Factor)
        if (file.exists(File)){
            data=fread(File)
            data = data[Category=="GO: Biological Process" & order(`Hit Count in Query List`, decreasing = TRUE)][1:10, Name]
            TFCountTable=rbind(TFCountTable, data.table(rep(Factor, times=length(data)), data))
        } else {
            message(sprintf("No data found for: [%s]", Factor))
        }
    }
}
#generic, and all terms:
TFCountTable=data.table()
Files = list.files("Data/ToppFun/TFanalysis", ".txt")
for (File in Files) {
    Factor=gsub("181217_TF-GO-", "", File)
    Factor=gsub("\\.txt", "", Factor)
    data=fread(sprintf("Data/ToppFun/TFanalysis/%s", File))
    data = data[Category=="GO: Biological Process" & order(`Hit Count in Query List`, decreasing = TRUE)][, .(Name, `Hit Count in Query List`)]
    TFCountTable=rbind(TFCountTable, data.table(rep(Factor, nrow(data)), data))
}    
TFCountTable=TFCountTable[V1!="TP53(2)"]
TFCountTable.S = TFCountTable[order(`Hit Count in Query List`, decreasing = T),.SD[1:10], V1] #gets top 10 hits for each TF
TFCountTable.S2 = TFCountTable.S[,.N, Name][order(N, decreasing = T)][1:25] #sums top hits, gets top 25 of those
#ggplot(TFCountTable.S2, aes(x=reorder(Name, -N), y=N))+geom_col()+theme(axis.text.x = element_text(angle=90, vjust = 0.2, hjust = 1)) #this is across ALL TFs

#Overlap histogram of fig2 terms
Fig2TermTFs = TFCountTable[Name %in% unique(BioData[,Name]) & V1 %in% Table3[,`Upstream Regulator`]]
Fig2TermTFs = dcast(Fig2TermTFs, Name~., value.var = "V1", fun.aggregate = paste, collapse=",")
setnames(Fig2TermTFs, ".", "TFsAffectingPathway")
Fig2TermTFs[,NTFs := str_count(TFsAffectingPathway, ",")+1]
#figure 4 - leafcutter ====
tmsg("Processing leafcutter data (beta)...")
#MOCKUP, as leafcutter may be malfunctioning
Fig4A = ggplot(Fig4AData, aes(x=Comparison, y=N))+geom_col(width = 0.75)+ylab("Alternative splicing events")+
    theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5))
#FigS4AData = Fin[status!="Success",.N, .(Comparison,status)]

#volcano plots====
DESeq2Results.Spline = GrabResult(SplineModel, "ns.ImputedAgeDays..df...3.3", "SplineFull", 0.05, FALSE) 
DESeq2Results.Spline[,pValLog10:=-log10(pValue)]
Volcano1 = ggplot(DESeq2Results.Spline, aes(x=log2FC, y=pValLog10, shape=pValAdjust<=0.05, color=pValAdjust<=0.05))+
    geom_point()+theme(legend.position = "none")+
    labs(x=expression(paste("log"[2]," fold change")), y=expression(paste("log"[10]," adjusted p-value")))

DESeq2Results.Group = GetAllResultsForVariable(GroupModel, "Timepoint", 0.05, filterNonPass=FALSE)
DESeq2Results.Group[,pValLog10:=-log10(pValue)]
Volcano2 = ggplot(DESeq2Results.Group[Comparison %in% c("Timepoint|1M|1W", "Timepoint|3M|1M", "Timepoint|AD|3M", "Timepoint|AD|1W")], 
                  aes(x=log2FC, y=pValLog10, shape=pValAdjust<=0.05, color=pValAdjust<=0.05))+
    geom_point()+facet_wrap(~Comparison)+theme(legend.position="none")+
    labs(x=expression(paste("log"[2]," fold change")), y=expression(paste("log"[10]," adjusted p-value")))


#figure 8: LB development / GWPM C v HF Overlap====
LBGenes = DEGs.ERCC.1WAD[,.(geneID, log2FC, Comparison)]
GWPMGroupModel = readRDS("Data/RNAseq/Output/GroupGWPM.2excl.maxit5000.BatchGroup.Wald.RDS")
GWPM.CvHFModel.DEGs = GrabResult(GWPMGroupModel, c("Group", "HF", "C"), comparisonStr = "HF v C", Alpha = 0.05, Filter = TRUE)[,.(geneID, log2FC, Comparison)]

sprintf("LB's spline.df3 models returns %d DEGs. GM's C v HF model returns %d DEGs. The models share %d DEGs (%2.2f %% overlap w/ LB data, %2.2f %% overlap with GWPM data).", 
        nrow(LBGenes), nrow(GWPM.CvHFModel.DEGs), length(intersect(LBGenes[,geneID], GWPM.CvHFModel.DEGs[,geneID])), 
        ( length(intersect(LBGenes[,geneID], GWPM.CvHFModel.DEGs[,geneID])) / nrow(LBGenes) )*100,
        ( length(intersect(LBGenes[,geneID], GWPM.CvHFModel.DEGs[,geneID])) / nrow(GWPM.CvHFModel.DEGs) )*100
)

setkey(LBGenes, "geneID")
DEGsBothModels = GWPM.CvHFModel.DEGs[LBGenes][!is.na(log2FC) & !is.na(i.log2FC)]
sprintf("%d of the %d genes (%2.2f%%) that are significant in both models show significant differential expression opposite their direction in normal development in heart failure", 
        DEGsBothModels[sign(i.log2FC)!=sign(log2FC),.N], nrow(DEGsBothModels), (DEGsBothModels[sign(i.log2FC)!=sign(log2FC),.N] / nrow(DEGsBothModels))*100)

DEGsBothModels[,Code:=paste0(sign(log2FC), sign(i.log2FC))]
DEGsBothModels[,Code:=factor(Code, levels=c("11", "-1-1", "1-1", "-11"), labels=c("+ in HF and development ",  "- in HF and development", 
                                                                                  "+ in HF, - in development", "- in HF, + in development"))]
Fig3.8B = ggplot(DEGsBothModels[,.N,Code], aes(x=Code, y=N, fill=Code))+geom_col()+xlab("")+ylab("Number of shared DEGs")+theme(legend.position = "none", plot.margin = margin(0, 15, 0, 0))+
        scale_y_continuous(breaks=seq(0, 1000, 250), limits = c(0, 1000), expand = c(0,0))+coord_flip()
setkey(DEGsBothModels, geneID)
setkey(DictOrtholog, "Gene stable ID")
#ffwrite(DictOrtholog[DEGsBothModels], "DEGs.GWPMCvHF-LB1WvAD-HumanHomologsOne2One")

Fig3.8C = ggplot(DEGsBothModels, aes(x=log2FC, y=i.log2FC, color=Code, group=Code))+geom_point()+
    geom_smooth(method="lm", mapping = aes(group=NA), color="black", se = F)+theme(legend.position="none", plot.margin = margin(2,2,5,0))+
    labs(x=expression(paste("log"[2]," fold change in Development, 1W to Adult")), y=expression(paste("log"[2]," fold change in progression to HF")))+
    geom_hline(yintercept = 0)+geom_vline(xintercept = 0)
    
Output3.1 = lm(log2FC~i.log2FC, data=DEGsBothModels); summary(Output3.1)

#Figure 9 - those which do not recover to control levels after HF====
#Define "recover to control levels " - can't be any 1W 1M 3M AD levels as those are LEFT atrium, C HF R are RIGHT atrium, so it must be internal to GWPM's data

#Let's grab the counts for all genes that change when HF occurs:
GWPM.Counts = as.data.table(DESeq2::counts(GWPMGroupModel), keep.rownames = "geneID")[geneID %in% GWPM.CvHFModel.DEGs[,geneID]]
GWPM.Counts = melt(GWPM.Counts, id.vars="geneID")
GWPM.Counts[,group:= str_extract(variable, "C|R|HF")]
GWPM.Counts[,group:= factor(group, levels=c("C", "HF", "R"))]
#And calculate their average expression across our samples...
GWPM.Counts.Avgs = GWPM.Counts[,.(meanCounts=mean(value), sdCounts=sd(value)), .(group, geneID)]
rm(GWPM.Counts)
GWPM.Counts.Avgs[,upperBound := meanCounts + sdCounts]
GWPM.Counts.Avgs[,lowerBound := meanCounts - sdCounts]

GWPM.Counts.AvgsPerGene = dcast(GWPM.Counts.Avgs, geneID~group, value.var = "meanCounts")
#Using control as reference, of course, let's calculate percent change...
GWPM.Counts.AvgsPerGene[,CtoHF := ((HF-C)/C)*100]
GWPM.Counts.AvgsPerGene[,CtoR := ((R-C)/C)*100]
GWPM.Counts.AvgsPerGene2 = GWPM.Counts.AvgsPerGene[C > 50] #must have at least 50 counts in control as a lower bound

ggplot(GWPM.Counts.AvgsPerGene2, aes(x=CtoR))+geom_histogram(binwidth = 1)+xlab("% difference between control and recovery expression")+
        ylab("Number of DEGs within this percentage")+theme(plot.margin = margin(0,0,5,5))+ geom_density(aes(y=..density..*nrow(GWPM.Counts.AvgsPerGene)), size=1)
GWPM.Counts.AvgsPerGene2[,summary(CtoR)]

ggplot(GWPM.Counts.AvgsPerGene, aes(x=CtoR))+geom_histogram(binwidth = 1)+xlab("% difference between control and recovery expression")+
    ylab("Number of DEGs within this percentage")+theme(plot.margin = margin(0,0,5,5))+ geom_density(aes(y=..density..*nrow(GWPM.Counts.AvgsPerGene)), size=1)

GWPM.Counts.AvgsPerGene2[CtoHF>0,`Regulation in Heart Failure`:= "Upregulated in HF"]
GWPM.Counts.AvgsPerGene2[CtoHF<0,`Regulation in Heart Failure`:= "Downregulated in HF"]
ggplot(GWPM.Counts.AvgsPerGene2, aes(x=CtoR, color=`Regulation in Heart Failure`))+geom_histogram(binwidth = 1)+xlab("% difference between control and recovery expression")+
    ylab("Number of DEGs within this percentage")+theme(plot.margin = margin(0,0,5,5))+ geom_density(aes(y=..density..*nrow(GWPM.Counts.AvgsPerGene)), size=1)+
    scale_color_brewer(type = "qual")+facet_grid(`Regulation in Heart Failure`~.)+theme(legend.position = "none")+geom_vline(xintercept = 0, size=1)

GWPM.Counts.AvgsPerGene[,summary(CtoR)]

NWithin20Percent = GWPM.Counts.AvgsPerGene[CtoR < 20 & CtoR > -20, .N] 
sprintf("Of %d genes, %d (%f%%) return to within +/-20%% of mean control expression values", 
        nrow(GWPM.CvHFModel.DEGs), NWithin20Percent, round((NWithin20Percent/nrow(GWPM.CvHFModel.DEGs))*100, 2))


GWPM.Counts.AvgsPerGene[,mean(CtoR), CtoHF>0]

cor.test(GWPM.Counts.AvgsPerGene[,CtoR], GWPM.Counts.AvgsPerGene[,CtoHF])

#export====
tmsg("all tasks complete, removing all non-Figure objects")
#rm(list = (ls()[!grepl("Fig|Table", ls())]))

#Richards comparison----
source("Code/Desktop/functions.GO.r")
#Energy metabolism
GOData.Spline = fread("Data/PANTHER/181218_ERCC.ImputeB4Exclude.batchCutRINOX3nsImputedAgeDaysdf3.maxit5000.LRT.GOFullAndSlim.One2OneOrthosOnly.txt")
DictOrtholog.Human = unique(DictOrtholog[,.(`Gene stable ID`, `Human gene stable ID`, `Human homology type`)][`Human homology type`=="ortholog_one2one"])
Spline.HH = DictOrtholog.Human[DEGs.ERCC.SplineDF3]
setkey(GOData.Spline, "MappedIDs")
setkey(Spline.HH, "Human gene stable ID")
GOData.Spline = Spline.HH[GOData.Spline]

#using spline data
GOAnalysis.Spline = DoTheWholePANTHER(GOData.Spline, foldChangeTerm = "log2FC", geneTerm = "Gene stable ID", 
                                      colIDs = c("MF"="GO.MF.Complete", "CC"="GO.CC.Complete", "BP"="GO.BP.Complete"))
FattyData1 = GOAnalysis.Spline$Biological.DT[GO.BP.Complete == "fatty acid beta-oxidation"]
GenerateGeneSummaryPlots(FattyData1[,unique(`Gene stable ID`)])$Plot

FattyData2 = GOAnalysis.Spline$Biological.DT[GO.BP.Complete %like% "fatty acid meta"]
GenerateGeneSummaryPlots(FattyData2[,unique(`Gene stable ID`)])$Plot

Fig3.7 = GenerateGeneSummaryPlots(c("ENSOARG00000019427", "ENSOARG00000005581"))$Plot+facet_wrap(~geneName) 
#Figure3.7v2 = GenerateGeneSummaryPlots(Dictionary[grepl("^PPAR", `Gene name`), `Gene stable ID`], Dataset="ERCC+GWPM", Model = "1WAD+CvHF")$Plot
#PPAR alpha... totally falls off a cliff, when it should be constant or increase?? the other one also decreases. FML.
#The only literature source is in a VERY sketchy journal. 

Table3.8 = dcast(GOAnalysis.Spline$Biological.DT[grepl("fatty acid", GO.BP.Complete) & grepl("oxidation", GO.BP.Complete), .N, 
                                      .(IsUpregulated=log2FC>0, GO.BP.Complete)], GO.BP.Complete~IsUpregulated, value.var="N")

#GenerateGeneSummaryPlots(GOAnalysis.Spline$Biological.DT[grepl("fatty acid", GO.BP.Complete) & grepl("oxidation", GO.BP.Complete), unique(`Gene stable ID`)])$Plot

#using 1WAD model
GOData.1WAD = fread("Data/PANTHER/190206_ERCC.ImputeB4Exclude.batchCutRINOX3Timepoint.maxit5000.Wald1WAD.GOFullAndSlim.One2OneOrthosOnly.txt")
`1WAD.HH` = DictOrtholog.Human[DEGs.ERCC.1WAD]
setkey(GOData.1WAD, "MappedIDs")
setkey(`1WAD.HH`, "Human gene stable ID")
GOData.1WAD = `1WAD.HH`[GOData.1WAD]
GOAnalysis.1WAD = DoTheWholePANTHER(GOData.1WAD, foldChangeTerm = "log2FC", geneTerm = "Gene stable ID", 
                                    colIDs = c("MF"="GO.MF.Complete", "CC"="GO.CC.Complete", "BP"="GO.BP.Complete"))
GOAnalysis.1WAD$Biological.DT[GO.BP.Complete == "fatty acid oxidation"] #same 2 genes

dcast(GOAnalysis.1WAD$Biological.DT[GO.BP.Complete %like% "fatty acid", .N, .(IsUpregulated=log2FC>0, GO.BP.Complete)], GO.BP.Complete~IsUpregulated, value.var="N")

       
#Similarity to human studies====
set1=GenerateGeneSummaryPlots(c("GYS1", "GYS2", "NPPA", "MYH6", "MYH7", "GLUT1"="SLC2A1", "GLUT4"="SLC2A4", "PDK2", "PDK4", "ACADM", "CS"), Dataset="ERCC+GWPM", Model = "1WAD+CvHF", forceUserNames = TRUE)$Plot
set2=GenerateGeneSummaryPlots("ACTA1", Dataset="ERCC+GWPM", Model = "1WAD+CvHF")$Plot
set3=GenerateGeneSummaryPlots(c("MYH6"="ENSOARG00000002058"), Dataset="ERCC+GWPM", Model = "1WAD+CvHF")$Plot

GenerateGeneSummaryPlots(Dictionary[grepl("^ACTA\\d{1,2}", `Gene name`), `Gene stable ID`], Dataset="ERCC+GWPM", Model = "1WAD+CvHF")$Plot

GenerateGeneSummaryPlots("FNDC1", Dataset="ERCC+GWPM", Model = "1WAD+CvHF")$Plot

set4=GenerateGeneSummaryPlots(c("NPPA", "NPPB"), Dataset="ERCC+GWPM", Model = "1WAD+CvHF")$Plot
set5=GenerateGeneSummaryPlots(Dictionary[grepl("^TNNT\\d{1,2}", `Gene name`), `Gene stable ID`], Dataset="ERCC+GWPM", Model = "1WAD+CvHF")$Plot
set6=GenerateGeneSummaryPlots(Dictionary[grepl("^TNNI\\d{1,2}", `Gene name`), `Gene stable ID`], Dataset="ERCC+GWPM", Model = "1WAD+CvHF")$Plot
