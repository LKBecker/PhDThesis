library(data.table)
library(ggplot2)
library(stringr)
library(cowplot)
source("Code/Desktop/functions.DESeq2.r")

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
#functions====
RunEntireQC <- function(Path){
    #Path = "Data/FastQC/RawData/"
    SeqQCFiles <- list.files(Path, pattern = "*.txt", full.names = T)
    SeqQCData							  <- list()
    SeqQCData$Summary        <- data.table()
    SeqQCData$Adapters       <- data.table()
    SeqQCData$BasicStats     <- data.table()
    SeqQCData$KmerContent    <- data.table()
    SeqQCData$OverrepSeqs    <- data.table()
    SeqQCData$PerBaseNContent<- data.table()
    SeqQCData$PerBaseSeqCont <- data.table()
    SeqQCData$PerBaseSeqQual <- data.table()
    SeqQCData$PerSeqGCCont   <- data.table()
    SeqQCData$PerSeqQualScore<- data.table()
    SeqQCData$PerTileSeqQual <- data.table()
    SeqQCData$SeqDuplicLevels<- data.table()
    SeqQCData$SeqLengths     <- data.table()
    
    for(file in SeqQCFiles){
        FileType = str_match(file, "WTCHG_\\d+_\\d+_\\d_fastqc_(.*)\\.txt")[,2]
        FileName = str_match(file, "(WTCHG_\\d+_\\d+_\\d)_fastqc_.*\\.txt")[,2]
        tryCatch(expr = { QCData = fread(file)}, error = function (e) { QCData = data.table()  })
        QCData[,Dataset:=FileName]
        switch(FileType,
               `Summary`                       ={ SeqQCData$Summary = rbind(SeqQCData$Summary, QCData, fill=T) },
               `Adapter Content`               ={ SeqQCData$Adapters = rbind(SeqQCData$Adapters, QCData, fill=T) },
               `Basic Statistics`              ={ SeqQCData$BasicStats = rbind(SeqQCData$BasicStats, QCData, fill=T) },
               `Kmer Content`                  ={ SeqQCData$KmerContent = rbind(SeqQCData$KmerContent, QCData, fill=T) },
               `Overrepresented sequences`     ={ SeqQCData$OverrepSeqs = rbind(SeqQCData$OverrepSeqs, QCData, fill=T) },
               `Per base N content`            ={ SeqQCData$PerBaseNContent = rbind(SeqQCData$PerBaseNContent, QCData, fill=T) },
               `Per base sequence content`     ={ SeqQCData$PerBaseSeqCont = rbind(SeqQCData$PerBaseSeqCont, QCData, fill=T) },
               `Per base sequence quality`     ={ SeqQCData$PerBaseSeqQual = rbind(SeqQCData$PerBaseSeqQual, QCData, fill=T) },
               `Per sequence GC content`       ={ SeqQCData$PerSeqGCCont = rbind(SeqQCData$PerSeqGCCont, QCData, fill=T) },
               `Per sequence quality scores`   ={ SeqQCData$PerSeqQualScore = rbind(SeqQCData$PerSeqQualScore, QCData, fill=T) },
               `Per tile sequence quality`     ={ SeqQCData$PerTileSeqQual = rbind(SeqQCData$PerTileSeqQual, QCData, fill=T) },
               `Sequence Duplication Levels`   ={ SeqQCData$SeqDuplicLevels = rbind(SeqQCData$SeqDuplicLevels, QCData, fill=T) },
               `Sequence Length Distribution`  ={ SeqQCData$SeqLengths = rbind(SeqQCData$SeqLengths, QCData, fill=T) }
        )
    } #for file in SeqQCFiles
    return(SeqQCData)
}
GraphCountStats <- function(Directory, htseqMode, GTFSource, RefGenomeVersion, description){
    if(!substr(Directory, str_length(Directory), str_length(Directory)) == "/") { Directory = paste0(Directory, "/") }
    FileNames2 				<- dir(path = Directory, pattern = "*.counts")
    DT <- data.table()
    i <- 1
    for (SampleName in FileNames2) {
        DT2 <- fread(paste0(Directory, SampleName))
        DT2[,ID:=i]
        DT2[,Sample:=str_match(SampleName, "Sample-(\\d{6}-\\d{2}_\\w{2})_")[,2]]
        DT2[,Group:=str_match(SampleName, "Sample-\\d{6}-\\d{2}_(\\w{2})_")[,2]]
        i<-i+1
        DT <- rbind(DT, DT2)
    }
    rm(DT2, i, SampleName)
    colnames(DT)<-c("Type", "Count", "ID", "Sample", "Group")
    Totals 	<- DT[!(Type %like% "__"),.(Count=sum(Count)), .(ID, Sample, Group)]
    Totals[,Type:="aligned"]
    Losses 	<- DT[Type %like% "__"]
    
    Stats <- rbind(Totals, Losses)
    Stats[,Type := factor(Type, levels = c("__alignment_not_unique", "__ambiguous", "__no_feature", "__not_aligned", "__too_low_aQual", "aligned"),
                          labels = c("Alignment not unique", "Ambiguous alignment", "No GTF feature", "Not aligned", "Quality too low", "Successfully aligned"))]
    Stats[,Description:=description]
    Stats[,RefGenome:=RefGenomeVersion]
    Stats[,GTF:=GTFSource]
    return(Stats)
}
#FastQC====
# SeqQCData = RunEntireQC("Data/FastQC/RawData/")
# Summary        <- SeqQCData$Summary
# Adapters       <- SeqQCData$Adapters
# BasicStats     <- SeqQCData$BasicStats

#HISAT2 Graph====
#load LAZYHISAT.py results
SampleReadData = fread("../../../Data/results/LKB_RNA/180725_HisatRunOvAr3.1.txt")
SampleReadData[,group:=factor(str_extract(Sample, "1W|(1|3)M|AD"), levels=c("1W", "1M", "3M", "AD"))]
SampleReadData[,Sample:=factor(Sample, levels=SampleReadData[order(group), Sample])]

#SampleReadData[,OverallAlignTest := ((AlignedConcord1Time + AlignedConcordManyTime + Concord0Discord1 + (NAP.1AlignMates/2) + (NAP.ManyAlignMates/2)) / Reads) * 100]

SampleReadData = SampleReadData[,.(Reads, Sample, group, AlignedConcord1Time, AlignedConcordManyTime, Concord0Discord1,
                                   NAP.NAMates=NAP.NAMates/2, NAP.1AlignMates=NAP.1AlignMates/2, NAP.ManyAlignMates=NAP.ManyAlignMates/2)]
setnames(SampleReadData, c("AlignedConcord1Time", "AlignedConcordManyTime", "Concord0Discord1", "NAP.NAMates", "NAP.1AlignMates", "NAP.ManyAlignMates"),
         c("Once, concordantly", "Many, concordantly", "Once, discordantly",
           "Never", "Once, outside of mate", "Many, outside of mate"))
SampleReadData.M = melt(SampleReadData, id.vars = c("Sample", "group"))
setnames(SampleReadData.M, "variable", "Alignment")

Fig2.2 = ggplot(SampleReadData.M[Alignment!="Reads"], aes(x=Sample, y=value, fill=Alignment))+geom_col()+ylab("Reads")+
    theme(axis.text.x = element_text(angle=90, vjust = 0.2), legend.position = "bottom", axis.title.x = element_text(margin = margin(5,5,-10,5)),
          legend.margin = margin(1,1,20,1))+
    scale_fill_discrete(guide = guide_legend(nrow = 2, keywidth=0.8, keyheight=0.8, title.hjust = 1, 
                                             label.theme = element_text(margin = margin(1,1,5,1), size = geom.text.size)))

sprintf("An average of %0.2f%% of reads were successfully feature counted", (SampleReadData.M[,sum(value)]/SampleReadData.M[,sum(value)])*100)

unique(SampleReadData.M[,.SD[Alignment=="Once, concordantly", value]/.SD[Alignment=="Reads", value], Sample])[,mean(V1)]

#htseq-count graph====
OAv3.1 = GraphCountStats("../../../Data/results/LKB_RNA/Counts.RTFM.Attempt4/", "union", "ENSEMBL", "3.1", "ENSEMBL Ovis Aries v3.1")
OAv3.1[,Group:=factor(Group, levels=c("1W", "1M", "3M", "AD"))]
OAv3.1[,Sample:=factor(Sample, levels = unique(OAv3.1[order(Group), Sample]))]
OAv3.1[Type=="Successfully aligned", Feature:="Successful"]
OAv3.1[Type=="Ambiguous alignment", Feature:="Ambiguous alignment"]
OAv3.1[Type=="No GTF feature", Feature:="No GTF Feature"]
OAv3.1[Type=="Not aligned", Feature:="Failed"]
OAv3.1[Type=="Quality too low", Feature:="Quality Too Low"]
OAv3.1[Type=="Alignment not unique", Feature:="Alignment not Unique"]
OAv3.1[,Feature:=factor(Feature, levels=c("Alignment not Unique", "Ambiguous alignment", "No GTF Feature", "Failed", "Quality Too Low", "Successful"))]

sprintf("An average of %0.2f%% of reads were successfully feature counted", (OAv3.1[,sum(Count)]/OAv3.1[,sum(Count)])*100)

unique(OAv3.1[,.(.SD[Feature=="Successful", Count]/sum(Count), Group), Sample])[order(Group)]

Fig2.3 = ggplot(OAv3.1, aes(x=Sample, y=Count, fill=Feature, order=Group))+geom_col()+
    theme(axis.text.x = element_text(angle=90, vjust = 0.2, margin = margin(8, 1, 0, 1)), legend.position = "bottom", plot.margin = margin(5,5,5,0),
          axis.title.x = element_text(margin = margin(5,5,-10,5)))+ylab("Counts")+
    scale_x_discrete(expand = c(0,0))+scale_y_continuous(expand=c(0,0))+
    scale_fill_discrete(guide = guide_legend(nrow = 1, keywidth=0.8, keyheight=0.8, title.hjust = 1, 
                                             label.theme = element_text(margin = margin(1,5,1,5), size = geom.text.size)))

Fig2.3

#Volcano plots====
SplineModel =readRDS("Data/RNAseq/Output/SplineModel.ERCC.ImputeB4Exclude.2excl.maxit5000.batchCutRINOXnsImputedAgeDaysdf3.LRT.RDS")
DESeq2Results.Spline = GrabResult(SplineModel, "ns.ImputedAgeDays..df...3.3", "LRT", 0.05, FALSE)
DESeq2Results.Spline[,pValNergLog10:=-log10(pValAdjust)]
DESeq2Results.Spline[pValAdjust<0.05,Signif:=pValNergLog10]
DESeq2Results.Spline[pValAdjust>=0.05,notSignif:=pValNergLog10]

VolcanoSpline = ggplot(DESeq2Results.Spline, aes(x=log2FC, y= -log10(pValAdjust), shape=pValAdjust<=0.05, color=pValAdjust<=0.05))+
    geom_point()+theme(legend.position = "none")+geom_vline(xintercept = c(-2, 2))+geom_hline(yintercept = -log(0.05, 10))+
    labs(x=expression(paste("log"[2]," fold change")), y=expression(paste("-log"[10]," adjusted p-value")))
    
VolcanoSpline

