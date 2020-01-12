#setwd("H:/Data/PhD/Data/results/RNA QC/Nanodrop/")
library("data.table")
library("ggplot2")
GraphData <- fread("./2016-11-22 - ErroneousMeasurementsCut.ndv", header = T, skip = 4)
set(GraphData, j = GraphData[, which(colnames(GraphData) %in% c("User ID", "Measurement Type", "Serial #", "Config.", "Date", "Time"))], value=NULL)
GraphData.M <- melt.data.table(GraphData, id.vars = "Sample ID")
GraphData.M[`Sample ID` %in% c("161121-01","161121-03"), SampleType:="TRIzol + Column"]
GraphData.M[`Sample ID` %in% c("161121-02","161121-04"), SampleType:="TRIzol"]
GraphData.M[`Sample ID` %in% c("161121-01","161121-02"), Tissue:="Sample 1"]
GraphData.M[`Sample ID` %in% c("161121-03","161121-04"), Tissue:="Sample 2"]

#Shows that column treatment reduces phenol contamination (peaks at 230)
GraphData.ND <- GraphData.M[variable %in% seq(220, 350)]
ggplot(GraphData.ND, aes(x=variable, y=value, color=`Tissue`, group=`Sample ID`))+geom_line(size=1.1)+facet_grid(.~SampleType)+
	scale_x_discrete(breaks=seq(220, 350, 10))+theme(panel.spacing = unit(2, "lines"))+xlab("Wavelength (nm)")+ylab("10mm Absorbance")

#this one is kinda crap
GraphData.ND2 <- GraphData.M[variable %like% "260/"]
ggplot(GraphData.ND2, aes(x=SampleType, y=value, color=`Tissue`))+geom_crossbar(aes(ymin=value, ymax=value))+
	facet_grid(variable~., scales="free_x")+ylab("Absorbance ratio")+xlab("Extraction protocol")

RINdata <- fread("../RINData.txt", na.strings = c("N/A", "#N/A", ""))
RINdata[,Date:=as.Date.character(Date, "%d/%m/%Y")]
RINdata[,UsesBeatBeating:=(Date>as.Date.character("21/11/2016", "%d/%m/%Y"))]
RINdata[!(is.na(RIN)),.(.N, mean(RIN)), UsesBeatBeating]

RINdata.Rep <- RINdata[,.(Date, Method, Result_ID, RIN, `260/230`, `260/280`, `ng/ul`)]
#write.table(RINdata.Rep, "clipboard", sep = "\t", row.names = F, quote = F)

# RINdata.M <- melt.data.table(RINdata.Rep, id.vars = c("Date", "Method", "Result_ID"))
# ggplot(RINdata.M[variable=="ng/ul"], aes(x=Method, y=value))+geom_point()

#Another Nanodrop Graph
GraphData2 <- fread("./2016-12-06_Batch1.ndv", header = T, skip = 4)
set(GraphData2, j = GraphData2[, which(colnames(GraphData2) %in% c("User ID", "Measurement Type", "Serial #", "Config.", "Date", "Time"))], value=NULL)
GraphData.M2 <- melt.data.table(GraphData2, id.vars = "Sample ID")
GraphData.ND2 <- GraphData.M2[variable %in% seq(220, 350)]
ggplot(GraphData.ND2, aes(x=variable, y=value, color=`Sample ID`, group=`Sample ID`))+geom_line(size=1.1)+
	scale_x_discrete(breaks=seq(220, 350, 10))+theme(panel.spacing = unit(2, "lines"))+xlab("Wavelength (nm)")+ylab("10mm Absorbance")
