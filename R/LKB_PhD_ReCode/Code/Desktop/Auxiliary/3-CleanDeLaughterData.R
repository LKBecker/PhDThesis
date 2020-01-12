library(data.table); library(stringr)
source("Code/Desktop/functions.r")
rm(SplitDataTableWithMultiRows, AskYN, completeDT, DoHeatmap, GenerateGeneSummaryPlots, give.n, LeafCutterContainerToResultsObj, se, tmsg, tstamp)
#Old way - clean FPKM matrix====
# source("Code/Desktop/functions.r")
# 
DeLaughterData = fread("Data/MuSiC/DeLaughter/DevCell_FPKMandMeta.csv")                         # Load whole file
#DeLaughterData = DeLaughterData[, c(1, which(DeLaughterData[1,]=="myocardial")), with=FALSE]    # with=FALSE deactivates "expression mode" and causes the containing numbers to become column IDs
DeLaughterData = DeLaughterData[, c(1, which(DeLaughterData[3,]=="p21")), with=FALSE]    # with=FALSE deactivates "expression mode" and causes the containing numbers to become column IDs
DeLaughterMeta = transpose(DeLaughterData[1:6, 1:ncol(DeLaughterData)])                         # first six rows are metadata, grab those
colnames(DeLaughterMeta)=unlist(DeLaughterMeta[1,])                                             # first row of big data matrix is actually colnames
DeLaughterMeta = cbind(DeLaughterMeta, data.table("Library_Name"=colnames(DeLaughterData)))     # matrix headers are data too (library names)
DeLaughterMeta = DeLaughterMeta[Cell_Type!="Cell_Type"]                                         # now that cols are named, take out row 1

DeLaughterData = DeLaughterData[7:nrow(DeLaughterData)]                                         # remove first six rows (metadata) and first column which have been extracted
#DeLaughterGeneNames = unlist(DeLaughterData[,1])
#DeLaughterData[,Library_Name:=NULL]
#DeLaughterMatrix = as.matrix(DeLaughterData)

#ffwrite(DeLaughterMeta, "DeLaughter.Metadata", Folder="Data/MuSiC/DeLaughter/")

#New way - DIY from FASTQ files====
RNAFiles = fread("Data/MuSiC/DeLaughter/DeLaughterFastqFiles-2.txt", header = FALSE)
setnames(RNAFiles, "V1", "FileName")
RNAFiles[,Library_Name:=str_match(FileName, "(LIB008986-TRA00020\\d{3})_S1")[,2]]
RNAFiles[,Library_Name:=gsub(x = Library_Name, pattern = "-", replacement = "_")]
RNAFiles[,Orientation:=str_match(FileName, "LIB008986-TRA00020\\d{3}_S1_L001_(R1|R2)")[,2]]

DeLaughterMeta = DeLaughterMeta[Timepoint == "p21"]
RNAFiles.Incl = RNAFiles[Library_Name %in% DeLaughterMeta[,Library_Name]]
RNAFiles.Incl = merge(RNAFiles.Incl, DeLaughterMeta, by="Library_Name", all.x=T)
#RNAFiles.Incl = RNAFiles.Incl[Cell_Type=="myocardial"] #FILTERING HERE
RNAFiles.Incl[,FileName:=gsub("_R(1|2)\\.fastq\\.gz", "", FileName)]
RNAFiles.Incl[,Library_Name:=gsub(x = Library_Name, pattern = "_", replacement = ".")]

RNAFiles.Incl.CM = unique(RNAFiles.Incl[Cell_Type=="myocardial",.(Library_Name, FileName, Cell_Name, Cell_Type)])
RNAFiles.Incl.CM[,Cell_Number:=as.numeric(str_extract(Cell_Name, "\\d{1,2}$"))]
RNAFiles.Incl.CM[,NewFileName:=paste0("p21_", sprintf("S%d_", Cell_Number), "L001_"), by=Cell_Number]

RNAFiles.Incl.FB = unique(RNAFiles.Incl[Cell_Type=="fibroblast",.(Library_Name, FileName, Cell_Name, Cell_Type)])
RNAFiles.Incl.FB[,Cell_Number:=as.numeric(str_extract(Cell_Name, "\\d{1,2}$"))]
RNAFiles.Incl.FB[,NewFileName:=paste0("p21_", sprintf("S%d_", Cell_Number), "L001_"), by=Cell_Number]

RNAFiles.Incl.EC = unique(RNAFiles.Incl[Cell_Type=="endothelial",.(Library_Name, FileName, Cell_Name, Cell_Type)])
RNAFiles.Incl.EC[,Cell_Number:=as.numeric(str_extract(Cell_Name, "\\d{1,2}$"))]
RNAFiles.Incl.EC[,NewFileName:=paste0("p21_", sprintf("S%d_", Cell_Number), "L001_"), by=Cell_Number]

#RNAFiles.Excl = RNAFiles[!(Library_Name %in% DeLaughterMeta[,Library_Name])] # 26 files (13 cells) which aren't characterised.
rm(RNAFiles, RNAFiles.Incl)

#In lieu of replicating the entirety of cellranger, we used STAR for read counting...
DeLaughter.STARCountFiles = data.table(FileName=list.files("Data/MuSiC/DeLaughter/DeLaughter.v1/", ".tab", full.names = T))
DeLaughter.STARCountFiles[,Library_Name:=gsub("-", "_", str_extract(FileName, "LIB0089\\d{2}-TRA00020\\d{3}"))]

DeLaughterHuegMatrix = data.table()
for (i in 1:nrow(DeLaughter.STARCountFiles)) {
    tempdata = fread(DeLaughter.STARCountFiles[[i, 1]])
    tempdata[,Library:=DeLaughter.STARCountFiles[[i, 2]]]
    DeLaughterHuegMatrix = rbind(DeLaughterHuegMatrix, tempdata, fill=T)
}
DeLaughterHuegMatrix[,V3:=NULL]
DeLaughterHuegMatrix[,V4:=NULL]
colnames(DeLaughterHuegMatrix) = c("gene", "counts", "library")
DeLaughterHuegMatrix = DeLaughterHuegMatrix[which(substr(DeLaughterHuegMatrix[,gene],1,2)!="N_")]
DeLaughterHuegMatrix = dcast(DeLaughterHuegMatrix, gene~library, value.var = "counts")
DeLRows=DeLaughterHuegMatrix[,gene]
DeLaughterHuegMatrix[,gene:=NULL]
DeLaughterHuegMatrix = as.matrix(DeLaughterHuegMatrix)
rownames(DeLaughterHuegMatrix)=DeLRows
rm(DeLRows, DeLaughter.STARCountFiles, i, tempdata)

saveRDS(DeLaughterHuegMatrix, "Data/MuSiC/DeLaughter/DevCell_LKB-STAR_p21-CM-FB-EC.rds")
#Comparison: ENSMUSG00000035885 or Cox8a

Comp.LB = DeLaughterHuegMatrix[which(rownames(DeLaughterHuegMatrix)=="ENSMUSG00000035885"),]
Comp.Orig = DeLaughterData[Library_Name=="Cox8a"]
Comp.Orig[,Library_Name:=NULL]
cor.test(unname(Comp.LB), as.numeric(unname(Comp.Orig)))
