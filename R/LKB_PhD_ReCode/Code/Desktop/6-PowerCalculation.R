#base coverage = (number of reads mapped to exons * average read length) / total length of all exons
#depth of sequencing = (total number of reads * average read length) / total length of all exons
library(data.table); library(RNASeqPower); library(DESeq2)

Exons = fread("Data/biomart/parsedOvAr3.1GTF.tsv")
Exons = Exons[feature=="exon"]
TotalExonLength = Exons[,sum(end-start)]

MappedCountFiles = list.files("Data/RNAseq/Input/CountFiles.ERCC/", "*.count", full.names = T)
ExonMappedReads = data.table()
for (File in MappedCountFiles) {
    temp = fread(File)
    ExonMappedReads = rbind(ExonMappedReads, temp[!grepl("^__|^ERCC-", V1), .(FileName=File, nExonReads=sum(V2))])
}
ExonMappedReads[,FileName:=gsub("Data/RNAseq/Input/CountFiles.ERCC/", "", FileName)]
ExonMappedReads[,FileName:=gsub("_OAv3\\.1\\+ERCC_aligned_ENSEMBL_Union_RevStranded\\.counts", "", FileName)]
ExonMappedReads[,FileName:=gsub("Sample", "Sample-", FileName)]
ExonMappedReads[,TotalExonLength:=TotalExonLength]
ExonMappedReads[,AverageReadLength:= 75]

ExonMappedReads[,EstimatedCoverage:= (nExonReads * AverageReadLength) / TotalExonLength] #656??????? IMPLAUSIBLE
ExonMappedReads[,mean(EstimatedCoverage)] #21.87154...? i guess??

#The dispersion is essentially the squared coefficient of variation: 
#Var = mu + dispersion * mu^2 
#Var / mu^2 = 1/mu + dispersion
#CV^2 = 1/mu + dispersion
#When the counts are large, e.g. 100, then 1/100 is small compared to a typical dispersion value, so you have CV ~= sqrt(dispersion). 
#The final dispersion values can be used: dispersions(dds)
GroupModel=readRDS("Data/RNAseq/Output/GroupModel.ERCC.ImputeB4Exclude.2excl.maxit5000.batchCutRINOXTimepoint.Wald.RDS")
a = sqrt(DESeq2::dispersions(GroupModel)) #Gene-Wise dispersion

Groups = c("1W", "AD")
GroupLen = c()
for (GRP in Groups) { GroupLen = c(GroupLen, GRP=ExonMappedReads[FileName %like% GRP, .N]) }
names(GroupLen)=Groups

RNASeqPower::rnapower(depth = ExonMappedReads[,mean(EstimatedCoverage)], n=GroupLen[1], n2=GroupLen[2], cv=median(a), effect=seq(0,5,0.5), 
                      alpha = 0.05)
