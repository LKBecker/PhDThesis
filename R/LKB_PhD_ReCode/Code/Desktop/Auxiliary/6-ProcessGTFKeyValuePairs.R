library(data.table)
exons_file = fread("Data/biomart/Ovis_aries.Oar_v3.1.95.gtf")

#exons_file = exons_file[V3 == "exon"]

AllHeadings = unlist(strsplit(exons_file[,V9], "; "))
AllHeadings = unique(sapply(strsplit(AllHeadings," "), `[`, 1))

sapply(AllHeadings, function(x) { exons_file[, (x) := str_match(V9, sprintf('%s "(.*?)";', x))[,2]] } )
#yeah this is a lot faster than a for loop, who'd have thought

exons_file[,V9:=NULL]
colnames(exons_file)=c("chr", "source", "feature", "start", "end", "score", "strand", "frame", AllHeadings)

fwrite(exons_file, "Data/biomart/parsedOvAr3.1GTF.tsv", row.names = F)
