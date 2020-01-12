#Libraries----
library(data.table)
#library(stringr)
if (!require(fgsea)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("fgsea")
  if (!require("Rcpp")) { install.packages("Rcpp"); require("Rcpp") }
  if (!require("tibble")) { install.packages("tibble"); require("tibble") }
  library(fgsea)
}
#script start----
setwd("/mnt/iusers01/bk01-icvs/mqbprlb2/PhD/tools/jobgen/jobs/LKB_RNA_Analysis/R/FGSEAv1/attempt5.ltcc")
CONF_MIN_ARGS = 2

StartTime = Sys.getenv("StartTime", unset=format(Sys.time(), "%y%m%d_%H%M_"))
CONF_N_CORES = as.numeric(Sys.getenv("NSLOTS", unset=1))

args <- commandArgs(TRUE)
#args = c("181018-113415_DESeq2.2rm.batchCutRinOx3nsImputedAgeDaysDF3.tsv", "LTCC_Custom.gmt", "1e+7", "geneID")
if (length(args) < CONF_MIN_ARGS) { stop(paste0("Fatal error: Script requires at least ", CONF_MIN_ARGS, " arguments to run.")) }

File = args[1] #file containing t-stat, gained NOT from DESeq2::results() but directly out of the DESeq2 object
if(!(class(File) == "character")) { stop(paste0("Error: ", File, "is not a valid string for a filename")) }
if(!(file.exists(File))) { stop(paste0("Error: ", File, "does not exist in directory '", getwd(), "'")) }
tstats = fread(File)

#setting up variables----
PathwayFile = as.character(args[2])
pathways_gmt 	 = gmtPathways(PathwayFile)

CONF_N_REPS = as.numeric(args[3])
if(is.na(CONF_N_REPS)) { stop(paste0("Error: ", args[3], " cannot be converted to a number and cannot be used as number of iterations.")) }

CONF_ID_Arg = as.character(args[4])
if(is.na(CONF_ID_Arg)) { CONF_ID_Arg = "Gene name" }

#More than one geneID may map to a name; thus, we select the one with the largest absolute t_stat as the relevant one.
#Assumes input has Gene name, t_stat
if (!("t_stat" %in% colnames(tstats)) & !("stat" %in% colnames(tstats))) { stop("Input file must have \"t_stat\" or \"stat\" column" ) }
if ("t_stat" %in% colnames(tstats)) { setnames(tstats, "t_stat", "stat") }
if (!(CONF_ID_Arg %in% colnames(tstats))) { stop(paste0("Input file must have '", CONF_ID_Arg, "' column")) }

#DESeq2 'raw' output
setorder(tstats, stat)
stats <- tstats[,stat] #must be ranked for fGSEA; compare data(exampleRanks)

names(stats) <- tstats[, get(CONF_ID_Arg)]

res <- fgsea(pathways_gmt, stats, nperm=CONF_N_REPS, nproc=CONF_N_CORES) #~30 minutes at 1e6, do not run off the cluster at 1e9 unless you enjoy your system slowed to a crawl due to 100% CPU utilisation for about 16 hours
setorder(res, padj)

OutFolder <- paste0(getwd(), "/", StartTime, "FGSEA")
if(!dir.exists(OutFolder)) { dir.create(OutFolder) }

FileName <- paste0(OutFolder, "/FGSEA-", CONF_N_REPS, "-", PathwayFile, "-Results")
saveRDS(res, paste0(FileName, ".RDS"))
message("Saved to RDS.")
