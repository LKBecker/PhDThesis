#DESeq2 (main.R, )====
CONF = list()
CONF[["VERSION"]]		        = "1.3-ForFigures"
CONF[["FILE_EXCLUDE_LIST"]]	    = c("161206-07_3M", "161207-05_3M")
CONF[["MaxCores"]] 				= 4
CONF[["MIN_FDR"]]				= 0.05
CONF[["MIN_LFC"]]				= 0.0
CONF[["MIN_GENES_WITH_MIN_LFC"]]= 0025
CONF[["MAX_WALD_ITERATIONS"]]   = 20000 #do not set to 999999, it takes four fucking hours and still doesn't converge

#Expression Patterns (Graphs.General.R)====
CONF[["USE_LEGACY_ALGO"]]       = T
CONF[["EXPRESSED_PCT_LIMIT"]] 	= 0.50  #GTex uses 20%, but at 7-8 samples, 50% is 4/7
CONF[["EXPRESSED_MIN_COUNTS"]]  = 6    #GTex uses 6
CONF[["NORMALISE_COUNTS"]] 		= TRUE
CONF[["TPM_LIMIT"]]             = 0.1   #GTex 0.1

#GGplot2 shenanigans====
#CONF[["GG_DEFAULT_DPI"]]= 72
CONF[["GG_TARGET_DPI"]]         = 96
#CONF_GG_SCALE_FACTOR= CONF_GG_TARGET_DPI / CONF_GG_DEFAULT_DPI

CONF[["COR_CUTOFF"]]         = 0.80