if (!require(GO.db)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite("GO.db")
  require(GO.db)
}
require(data.table)

MFROOT <- "GO:0003674"
BPROOT <- "GO:0008150"
CCROOT <- "GO:0005575"

BuildGORankTable <- function(nodeMode){
  if(!(nodeMode %in% c("BP", "CC", "MF"))){ stop("mode must be one of BP, CC or MF") }
  NodeTable = NULL;
  switch(nodeMode,
    "BP"={ ChildTable = GOBPCHILDREN; NodeTable = data.table(Node=GOBPOFFSPRING$"GO:0008150"); GOStack = c(BPROOT)},
    "CC"={ ChildTable = GOCCCHILDREN; NodeTable = data.table(Node=GOCCOFFSPRING$"GO:0005575"); GOStack = c(CCROOT)},
    "MF"={ ChildTable = GOMFCHILDREN; NodeTable = data.table(Node=GOMFOFFSPRING$"GO:0003674"); GOStack = c(MFROOT)}
  )
  NodeTable[,Rank:=Inf]
  NodeTable = rbind(NodeTable, data.table(Node=GOStack[1], Rank=Inf))
  setkey(NodeTable, Node)
  CurrentRank = -1
  while(length(GOStack)!=0){
    CurrentRank = CurrentRank + 1

    GOStack.Old = intersect(GOStack, NodeTable[!is.infinite(Rank), Node]) #Children already ranked 
    NodeTable[Node %in% GOStack.Old, Rank:=min(CurrentRank, Rank), Node]  #Are given smallest possible rank
    #Must be done by=Node, else you'll get an accidental bug where min is a GLOBAL min of all members of GOStack.Old!
    
    GOStack.New = setdiff(GOStack, NodeTable[!is.infinite(Rank), Node]) #Children not yet ranked
    NodeTable[Node %in% GOStack.New, Rank:=CurrentRank]  #All current children are given the current rank; no need to do by=Node

    GOStack = unlist(lapply(GOStack.New, function (x) { get(x, ChildTable) })) #[which(names(get(x, ChildTable)) == "is_a")]
    #new stack consists of all the children~~, where the relationship is "is_a"~~
    #lapplying to GOStack uses unnecessary cycles to revisit nodes that can't get a lower rank, and leads to a repetition bug if not by=Node is used.
  }
  NodeTable
}

NodeMode = c("BP", "CC", "MF")[Sys.getenv("SGE_TASK_ID", 0)]
Node = BuildGORankTable(NodeMode)
saveRDS(Node, paste0(format(Sys.time(), "%y%m%d-%H%M_"), "GOTravellingSalesmanRank.", NodeMode, ".RDS"))
