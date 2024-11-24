library(SINCERA)
library(cluster)

dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")
source("RIO.R")

input <- function(inputfile) {
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
    pfix = prefix()
  if (length(pfix) != 0) {
     pfix <<- paste(pfix, "/", sep="")
  }
}

run <- function() {}

output <- function(outputfile) {
sc <- readRDS(paste(pfix, parameters["siggenes", 2], sep="/"))
	# get the signature genes as a data frame
# Use those genes for functional enrichment
siggenes <- getSigGenes(sc)

groups <- readSequential(paste(pfix, parameters["groups", 2], sep="/"))
types <- readSequential(paste(pfix, parameters["type", 2], sep="/"))
markertypes <- readSequential(paste(pfix, parameters["markertype", 2], sep="/"))
markersymbol <- readSequential(paste(pfix, parameters["markersymbol", 2], sep="/"))


# map cell clusters to cell types
sc <- setCellType(sc, do.reset = T)
#sc <- setCellType(sc, groups=c("1", "2", "3","5"), types=c("AT2", "Basal", "Indeterminate","Goblet"))
sc <- setCellType(sc, groups=groups, types=types)


# add cell type marker information to sincera
markers <- data.frame(TYPE=markertypes,
                      SYMBOL=markersymbol)
#markers <- data.frame(TYPE=c(rep("AT2", 4), rep("Basal", 3), rep("Goblet", 4)),
#                      SYMBOL=c("SFTPB","ABCA3","SLC34A2","LPCAT1","KRT5","KRT14","TP63","SPDEF","MUC5AC","MUC5B","SCGB3A2"))
sc <- setCellTypeMarkers(sc, markers)

head(getCellTypeMarkers(sc))

# perform rank-aggregation based cell type validation using marker expression
sc <- celltype.validation(sc)

saveRDS(sc, outputfile)
}
