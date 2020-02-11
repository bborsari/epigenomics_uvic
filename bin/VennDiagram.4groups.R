#*****************
# OPTION PARSING *
#*****************

suppressPackageStartupMessages(library("optparse"))

option_list <- list (
  
  make_option ( c("--setA"),
                help = "stomach H3K4me3" ),
  
  make_option ( c("--setB"), 
                help = "stomach POLR2A" ),
  
  make_option ( c("--setC"), 
                help = "sigmoid_colon H3K4me3" ),
  
  make_option ( c("--setD"),
                help = "sigmoid_colon POLR2A" ),
  
  make_option ( c("--output"), help = "Output file name") 
  
)


parser <- OptionParser(
  usage = "%prog [options]", 
  option_list=option_list,
  description = "\nPlots a 4-group Venn diagram."
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options




#************
# LIBRARIES *
#************

library(VennDiagram)
library(ggplot2)



#********
# BEGIN *
#********

# 1. import dataframes
setA <- read.table( file = opt$setA, header = F, quote = NULL, sep="\t", stringsAsFactors = F )
setB <- read.table( file = opt$setB, header = F, quote = NULL, sep="\t", stringsAsFactors = F )
setC <- read.table( file = opt$setC, header = F, quote = NULL, sep="\t", stringsAsFactors = F )
setD <- read.table( file = opt$setD, header = F, quote = NULL, sep="\t", stringsAsFactors = F )

# 2. retrieve genes lists
setA <- setA[, 1]
setB <- setB[, 1]
setC <- setC[, 1]
setD <- setD[, 1]



venn.diagram(
  x = list(setA, setB, setC, setD),
  category.names = c("stomach H3K4me3" , "stomach POLR2A" , "sigmoid_colon H3K4me3", "sigmoid_colon POLR2A"),
  filename = opt$output,
  output=TRUE, 
  col = "black",
  fill = c(alpha("#1f78b4", 1), alpha('#cc4c02', 1), alpha("#1f78b4", 0.4), alpha('#cc4c02', 0.4)),
  width = 4000
)
