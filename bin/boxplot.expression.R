#*****************
# OPTION PARSING *
#*****************

suppressPackageStartupMessages(library("optparse"))

option_list <- list (
  
  make_option ( c("--expression"),
                help = "expression matrix" ),
  
  make_option ( c("--marked_both_tissues"),
                help = "list of genes marked in both tissues" ),
  
  make_option ( c("--stomach_specific"), 
                help = "list of genes with stomach-specific marking" ),
  
  make_option ( c("--sigmoid_colon_specific"), 
                help = "list of genes with sigmoid_colon-specific marking" ),
  
  make_option ( c("--not_marked"),
                help = "list of genes not marked in any tissue" ),
  
  make_option ( c("--output"), help = "Output file name") 
  
)


parser <- OptionParser(
  usage = "%prog [options]", 
  option_list=option_list,
  description = "\nPlots distributions of expression values."
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options




#************
# LIBRARIES *
#************

library(ggplot2)
library(reshape2)


#********
# BEGIN *
#********

# debugging options
# expression.matrix <- read.table( file = "~/Documents/learn/courses/epigenomics_uvic/ChIP-seq/downstream.analyses/build.matrix/expression.matrix.tsv", header = T, quote = NULL, sep="\t" )
# marked.both.tissues <- read.table( file = "~/Documents/learn/courses/epigenomics_uvic/ChIP-seq/downstream.analyses/peaks.analysis/genes.marked.both.tissues.txt", h=F, quote = NULL, sep="\t",
#                                    stringsAsFactors = F)
# stomach.specific <- read.table( file = "~/Documents/learn/courses/epigenomics_uvic/ChIP-seq/downstream.analyses/peaks.analysis/genes.with.stomach.specific.peaks.txt", header = F, quote = NULL, sep="\t",
#                                 stringsAsFactors = F)
# sigmoid_colon.specific <- read.table( file = "~/Documents/learn/courses/epigenomics_uvic/ChIP-seq/downstream.analyses/peaks.analysis/genes.with.sigmoid_colon.specific.peaks.txt", header = F, quote = NULL, sep="\t",
#                                       stringsAsFactors = F)
# not.marked <- read.table( file = "~/Documents/learn/courses/epigenomics_uvic/ChIP-seq/downstream.analyses/peaks.analysis/genes.not.marked.txt", header = F, quote = NULL, sep="\t",
#                           stringsAsFactors = F)

# 1. read expression matrix
expression.matrix <- read.table( file = opt$expression, header = T, quote = NULL, sep="\t",
                                 stringsAsFactors = F)


# 2. remove "." in the gene_id
rownames(expression.matrix) <- gsub("\\..*", "", rownames(expression.matrix))


# 3. read list of genes marked in both tissues
marked.both.tissues <- read.table( file = opt$marked_both_tissues, header = F, quote = NULL, sep="\t",
                                   stringsAsFactors = F)
marked.both.tissues <- marked.both.tissues$V1

# 4. read list of genes with stomach-specific marking
stomach.specific <- read.table( file = opt$stomach_specific, header = F, quote = NULL, sep="\t",
                                stringsAsFactors = F)
stomach.specific <- stomach.specific$V1

# 5. read list of genes with sigmoid_colon-specific marking
sigmoid_colon.specific <- read.table( file = opt$sigmoid_colon_specific, header = F, quote = NULL, sep="\t",
                                      stringsAsFactors = F)
sigmoid_colon.specific <- sigmoid_colon.specific$V1

# 6. read list of genes not marked in either of the tissues
not.marked <- read.table( file = opt$not_marked, header = F, quote = NULL, sep="\t",
                          stringsAsFactors = F)
not.marked <- not.marked$V1


# 7. retrieve expression for 

# 7.1. genes marked in both tissues
expression.marked.both.tissues <- expression.matrix[rownames(expression.matrix) %in% marked.both.tissues, ]
expression.marked.both.tissues$group <- "marked both tissues"

# 7.2. genes w/ stomach-specific marking
expression.stomach.specific <- expression.matrix[rownames(expression.matrix) %in% stomach.specific, ]
expression.stomach.specific$group <- "stomach-specific marking"

# 7.3. genes w/ sigmoid_colon-specific marking
expression.sigmoid_colon.specific <- expression.matrix[rownames(expression.matrix) %in% sigmoid_colon.specific, ]
expression.sigmoid_colon.specific$group <- "sigmoid_colon-specific marking"

# 7.4. genes not marked in either of the tissues
expression.not.marked <- expression.matrix[rownames(expression.matrix) %in% not.marked, ]
expression.not.marked$group <- "not marked"


# 8. prepare dataframe for plot
df <- rbind(expression.marked.both.tissues,
            expression.stomach.specific,
            expression.sigmoid_colon.specific,
            expression.not.marked)

df.plot <- melt(df)


# 9. define palette
palette <- c("stomach" = "#dfc27d",
             "sigmoid_colon" = "#8c510a")


# 10. make plot
pdf(opt$output, # define the output filename
    width = 12, # specify width and height of the plot
    height=4) 
ggplot(df.plot, aes(x=variable, y=log2(value+1), fill=variable)) +
  geom_boxplot() +
  guides(fill=F) +
  scale_fill_manual(values = palette) +
  facet_wrap(~group, scales = "free_y", ncol = 4) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        strip.text.x = element_text(size=12),
        axis.text.x = element_text(size=12)) +
  ylab("log2(TPMs+1)")
dev.off()