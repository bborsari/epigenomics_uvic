#*****************
# OPTION PARSING *
#*****************

suppressPackageStartupMessages(library("optparse"))

option_list <- list (
  
  make_option ( c("--expression"),
                help = "expression matrix" ),
  
  make_option ( c("--mark"),
                help = "H3K4me3 matrix" ),
  
  make_option ( c("--tissue"), 
                help = "tissue"),
  
  make_option ( c("--output"), help = "Output file name") 
  
)


parser <- OptionParser(
  usage = "%prog [options]", 
  option_list=option_list,
  description = "\nScatterplot of expression and mark levels."
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options




#************
# LIBRARIES *
#************

library(ggplot2)



#********
# BEGIN *
#********


# 1. read expression matrix
expression.matrix <- read.table( file = opt$expression, header = T, quote = NULL, sep="\t" )

# 2. read mark matrix
mark.matrix <- read.table( file = opt$mark, header = T, quote = NULL, sep="\t" )

# 3. check order of rownames is the same
stopifnot(identical(rownames(expression.matrix),
                    rownames(mark.matrix)))

# 4. prepare dataframe for plot
df.plot <- data.frame(x = expression.matrix[, opt$tissue],
                      y = mark.matrix[, opt$tissue])

# 5. compute correlation coefficients
# 5.1. Pearson
pearson.cor <- round(cor(log2(df.plot$x+1), df.plot$y, method = "pearson"), 2)
# 5.2. Spearman
spearman.cor <- round(cor(log2(df.plot$x+1), df.plot$y, method = "spearman"), 2)

# 6. make plot
pdf(opt$output, # define the output filename
    width = 4, # specify width and height of the plot
    height=4) 
ggplot(df.plot, aes(x=log2(df.plot$x +1),
                    y=df.plot$y)) +
  geom_point(shape=21) +
  xlab("expression - log2(TPM+1)") +
  ylab("H3K4me3") +
  geom_smooth(method="lm") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = opt$tissue) +
  annotate(geom="text", x=2.5, y=80, 
           label=paste0("Pearson cc: ", pearson.cor),
           color="black") +
  annotate(geom="text", x=2.5, y=70, 
           label=paste0("Spearman cc: ", spearman.cor),
           color="black")
dev.off()


