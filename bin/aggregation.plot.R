#*****************
# OPTION PARSING *
#*****************

suppressPackageStartupMessages(library("optparse"))

option_list <- list (
  
  make_option ( c("--most"),
                help = ".most.expressed.genes.aggregate.tsv" ),
  
  make_option ( c("--least"),
                help = ".least.expressed.genes.aggregate.tsv" ),
  
  make_option ( c("--tissue"), 
                help = "tissue"),
  
  make_option ( c("--output"), help = "Output file name") 
  
)


parser <- OptionParser(
  usage = "%prog [options]", 
  option_list=option_list,
  description = "\nPlots aggregation signal."
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


# 1. read dataframe of bwtool aggregate - most expressed genes
most.expressed <- read.table( file = opt$most, header = F, quote = NULL, sep="\t" )

# 2. read dataframe of bwtool aggregate - least expressed genes
least.expressed <- read.table( file = opt$least, header = F, quote = NULL, sep="\t" )

# 3. add colnames
colnames(most.expressed) <- c("position", "signal")
colnames(least.expressed) <- c("position", "signal")

# 4. add column with group specification, it will be needed to define the colour of the plot
most.expressed$group <- "1000 most expressed genes"
least.expressed$group <- "1000 least expressed genes"

# 5. merge the 2 dataframes into a single one
df.plot <- rbind(most.expressed, least.expressed)

# 6. define color palette
palette <- c("1000 least expressed genes" = "#00BFC4",
             "1000 most expressed genes" = "#F8766D")

# 7. make plot
pdf(opt$output, # define the output filename
    width = 6, # specify width and height of the plot
    height=3) 
ggplot(df.plot, # define the dataframe to be used for the plot
       aes(x=position, # the column that corresponds to the x axis
           y=signal, # the column that corresponds to the y axis
           group=group, # additional parameter required for line plots
           color=group)) + # the column that contains info for the colour
  geom_line() + # choose the type of plot you want to make: in this case, a line plot
  scale_color_manual(values = palette) + # define your own palette colour
  labs(title = opt$tissue) + # add as plot title the tissue name
  theme(plot.title = element_text(hjust = 0.5)) + # move plot title to the centre
  xlab("position from TSS (bp)") + # add axis x title
  ylab("fold-change signal") # add axis y title
dev.off()



