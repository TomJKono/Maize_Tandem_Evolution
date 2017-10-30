# The same plot as the genomewide plot, except it is only for chromosome 2,
# which has a nice representation of gene density, TE density, tandem
# duplicates, and subgenome assignments

# Poor decision, but use ggplot2
library(ggplot2)
library(reshape2)

# Define colors
m1_col <- "#b2df8a"
m2_col <- "#a6cee3"
dna_te <- "#aaaaaa"
rna_te <- "#666666"
genes <- "#000000"
tandem <- "#984ea3"

# Multiplot function from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#   Read data files
b73_gene_te <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Data/B73_GeneDensity_TEDensity.txt", header=TRUE)
b73_syn_blocks <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Data/ABB_Synteny/b73-syntenic-blocks-merged.txt", header=FALSE)
b73_tandem <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Data/References/B73_True_Tandem_Positions.gff", header=FALSE, sep="\t")
ph207_gene <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Data/PH207_GeneDensity.txt", header=TRUE)
ph207_syn_blocks <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Data/ABB_Synteny/ph207-syntenic-blocks-merged.txt", header=FALSE)
ph207_tandem <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Data/References/PH207_True_Tandem_Positions.gff", header=FALSE)

#   Trim down the tandem positions to just the chromosome and start
b73_tandem <- b73_tandem[,c("V1", "V4")]
names(b73_tandem) <- c("Chromosome", "Start")
# Select only chromosome 2
b73_tandem <- b73_tandem[b73_tandem$Chromosome == 2,]

ph207_tandem <- ph207_tandem[, c("V1", "V4")]
names(ph207_tandem) <- c("Chromosome", "Start")
ph207_tandem <- ph207_tandem[ph207_tandem$Chromosome == 2,]

#   Get maize1 and maize2 blocks
b73_m1 <- b73_syn_blocks[b73_syn_blocks$V5 == "maize1", c(1, 2, 3)]
b73_m2 <- b73_syn_blocks[b73_syn_blocks$V5 == "maize2", c(1, 2, 3)]
ph207_m1 <- ph207_syn_blocks[ph207_syn_blocks$V5 == "maize1", c(1, 2, 3)]
ph207_m2 <- ph207_syn_blocks[ph207_syn_blocks$V5 == "maize2", c(1, 2, 3)]


#   Set the names
names(b73_m1) <- c("Chromosome", "Start", "End")
names(b73_m2) <- c("Chromosome", "Start", "End")
names(ph207_m1) <- c("Chromosome", "Start", "End")
names(ph207_m2) <- c("Chromosome", "Start", "End")

# Get only chromosome 2
b73_m1 <- b73_m1[b73_m1$Chromosome == 2,]
b73_m2 <- b73_m2[b73_m2$Chromosome == 2,]
ph207_m1 <- ph207_m1[ph207_m1$Chromosome == 2,]
ph207_m2 <- ph207_m2[ph207_m2$Chromosome == 2,]


# Select only chromosome 2
b73_gene_te <- b73_gene_te[b73_gene_te$Chromosome == 2,]
ph207_gene <- ph207_gene[ph207_gene$Chromosome == 2,]

#   Make a plot
pdf(file="B73-PH207_Chr2_Tandem.pdf", width=3, height=3)
# Set margins. The order is bottom, left, top, right
par(mar=c(4, 4, 0.1, 0.1), mgp=c(2, 1, 0))
# Plot for B73
p1 <- ggplot(NULL) +
    geom_rect(data=b73_m1, aes(xmin=Start/1000000, xmax=End/1000000, ymin=0, ymax=100), fill=m1_col) +
    geom_rect(data=b73_m2, aes(xmin=Start/1000000, xmax=End/1000000, ymin=0, ymax=100), fill=m2_col) +
    geom_line(data=b73_gene_te, aes(x=Midpoint/1000000, y=GenesPerMb), color=genes, size=0.5) +
    geom_line(data=b73_gene_te, aes(x=Midpoint/1000000, y=DNATEsPerMb), color=dna_te, size=0.5) +
    geom_line(data=b73_gene_te, aes(x=Midpoint/1000000, y=RNATEsPerMb), color=rna_te, size=0.5) +
    geom_segment(data=b73_tandem, mapping=aes(x=Start/1000000, xend=Start/1000000, y=100, yend=80), color=tandem, size=0.25) +
    scale_y_continuous(name="", breaks=seq(0, 90, 30), labels=seq(0, 90, 30)) +
    scale_x_continuous(name="Physical Position (Mb)", breaks=seq(0, 300, 50), labels=seq(0, 300, 50)) +
    theme_bw() +
    theme(
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        strip.background=element_blank()
        )

# Plot PH207, without the RNA and DNA TE lines.
p2 <- ggplot(NULL) +
    geom_rect(data=ph207_m1, aes(xmin=Start/1000000, xmax=End/1000000, ymin=0, ymax=100), fill=m1_col) +
    geom_rect(data=ph207_m2, aes(xmin=Start/1000000, xmax=End/1000000, ymin=0, ymax=100), fill=m2_col) +
    geom_line(data=ph207_gene, aes(x=Midpoint/1000000, y=GenesPerMb), color=genes, size=0.5) +
    geom_segment(data=ph207_tandem, mapping=aes(x=Start/1000000, xend=Start/1000000, y=100, yend=80), color=tandem, size=0.25) +
    scale_y_continuous(name="", breaks=seq(0, 90, 30), labels=seq(0, 90, 30)) +
    scale_x_continuous(name="Physical Position (Mb)", breaks=seq(0, 300, 50), labels=seq(0, 300, 50)) +
    theme_bw() +
    theme(
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        strip.background=element_blank()
        )
multiplot(p1, p2, cols=1)
dev.off()
