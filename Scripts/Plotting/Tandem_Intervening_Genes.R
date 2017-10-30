# Plot the distribution of intervening genes between tandem duplicates. The
# script that was used to calculate cluster "widths" includes the tandem genes
# in the count, so we subtract the number of tandem duplciates from the width

# Read in the data
b73_widths <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Filtering/B73_Cluster_Widths.txt", header=TRUE)
ph207_widths <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Filtering/PH207_Cluster_Widths.txt", header=TRUE)

# Define colors for B73 and PH207 for consistency
#   This is a nice blue that is print-friendly
b73_color <- rgb(2/256, 112/256, 189/256)
#   And a nice red that is print-friendly
ph207_color <- rgb(237/256, 28/256, 36/256)

b_w <- as.numeric(unlist(strsplit(as.character(b73_widths$Width), ",")))
p_w <- as.numeric(unlist(strsplit(as.character(ph207_widths$Width), ",")))

summary(b_w)
summary(p_w)

# Make a matrix to plot
toplot <- matrix(0, ncol=2, nrow=16)
# Count up how many intervening genes there are between 0 and 15
for(s in seq(0, 15)) {
    # R is 1-indexed, so to actually fill the matrix, we have to add 1 to the
    # number of genes considered
    r <- s+1
    b_c <- sum(b_w == s)
    p_c <- sum(p_w == s)
    toplot[r, 1] <- b_c
    toplot[r, 2] <- p_c
}

pdf(file="Cluster_Intervening_Genes.pdf", 6, 6)
at <- barplot(
    t(toplot),
    beside=TRUE,
    col=c(b73_color, ph207_color),
    xlab="Number of Intervening Genes",
    ylab="Count",
    main="Number of Intervening Genes in Tandem Clusters",
    axes=F)
axis(side=2)
axis(side=1, at=apply(at, 2, mean)[c(TRUE, FALSE)], labels=as.character(0:15)[c(TRUE, FALSE)])
legend("topright", c("B73", "PH207"), fill=c(b73_color, ph207_color))
dev.off()

pdf(file="Cluster_Intervening_Genes_Pub.pdf", 3, 3)
# Set margins. The order is bottom, left, top, right
par(mar=c(4, 4, 0.1, 0.1), mgp=c(2, 1, 0))
at <- barplot(
    t(toplot),
    beside=TRUE,
    col=c(b73_color, ph207_color),
    xlab="Number of Intervening Genes",
    ylab="Count",
    main="",
    axes=F)
axis(side=2, at=c(0, 500, 1000, 1500), labels=c("0", "500", "1000", "1500"))
axis(side=1, at=apply(at, 2, mean)[c(TRUE, FALSE)], labels=NA)
mtext(as.character(0:15)[c(TRUE, FALSE)], side=1, at=apply(at, 2, mean)[c(TRUE, FALSE)], padj=1)
legend("topright", c("B73", "PH207"), fill=c(b73_color, ph207_color))
box(which="plot")
dev.off()
