# Script to plot the distributions of cluster sizes for B73 and PH207

# Read in data files
b73_tandem <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Filtering/B73_True_Tandem_Clusters.txt", header=FALSE)
ph207_tandem <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Filtering/PH207_True_Tandem_Clusters.txt", header=FALSE)

# Define colors for B73 and PH207 for consistency
#   This is a nice blue that is print-friendly
b73_color <- rgb(2/256, 112/256, 189/256)
#   And a nice red that is print-friendly
ph207_color <- rgb(237/256, 28/256, 36/256)

# Define a function to return the number of genes in the comma-separated list
count_genes <- function(glist) {
    # Cast to character
    c_glist <- as.character(glist)
    # Split on commas, cast to vector
    genes <- unlist(strsplit(c_glist, ","))
    # Return the length
    return(length(genes))
}

# Apply the counting function over the second column of the data files
b_counts <- sapply(b73_tandem$V2, count_genes)
p_counts <- sapply(ph207_tandem$V2, count_genes)

# Count up the numbers. We use a custom function here because we want to have
# counts of 0 for some cluster sizes
b_tab <- sapply(
    seq(2, 20),
    function(x) {
        return(sum(b_counts == x))
        })
p_tab <- sapply(
    seq(2, 20),
    function(x) {
        return(sum(p_counts == x))
        })

# Put them into a matrix for plotting
to_plot <- matrix(
    c(b_tab, p_tab),
    ncol=2,
    byrow=FALSE)
# And plot it
pdf(file="Tandem_Cluster_Sizes.pdf", 6, 6)
at <- barplot(
    t(to_plot),
    beside=TRUE,
    col=c(b73_color, ph207_color),
    axes=FALSE,
    ylab="Number of Clusters",
    xlab="Number of Genes in Cluster",
    main="Distribution of Cluster Sizes")
axis(side=2)
axis(side=1, labels=seq(2, 20)[c(TRUE, FALSE)], at=apply(at, 2, mean)[c(TRUE, FALSE)])
legend("topright", c("B73 Tandem", "PH207 Tandem"), fill=c(b73_color, ph207_color))
dev.off()

pdf(file="Tandem_Cluster_Sizes_Pub.pdf", 3, 3)
# Set margins. The order is bottom, left, top, right
par(mar=c(4, 4, 0.1, 0.1), mgp=c(2, 1, 0))
at <- barplot(
    t(to_plot),
    beside=TRUE,
    col=c(b73_color, ph207_color),
    axes=FALSE,
    ylab="Number of Clusters",
    xlab="Number of Genes in Cluster",
    main="")
axis(side=2, at=c(0, 250, 500, 750, 1000, 1250), labels=c(0, 250, 500, 750, 1000, ""))
axis(side=1, at=apply(at, 2, mean)[c(TRUE, FALSE, FALSE)], labels=NA)
legend("topright", c("B73", "PH207"), fill=c(b73_color, ph207_color))
mtext(seq(2, 20)[c(TRUE, FALSE, FALSE)], side=1, at=apply(at, 2, mean)[c(TRUE, FALSE, FALSE)], padj=1)
box(which="plot", col="black")
dev.off()
