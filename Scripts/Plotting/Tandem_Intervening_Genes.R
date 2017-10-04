# Plot the distribution of intervening genes between tandem duplicates. The
# script that was used to calculate cluster "widths" includes the tandem genes
# in the count, so we subtract the number of tandem duplciates from the width

# Read in the data
b73_widths <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Filtering/B73_Cluster_Widths.txt", header=TRUE)
ph207_widths <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Filtering/PH207_Cluster_Widths.txt", header=TRUE)

# Define a function to subtract the number of tandem duplicate genes from the
# cluster "width"
calc_width <- function(x) {
    ngenes <- length(unlist(strsplit(as.character(x["Genes"]), ",")))
    w <- as.numeric(x["Width"]) - ngenes
    return(w)
}

b_w <- apply(b73_widths, 1, calc_width)
p_w <- apply(ph207_widths, 1, calc_width)

# Make a matrix to plot
toplot <- matrix(0, ncol=2, nrow=17)
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
# Then, add the final row, for greater than 15
toplot[17, 1] <- sum(b_w > 15)
toplot[17, 2] <- sum(p_w > 15)

pdf(file="Cluster_Intervening_Genes.pdf", 6, 6)
at <- barplot(
    t(toplot),
    beside=TRUE,
    col=c("black", "grey"),
    xlab="Number of Intervening Genes",
    ylab="Count",
    main="Number of Intervening Genes in Tandem Clusters",
    axes=F)
axis(side=2)
axis(side=1, at=apply(at, 2, mean)[c(TRUE, FALSE)], labels=c(as.character(0:15), ">15")[c(TRUE, FALSE)])
legend("topright", c("B73", "PH207"), fill=c("black", "grey"))
dev.off()
