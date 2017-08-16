# Make a heatmap plot of the homologous tandem cluster sizes

# Read the data files
homology <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Homology/B_P_Homologous_Tandems.txt", header=TRUE)

# We want to make a vector of counts for each row. Define a function to do this
genecount <- function(x) {
    y <- as.character(x)
    z <- unlist(strsplit(y, ","))
    return(length(z))
}

b_counts <- sapply(homology$B73_Genes, genecount)
p_counts <- sapply(homology$PH207_Genes, genecount)

# cbind them to process later
cts <- cbind(b_counts, p_counts)

# Define a matrix to hold the heatmap data
heat_mat <- matrix(0, ncol=20, nrow=20)
# Iterate through the gene counts and set the matrix values
for(value in 1:length(b_counts)) {
    row <- cts[value, 1]
    col <- cts[value, 2]
    heat_mat[[row, col]] <- heat_mat[[row, col]] + 1
    }

pdf(file="B73_PH207_Homologous_Sizes.pdf", 3, 3)
# Set margins. The order is bottom, left, top, right
par(mar=c(4, 4, 0.1, 0.1))
image(
    log(1+heat_mat),
    x=seq(1:20),
    y=seq(1:20),
    col=heat.colors(12),
    xlab="B73 Cluster Size",
    ylab="PH207 Cluster Size",
    cex.axis=0.8
    )
abline(c(0, 0), 1, lwd=2, col="black")
dev.off()
