# Script to plot the estiamted duplication dates of tandem clusters that are
# in cassettes

# Read the data files
b_cassettes <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Dating/B73_Two_Cluster_Cassette_Dates.txt", header=FALSE)
p_cassettes <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Dating/PH207_Two_Cluster_Cassette_Dates.txt", header=FALSE)


# Define colors for B73 and PH207 for consistency
#   This is a nice blue that is print-friendly
b73_color <- rgb(2/256, 112/256, 189/256)
#   And a nice red that is print-friendly
ph207_color <- rgb(237/256, 28/256, 36/256)

# Generate vectors of numbers for the X and Y coordinates
b_ages_x <- sapply(b_cassettes$V1, function(x) {
    return(unlist(strsplit(as.character(x), ","))[1])})
b_ages_y <- sapply(b_cassettes$V1, function(x) {
    return(unlist(strsplit(as.character(x), ","))[2])})
p_ages_x <- sapply(p_cassettes$V1, function(x) {
    return(unlist(strsplit(as.character(x), ","))[1])})
p_ages_y <- sapply(p_cassettes$V1, function(x) {
    return(unlist(strsplit(as.character(x), ","))[2])})


# Plot them
pdf(file="Cassette_Ages.pdf", 3, 3)
# Set margins. The order is bottom, left, top, right
par(mar=c(4, 4, 0.1, 0.1), mgp=c(2, 1, 0))
plot(b_ages_y ~ b_ages_x,
    pch=19,
    xlab="Est. Age of Cluster 1 (MYA)",
    ylab="Est. Age of Cluster 2 (MYA)",
    main="",
    col=b73_color,
    axes=FALSE)
points(p_ages_y ~ p_ages_x, pch=19, col=ph207_color)
abline(c(0, 0), c(1, 1), lwd=2, col="black", lty=2)
axis(side=1, at=c(0, 2, 4, 6, 8, 10, 11.899), labels=c(0, 2, 4, 6, 8, 10, 12))
axis(side=2, at=c(0, 2, 4, 6, 8, 10, 11.899), labels=c(0, 2, 4, 6, 8, 10, 12))
box()
dev.off()
