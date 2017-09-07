# Make a heatmap plot of the homologous tandem cluster sizes

# Read the data files
syntenic <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Homology/Syntenic_Homologous_Sizes.txt", header=TRUE)
nonsyntenic <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Homology/Nonsyntenic_Homologous_Sizes.txt", header=TRUE)

# Define a matrix to hold the heatmap data
syn_hmp <- matrix(0, ncol=21, nrow=21)
# Iterate through the gene counts and set the matrix values
for(value in 1:nrow(syntenic)) {
    row <- syntenic[value, 1] + 1
    col <- syntenic[value, 2] + 1
    syn_hmp[[row, col]] <- syn_hmp[[row, col]] + 1
    }

nonsyn_hmp <- matrix(0, ncol=21, nrow=21)
# Iterate through the gene counts and set the matrix values
for(value in 1:nrow(nonsyntenic)) {
    row <- nonsyntenic[value, 1] + 1
    col <- nonsyntenic[value, 2] + 1
    nonsyn_hmp[[row, col]] <- nonsyn_hmp[[row, col]] + 1
    }

cor(syntenic$B_Size, syntenic$P_Size)
cor(nonsyntenic$B_Size, nonsyntenic$P_Size)

pdf(file="Syntenic_Homologous_Sizes_Pub.pdf", 3, 3)
# Set margins. The order is bottom, left, top, right
par(mar=c(4, 4, 0.1, 0.1))
image(
    log(1+syn_hmp),
    x=seq(0:20),
    y=seq(0:20),
    col=heat.colors(12),
    xlab="B73 Cluster Size",
    ylab="PH207 Cluster Size",
    cex.axis=0.8
    )
abline(c(0, 0), 1, lwd=2, col="black")
dev.off()

pdf(file="Nonsyntenic_Homologous_Sizes_Pub.pdf", 3, 3)
# Set margins. The order is bottom, left, top, right
par(mar=c(4, 4, 0.1, 0.1))
image(
    log(1+nonsyn_hmp),
    x=seq(0:20),
    y=seq(0:20),
    col=heat.colors(12),
    xlab="B73 Cluster Size",
    ylab="PH207 Cluster Size",
    cex.axis=0.8
    )
abline(c(0, 0), 1, lwd=2, col="black")
dev.off()

# Make plots for the markdown doc
pdf(file="Syntenic_Homologous_Sizes.pdf", 6, 6)
image(
    log(1+syn_hmp),
    x=seq(0:20),
    y=seq(0:20),
    col=heat.colors(12),
    xlab="B73 Cluster Size",
    ylab="PH207 Cluster Size",
    main="Homologous Syntenic Cluster Sizes",
    cex.axis=0.8
    )
abline(c(0, 0), 1, lwd=2, col="black")
dev.off()

pdf(file="Nonsyntenic_Homologous_Sizes.pdf", 6, 6)
image(
    log(1+nonsyn_hmp),
    x=seq(0:20),
    y=seq(0:20),
    col=heat.colors(12),
    xlab="B73 Cluster Size",
    ylab="PH207 Cluster Size",
    main="Homologous Nonsyntenic Cluster Sizes",
    cex.axis=0.8
    )
abline(c(0, 0), 1, lwd=2, col="black")
dev.off()
