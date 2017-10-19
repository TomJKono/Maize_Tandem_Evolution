# Make plots of the number of exons and the length of the gene in bp for both
# tandem duplicates and genomewide genes

# Read the data files
b_genome <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/TEs/B73_FullGene_Genomewide_TE_NumExons.txt", header=TRUE)
p_genome <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/TEs/PH207_FullGene_Genomewide_TE_NumExons.txt", header=TRUE)

# Read the clusters. We are only interested in the gene IDs
b_tandem <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Filtering/B73_True_Tandem_Clusters.txt", header=FALSE)$V2
p_tandem <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Filtering/PH207_True_Tandem_Clusters.txt", header=FALSE)$V2

# Convert the comma-separated values into a character vector
b_tandem <- unlist(strsplit(as.character(b_tandem), ","))
p_tandem <- unlist(strsplit(as.character(p_tandem), ","))

# Then, separate the values for tandems
b_tand_dat <- b_genome[b_genome$GenesInCluster %in% b_tandem,]
p_tand_dat <- p_genome[p_genome$GenesInCluster %in% p_tandem,]

# And the non-tandems
b_nontand_dat <- b_genome[!(b_genome$GenesInCluster %in% b_tandem),]
p_nontand_dat <- p_genome[!(p_genome$GenesInCluster %in% p_tandem),]

# Then, make a plot of length v. number of exons for both tandem and non-tandem
pdf(file="B73_NumExons_Length.pdf", 3, 3)
# Set margins. The order is bottom, left, top, right
par(mar=c(4, 4, 0.1, 0.1), mgp=c(2, 1, 0))
plot(
    b_nontand_dat$Lengths/1000 ~ b_nontand_dat$NumExons,
    col="grey",
    pch=19,
    cex=0.5,
    main="",
    xlab="Number of Exons",
    ylab="Length (kb)")
points(
    b_tand_dat$Lengths/1000 ~ b_tand_dat$NumExons,
    col="black",
    pch=19,
    cex=0.5)
legend("topright", c("Non-Tandem", "Tandem"), col=c("grey", "black"), pch=19)
dev.off()

pdf(file="PH207_NumExons_Length.pdf", 3, 3)
# Set margins. The order is bottom, left, top, right
par(mar=c(4, 4, 0.1, 0.1), mgp=c(2, 1, 0))
plot(
    p_nontand_dat$Lengths/1000 ~ p_nontand_dat$NumExons,
    col="grey",
    pch=19,
    cex=0.5,
    main="",
    xlab="Number of Exons",
    ylab="Length (kb)")
points(
    p_tand_dat$Lengths/1000 ~ p_tand_dat$NumExons,
    col="black",
    pch=19,
    cex=0.5)
legend("topright", c("Non-Tandem", "Tandem"), col=c("grey", "black"), pch=19)
dev.off()
