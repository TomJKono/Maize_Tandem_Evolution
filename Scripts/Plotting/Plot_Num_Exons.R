# Plot the number of exons in tandems with respect to the genome wide
# distribution of representative transcripts.

# Read the data files. The tandem data have headers, and the genomewide data
# do not.
b_tandem <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/TEs/B73_TE_NumExons.txt", header=TRUE)
p_tandem <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/TEs/PH207_TE_NumExons.txt", header=TRUE)
b_genome <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/TEs/B73_Genomewide_NumExons.txt", header=FALSE)$V1
p_genome <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/TEs/PH207_Genomewide_NumExons.txt", header=FALSE)$V1


# Define colors for consistency
tandem_col <- "#984ea3"
gw_col <- "#999999"

# Get the number of exons in tandems, which is stored as a vector of
# comma-separated values.
b_tandem <- as.numeric(unlist(strsplit(as.character(b_tandem$NumExons), ",")))
p_tandem <- as.numeric(unlist(strsplit(as.character(p_tandem$NumExons), ",")))

print(c(mean(b_genome), median(b_genome)))
print(c(mean(b_tandem), median(b_tandem)))
print(c(mean(p_genome), median(p_genome)))
print(c(mean(p_tandem), median(p_tandem)))

t.test(b_tandem, b_genome)
t.test(p_tandem, p_genome)

# Trim down the data range to just 10 exons
b_genome[b_genome > 10] <- 11
p_genome[p_genome > 10] <- 11
b_tandem[b_tandem > 10] <- 11
p_tandem[p_tandem > 10] <- 11

# Count them up and transform them to proportions
b_genome_prop <- table(b_genome)/length(b_genome)
b_tandem_prop <- table(b_tandem)/length(b_tandem)
p_genome_prop <- table(p_genome)/length(p_genome)
p_tandem_prop <- table(p_tandem)/length(p_tandem)

# Combine them for plotting
b_plot <- matrix(c(b_genome_prop, b_tandem_prop), ncol=2)
p_plot <- matrix(c(p_genome_prop, p_tandem_prop), ncol=2)

# Plot them
pdf(file="B73_Tandem_NumExons.pdf", 6, 6)
at <- barplot(t(b_plot),
    xlab="Number of Exons",
    ylab="Proportion",
    main="Number of Exons in B73 Representative Transcripts",
    col=c(gw_col, tandem_col),
    beside=TRUE,
    axes=FALSE)
axis(side=2)
axis(side=1, at=apply(at, 2, mean), labels=c(as.character(1:10), ">10"))
legend("topright", c("Genomewide", "Tandem"), fill=c(gw_col, tandem_col))
dev.off()

pdf(file="PH207_Tandem_NumExons.pdf", 6, 6)
at <- barplot(t(p_plot),
    xlab="Number of Exons",
    ylab="Proportion",
    main="Number of Exons in PH207 Representative Transcripts",
    col=c(gw_col, tandem_col),
    beside=TRUE,
    axes=FALSE)
axis(side=2)
axis(side=1, at=apply(at, 2, mean), labels=c(as.character(1:10), ">10"))
legend("topright", c("Genomewide", "Tandem"), fill=c(gw_col, tandem_col))
dev.off()

# Plot them
pdf(file="B73_Tandem_NumExons_Pub.pdf", 3, 3)
# Set margins. The order is bottom, left, top, right
par(mar=c(4, 4, 0.1, 0.1), mgp=c(2, 1, 0))
at <- barplot(t(b_plot),
    xlab="Number of Exons",
    ylab="Proportion of Genes",
    main="",
    col=c(gw_col, tandem_col),
    beside=TRUE,
    axes=FALSE)
axis(side=2)
axis(side=1, at=apply(at, 2, mean), labels=c(as.character(1:10), ">10"))
legend("topright", c("Genomewide", "Tandem"), fill=c(gw_col, tandem_col))
box()
dev.off()

pdf(file="PH207_Tandem_NumExons_Pub.pdf", 3, 3)
# Set margins. The order is bottom, left, top, right
par(mar=c(4, 4, 0.1, 0.1), mgp=c(2, 1, 0))
at <- barplot(t(p_plot),
    xlab="Number of Exons",
    ylab="Proportion of Genes",
    main="",
    col=c(gw_col, tandem_col),
    beside=TRUE,
    axes=FALSE)
axis(side=2)
axis(side=1, at=apply(at, 2, mean), labels=c(as.character(1:10), ">10"))
legend("topright", c("Genomewide", "Tandem"), fill=c(gw_col, tandem_col))
box()
dev.off()
