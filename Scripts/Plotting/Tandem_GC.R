# Script to plot GC content distributions for various partitions of tandem
# duplicates.

# Read the data tables
b_genome <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Data/References/B73_Genes_GC_Prop.txt", header=F)$V2
p_genome <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Data/References/PH207_Genes_GC_Prop.txt", header=F)$V2
tandem <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Dating/Tandem_Dates_with_GC.txt", header=T)

# Make a scatterplot of age v. GC content (mean of two genes)
tandem$Mean_GC <- (tandem$Gene1_GC + tandem$Gene2_GC)/2
# Separate syntenic and nonsyntenic
syn <- tandem[tandem$Class == "Syntenic",]
nonsyn <- tandem[tandem$Class == "Nonsyntenic",]

pdf(file="Tandem_GC_Ages.pdf", 6, 6)
plot(
    syn$Age ~ syn$Mean_GC,
    pch=19,
    col="black",
    cex=0.25,
    xlab="Mean Proportion GC of Tandem Duplicates",
    ylab="Estimated Age (MYA)",
    main="Age v. GC Content of Tandem Duplicates",
    ylim=c(0, 15))
points(
    nonsyn$Age ~ nonsyn$Mean_GC,
    pch=19,
    col="red",
    cex=0.25)
legend(
    "topright",
    c("Syntenic", "Nonsyntenic"),
    col=c("black", "red"),
    pch=19)
dev.off()


# Plot distributions of the GC content for all genes and tandem genes
pdf(file="All_Genes_GC.pdf", 6, 6)
plot(
    density(b_genome),
    col="blue",
    lwd=2,
    xlab="Proportion GC in Longest Transcript",
    ylab="Density",
    main="GC Content Distributions",
    xlim=c(0.2, 0.9),
    ylim=c(0, 6.5))
lines(density(p_genome), col="red", lwd=2)
lines(density(c(syn$Gene1_GC, syn$Gene2_GC)), col="black", lwd=2)
lines(density(c(nonsyn$Gene1_GC, nonsyn$Gene2_GC)), col="grey", lwd=2)
legend(
    "topleft",
    c("B73 Genes", "PH207 Genes", "Syntenic Tandems", "Nonsyntenic Tandems"),
    col=c("blue", "red", "black", "grey"),
    lwd=2)
dev.off()

# Perhaps it is more informative to look at "old" and "new" tandems instead
old <- tandem[tandem$Age >= 10,]
new <- tandem[tandem$Age <= 2,]
pdf(file="Binary_Age_GC.pdf", 6, 6)
plot(
    density(c(old$Gene1_GC, old$Gene2_GC)),
    col="black",
    lwd=2,
    xlab="Proportion GC in Longest Transcript",
    ylab="Density",
    main="GC Content Distribuions for Ancient and Recent Tandems")
lines(density(c(new$Gene1_GC, new$Gene2_GC)), col="red", lwd=2)
legend(
    "topleft",
    c(">= 10MYA (N=2,107)", "<= 2MYA (N=1,039)"),
    col=c("black", "red"),
    lwd=2)
dev.off()
