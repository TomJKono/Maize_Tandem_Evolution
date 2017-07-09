# Make plots of the mean number of gene sequences per species for tandems that
# are split or not-split across orthogroups.

setwd("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Orthofinder")

# Read in the data
b_false_single_genes <- read.table("B73_False_Single_OGGenes.txt", header=FALSE)
b_true_single_genes <- read.table("B73_True_Single_OGGenes.txt", header=FALSE)
b_false_multi_genes <- read.table("B73_False_Multi_OGGenes.txt", header=FALSE)
b_true_multi_genes <- read.table("B73_True_Multi_OGGenes.txt", header=FALSE)


p_false_single_genes <- read.table("PH207_False_Single_OGGenes.txt", header=FALSE)
p_true_single_genes <- read.table("PH207_True_Single_OGGenes.txt", header=FALSE)
p_false_multi_genes <- read.table("PH207_False_Multi_OGGenes.txt", header=FALSE)
p_true_multi_genes <- read.table("PH207_True_Multi_OGGenes.txt", header=FALSE)

# Define a function to return a numeric vector of means
comma_mean <- function(string) {
    string <- as.character(string)
    counts <- unlist(strsplit(string, ","))
    counts <- as.numeric(counts)
    return(mean(counts))
}

# Get means for each of the partitions
b_false_single_means <- sapply(b_false_single_genes$V3, comma_mean)
b_true_single_means <- sapply(b_true_single_genes$V3, comma_mean)
b_false_multi_means <- sapply(b_false_multi_genes$V3, comma_mean)
b_true_multi_means <- sapply(b_true_multi_genes$V3, comma_mean)

p_false_single_means <- sapply(p_false_single_genes$V3, comma_mean)
p_true_single_means <- sapply(p_true_single_genes$V3, comma_mean)
p_false_multi_means <- sapply(p_false_multi_genes$V3, comma_mean)
p_true_multi_means <- sapply(p_true_multi_genes$V3, comma_mean)

# Make some plots
pdf(file="B73_GenesPerSp_SplitOGs.pdf", 6, 6)
single <- c(b_false_single_means, b_true_single_means)
multi <- c(b_false_multi_means, b_true_multi_means)
hist(
    log(single),
    breaks=20,
    col=rgb(0, 0, 1, 0.5),
    xlab="log(Mean Number of Genes per Species in Orthogroup)",
    ylab="Count",
    main="Distribution of Genes per Species in\nOrthogroups With B73 Tandem Duplicates",
    ylim=c(0, 250))
hist(
    log(multi),
    breaks=20,
    col=rgb(1, 0, 0, 0.5),
    add=TRUE)
legend(
    "topright",
    c("Not Split", "Split"),
    fill=c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5)))
dev.off()

pdf(file="B73_GenesPerSp_TrueFalse.pdf", 6, 6)
truetand <- c(b_true_multi_means, b_true_single_means)
falsetand <- c(b_false_multi_means, b_false_single_means)
hist(
    log(truetand),
    breaks=20,
    col=rgb(0, 0, 1, 0.5),
    xlab="log(Mean Number of Genes per Species in Orthogroup)",
    ylab="Count",
    main="Distribution of Genes per Species in\nOrthogroups With B73 Tandem Duplicates",
    ylim=c(0, 250))
hist(
    log(falsetand),
    breaks=20,
    col=rgb(1, 0, 0, 0.5),
    add=TRUE)
legend(
    "topright",
    c("True Tandem", "False Tandem"),
    fill=c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5)))
dev.off()

pdf(file="PH207_GenesPerSp_SplitOGs.pdf", 6, 6)
single <- c(p_false_single_means, p_true_single_means)
multi <- c(p_false_multi_means, p_true_multi_means)
hist(
    log(single),
    breaks=20,
    col=rgb(0, 0, 1, 0.5),
    xlab="log(Mean Number of Genes per Species in Orthogroup)",
    ylab="Count",
    main="Distribution of Genes per Species in\nOrthogroups With PH207 Tandem Duplicates",
    ylim=c(0, 400))
hist(
    log(multi),
    breaks=20,
    col=rgb(1, 0, 0, 0.5),
    add=TRUE)
legend(
    "topright",
    c("Not Split", "Split"),
    fill=c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5)))
dev.off()

pdf(file="PH207_GenesPerSp_TrueFalse.pdf", 6, 6)
truetand <- c(p_true_multi_means, p_true_single_means)
falsetand <- c(p_false_multi_means, p_false_single_means)
hist(
    log(truetand),
    breaks=20,
    col=rgb(0, 0, 1, 0.5),
    xlab="log(Mean Number of Genes per Species in Orthogroup)",
    ylab="Count",
    main="Distribution of Genes per Species in\nOrthogroups With PH207 Tandem Duplicates",
    ylim=c(0, 400))
hist(
    log(falsetand),
    breaks=20,
    col=rgb(1, 0, 0, 0.5),
    add=TRUE)
legend(
    "topright",
    c("True Tandem", "False Tandem"),
    fill=c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5)))
dev.off()
