# Script to generate plots of the distribution of mean gene spacing for tandems
# that are split versus those that are not split

setwd("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Orthofinder")

# Read in the data files
b_int <- read.table("../../Data/B73_Cluster_Intervals.txt", header=TRUE)
p_int <- read.table("../../Data/PH207_Cluster_Intervals.txt", header=TRUE)
b_split <- as.character(read.table("B73_Split_Clusters.txt", header=F)$V1)
b_nosplit <- as.character(read.table("B73_Nonsplit_Clusters.txt", header=F)$V1)
p_split <- as.character(read.table("PH207_Split_Clusters.txt", header=F)$V1)
p_nosplit <- as.character(read.table("PH207_Nonsplit_Clusters.txt", header=F)$V1)

# Then, get the intervals for those that were spit and those that were not
b_s_int <- b_int[b_int$ClusterID %in% b_split,]
b_n_int <- b_int[b_int$ClusterID %in% b_nosplit,]
p_s_int <- p_int[p_int$ClusterID %in% p_split,]
p_n_int <- p_int[p_int$ClusterID %in% p_nosplit,]

# Define a function to calculate the mean interval size
intmean <- function(intervals) {
    ints <- unlist(strsplit(as.character(intervals), ","))
    return(mean(as.numeric(ints)))
}

# Use it to split the comma-separated values and calculate the mean
b_split_meanint <- sapply(b_s_int$Intervals, intmean)
b_nosplit_meanint <- sapply(b_n_int$Intervals, intmean)
p_split_meanint <- sapply(p_s_int$Intervals, intmean)
p_nosplit_meanint <- sapply(p_n_int$Intervals, intmean)

# Make plots
pdf(file="B73_Tandem_Intervals_Split.pdf", 6, 6)
hist(
    b_nosplit_meanint,
    breaks=15,
    col=rgb(0, 0, 1, 0.5),
    xlab="Mean Number of Genes Between Tandem Duplicates",
    ylab="Count",
    main="Distribution of Mean Number of Intervening Genes\nfor B73 Tandem Duplicates",
    ylim=c(0, 1200))
hist(
    b_split_meanint,
    breaks=15,
    col=rgb(1, 0, 0, 0.5),
    add=TRUE)
legend(
    "topright",
    c("Not Split", "Split"),
    fill=c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5)))
dev.off()

pdf(file="PH207_Tandem_Intervals_Split.pdf", 6, 6)
hist(
    p_nosplit_meanint,
    breaks=20,
    col=rgb(0, 0, 1, 0.5),
    xlab="Mean Number of Genes Between Tandem Duplicates",
    ylab="Count",
    main="Distribution of Mean Number of Intervening Genes\nfor PH207 Tandem Duplicates",
    ylim=c(0, 2000))
hist(
    p_split_meanint,
    breaks=20,
    col=rgb(1, 0, 0, 0.5),
    add=TRUE)
legend(
    "topright",
    c("Not Split", "Split"),
    fill=c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5)))
dev.off()
