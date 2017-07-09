# Script to generate distributions of pairwise similarity between tandems that
# are split across orthogroups or not

# Read the data files
b_div <- read.table("/Users/tomkono/Dropbox/Projects/Fractionation/Data/Tandem_Divergence/B73_MzOnly_Summary.txt", header=T)
p_div <- read.table("/Users/tomkono/Dropbox/Projects/Fractionation/Data/Tandem_Divergence/PH207_MzOnly_Summary.txt", header=T)
ogs <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Orthofinder/Cluster_Orthogroup_Numbers.txt", header=F)

# Remove the ".fasta" from the "locus" names
b_div$locus <- sapply(b_div$locus, function(x) {return(unlist(strsplit(as.character(x), ".", fixed=TRUE))[1])})
p_div$locus <- sapply(p_div$locus, function(x) {return(unlist(strsplit(as.character(x), ".", fixed=TRUE))[1])})

# Get the B73 and the PH207 clusters
b_ogs <- ogs[ogs$V1 == "B73",]
p_ogs <- ogs[ogs$V1 == "PH207",]

# Then get which clusters were split, and which were not
b_split <- as.character(b_ogs$V2[b_ogs$V4 > 1])
b_nosplit <- as.character(b_ogs$V2[b_ogs$V4 == 1])
p_split <- as.character(p_ogs$V2[p_ogs$V4 > 1])
p_nosplit <- as.character(p_ogs$V2[p_ogs$V4 == 1])

# Plot them
pdf(file="B73_Tandem_Splitting_lt50.pdf", height=6, width=6)
hist(
    b_div$ThetaPi[b_div$locus %in% b_split & b_div$ThetaPi <= 0.5],
    breaks=50,
    col=rgb(0, 0, 1, 0.5),
    main="Distribution of Pairwise Differences\nof B73 Tandem Duplicates",
    xlab="Pairwise Differences",
    ylab="Count",
    ylim=c(0, 150))
hist(
    b_div$ThetaPi[b_div$locus %in% b_nosplit & b_div$ThetaPi <= 0.5],
    breaks=50,
    add=TRUE,
    col=rgb(1, 0, 0, 0.5))
abline(v=0.25, col="green", lwd=2)
legend("topright",
    c("Split Across OGs", "Not Split"),
    fill=c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5))
    )
dev.off()

pdf(file="PH207_Tandem_Splitting_lt50.pdf", height=6, width=6)
hist(
    p_div$ThetaPi[p_div$locus %in% b_split & p_div$ThetaPi <= 0.5],
    breaks=50,
    col=rgb(0, 0, 1, 0.5),
    main="Distribution of Pairwise Differences\nof PH207 Tandem Duplicates",
    xlab="Pairwise Differences",
    ylab="Count",
    ylim=c(0, 30))
hist(
    p_div$ThetaPi[p_div$locus %in% b_nosplit & p_div$ThetaPi <= 0.5],
    breaks=50,
    add=TRUE,
    col=rgb(1, 0, 0, 0.5))
abline(v=0.25, col="green", lwd=2)
legend("topright",
    c("Split Across OGs", "Not Split"),
    fill=c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5))
    )
dev.off()
