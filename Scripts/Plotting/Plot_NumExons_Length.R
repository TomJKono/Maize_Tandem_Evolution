# Make plots of the number of exons and the length of the gene in bp for both
# tandem duplicates and genomewide genes

# Read the data files
b_genome <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/TEs/B73_FullGene_Genomewide_TE_NumExons.txt", header=TRUE)
p_genome <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/TEs/PH207_FullGene_Genomewide_TE_NumExons.txt", header=TRUE)

# Read the clusters. We are only interested in the gene IDs
b_tandem <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Filtering/B73_True_Tandem_Clusters.txt", header=FALSE)$V2
p_tandem <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Filtering/PH207_True_Tandem_Clusters.txt", header=FALSE)$V2

# Define colors for B73 and PH207 for consistency
#   This is a nice blue that is print-friendly
b73_nontandem <- "#999999"
b73_tandem <- rgb(2/256, 112/256, 189/256)
#   And a nice red that is print-friendly
ph207_nontandem <- "#999999"
ph207_tandem <- rgb(237/256, 28/256, 36/256)

# Convert the comma-separated values into a character vector
b_tandem <- unlist(strsplit(as.character(b_tandem), ","))
p_tandem <- unlist(strsplit(as.character(p_tandem), ","))

# Then, separate the values for tandems
b_tand <- b_genome[b_genome$GenesInCluster %in% b_tandem,]
p_tand <- p_genome[p_genome$GenesInCluster %in% p_tandem,]

# And the non-tandems
b_nontand <- b_genome[!(b_genome$GenesInCluster %in% b_tandem),]
p_nontand <- p_genome[!(p_genome$GenesInCluster %in% p_tandem),]

# Finally, separate by multi- and single- exon genes
b_tand_se <- b_tand[b_tand$NumExons == 1,]
b_tand_me <- b_tand[b_tand$NumExons > 1,]
p_tand_se <- p_tand[p_tand$NumExons == 1,]
p_tand_me <- p_tand[p_tand$NumExons > 1,]

b_nontand_se <- b_nontand[b_nontand$NumExons == 1,]
b_nontand_me <- b_nontand[b_nontand$NumExons > 1,]
p_nontand_se <- p_nontand[p_nontand$NumExons == 1,]
p_nontand_me <- p_nontand[p_nontand$NumExons > 1,]

# Put them into a data frame to plot
toplot <- data.frame(
    Value=c(b_nontand_se$Lengths, b_tand_se$Lengths, b_nontand_me$Lengths, b_tand_me$Lengths,
            p_nontand_se$Lengths, p_tand_se$Lengths, p_nontand_me$Lengths, p_tand_me$Lengths),
    Label=c(
        rep("B NT S", nrow(b_nontand_se)),
        rep("B T S", nrow(b_tand_se)),
        rep("B NT M", nrow(b_nontand_me)),
        rep("B T M", nrow(b_tand_me)),
        rep("P NT S", nrow(p_nontand_se)),
        rep("P T S", nrow(p_tand_se)),
        rep("P NT M", nrow(p_nontand_me)),
        rep("P T M", nrow(p_tand_me))
    ))

pdf(file="NumExons_v_Length.pdf", height=3, width=6.5)
# Set margins. The order is bottom, left, top, right
par(mar=c(4, 4, 0.1, 0.1), mgp=c(2, 1, 0))
boxplot(
    toplot$Value/1000 ~ toplot$Label,
    log="y",
    lwd=2,
    border=c(
        b73_nontandem, b73_tandem,
        b73_nontandem, b73_tandem,
        ph207_nontandem, ph207_tandem,
        ph207_nontandem, ph207_tandem),
    at=c(1, 2, 3, 4, 6, 7, 8, 9),
    pch=19,
    cex=0.5,
    las=2,
    ylab="Gene Length",
    axes=F,
    main=""
    )
axis(side=2, at=c(0.1, 1, 10, 100), labels=c("100bp", "1kb", "10kb", "100kb"), cex.axis=0.95)
axis(side=1, at=c(1.5, 3.5, 6.5, 8.5), 
    labels=c("Single-Exon", "Multi-Exon", "Single-Exon", "Multi-Exon"))
mtext(text=c("B73", "PH207"), side=1, line=3, at=c(2.5, 7.5))
box()
dev.off()
