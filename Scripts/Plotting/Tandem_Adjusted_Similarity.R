# Plot an "adjusted" similarity of tandem duplicates. Basically, it scales the
# pairwise similarity value by the proportion of the alignment that is gaps

b_genome <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Divergence/B73_Genomewide_Stats.txt", header=TRUE)
p_genome <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Divergence/PH207_Genomewide_Stats.txt", header=TRUE)
b_tandem <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Divergence/B73_Tandem_Alignment_Stats.txt", header=TRUE)
p_tandem <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Divergence/PH207_Tandem_Alignment_Stats.txt", header=TRUE)

# Calculate the adjusted similarities
bg_sim <- 1 - b_genome$ThetaPi
bg_ug <- b_genome$nsites_ug/b_genome$nsites
bg_adj <- bg_sim*bg_ug

pg_sim <- 1 - p_genome$ThetaPi
pg_ug <- p_genome$nsites_ug/p_genome$nsites
pg_adj <- pg_sim * pg_ug

bt_sim <- 1 - b_tandem$ThetaPi
bt_ug <- b_tandem$nsites_ug/b_tandem$nsites
bt_adj <- bt_sim*bt_ug

pt_sim <- 1 - p_tandem$ThetaPi
pt_ug <- p_tandem$nsites_ug/p_tandem$nsites
pt_adj <- pt_sim*pt_ug

# And plot the densities of them
pdf(file="Tandem_Adjusted_Similarity.pdf", 6, 6)
# Start with the genomewide plot
plot(
    density(bg_adj),
    col="darkblue",
    lwd=2,
    main="Adjusted Pairwise Similarity",
    xlab="Pairwise Similarity * Proportion Ungapped",
    ylab="Density")
lines(
    density(pg_adj),
    col="darkred",
    lwd=2)
lines(
    density(bt_adj),
    col="blue",
    lwd=2,
    lty=2)
lines(
    density(pt_adj),
    col="red",
    lwd=2,
    lty=2)
abline(
    v=0.30,
    col="darkgreen",
    lwd=2)
legend(
    "topright",
    c("B73 Genome-wide", "PH207 Genome-wide", "B73 Tandem", "PH207 Tandem"),
    col=c("darkblue", "darkred", "blue", "red"),
    lwd=2,
    lty=c(1, 1, 2, 2))
dev.off()

# And plot the densities of them, above 30%
pdf(file="Tandem_Adjusted_Similarity_Gt30perc.pdf", 6, 6)
# Start with the genomewide plot
plot(
    density(bg_adj[bg_adj >= 0.3]),
    col="darkblue",
    lwd=2,
    main="Adjusted Pairwise Similarity",
    xlab="Pairwise Similarity * Proportion Ungapped",
    ylab="Density",
    xlim=c(0, 1),
    ylim=c(0, 3))
lines(
    density(pg_adj[pg_adj >= 0.3]),
    col="darkred",
    lwd=2)
lines(
    density(bt_adj[bt_adj >= 0.3]),
    col="blue",
    lwd=2,
    lty=2)
lines(
    density(pt_adj[pt_adj >= 0.3]),
    col="red",
    lwd=2,
    lty=2)
legend(
    "topright",
    c("B73 Genome-wide", "PH207 Genome-wide", "B73 Tandem", "PH207 Tandem"),
    col=c("darkblue", "darkred", "blue", "red"),
    lwd=2,
    lty=c(1, 1, 2, 2))
dev.off()
# Plot the the relationships between adjusted similarity and just the estimated
# similarity
# pdf(file="B73_AdjustedSim_v_PairwiseSim.pdf", 6, 6)
# plot(
#     bt_sim ~ bt_adj,
#     pch=19,
#     cex=0.5,
#     col=rgb(0, 0, 1, 0.5),
#     xlim=c(0, 1),
#     ylim=c(0, 1),
#     main="Relationship Between Adjusted and\nRaw Pairwise Similarity for B73 Tandems",
#     xlab="Pairwise Similarity * Proportion Ungapped",
#     ylab="Pairwise Similarity"
#     )
# dev.off()

# pdf(file="PH207_AdjustedSim_v_PairwiseSim.pdf", 6, 6)
# plot(
#     pt_sim ~ pt_adj,
#     pch=19,
#     cex=0.5,
#     col=rgb(1, 0, 0, 0.5),
#     xlim=c(0, 1),
#     ylim=c(0, 1),
#     main="Relationship Between Adjusted and\nRaw Pairwise Similarity for PH207 Tandems",
#     xlab="Pairwise Similarity * Proportion Ungapped",
#     ylab="Pairwise Similarity"
#     )
# dev.off()

# # Plot the pairwise sim and ungapped against each other
# pdf(file="B73_PropUngap_AdjSimilarity.pdf", 6, 6)
# plot(
#     bt_ug ~ bt_adj,
#     pch=19,
#     cex=0.5,
#     col=rgb(0, 0, 1, 0.5),
#     xlim=c(0, 1),
#     ylim=c(0, 1),
#     main="Relationship Between Proportion of Ungapped Sites\nand Similarity for B73 Tandems",
#     xlab="Pairwise Similarity * Proportion Ungapped",
#     ylab="Proportion Ungapped Sites")
# dev.off()

# pdf(file="PH207_PropUngap_AdjSimilarity.pdf", 6, 6)
# plot(
#     pt_ug ~ pt_adj,
#     pch=19,
#     cex=0.5,
#     col=rgb(1, 0, 0, 0.5),
#     xlim=c(0, 1),
#     ylim=c(0, 1),
#     main="Relationship Between Proportion of Ungapped Sites\nand Similarity for PH207 Tandems",
#     xlab="Pairwise Similarity * Proportion Ungapped",
#     ylab="Proportion Ungapped Sites")
# dev.off()
