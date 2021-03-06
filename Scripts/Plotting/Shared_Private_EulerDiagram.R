# Generate an Euler diagram that shows the shared/private relationships for
# both genotype and subgenome

library(eulerr)

#   This is a nice blue that is print-friendly
b73_color <- rgb(2/256, 112/256, 189/256)
#   And a nice red that is print-friendly
ph207_color <- rgb(237/256, 28/256, 36/256)

# These are counts that are derived from external Python scripts. We hard-code
# the values here. 
b1_p1_shared <- 534
b2_p2_shared <- 217
b1_priv <- 387
b2_priv <- 190
p1_priv <- 145
p2_priv <- 84
b_ns <- 115
p_ns <- 90
shared_ns <- 134
b_nohom <- 181
p_nohom <- 263

# In this vector, "A" refers to B73, and "B" is PH207. The private counts are
# just the letter, and "A&B" denotes sharing.
dat <- c(
    "A"=b_ns+b1_priv+b2_priv+b_nohom,
    "B"=p_ns+p1_priv+p2_priv+p_nohom,
    "A&B"=b1_p1_shared+b2_p2_shared+shared_ns)

# Plot it
eu <- euler(dat)
pdf(file="Genotype_Sharing.pdf", width=3, height=2.5)
# Set margins. The order is bottom, left, top, right
par(mar=c(0.1, 0.1, 0.1, 0.1))
plot(eu,
    labels=c("B73", "PH207"),
    main="",
    col=c(b73_color, ph207_color),
    border=c(b73_color, ph207_color),
    fill=c(b73_color, ph207_color),
    fill_opacity=0.1,
    lwd=2,
    counts=TRUE)
dev.off()
