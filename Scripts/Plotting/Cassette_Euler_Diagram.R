# Make an Euler diagram for the sharing of cassettes between B73 and PH207.
# The values were calculated from Shared_Private_Cassettes.py, and were
# hard-coded into the script.

library(eulerr)

#   This is a nice blue that is print-friendly
b73_color <- rgb(2/256, 112/256, 189/256)
#   And a nice red that is print-friendly
ph207_color <- rgb(237/256, 28/256, 36/256)

# These are counts that are derived from external Python scripts. We hard-code
# the values here. 
shared <- 14
b_priv <- 44
p_priv <- 46

# In this vector, "A" refers to B73, and "B" is PH207. The private counts are
# just the letter, and "A&B" denotes sharing.
dat <- c(
    "A"=b_priv,
    "B"=p_priv,
    "A&B"=shared)

# Plot it
eu <- euler(dat)
pdf(file="Cluster_Sharing.pdf", width=3, height=2.5)
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
