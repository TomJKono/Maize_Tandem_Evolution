# Script to plot the estiamted duplication dates of tandem clusters that are
# in cassettes

# Read the data files
b_cassettes <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Dating/B73_Cassette_Dates.txt", header=T)
p_cassettes <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Dating/PH207_Cassette_Dates.txt", header=T)

# Generate the ages to plot
b_age_range <- sapply(
    b_cassettes$Ages,
    function (x) {
    ages <- as.numeric(unlist(strsplit(as.character(x), ",")))
    ages <- ages[!is.na(ages)]
    return(max(ages) - min(ages))
    })
# Remove inf values
b_age_range <- b_age_range[is.finite(b_age_range)]

p_age_range <- sapply(
    p_cassettes$Ages,
    function (x) {
    ages <- as.numeric(unlist(strsplit(as.character(x), ",")))
    ages <- ages[!is.na(ages)]
    return(max(ages) - min(ages))
    })
# Remove inf values
p_age_range <- p_age_range[is.finite(p_age_range)]

# Plot them
pdf(file="Cassette_Ages.pdf", 3, 3)
# Set margins. The order is bottom, left, top, right
par(mar=c(4, 4, 0.1, 0.1), mgp=c(2, 1, 0))
plot(
    density(b_age_range),
    col="black",
    lwd=2,
    xlim=c(-2, 15),
    ylim=c(0, 0.11),
    xlab="Age Range (MYA)",
    ylab="Density",
    main="")
lines(density(p_age_range), col="black", lwd=2, lty=3)
legend(
    "topright",
    c("B73 Cassettes", "PH207 Cassettes"),
    col="black",
    lwd=2,
    lty=c(1, 3),
    cex=0.7)
dev.off()
