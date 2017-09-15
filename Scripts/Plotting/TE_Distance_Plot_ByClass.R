#   Plot the distributions of the distances to TEs of each of the major
#   annotated classes. LTRs and solo_LTRs are combined into the same plot.

# Read the genomewide distances
gw_ltr <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/TEs/B73_Genomewide_LTR_Proximity.txt", header=TRUE)
gw_line <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/TEs/B73_Genomewide_LINE_Proximity.txt", header=TRUE)
gw_sine <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/TEs/B73_Genomewide_SINE_Proximity.txt", header=TRUE)
gw_tir <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/TEs/B73_Genomewide_TIR_Proximity.txt", header=TRUE)
gw_heli <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/TEs/B73_Genomewide_Helitron_Proximity.txt", header=TRUE)

# Read the tandem distances
tandem_ltr <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/TEs/B73_LTR_Proximity.txt", header=TRUE)
tandem_line <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/TEs/B73_LINE_Proximity.txt", header=TRUE)
tandem_sine <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/TEs/B73_SINE_Proximity.txt", header=TRUE)
tandem_tir <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/TEs/B73_TIR_Proximity.txt", header=TRUE)
tandem_heli <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/TEs/B73_Helitron_Proximity.txt", header=TRUE)

# Upstream+Downstream is the total distance
gw_ltr_win <- gw_ltr$UpstreamBP + gw_ltr$DownstreamBP
gw_line_win <- gw_line$UpstreamBP + gw_line$DownstreamBP
gw_sine_win <- gw_sine$UpstreamBP + gw_sine$DownstreamBP
gw_tir_win <- gw_tir$UpstreamBP + gw_tir$DownstreamBP
gw_heli_win <- gw_heli$UpstreamBP + gw_heli$DownstreamBP

tandem_ltr_win <- tandem_ltr$UpstreamBP + tandem_ltr$DownstreamBP
tandem_line_win <- tandem_line$UpstreamBP + tandem_line$DownstreamBP
tandem_sine_win <- tandem_sine$UpstreamBP + tandem_sine$DownstreamBP
tandem_tir_win <- tandem_tir$UpstreamBP + tandem_tir$DownstreamBP
tandem_heli_win <- tandem_heli$UpstreamBP + tandem_heli$DownstreamBP

# Remove NA for genes on scaffolds or at the ends of chromosomes
gw_ltr_win <- gw_ltr_win[!is.na(gw_ltr_win)]
gw_line_win <- gw_line_win[!is.na(gw_line_win)]
gw_sine_win <- gw_sine_win[!is.na(gw_sine_win)]
gw_tir_win <- gw_tir_win[!is.na(gw_tir_win)]
gw_heli_win <- gw_heli_win[!is.na(gw_heli_win)]

tandem_ltr_win <- tandem_ltr_win[!is.na(tandem_ltr_win)]
tandem_line_win <- tandem_line_win[!is.na(tandem_line_win)]
tandem_sine_win <- tandem_sine_win[!is.na(tandem_sine_win)]
tandem_tir_win <- tandem_tir_win[!is.na(tandem_tir_win)]
tandem_heli_win <- tandem_heli_win[!is.na(tandem_heli_win)]

# Report the means
print(c(mean(gw_ltr_win), median(gw_ltr_win), mean(tandem_ltr_win), median(tandem_ltr_win)))
print(c(mean(gw_line_win), median(gw_line_win), mean(tandem_line_win), median(tandem_line_win)))
print(c(mean(gw_sine_win), median(gw_sine_win), mean(tandem_sine_win), median(tandem_sine_win)))
print(c(mean(gw_tir_win), median(gw_tir_win), mean(tandem_tir_win), median(tandem_tir_win)))
print(c(mean(gw_heli_win), median(gw_heli_win), mean(tandem_heli_win), median(tandem_heli_win)))

#   Put them into a data frame to plot
toplot <- data.frame(
    Label=c(
        rep("G LTR", length(gw_ltr_win)),
        rep("T LTR", length(tandem_ltr_win)),
        rep("G LINE", length(gw_line_win)),
        rep("T LINE", length(tandem_line_win)),
        rep("G SINE", length(gw_sine_win)),
        rep("T SINE", length(tandem_sine_win)),
        rep("G TIR", length(gw_tir_win)),
        rep("T TIR", length(tandem_tir_win)),
        rep("G Helitron", length(gw_heli_win)),
        rep("T Helitron", length(tandem_heli_win))),
    Value=c(
        gw_ltr_win,
        tandem_ltr_win,
        gw_line_win,
        tandem_line_win,
        gw_sine_win,
        tandem_sine_win,
        gw_tir_win,
        tandem_tir_win,
        gw_heli_win,
        tandem_heli_win)
)
#   Set the ordering
toplot$Label <- factor(toplot$Label, levels=c(
    "G LTR", "T LTR",
    "G LINE", "T LINE",
    "G SINE", "T SINE",
    "G TIR", "T TIR",
    "G Helitron", "T Helitron"))
#   Plot it
pdf(file="B73_True_Tandems_TE_Distances.pdf", width=6, height=6)
boxplot(
    toplot$Value/1000 ~ toplot$Label,
    log="y",
    lwd=1,
    border=c("black", "blue"),
    at=c(1, 2, 4, 5, 7, 8, 10, 11, 13, 14),
    pch=19,
    cex=0.5,
    las=2,
    ylab="Flanking Distance to Nearest TEs (kb)",
    axes=F,
    main="Distances From B73 Tandem Duplicates to Flanking TEs"
    )
axis(side=2, at=c(0.1, 2, 10, 100, 1000, 10000, 100000), labels=c("0.1", "2", "10", "100", "1000", "10000", "100000"), las=2)
axis(side=1, at=c(1.5, 4.5, 7.5, 10.5, 13.5), labels=c("LTR", "LINE", "SINE", "TIR", "Helitron"))
legend("topright", c("Genomewide", "Tandem Duplicates"), col=c("black", "blue"), lwd=1)
dev.off()

# Let's do the same, but with distance to closest TE
closest_te <- function(x) {
    up <- as.numeric(x["UpstreamBP"])
    down <- as.numeric(x["DownstreamBP"])
    closest <- min(c(up, down))
    return(closest)
}

gw_ltr$ClosestTE <- apply(gw_ltr, 1, closest_te)
gw_line$ClosestTE <- apply(gw_line, 1, closest_te)
gw_sine$ClosestTE <- apply(gw_sine, 1, closest_te)
gw_tir$ClosestTE <- apply(gw_tir, 1, closest_te)
gw_heli$ClosestTE <- apply(gw_heli, 1, closest_te)

tandem_ltr$ClosestTE <- apply(tandem_ltr, 1, closest_te)
tandem_line$ClosestTE <- apply(tandem_line, 1, closest_te)
tandem_sine$ClosestTE <- apply(tandem_sine, 1, closest_te)
tandem_tir$ClosestTE <- apply(tandem_tir, 1, closest_te)
tandem_heli$ClosestTE <- apply(tandem_heli, 1, closest_te)

gw_ltr_closest <- gw_ltr$ClosestTE[!is.na(gw_ltr$ClosestTE)]
gw_line_closest <- gw_line$ClosestTE[!is.na(gw_line$ClosestTE)]
gw_sine_closest <- gw_sine$ClosestTE[!is.na(gw_sine$ClosestTE)]
gw_tir_closest <- gw_tir$ClosestTE[!is.na(gw_tir$ClosestTE)]
gw_heli_closest <- gw_heli$ClosestTE[!is.na(gw_heli$ClosestTE)]
tandem_ltr_closest <- tandem_ltr$ClosestTE[!is.na(tandem_ltr$ClosestTE)]
tandem_line_closest <- tandem_line$ClosestTE[!is.na(tandem_line$ClosestTE)]
tandem_sine_closest <- tandem_sine$ClosestTE[!is.na(tandem_sine$ClosestTE)]
tandem_tir_closest <- tandem_tir$ClosestTE[!is.na(tandem_tir$ClosestTE)]
tandem_heli_closest <- tandem_heli$ClosestTE[!is.na(tandem_heli$ClosestTE)]

print(c(mean(gw_ltr_closest), median(gw_ltr_closest), mean(tandem_ltr_closest), median(tandem_ltr_closest)))
print(c(mean(gw_line_closest), median(gw_line_closest), mean(tandem_line_closest), median(tandem_line_closest)))
print(c(mean(gw_sine_closest), median(gw_sine_closest), mean(tandem_sine_closest), median(tandem_sine_closest)))
print(c(mean(gw_tir_closest), median(gw_tir_closest), mean(tandem_tir_closest), median(tandem_tir_closest)))
print(c(mean(gw_heli_closest), median(gw_heli_closest), mean(tandem_heli_closest), median(tandem_heli_closest)))

toplot <- data.frame(
    Label=c(
        rep("G LTR", length(gw_ltr_closest)),
        rep("T LTR", length(tandem_ltr_closest)),
        rep("G LINE", length(gw_line_closest)),
        rep("T LINE", length(tandem_line_closest)),
        rep("G SINE", length(gw_sine_closest)),
        rep("T SINE", length(tandem_sine_closest)),
        rep("G TIR", length(gw_tir_closest)),
        rep("T TIR", length(tandem_tir_closest)),
        rep("G Helitron", length(gw_heli_closest)),
        rep("T Helitron", length(tandem_heli_closest))),
    Value=c(
        gw_ltr_closest,
        tandem_ltr_closest,
        gw_line_closest,
        tandem_line_closest,
        gw_sine_closest,
        tandem_sine_closest,
        gw_tir_closest,
        tandem_tir_closest,
        gw_heli_closest,
        tandem_heli_closest)
    )
#   Set the ordering
toplot$Label <- factor(toplot$Label, levels=c(
    "G LTR", "T LTR",
    "G LINE", "T LINE",
    "G SINE", "T SINE",
    "G TIR", "T TIR",
    "G Helitron", "T Helitron"))
#   Plot it
pdf(file="B73_True_Tandems_Closest_TE_Distances.pdf", width=6, height=6)
boxplot(
    toplot$Value/1000 ~ toplot$Label,
    log="y",
    lwd=1,
    border=c("black", "blue"),
    at=c(1, 2, 4, 5, 7, 8, 10, 11, 13, 14),
    pch=19,
    cex=0.5,
    las=2,
    ylab="Distance to nearest TE (kb)",
    axes=F,
    main="Distances From B73 Tandem Duplicates to Nearest TE"
    )
axis(side=2, at=c(0.001, 0.1, 2, 10, 100, 1000, 10000, 50000), labels=c("0.001", "0.1", "2", "10", "100", "1000", "10000", "50000"), las=2)
axis(side=1, at=c(1.5, 4.5, 7.5, 10.5, 13.5), labels=c("LTR", "LINE", "SINE", "TIR", "Helitron"))
legend("topright", c("Genomewide", "Tandem Duplicates"), col=c("black", "blue"), lwd=1)
dev.off()
