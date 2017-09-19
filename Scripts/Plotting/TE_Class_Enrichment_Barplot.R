# Make a barplot that shows the proportion of TE classes present in the full
# annotated gene sequences of tandem and non-tandem duplicated genes

# Read the data tables
b_gw_tes <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/TEs/B73_FullGene_Genomewide_TE_Summary.txt", header=TRUE)
b_tand_tes <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/TEs/B73_FullGene_Tandem_TE_Summary.txt", header=TRUE)
p_gw_tes <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/TEs/PH207_FullGene_Genomewide_TE_Summary.txt", header=TRUE)
p_tand_tes <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/TEs/PH207_FullGene_Tandem_TE_Summary.txt", header=TRUE)

# A bit of an ugly hard-coded step: B73 tandems do not have any LINE elements
# in their full sequences
toadd <- data.frame(
    TE_Type="LINE_element",
    Count=0)
b_tand_tes <- rbind(b_tand_tes, toadd)

# Also, no SINEs for some reason
nosine <- data.frame(
    TE_Type="SINE_element",
    Count=0)
b_gw_tes <- rbind(b_gw_tes, nosine)
b_tand_tes <- rbind(b_tand_tes, nosine)
p_gw_tes <- rbind(p_gw_tes, nosine)
p_tand_tes <- rbind(p_tand_tes, nosine)

# Next, calculate proportions
b_gw_tes$Prop <- b_gw_tes$Count / sum(b_gw_tes$Count)
b_tand_tes$Prop <- b_tand_tes$Count / sum(b_tand_tes$Count)
p_gw_tes$Prop <- p_gw_tes$Count / sum(p_gw_tes$Count)
p_tand_tes$Prop <- p_tand_tes$Count / sum(p_tand_tes$Count)

# And make 0 values infinitesimal for the log scale
b_gw_tes$Prop[b_gw_tes$Prop == 0] <- 1e-16
b_tand_tes$Prop[b_tand_tes$Prop == 0] <- 1e-16
p_gw_tes$Prop[p_gw_tes$Prop == 0] <- 1e-16
p_tand_tes$Prop[p_tand_tes$Prop == 0] <- 1e-16

# Put them in a matrix to plot. Rows will be TE types and columns will be
# partitions
toplot <- matrix(c(
    b_gw_tes$Prop[b_gw_tes$TE_Type == "terminal_inverted_repeat_element"],
    b_tand_tes$Prop[b_tand_tes$TE_Type == "terminal_inverted_repeat_element"],
    p_gw_tes$Prop[p_gw_tes$TE_Type == "terminal_inverted_repeat_element"],
    p_tand_tes$Prop[p_tand_tes$TE_Type == "terminal_inverted_repeat_element"],

    b_gw_tes$Prop[b_gw_tes$TE_Type == "helitron"],
    b_tand_tes$Prop[b_tand_tes$TE_Type == "helitron"],
    p_gw_tes$Prop[p_gw_tes$TE_Type == "helitron"],
    p_tand_tes$Prop[p_tand_tes$TE_Type == "helitron"],

    b_gw_tes$Prop[b_gw_tes$TE_Type == "LTR_retrotransposon"] + b_gw_tes$Prop[b_gw_tes$TE_Type == "solo_LTR"],
    b_tand_tes$Prop[b_tand_tes$TE_Type == "LTR_retrotransposon"] + b_tand_tes$Prop[b_tand_tes$TE_Type == "solo_LTR"],
    p_gw_tes$Prop[p_gw_tes$TE_Type == "LTR_retrotransposon"] + p_gw_tes$Prop[p_gw_tes$TE_Type == "solo_LTR"],
    p_tand_tes$Prop[p_tand_tes$TE_Type == "LTR_retrotransposon"] + p_tand_tes$Prop[p_tand_tes$TE_Type == "solo_LTR"],

    b_gw_tes$Prop[b_gw_tes$TE_Type == "LINE_element"],
    b_tand_tes$Prop[b_tand_tes$TE_Type == "LINE_element"],
    p_gw_tes$Prop[p_gw_tes$TE_Type == "LINE_element"],
    p_tand_tes$Prop[p_tand_tes$TE_Type == "LINE_element"],

    b_gw_tes$Prop[b_gw_tes$TE_Type == "SINE_element"],
    b_tand_tes$Prop[b_tand_tes$TE_Type == "SINE_element"],
    p_gw_tes$Prop[p_gw_tes$TE_Type == "SINE_element"],
    p_tand_tes$Prop[p_tand_tes$TE_Type == "SINE_element"]),
    ncol=5)

pdf(file="TE_Class_Enrichment.pdf", 6, 6)
at <- barplot(
    toplot,
    axes=FALSE,
    beside=TRUE,
    log="y",
    ylim=c(1e-16, 1),
    xlab="TE Class",
    ylab="Proportion of All TEs Present",
    main="TE Class Enrichment in Full Gene Sequences",
    col=c("darkblue", "lightblue", "darkred", "red"))
axis(side=2)
axis(side=1, at=apply(at, 2, mean), labels=c("TIR", "Helitron", "LTR", "LINE", "SINE"))
legend("topright", c("B73 Genomewide", "B73 Tandem", "PH207 Genomewide", "PH207 Tandem"), fill=c("darkblue", "lightblue", "darkred", "red"), cex=0.9)
dev.off()
