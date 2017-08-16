#   Make a big plot of the genome showing gene density, TE density, synteny
#   with Sorghum/rice. Pericentromeres will come later, when annotated. This
#   is all for B73v4.

#   Make a poor decision, use ggplot2.
library(ggplot2)
library(reshape2)

#   Read data files
gene_te <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Data/B73_GeneDensity_TEDensity.txt", header=TRUE)
syntenic_blocks <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Data/B73vSbandOs_merged_blocks_final.txt", header=FALSE)
tandem_pos <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Data/B73_True_Tandem_Positions.gff", header=FALSE, sep="\t")

#   Trim down the tandem positions to just the chromosome and start
tandem_pos <- tandem_pos[,c("V1", "V4")]
names(tandem_pos) <- c("Chromosome", "Start")
#   Remove those not on chromosomes 1:10
tandem_pos <- tandem_pos[tandem_pos$Chromosome %in% as.character(1:10),]
#   Get maize1 and maize2 blocks
m1_blocks <- syntenic_blocks[syntenic_blocks$V4 == "maize1", c(1, 2, 3)]
m2_blocks <- syntenic_blocks[syntenic_blocks$V4 == "maize2", c(1, 2, 3)]

#   Set the names
names(m1_blocks) <- c("Chromosome", "Start", "End")
names(m2_blocks) <- c("Chromosome", "Start", "End")

m1_blocks$Chromosome <- factor(m1_blocks$Chromosome, levels=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
m2_blocks$Chromosome <- factor(m2_blocks$Chromosome, levels=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
gene_te$Chromosome <- factor(gene_te$Chromosome, levels=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
tandem_pos$Chromosome <- factor(tandem_pos$Chromosome, levels=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))

#   Make a plot
pdf(file="B73_True_Tandem_Dist.pdf", width=6.5, height=8)
p <- ggplot(NULL) +
    geom_rect(data=m1_blocks, aes(xmin=Start/1000000, xmax=End/1000000, ymin=0, ymax=100), fill="#b2df8a") +
    geom_rect(data=m2_blocks, aes(xmin=Start/1000000, xmax=End/1000000, ymin=0, ymax=100), fill="#a6cee3") +
    geom_line(data=gene_te, aes(x=Midpoint/1000000, y=GenesPerMb), color="black", size=0.5) +
    geom_line(data=gene_te, aes(x=Midpoint/1000000, y=DNATEsPerMb), color="#1f78b4", size=0.5) +
    geom_line(data=gene_te, aes(x=Midpoint/1000000, y=RNATEsPerMb), color="#33a02c", size=0.5) +
    geom_segment(data=tandem_pos, mapping=aes(x=Start/1000000, xend=Start/1000000, y=100, yend=80), color="#984ea3", size=0.25) +
    scale_y_continuous(name="Features Per Mb", breaks=seq(0, 90, 15), labels=seq(0, 90, 15)) +
    scale_x_continuous(name="Physical Position (Mb)", breaks=seq(0, 350, 50), labels=seq(0, 350, 50)) +
    theme_bw() +
    theme(
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        strip.background=element_blank(),
        strip.text.x=element_text(size=14)
        ) +
    facet_wrap(~Chromosome, nrow=5)
p
dev.off()
