# Make barplots of the distribution of cluster sizes in the two subgenomes.
# Combine data from B73 and PH207, because they are very similar.

c_sizes <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Filtering/Syntenic_Cluster_Assignments.txt", header=TRUE)

m1_sizes <- function(x) {
    b <- as.character(x["B1_Cluster"])
    p <- as.character(x["P1_Cluster"])
    b_genes <- as.character(x["B1_Genes"])
    p_genes <- as.character(x["P1_Genes"])
    if(is.na(b) && is.na(p)) {
        return(NA)
    }
    else {
        b_size <- length(unlist(strsplit(b_genes, ",")))
        p_size <- length(unlist(strsplit(p_genes, ",")))
        if(b_size == 1) {
            bdup <- NA
        }
        else {
            bdup <- b_size
        }
        if(p_size == 1) {
            pdup <- NA
        }
        else {
            pdup <- p_size
        }
        return(c(bdup, pdup))
    }
}

m2_sizes <- function(x) {
    b <- as.character(x["B2_Cluster"])
    p <- as.character(x["P2_Cluster"])
    b_genes <- as.character(x["B2_Genes"])
    p_genes <- as.character(x["P2_Genes"])
    if(is.na(b) && is.na(p)) {
        return(NA)
    }
    else {
        b_size <- length(unlist(strsplit(b_genes, ",")))
        p_size <- length(unlist(strsplit(p_genes, ",")))
        if(b_size == 1) {
            bdup <- NA
        }
        else {
            bdup <- b_size
        }
        if(p_size == 1) {
            pdup <- NA
        }
        else {
            pdup <- p_size
        }
        return(c(bdup, pdup))
    }
}

m1_dups <- unlist(apply(c_sizes, 1, m1_sizes))
m2_dups <- unlist(apply(c_sizes, 1, m2_sizes))

# Remove NAs
m1_dups <- m1_dups[!is.na(m1_dups)]
m2_dups <- m2_dups[!is.na(m2_dups)]

# Make a matrix for the plot
toplot <- matrix(0, ncol=2, nrow=10)
for(i in 2:10) {
    r <- i-1
    toplot[r, 1] <- sum(m1_dups == i)
    toplot[r, 2] <- sum(m2_dups == i)
}
toplot[10, 1] <- sum(m1_dups > 10)
toplot[10, 2] <- sum(m2_dups > 10)

# Then plot them
pdf(file="Syntenic_Cluster_Sizes.pdf", 6, 6)
at <- barplot(
    t(toplot),
    beside=TRUE,
    col=c("black", "grey"),
    xlab="N. Genes",
    ylab="Count",
    main="Size of Tandem Duplicate Clusters",
    axes=F)
axis(side=2)
axis(side=1, at=apply(at, 2, mean), labels=c(as.character(2:10), ">10"))
legend("topright", c("Maize1", "Maize2"), fill=c("black", "grey"))
dev.off()
