# Plot the distributions of syntenic tandem duplication ages.

syn_ages <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Dating/Syntenic_Duplicate_Ages.txt", header=TRUE)
nonsyn_ages <- read.table("/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Dating/Nonsyntenic_Duplicate_Ages.txt", header=TRUE)

# We want to be careful with how we identify shared/private duplications.
# For shared duplications, if tandems in B1/P1 or B2/P2 have the same ages,
# then we will assume they arose from the same duplication event. Because gene
# duplications are not expected to follow the infinite sites mutation model,
# we will return other duplications as private to their partitions. This will
# hopefully better reflect the age distributions of the tandem duplicates. We
# have to treat nonsyntenic duplications slightly differently from syntenic
# duplications, because there is no M1/M2 for nonsyntenic duplications.
parse_syntenic <- function(x) {
    # Get the ages of B1, B2, P1, P2
    b1_age <- as.numeric(unlist(strsplit(as.character(x["B1_Ages"]), ",")))
    b2_age <- as.numeric(unlist(strsplit(as.character(x["B2_Ages"]), ",")))
    p1_age <- as.numeric(unlist(strsplit(as.character(x["P1_Ages"]), ",")))
    p2_age <- as.numeric(unlist(strsplit(as.character(x["P2_Ages"]), ",")))
    # Remove NAs
    b1_age <- b1_age[!is.na(b1_age)]
    b2_age <- b2_age[!is.na(b2_age)]
    p1_age <- p1_age[!is.na(p1_age)]
    p2_age <- p2_age[!is.na(p2_age)]
    # Identify any shared ages among them
    m1_shared <- intersect(b1_age, p1_age)
    m2_shared <- intersect(b2_age, p2_age)
    all_shared <- intersect(m1_shared, m2_shared)
    # Count how many shared duplications there are.
    n_all <- unlist(sapply(all_shared, function(i) {
        b1 <- sum(b1_age == i)
        b2 <- sum(b2_age == i)
        p1 <- sum(p2_age == i)
        p2 <- sum(p2_age == i)
        minshare <- min(c(b1, b2, p1, p2))
        return(minshare)}))
    n_m1 <- unlist(sapply(m1_shared, function(i) {
        b <- sum(b1_age == i)
        p <- sum(p1_age == i)
        minshare <- min(b, p)
        return(minshare)}))
    n_m2 <- unlist(sapply(m2_shared, function(i) {
        b <- sum(b2_age == i)
        p <- sum(p2_age == i)
        minshare <- min(b, p)
        return(minshare)}))
    # Identify ages that are private to B1, B2, P1, P2
    b1_priv <- b1_age[!(b1_age %in% m1_shared)]
    p1_priv <- p1_age[!(p1_age %in% m1_shared)]
    b2_priv <- b2_age[!(b2_age %in% m2_shared)]
    p2_priv <- p2_age[!(p2_age %in% m2_shared)]
    # Return these seven categories
    ret <- list(
        All_Shared=rep(all_shared, n_all),
        M1_Shared=rep(m1_shared, n_m1),
        M2_Shared=rep(m2_shared, n_m2),
        B1_Private=b1_priv,
        B2_Private=b2_priv,
        P1_Private=p1_priv,
        P2_Private=p2_priv
        )
    return(ret)
}

# Define a function for similar counting of the nonsyntenic clusters. We only
# have three categories though: Private to B73, Private to PH207, shared.
parse_nonsyntenic <- function(x) {
    b_age <- as.numeric(unlist(strsplit(as.character(x["B_Ages"]), ",")))
    p_age <- as.numeric(unlist(strsplit(as.character(x["P_Ages"]), ",")))
    # Remove NA
    b_age <- b_age[!is.na(b_age)]
    p_age <- p_age[!is.na(p_age)]
    # Identify shared duplications
    shared <- intersect(b_age, p_age)
    # Identify private ages
    b_priv <- b_age[!(b_age %in% shared)]
    p_priv <- p_age[!(p_age %in% shared)]
    # Return the three categories
    ret <- list(
        Shared=shared,
        B_Private=b_priv,
        P_Private=p_priv)
    return(ret)
}

s_ages <- apply(syn_ages, 1, parse_syntenic)
n_ages <- apply(nonsyn_ages, 1, parse_nonsyntenic)

# Get the values for each of the seven categories
all_shared <- c()
m1_shared <- c()
m2_shared <- c()
b1_priv <- c()
b2_priv <- c()
p1_priv <- c()
p2_priv <- c()
ns_shared <- c()
ns_b_priv <- c()
ns_p_priv <- c()
for(i in 1:length(s_ages)) {
    all_shared <- c(all_shared, s_ages[[i]]$All_Shared)
    m1_shared <- c(m1_shared, s_ages[[i]]$M1_Shared)
    m2_shared <- c(m2_shared, s_ages[[i]]$M2_Shared)
    b1_priv <- c(b1_priv, s_ages[[i]]$B1_Private)
    b2_priv <- c(b2_priv, s_ages[[i]]$B2_Private)
    p1_priv <- c(p1_priv, s_ages[[i]]$P1_Private)
    p2_priv <- c(p2_priv, s_ages[[i]]$P2_Private)
}
for(i in 1:length(n_ages)) {
    ns_shared <- c(ns_shared, n_ages[[i]]$Shared)
    ns_b_priv <- c(ns_b_priv, n_ages[[i]]$B_Private)
    ns_p_priv <- c(ns_p_priv, n_ages[[i]]$P_Private)
}

length(all_shared)
length(m1_shared)
length(m2_shared)
length(b1_priv)
length(b2_priv)
length(p1_priv)
length(p2_priv)
length(ns_shared)
length(ns_b_priv)
length(ns_p_priv)


# Plot all 10 lines. A big, ugly plot.
pdf(file="Tandem_Ages.pdf", height=6, width=6)
plot(
    density(all_shared),
    xlim=c(-2, 15),
    ylim=c(0, 0.45),
    main="Distribution of Estimated Tandem Duplicate Ages",
    xlab="Estimated Age (MYA)",
    ylab="Density",
    col="#984ea3",
    lwd=2)
lines(density(m1_shared), col="#4daf4a", lwd=2)
lines(density(m2_shared), col="#9fff9c", lwd=2)
lines(density(b1_priv), col="#377eb8", lwd=2)
lines(density(b2_priv), col="#a9affa", lwd=2)
lines(density(p1_priv), col="#e41a1c", lwd=2)
lines(density(p2_priv), col="#ff9c9f", lwd=2)
lines(density(ns_shared), col="#1d6f1a", lwd=2, lty=2)
lines(density(ns_b_priv), col="#173e68", lwd=2, lty=2)
lines(density(ns_p_priv), col="#640000", lwd=2, lty=2)
legend(
    "topleft",
    c("All Shared", "Subgenome 1 Shared", "Subgenome 2 Shared", "B1 Private", "B2 Private", "P1 Private", "P2 Private", "Shared Nonsyntenic", "B Priv. Nonsyntenic", "P Priv. Nonsyntenic"),
    col=c("#984ea3", "#4daf4a", "#9fff9c", "#377eb8", "#a9affa", "#e41a1c", "#ff9c9f", "#1d6f1a", "#173e68", "#640000"),
    lwd=2,
    lty=c(1, 1, 1, 1, 1, 1, 1, 2, 2, 2),
    cex=0.75)
dev.off()

# Make a nice plot with just five lines. This will be in the paper:
#   Subgenome 1, shared B/P
#   Subgenome 2, shared B/P
#   Subgenome 1, private B or P
#   Subgenome 2, private B or P
#   Nonsyntenic
subg1_priv <- c(b1_priv, p1_priv)
subg2_priv <- c(b2_priv, p2_priv)
ns <- c(ns_b_priv, ns_p_priv, ns_shared)
pdf(file="Tandem_Ages_Clean_Pub.pdf", height=3, width=3)
# Set margins. The order is bottom, left, top, right
par(mar=c(4, 4, 0.1, 0.1), mgp=c(2, 1, 0))
plot(
    density(m1_shared),
    xlim=c(-2, 15),
    ylim=c(0, 0.35),
    main="",
    xlab="Estimated Age (MYA)",
    ylab="Density",
    col="black",
    lwd=2)
lines(density(m2_shared), col="darkgrey", lwd=2)
lines(density(subg1_priv), col="black", lwd=2, lty=3)
lines(density(subg2_priv), col="darkgrey", lwd=2, lty=3)
lines(density(ns), col="red", lwd=2)
legend(
    "topleft",
    c("Maize 1, Shared", "Maize 2, Shared", "Maize 1, Private", "Maize 2, Private", "Nonsyntenic"),
    col=c("black", "darkgrey", "black", "darkgrey", "red"),
    lwd=2,
    lty=c(1, 1, 3, 3, 1),
    cex=0.7
    )
dev.off()
