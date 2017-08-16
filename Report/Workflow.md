# Maize Tandem Duplicate Evolution
This document describes the detailed procedure that was taken for the various
maize tandem duplicate evolution analyses and data handling.

## Overview
The general categories of the analysis fall into four categories:

- Tandem duplicate identification
- Tandem duplicate summaries
- Divergence date estimation
- Evolutionary hypothesis tests

They will be detailed in the following sections. Commands and scripts will be
given.

## Tandem Duplicate Identification
Tandem duplicates were identified in two ways:

- Filtering SynMap tandem duplicate output
- Sequence similarity of adjacent genes

SynMap identifies tandem duplicates in a slightly strange way: different genes
within a window that have BLAST hits to the same gene in the other genome. It
does not explicitly check sequence similarity of genes within the same genome,
so there is no guarantee that they are true tandem duplicates. It also
identifies tandem duplicates before it clusters BLAST matches into syntenic
blocks, so there will be nonsyntenic duplications identifed this way.

To refine the SynMap clusters and identify additional tandem duplicates, we
calculated an adjusted pairwise similarity metric for each pair of adjacent
genes in B73v4 and PH207v1. Adjusted pairwise simiarity is pairwise similarity,
down-weighted for the number of gaps opened in alignment of the two sequences.
The procedure for calculating this was as follows:

For each pair of adjacent genes:

1. Translate the CDS of the longest transcripts into amino acids.
2. Align them with ClustalOmega, with 10 iterations of refinement.
3. Back-translate the amino acid alignment into nucleotides.
4. Calculate the proportion of gaps and pairwise similarity (1 - `ThetaPi`)
   with `compute` from the `analysis` package. `compute` can interpret shell
   globbing characters (need to be in single quotes):

    ```bash
    compute -i '*.fasta' > Genomewide_Stats.txt
    ```

5. Calculate adjusted pairwise similarity as `(1 - ThetaPi) * (nsites_ug/nsites)`

This procedure was performed for both B73 and PH207. We ran it on both all
adjacent genes genome-wide, and for pairs of genes within SynMap putative
tandem duplicate clusters.

The alignment is implemented in `Scripts/Analysis/Adjacent_Gene_Similarity.py`,
and the plotting is in `Scripts/Plotting/Tandem_Adjusted_Similarity.R`.

A distribution of adjusted pairwise similarity for adjacent B73 genes, adjacent
PH207 genes, B73 putative tandem genes, and PH207 putative tandem genes is
shown below.

![Distribution of Adjusted Pairwise Similarity](Tandem_Adjusted_Similarity.png)

Based on this plot, putative tandem duplicate genes from SynMap and adjacent
genes with at least 0.3 adjusted similarity were kept as true tandem duplicates.

## Tandem Duplicate Summaries
The first summary of tandem duplicates is a raw count:

| Genotype | Maize1 Clusters | Maize2 Clusters | Nonsyntenic Clusters |
|----------|-----------------|-----------------|----------------------|
| B73      | 938             | 420             | 1,297*               |
| PH207    | 691             | 316             | 1,005*               |
*: Some clusters overlap syntenic and nonsyntenic genes. This happens when a
tandem duplicate happens at the edge of a syntenic block, and the cluster
includes both syntenic and nonsyntenic genes.

## Divergence Date Estimation
## Evolutionary Hypothesis Tests
