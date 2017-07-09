#!/usr/bin/env python
"""Script to parse the pairwise comparisons of CoGe-identified tandem duplicate
gene clusters, and identify a set of tandem duplicate clusters. From inspection
of the distribution of adjacent gene similarities, we have determined that 30%
is an appropriate threshold above which to call a pair of genes true tandem
duplicates. Takes two arguments:
    1) Output from 'compute' run on all tandem pairwise alignments
    2) Same output from all genomewide adjacent genes
"""

import sys


def adjusted_similarity(pi, ungapped, total):
    """Calculate the adjusted pairwise similarity, based on the summary
    statistics from 'compute.' We calculate the adjusted pairwise similarity
    as a scaling of pairwise similarity: (1 - theta_pi) * (sites_ug / tot_sites)
    which down-weights by the proportion of gapped sites."""
    similarity = 1 - pi
    prop_ungapped = ungapped / total
    return similarity * prop_ungapped


def identify_links(summary):
    """Read the table from 'compute' and return tuples that give the
    relationships that pass the 30% adjusted similarity threshold."""
    links = []
    with open(summary, 'r') as f:
        for index, line in enumerate(f):
            if index <= 3:
                continue
            else:
                tmp = line.strip().split()
                # The first field has the gene names in it
                g1, g2 = tmp[0].split('.')[0].split('-')
                # Cast everything to float so we don't have truncation
                tot_sites = float(tmp[2])
                ug_sites = float(tmp[3])
                difference = float(tmp[12])
                # Calculate the adjusted similarity
                adj_sim = adjusted_similarity(difference, ug_sites, tot_sites)
                if adj_sim < 0.3:
                    continue
                else:
                    links.append((g1, g2))
    return links


def find_clusters(edges):
    """Using the set of all edges, identify clusters. The unique set of nodes
    can be built from the edges. A cluster is a group of genes with high
    similarity to each other, but low similarity to others."""

    # Define a helper function to return a list of connected genes
    def connected(g, e):
        """Given a gene and a list of edges, return a list of genes that are
        connected to gene g."""
        neighbors = []
        for pair in e:
            if g not in pair:
                continue
            else:
                # Convert it to a set and throw out the query gene
                p = set(pair)
                p.discard(g)
                # Get the element back from the set
                n = list(p)[0]
                neighbors.append(n)
        return set(neighbors)

    # Get a set of genes (nodes)
    all_nodes = set([a[0] for a in edges] + [a[1] for a in edges])
    # start a list of clusters
    clusters = []
    # Then iterate through the genes. For each gene, get the other genes that
    # it's connected to. For each of the connected genes, get a list of its
    # connected genes, excluding those we've already seen. When the list of
    # new connections is 0, we are done
    while all_nodes:
        # Grab a gene off the list of all genes
        curr_gene = all_nodes.pop()
        # Start a set for the current cluster
        curr_clust = {curr_gene}
        # Then, start a list of genes to search
        search = [curr_gene]
        while search:
            # Pop off a gene
            n = search.pop()
            # Find its connected genes
            c = connected(n, edges)
            # Remove the connections that we've already seen
            c.difference_update(curr_clust)
            # and remove these from the set of genes to query
            all_nodes.difference_update(c)
            # Add them to the current cluster
            curr_clust.update(c)
            # Then add any non-visited connectons back to the list of things
            # to search
            search.extend(c)
        # Add the cluster to the list of clusters
        clusters.append(curr_clust)
    return clusters


def main(summary, genomewide):
    """Main function."""
    l = identify_links(summary)
    g = identify_links(genomewide)
    # Combine the tandems with the genomewide
    candidates = list(set(l+g))
    clst = find_clusters(candidates)
    # Then, convert the clusters to lists and sort them
    clst_gsorted = [sorted(list(c)) for c in clst]
    # And sort by the leftmost gene
    clst_sorted = sorted(clst_gsorted, key=lambda x: x[0])
    # Then, print out the cluster number and a list of genes to go with it
    for index, genes in enumerate(clst_sorted):
        # Convert the index to a string, then pad it with zeroes. We have on the
        # order of thousands of clusters, so we pad to 4 digits.
        s_i = str(index)
        cluster_num = 'Cluster_' + s_i.zfill(4)
        print cluster_num + '\t' + ','.join(genes)
    return


if len(sys.argv) != 3:
    print """
Script to parse the pairwise comparisons of CoGe-identified tandem duplicate
gene clusters, and identify a set of tandem duplicate clusters. From inspection
of the distribution of adjacent gene similarities, we have determined that 30%
is an appropriate threshold above which to call a pair of genes true tandem
duplicates. Takes two arguments:
    1) Output from 'compute' run on all tandem pairwise alignments
    2) Same output from all genomewide adjacent genes"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2])
