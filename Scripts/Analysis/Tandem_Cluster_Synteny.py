#!/usr/bin/env python
"""Generate a table of homologous genes/tandem duplicate clusters using the
assignment table from ABB. Will not find non-syntenic clusters. Takes three
arguments:
    1) B73 clusters
    2) PH207 clusters
    3) Synteny table
"""

import sys


def parse_clusters(clst):
    """Parse the cluster assignment, and store as a dictionary."""
    clust = {}
    with open(clst, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            c_id = tmp[0]
            genes = tmp[1].split(',')
            for g in genes:
                clust[g] = c_id
    return clust


def parse_synteny(syn):
    """Parse the synteny file and return it as a list. we do it this way
    because we will have to look up by either of the B73 or PH207 columns, and
    a dictionary will not help us here."""
    s = []
    with open(syn, 'r') as f:
        for line in f:
            s.append(line.strip().split()[0:5])
    # sort on ancestral gene, just because
    s = sorted(s, key=lambda x: x[1])
    return s


def homologous_dups(bdup, pdup, syn):
    """Using the synteny table, identify clusters of genes that are homologous.
    Return them as a list of lists."""
    hom = []
    # Iterate over syntenic relationships
    c_num = 0
    for rel in syn:
        anc = rel[0]
        b1 = rel[1]
        b2 = rel[2]
        p1 = rel[3]
        p2 = rel[4]
        # Get the clusters from the dictionary. If they are no present, drop
        # in a NA instead
        b1_c = bdup.get(b1, 'NA')
        b2_c = bdup.get(b2, 'NA')
        p1_c = pdup.get(p1, 'NA')
        p2_c = pdup.get(p2, 'NA')
        # Append these to the homology list
        hom.append([anc, b1_c, b2_c, p1_c, p2_c, b1, b2, p1, p2])
    return hom


def main(b_clusters, p_clusters, synteny):
    """Main function."""
    b_c = parse_clusters(b_clusters)
    p_c = parse_clusters(p_clusters)
    syn = parse_synteny(synteny)
    homology = homologous_dups(b_c, p_c, syn)
    # Print a header
    print 'Ancestral\tB1_Cluster\tB2_Cluster\tP1_Cluster\tP2_Cluster\tB1_Gene\tB2_Gene\tP1_Gene\tP2_Gene'
    for row in homology:
        print '\t'.join(row)
    return


if len(sys.argv) != 4:
    print """Generate a table of homologous genes/tandem duplicate clusters using the
assignment table from ABB. Will not find non-syntenic clusters. Takes three
arguments:
    1) B73 clusters
    2) PH207 clusters
    3) Synteny table"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2], sys.argv[3])
