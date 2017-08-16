#!/usr/bin/env python
"""Identify the homologous tandem dupicates between B73 and PH207, using the
homologous assignments from ABB's work. Takes three arguments:
    1) "Working gene list" from ABB
    2) B73 tandem duplicate clusters
    3) PH207 tandem duplicate clusters
"""

import sys
import pprint


def parse_clusters(clst_file):
    """Parse the cluster file and return a dictionary of the form
    {
        gene_id: cluster_num,
        gene_id: cluster_num,
        ...
    }"""
    c = {}
    with open(clst_file, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            for gene in tmp[1].split(','):
                c[gene] = tmp[0]
    return c


def main(homologues, b_clust, p_clust):
    """Main function."""
    b_c = parse_clusters(b_clust)
    p_c = parse_clusters(p_clust)
    cluster_rel = []
    with open(homologues, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            b_gene = tmp[3]
            p_gene = tmp[7]
            if b_gene in b_c:
                b_clst = b_c[b_gene]
            else:
                b_clst = None
            if p_gene in p_c:
                p_clst = p_c[p_gene]
            else:
                p_clst = None
            cluster_rel.append((b_clst, p_clst, b_gene, p_gene))
    cluster_rel = sorted(list(set(cluster_rel)))
    # Iterate through the cluster relationships and print out the table
    print 'B73_Cluster\tB73_Genes\tPH207_Cluster\tPH207_Genes'
    for rel in cluster_rel:
        b_cluster = rel[0]
        p_cluster = rel[1]
        if b_cluster:
            b_genes = sorted([key for key, val in b_c.iteritems() if val == b_cluster])
        else:
            b_cluster = 'NA'
            b_genes = [rel[2]]
        if p_cluster:
            p_genes = sorted([key for key, val in p_c.iteritems() if val == p_cluster])
        else:
            p_cluster = 'NA'
            p_genes = [rel[3]]
        # Skip lines that are not clusters in both genotypes
        if b_cluster == 'NA' and p_cluster == 'NA':
            continue
        else:
            print '\t'.join([b_cluster, ','.join(b_genes), p_cluster, ','.join(p_genes)])
    return


if len(sys.argv) != 4:
    print """
Identify the homologous tandem dupicates between B73 and PH207, using the
homologous assignments from ABB's work. Takes three arguments:
    1) "Working gene list" from ABB
    2) B73 tandem duplicate clusters
    3) PH207 tandem duplicate clusters"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2], sys.argv[3])
