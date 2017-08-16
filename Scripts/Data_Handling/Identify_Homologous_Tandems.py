#!/usr/bin/env python
"""Identify homologous tandem duplicate clusters between B73 and PH207. Takes
three arguments:
    1) Orthofinder B73-PH207 homologues file
    2) B73 tandem clusters
    3) PH207 tandem clusters

NOTE: Obsolete. Use homologous assignments from ABB instead, which made use of
syntenic assignments.
"""

import sys
import csv


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


def parse_homology(hom_file):
    """Parse the homologues file from Orthofinder. Return a list of tuples
    of strings with the following form
    [
        ('b1,b2,b3,...', 'p1,p2,p3,...'),
        ('bx,by,bz,...', 'px,py,pz,...'),
        ...
    ]"""
    hom = []
    with open(hom_file, 'r') as f:
        reader = csv.reader(f)
        for index, line in enumerate(reader):
            if index == 0:
                continue
            else:
                b = ','.join([
                    g.strip().split('_')[0] for g in line[1].split(',')])
                p = ','.join([
                    g.strip().split('_')[0] for g in line[2].split(',')])
                hom.append((b, p))
    return hom

def main(homologues, b_clust, p_clust):
    """Main function."""
    b_c = parse_clusters(b_clust)
    p_c = parse_clusters(p_clust)
    homology = parse_homology(homologues)
    for h in homology:
        print h
    return


if len(sys.argv) != 4:
    print """
Identify homologous tandem duplicate clusters between B73 and PH207. Takes
three arguments:
    1) Orthofinder B73-PH207 homologues file
    2) B73 tandem clusters
    3) PH207 tandem clusters"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2], sys.argv[3])
