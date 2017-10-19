#!/usr/bin/env python
"""Script to calculate the widths of tandem duplicate clusters. The "width" is
defined as the number of genes included in the window from the first to the
last gene, including the tandem duplicates. Takes two arguments:
    1) Tandem Duplicate gene IDs
    2) GFF
"""


import sys


def parse_gff(gff):
    """Read the GFF and return an ordered list of each gene on each chromosome.
    Assumes the GFF is sorted by start position."""
    gff_data = {}
    with open(gff, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                tmp = line.strip().split('\t')
                chrom = tmp[0]
                geneid = tmp[-1].split(';')[0][3:]
                if chrom not in gff_data:
                    gff_data[chrom] = [geneid]
                else:
                    gff_data[chrom].append(geneid)
    return gff_data


def calc_width(cluster, gff):
    """Given a list of gene IDs in a cluster, calculate the width of the
    cluster."""
    #   Assume that the first gene is on the same chromosome as the others
    for c in gff:
        if cluster[0] in gff[c]:
            chromgenes = gff[c]
            break
    #   Caclulate the indices of the genes in the chromosome
    indices = [chromgenes.index(g) for g in cluster]
    #   Sort it
    indices.sort()
    # Calculate the number of intervening genes between tandems
    intervening = []
    for i in range(1, len(indices)):
        # The number of intervening genes is one less than the differences of
        # indices. E.g., if tandems are directly adjacent, the difference in
        # their indices is 1, and there is 0=1-1 intervening genes.
        intervening.append(indices[i] - indices[i-1] - 1)
    return intervening


def main(tandems, gff):
    """Main function."""
    g_data = parse_gff(gff)
    print 'ClusterID\tGenes\tWidth\t'
    with open(tandems, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            c_id = tmp[0]
            c_list = tmp[1].split(',')
            w = calc_width(c_list, g_data)
            print c_id + '\t' + tmp[1] + '\t' + ','.join([str(i) for i in w])
    return


main(sys.argv[1], sys.argv[2])
