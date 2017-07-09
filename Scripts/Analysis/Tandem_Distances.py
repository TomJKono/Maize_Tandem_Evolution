#!/usr/bin/env python
"""Script to calculate the interval spacing of genes in tandem duplicates.
Between each consecutive pair of tandem duplicates, it will report the number
of intervening genes. Takes two arguments:
    1) Tandem duplicates CSV
    2) Coordinate-sorted GFF
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
    # N intervening is (difference of indices) - 1
    widths = []
    for index, coord in enumerate(indices[1:]):
        widths.append(str(coord - indices[index] - 1))
    return widths


def main(tandems, gff):
    """Main function."""
    g_data = parse_gff(gff)
    print 'ClusterID\tGenes\tIntervals\t'
    with open(tandems, 'r') as f:
        for line in f:
            g_list = line.strip().split()[1]
            c_list = g_list.split(',')
            w = calc_width(c_list, g_data)
            print line.strip() + '\t' + ','.join(w)
    return


main(sys.argv[1], sys.argv[2])
