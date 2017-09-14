#!/usr/bin/env python
"""Calculates the distance between the start coordinates of each pair of
tandemly duplicated genes. For multiple-gene clusters, it will only calculate
the adjacent gene distances. Takes two arguments:
    1) File with duplicate IDs, separated by commas, one cluster per line
    2) GFF
Writes to stdout."""

import sys
import math


def parse_gff(g):
    """Read the GFF and return a dictionary of the form
    { gene_id: (chromosome, start),
      gene_id: (chromosome, start),
      ...
    }
    """
    gff_dat = {}
    with open(g, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                tmp = line.strip().split('\t')
                if tmp[2] == 'gene':
                    metadata = tmp[-1].split(';')
                    for m in metadata:
                        if m.startswith('ID'):
                            g_id = m[3:]
                            start = int(tmp[3])
                            chrom = tmp[0]
                            gff_dat[g_id] = (chrom, start)
                            break
    return gff_dat


def calculate_distance(gff, g1, g2):
    """Calculate the distance in base pairs between the start coordinates of
    the two supplied gene models. Returns a positive integer, or NA if the two
    gene models are on separate chromosomes (which shouldn't happen...)"""
    if g1 == g2:
        return 'NA'
    elif g1 not in gff or g2 not in gff:
        return 'NA'
    elif gff[g1][0] != gff[g2][0]:
        return 'NA'
    else:
        dist = gff[g2][1] - gff[g1][1]
        return math.fabs(dist)


def calculate_intervening(gff, g1, g2):
    """Calculate the number of genes between adjacent pairs of tandem
    duplicates, using the GFF data."""
    if g1 == g2:
        return 'NA'
    elif g1 not in gff or g2 not in gff:
        return 'NA'
    elif gff[g1][0] != gff[g2][0]:
        return 'NA'
    else:
        #   Find the endpoints of the search interval as the absolute left and
        #   right bounds
        first = min(gff[g1][1], gff[g2][1])
        end = max(gff[g1][1], gff[g2][1])
        #   count the number of genes with a comprehension
        in_interval = [
            x
            for x
            in gff.iteritems()
            if x[1][0] == gff[g1][0] and x[1][0] == gff[g2][0] and x[1][1] > first and x[1][1] < end]
        return len(in_interval)


def main(tandem_file, gff):
    """Main function."""
    gff_data = parse_gff(gff)
    #   Print a header
    print '\t'.join([
        'Gene1',
        'Gene2',
        'Dist_BP',
        'N_In_Between'])
    #   Iteate through the tandems file
    with open(tandem_file, 'r') as f:
        for line in f:
            #   Split on commas, and sort them
            tandems = sorted(line.strip().split(','))
            #   For each adjacent pair...
            for i in xrange(0, len(tandems)-1):
                d = calculate_distance(
                    gff_data,
                    tandems[i],
                    tandems[i+1])
                n = calculate_intervening(
                    gff_data,
                    tandems[i],
                    tandems[i+1])
                print '\t'.join([tandems[i], tandems[i+1], str(d), str(n)])


main(sys.argv[1], sys.argv[2])
