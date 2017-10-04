#!/usr/bin/env python
"""A very simple script to calculate the proportion of overlap between gene
regions and TE regions, from 'bedtools intersect' output. This script expects
that bedtools was run with the gene GFF as intervals A, and the TE GFF was run
as intervals B. Takes one argument:
    1) bedtools intersect output
"""

import sys


def main(bedtools_out):
    """Main function."""
    with open(bedtools_out, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            gene_len = int(tmp[4]) - int(tmp[3])
            # Use a float here to avoid truncation
            overlap = float(tmp[-1])
            prop_overlap = str(overlap/gene_len)
            print '\t'.join(tmp + [prop_overlap])
    return


if len(sys.argv) != 2:
    print """A very simple script to calculate the proportion of overlap between gene
regions and TE regions, from 'bedtools intersect' output. Takes one argument:
    1) bedtools intersect output"""
    exit(1)
else:
    main(sys.argv[1])
