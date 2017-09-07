#!/usr/bin/env python
"""Combine the results of 'bedtools coverage' commands on 1Mb windows and
various genomic features to make a data table for a linear model test. Takes
six arguments:
    1) Gene overlap
    2) Tandem overlap
    3) RNA TE overlap
    4) DNA TE overlap
    5) Maize1 overlap
    6) Maize2 overlap
"""

import sys


def parse_table(tab):
    """Parse one of the bedtools output tables and return it as a list of
    floating point values."""
    vals = []
    with open(tab, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            # The proportion of overlap is the final column
            vals.append(float(tmp[-1]))
    return vals


def assign_subgenome(m1, m2):
    """Read through the maize1 overlap and maize2 overlap lists and assign
    a majority state to each interval. If a window is neither m1 nor m2, then
    it is nonsyntenic."""
    assignments = []
    for m1_a, m2_a in zip(m1, m2):
        tmp = {
            'Maize1': m1_a,
            'Maize2': m2_a,
            'Nonsyntenic': 1 - (m1_a + m2_a)
            }
        # Get the max value, and return its key
        a = max(tmp.iterkeys(), key=(lambda key: tmp[key]))
        assignments.append(a)
    return assignments


def main(gene, tandem, rna, dna, m1, m2):
    """Main function."""
    gene_o = parse_table(gene)
    tandem_o = parse_table(tandem)
    rna_o = parse_table(rna)
    dna_o = parse_table(dna)
    m1_o = parse_table(m1)
    m2_o = parse_table(m2)
    subgenomes = assign_subgenome(m1_o, m2_o)
    # Print out a nice table
    print 'Prop_Genes\tProp_Tandem\tProp_RNATE\tProp_DNATE\tSubgenome'
    for g, t, r, d, s in zip(gene_o, tandem_o, rna_o, dna_o, subgenomes):
        print '\t'.join([str(g), str(t), str(r), str(d), s])
    return


if len(sys.argv) != 7:
    print """Combine the results of 'bedtools coverage' commands on 1Mb windows and
various genomic features to make a data table for a linear model test. Takes
six arguments:
    1) Gene overlap
    2) Tandem overlap
    3) RNA TE overlap
    4) DNA TE overlap
    5) Maize1 overlap
    6) Maize2 overlap"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
