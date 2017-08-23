#!/usr/bin/env python
"""Generate tables for the GC content of tandem duplicates, and the estimated
ages of the duplications. Takes four arguments:
    1) Syntenic tandem duplicates ages
    2) Nonsyntenic tandem duplicates ages
    3) B73 GC contents
    4) PH207 GC contents
"""

import sys
import pprint


def parse_syn(syntenic):
    """Parse the syntenic table and store as a dictionary that associates a
    tandem duplicate with its age."""
    syn_ages = {}
    with open(syntenic, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split('\t')
                # keep only the tandems
                b_genes = []
                p_genes = []
                if ',' in tmp[2]:
                    b_genes += tmp[2].split(',')
                if ',' in tmp[3]:
                    b_genes += tmp[3].split(',')
                if ',' in tmp[4]:
                    p_genes += tmp[4].split(',')
                if ',' in tmp[5]:
                    p_genes += tmp[5].split(',')
                b_ages = []
                p_ages = []
                if tmp[6] != 'NA':
                    b_ages += tmp[6].split(',')
                if tmp[7] != 'NA':
                    b_ages += tmp[7].split(',')
                if tmp[8] != 'NA':
                    p_ages += tmp[8].split(',')
                if tmp[9] !='NA':
                    p_ages += tmp[9].split(',')
                # If there are B73 tandem duplicates, then we iterate through
                # and build the dictionary
                if b_ages:
                    ct = 1
                    for age in b_ages:
                        b1 = b_genes[ct-1]
                        b2 = b_genes[ct]
                        syn_ages[(b1, b2)] = age
                        ct += 1
                if p_ages:
                    ct = 1
                    for age in p_ages:
                        p1 = p_genes[ct-1]
                        p2 = p_genes[ct]
                        syn_ages[(p1, p2)] = age
                        ct += 1
    return syn_ages


def parse_nonsyn(nonsyn):
    """Parse the nonsyntenic table and store as a dictionary in the same format
    as the syntenic tandem duplicates."""
    nonsyn_ages = {}
    with open(nonsyn, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split('\t')
                b_genes = []
                p_genes = []
                if ',' in tmp[1]:
                    b_genes += tmp[1].split(',')
                if ',' in tmp[2]:
                    p_genes += tmp[2].split(',')
                b_ages = []
                p_ages = []
                if tmp[3] != 'NA':
                    b_ages += tmp[3].split(',')
                if tmp[4] != 'NA':
                    p_ages += tmp[4].split(',')
                if b_ages:
                    ct = 1
                    for age in b_ages:
                        b1 = b_genes[ct-1]
                        b2 = b_genes[ct]
                        nonsyn_ages[(b1, b2)] = age
                        ct += 1
                if p_ages:
                    ct = 1
                    for age in p_ages:
                        p1 = p_genes[ct-1]
                        p2 = p_genes[ct]
                        nonsyn_ages[(p1, p2)] = age
                        ct += 1
    return nonsyn_ages


def parse_gc(gc):
    """Parse the gene GC contents file and store as a dictionary."""
    prop_gc = {}
    with open(gc, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            prop_gc[tmp[0]] = tmp[1]
    return prop_gc


def main(syn, nonsyn, b_gc, p_gc):
    """Main function."""
    b73_gc = parse_gc(b_gc)
    ph207_gc = parse_gc(p_gc)
    s_dups = parse_syn(syn)
    ns_dups = parse_nonsyn(nonsyn)
    # Start iterating through the duplications and associate ages with GC
    # contents. Print out a header
    print 'Gene1\tGene2\tAge\tGene1_GC\tGene2_GC\tClass'
    for pair in sorted(s_dups):
        g1 = pair[0]
        g2 = pair[1]
        age = s_dups[pair]
        # Find out if we are in B73 or PH207
        if g1.startswith('Zm00001d'):
            g1_gc = b73_gc[g1]
            g2_gc = b73_gc[g2]
        else:
            g1_gc = ph207_gc[g1]
            g2_gc = ph207_gc[g2]
        # Print out the table
        print '\t'.join([g1, g2, age, g1_gc, g2_gc, 'Syntenic'])
    # Do the same for nonsyntenic
    for pair in sorted(ns_dups):
        g1 = pair[0]
        g2 = pair[1]
        age = ns_dups[pair]
        # Find out if we are in B73 or PH207
        if g1.startswith('Zm00001d'):
            g1_gc = b73_gc[g1]
            g2_gc = b73_gc[g2]
        else:
            g1_gc = ph207_gc[g1]
            g2_gc = ph207_gc[g2]
        # Print out the table
        print '\t'.join([g1, g2, age, g1_gc, g2_gc, 'Nonsyntenic'])
    return


if len(sys.argv) != 5:
    print """Generate tables for the GC content of tandem duplicates, and the estimated
ages of the duplications. Takes four arguments:
    1) Syntenic tandem duplicates ages
    2) Nonsyntenic tandem duplicates ages
    3) B73 GC contents
    4) PH207 GC contents"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
