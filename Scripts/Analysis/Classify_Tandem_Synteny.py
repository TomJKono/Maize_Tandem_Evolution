#!/usr/bin/env python
"""Show some more detail as to which genes are syntenic and nonsyntenic. Takes
three arguments:
    1) ABB gene key
    2) Synteny file
    3) Tandem duplicate clusters"""


import sys


def parse_key(genekey):
    """Parse the gene key and return a list of genes in syntenic relationships
    between the two genotypes."""
    syn = []
    with open(genekey, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            status = tmp[8]
            if status == 'curated-syntenic' or status == 'curated-syntenic-fusion':
                syn.append(tmp[3])
                syn.append(tmp[7])
    return syn


def parse_synteny(syn):
    """Parse the syntenic file and return lists for m1 and m2."""
    m1 = []
    m2 = []
    with open(syn, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            b1 = tmp[1]
            b2 = tmp[2]
            p1 = tmp[3]
            p2 = tmp[4]
            if b1 != 'NA':
                m1.append(b1)
            if b2 != 'NA':
                m2.append(b2)
            if p1 != 'NA':
                m1.append(p1)
            if p2 != 'NA':
                m2.append(p2)
    return (m1, m2)


def main(genekey, syntenic, tandems):
    """Main function."""
    synteny = parse_key(genekey)
    m1, m2 = parse_synteny(syntenic)
    with open(tandems, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            status = ['Syntenic' if g in synteny else 'Nonsyntenic' for g in tmp[1].split(',')]
            if 'Syntenic' in status:
                final = 'Syntenic'
            else:
                final = 'Nonsyntenic'
            for g in tmp[1].split(','):
                if g in m1:
                    subg = 'Maize1'
                    break
                if g in m2:
                    subg = 'Maize2'
                    break
            else:
                subg = 'Nonsyntenic'
            print line.strip() + '\t' + ','.join(status) + '\t' + final + '\t' + subg
    return


if len(sys.argv) != 4:
    print """Show some more detail as to which genes are syntenic and nonsyntenic. Takes
three arguments:
    1) ABB gene key
    2) Synteny file
    3) Tandem duplicate clusters"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2], sys.argv[3])
