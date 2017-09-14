#!/usr/bin/env python
"""Reads the list of tandem duplicates and tries to determine if one of them
is likely a pseudogene based on sequence similarity to known TEs and the number
of exons. Takes one argument:
    1) GFF of representative transcripts"""

import sys


def parse_gff(g):
    """Read GFF data. For each gene, store its length and the number of exons
    in its representative transcript. This will return a data structure of the
    following form:
    {
        gene_id: {
            Length: INT,
            Chromosome: STR,
            Start: INT,
            Stop: INT,
            Exons: [EX1, EX2, ..., ]
            },
        gene_id: {
            Length: INT,
            Chromosome: STR,
            Start: INT,
            Stop: INT,
            Exons: [EX1, EX2, ..., ]
            },
        ...
    }"""
    gff_dat = {}
    with open(g, 'r') as f:
        for line in f:
            #   Skip header lines
            if line.startswith('#'):
                continue
            else:
                tmp = line.strip().split('\t')
                #   We want to save information from genes
                if tmp[2] == 'gene':
                    metadata = tmp[-1].split(';')
                    for m in metadata:
                        if m.startswith('ID'):
                            gid = m[3:]
                            break
                    if gid not in gff_dat:
                        gff_dat[gid] = {
                            'Length': int(tmp[4]) - int(tmp[3]),
                            'Chromosome': tmp[0],
                            'Exons': []}
                elif tmp[2] == 'exon':
                    metadata = tmp[-1].split(';')
                    for m in metadata:
                        if m.startswith('Parent'):
                            parent = m[7:].split('_')[0]
                        if m.startswith('ID') or m.startswith('exon_id'):
                            eid = m[3:]
                    if parent in gff_dat:
                        gff_dat[parent]['Exons'].append(eid)
                elif tmp[2] == 'CDS':
                    metadata = tmp[-1].split(';')
                    for m in metadata:
                        if m.startswith('Parent'):
                            parent = m[7:].split('_')[0]
    return gff_dat


def main(gff):
    """Main function."""
    gff_data = parse_gff(gff)
    for g in sorted(gff_data):
        print len(gff_data[g]['Exons'])
    return


if len(sys.argv) != 2:
    print """Reads the list of tandem duplicates and tries to determine if one of them
is likely a pseudogene based on sequence similarity to known TEs and the number
of exons. Takes one argument:
    1) GFF of representative transcripts"""
    exit(1)
else:
    main(sys.argv[1])
