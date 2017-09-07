#!/usr/bin/env python
"""Print a table that gives the B73 size and PH207 size for duplicated genes,
keeping syntenic and nonsyntenic regions separate. Prints syntenic info to
stdout, and nonsyntenic info to stderr. Takes two arguments:
    1) Syntenic assignment table
    2) Nonsyntenic assignment table
"""

import sys


def parse_syn(s):
    """Parse the syntenic table and return two lists: B73 cluster sizes and
    PH207 cluster sizes. Subgenome 1 and subgenome 2 will not be distinguished
    in the output."""
    b_sizes = []
    p_sizes = []
    with open(s, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split()
                b1_genes = set(tmp[6].split(','))
                b2_genes = set(tmp[7].split(','))
                p1_genes = set(tmp[8].split(','))
                p2_genes = set(tmp[9].split(','))
                # Remove NA genes
                b1_genes.discard('NA')
                b2_genes.discard('NA')
                p1_genes.discard('NA')
                p2_genes.discard('NA')
                # Then, append the sizes to the lists
                b_sizes.append(len(b1_genes))
                p_sizes.append(len(p1_genes))
                b_sizes.append(len(b2_genes))
                p_sizes.append(len(p2_genes))
    return (b_sizes, p_sizes)


def parse_nonsyn(ns):
    """Parse the nonsyntenic table and return two lists: B73 cluster sizes and
    PH207 cluster sizes."""
    b_sizes = []
    p_sizes = []
    with open(ns, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split()
                b_genes = set(tmp[3].split(','))
                p_genes = set(tmp[4].split(','))
                # Remove NA
                b_genes.discard('NA')
                p_genes.discard('NA')
                b_sizes.append(len(b_genes))
                p_sizes.append(len(p_genes))
    return (b_sizes, p_sizes)


def main(syn, nonsyn):
    """Main function."""
    b_s_sizes, p_s_sizes = parse_syn(syn)
    b_n_sizes, p_n_sizes = parse_nonsyn(nonsyn)
    # Then, print syntenic sizes to stdout
    sys.stdout.write('B_Size\tP_Size\n')
    for b, p in zip(b_s_sizes, p_s_sizes):
        towrite = str(b) + '\t' + str(p) + '\n'
        sys.stdout.write(towrite)
    # And nonsyntenic to stderr
    sys.stderr.write('B_Size\tP_Size\n')
    for b, p in zip(b_n_sizes, p_n_sizes):
        towrite = str(b) + '\t' + str(p) + '\n'
        sys.stderr.write(towrite)
    return


if len(sys.argv) != 3:
    print """Print a table that gives the B73 size and PH207 size for duplicated genes,
keeping syntenic and nonsyntenic regions separate. Prints syntenic info to
stdout, and nonsyntenic info to stderr. Takes two arguments:
    1) Syntenic assignment table
    2) Nonsyntenic assignment table"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2])
