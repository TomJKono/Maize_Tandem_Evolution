#!/usr/bin/env python
"""Select tandem duplicate clusters for analysis with PAML. The simplest case
of doing tandem duplicate tests is to find tandem duplicates that are private
to one genotype, one subgenome, and all other genes are single copy. We also
set a filter on the orthogroups. If all members of the tandem cluster are not
part of the same orthogroup, then skip it. Takes two arguments:
    1) Syntenic tandem duplicate assignment file
    2) Tandem orthogroups file
"""

import sys
import pprint

# Define a constant for the patterns of duplicates that we want to analyze with
# PAML tests. Basically - genotype and subgenome private duplicates and all
# homeologous genes retained
PAML_CLUST = [
    (2, 1, 1, 1),
    (1, 2, 1, 1),
    (1, 1, 2, 1),
    (1, 1, 1, 2)
    ]

def classify_tandem(syn_line):
    """Return a vector of states for B1, B2, P1, P2. The states are:
    0: lost (deleted)
    1: single copy (retained)
    2: duplicated."""
    tmp = syn_line.strip().split()
    if tmp[2] == 'NA' and tmp[6] == 'NA':
        b1 = 0
    elif tmp[2] == 'NA' and tmp[6].startswith('Z'):
        b1 = 1
    elif tmp[2] != 'NA':
        b1 = 2
    if tmp[3] == 'NA' and tmp[7] == 'NA':
        b2 = 0
    elif tmp[3] == 'NA' and tmp[7].startswith('Z'):
        b2 = 1
    elif tmp[3] != 'NA':
        b2 = 2
    if tmp[4] == 'NA' and tmp[8] == 'NA':
        p1 = 0
    elif tmp[4] == 'NA' and tmp[8].startswith('Z'):
        p1 = 1
    elif tmp[4] != 'NA':
        p1 = 2
    if tmp[5] == 'NA' and tmp[9] == 'NA':
        p2 = 0
    elif tmp[5] == 'NA' and tmp[9].startswith('Z'):
        p2 = 1
    elif tmp[5] != 'NA':
        p2 = 2
    return (b1, b2, p1, p2)


def parse_ortho(og):
    """Parse the orthogroup file and store the tandem->OG relationship in
    a nested dictionary:
    {
        B73: {tandem: OG, tandem:OG, ...},
        PH207: {tandem: OG, tandem:OG, ...}
    }"""
    tand_og = {'B73': {}, 'PH207': {}}
    with open(og, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split('\t')
                # If the gene is not duplicated, then we continue
                if tmp[2] == 'N' or tmp[4] == 'N':
                    continue
                else:
                    geno = tmp[1]
                    clust = tmp[3]
                    og_id = tmp[5]
                    # Check if the tandem has already been seen
                    if clust in tand_og[geno]:
                        # If the OG and tandem combination has been seen, then
                        # we continue
                        if og_id in tand_og[geno][clust]:
                            continue
                        else:
                            tand_og[geno][clust].append(og_id)
                    else:
                        tand_og[geno][clust] = [og_id]
    return tand_og



def main(syn, ogs):
    """Main function."""
    # Parse the tandem orthogroup file
    t_og = parse_ortho(ogs)
    # Then, for each syntenic tandem cluster:
    with open(syn, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                # Figure out if this is a duplicate that we are interested in
                # analyzing
                state = classify_tandem(line)
                if state in PAML_CLUST:
                    tmp = line.strip().split()
                    # Decide which genotype we are querying - it will be the
                    # only not-NA one
                    b1 = tmp[2]
                    b2 = tmp[3]
                    p1 = tmp[4]
                    p2 = tmp[5]
                    if b1 != 'NA':
                        clust = b1
                        geno = 'B73'
                    elif b2 != 'NA':
                        clust = b2
                        geno = 'B73'
                    elif p1 != 'NA':
                        clust = p1
                        geno = 'PH207'
                    elif p2 != 'NA':
                        clust = p2
                        geno = 'PH207'
                    # Find out which orthogroups this tandem hits. If it's more
                    # than 1, we skip it.
                    tandem_og = t_og[geno].get(clust, None)
                    if not tandem_og:
                        continue
                    elif len(tandem_og) != 1:
                        continue
                    else:
                        # Otherwise, we print the duplicate ID, orthogroup ID,
                        # and the clusters
                        toprint = '\t'.join([
                            tmp[0],
                            tandem_og[0],
                            b1,
                            b2,
                            p1,
                            p2])
                        print toprint
                else:
                    continue
    return


if len(sys.argv) != 3:
    print """Select tandem duplicate clusters for analysis with PAML. The simplest case
of doing tandem duplicate tests is to find tandem duplicates that are private
to one genotype, one subgenome, and all other genes are single copy. We also
set a filter on the orthogroups. If all members of the tandem cluster are not
part of the same orthogroup, then skip it. Takes two arguments:
    1) Syntenic tandem duplicate assignment file
    2) Tandem orthogroups file"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2])
