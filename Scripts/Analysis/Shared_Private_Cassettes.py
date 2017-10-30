#!/usr/bin/env python
"""Count the shared/private cassettes between B73 and PH207, using the
gene key file from ABB. Basically, what we will do so store the association
between B73 genes and PH207 genes, then ask if genes in B73 and PH207 cassettes
reciprocally point to each other in this association. Takes three argments:
    1) Gene key file
    2) B73 Cassettes
    3) PH207 cassettes
"""

import sys


def parse_key(k):
    """Parse the key and store it in a list of tuples. We use that instead of
    a dict because the gene mapping is not one-to-one."""
    gk = []
    with open(k, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split()
                b = tmp[3]
                p = tmp[7]
                gk.append((b, p))
    return gk


def store_cass(c):
    """Parse the cassettes file and store it in a list of lists. We want to do
    that because we want to preserve the associations within cassettes."""
    cassettes = []
    with open(c, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split()
                c_genes = tmp[0].split(',')
                cassettes.append(c_genes)
    return cassettes


def check_homology(gk, b, p):
    """Check the homology of the B73 cassette and PH207 cassette by searching
    for reciprocal homology of genes in the cassettes."""
    # First, find the PH207 homologues of the B73 cassette genes
    b_homologues = []
    for b_gene in b:
        for pair in gk:
            if b_gene in pair:
                b_homologues.append(pair[1])
    # Then, find the B73 homologues of the PH207 cassette genes
    p_homologues = []
    for p_gene in p:
        for pair in gk:
            if p_gene in pair:
                p_homologues.append(pair[0])
    # Then, if the homolgoues and the cassettes share genes, then the cassettes
    # are shared
    b_share = set(b) & set(p_homologues)
    p_share = set(p) & set(b_homologues)
    if b_share and p_share:
        return True
    else:
        return False


def main(genekey, b_cass, p_cass):
    """Main function."""
    gene_key = parse_key(genekey)
    b_cassettes = store_cass(b_cass)
    p_cassettes = store_cass(p_cass)
    # Then, we want to iterate through each pair of cassettes, and find ones
    # link up. Store the indices of the cassettes that are homologous. Try to
    # be sure that things do not get double-counted.
    sharing = []
    for b_i, b_c in enumerate(b_cassettes):
        for p_i, p_c in enumerate(p_cassettes):
            shared = check_homology(gene_key, b_c, p_c)
            # if they are shared,...
            if shared:
                # But we haven't already seen this combination of indices
                if (p_i, b_i) in sharing:
                    continue
                else:
                    # Append the indices into the sharing list
                    sharing.append((b_i, p_i))
    # Print the cassettes that are shared
    print 'B73_Cassette_Genes\tPH207_Cassette_Genes'
    for pair in sharing:
        b = ','.join(b_cassettes[pair[0]])
        p = ','.join(p_cassettes[pair[1]])
        print b + '\t' + p

if len(sys.argv) != 4:
    print """Count the shared/private cassettes between B73 and PH207, using the
gene key file from ABB. Basically, what we will do so store the association
between B73 genes and PH207 genes, then ask if genes in B73 and PH207 cassettes
reciprocally point to each other in this association. Takes three argments:
    1) Gene key file
    2) B73 Cassettes
    3) PH207 cassettes"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2], sys.argv[3])
