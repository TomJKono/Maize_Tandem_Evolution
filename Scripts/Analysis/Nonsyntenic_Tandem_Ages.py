#!/usr/bin/env python
"""Estimate the ages of the nonsyntenic tandem duplications using annotated
NEXUS trees from BEAST/TreeAnnotator. Be sure that the raw output from
TreeANnotator has the ampersands (&) removed, because they cause the Biopython
parser to choke. Assume that the distances between tandem duplicate genes
represents their age, and that substitutions are clock-like. Takes two
arguments:
    1) Nonsyntenic tandem duplicate assignment file
    2) Directory of fixed annotated trees
"""

import sys
import os
try:
    from Bio import Phylo
except ImportError:
    print 'This script requires Biopython.'
    exit(1)


def parse_assignment(nonsyntenic):
    """Parse the nonsyntenic file and return it as a list of lists. Keep track
    of the nonsyntenic duplication ID and the genes in each of the
    subgenomes."""
    nonsyn = []
    with open(nonsyntenic, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                nonsyn.append(line.strip().split())
    return nonsyn


def tree_distances(tree, genes):
    """Given a parsed tree object and a list of gene IDs, return the distance
    between them. If one gene ID is given, return NA. Calculate adjacent gene
    distances, and return as a list."""
    gene_ids = genes.split(',')
    tree_names = [n.name for n in tree.find_clades(terminal=True)]
    if len(gene_ids) == 1:
        return ['NA']
    # We assume the number of duplications is the number of genes -1. That is,
    # we assume that duplications affect only one gene at a time. It may be
    # a bad assumption, but we don't want to over-count old duplications
    # that are shared.
    gene_ids = sorted(gene_ids)
    distances = []
    for i in range(0, len(gene_ids)-1):
        j = i + 1
        g1 = [x for x in tree_names if x.startswith(gene_ids[i])][0]
        g2 = [x for x in tree_names if x.startswith(gene_ids[j])][0]
        distances.append(str(tree.distance(g1, g2)/2))
    return distances


def main(nonsyntenic, tree_dir):
    """Main function."""
    t_dir = os.path.abspath(os.path.expanduser(tree_dir))
    nonsyntenic_rel = parse_assignment(nonsyntenic)
    print 'Duplicate_ID\tB_Genes\tP_Genes\tB_Ages\tP_Ages'
    for r in nonsyntenic_rel:
        dup_id = r[0]
        b_genes = r[3]
        p_genes = r[4]
        # Parse the tree
        try:
            tfile = os.path.join(t_dir, dup_id + '_Ann.tree')
            tree = Phylo.read(tfile, 'nexus')
        except:
            continue
        # Get distances
        b_distances = tree_distances(tree, b_genes)
        p_distances = tree_distances(tree, p_genes)
        # Format for printing
        b_distances = ','.join(b_distances)
        p_distances = ','.join(p_distances)
        # Print it
        print '\t'.join([
            dup_id,
            ','.join(sorted(b_genes.split(','))),
            ','.join(sorted(p_genes.split(','))),
            b_distances,
            p_distances])
    return


if len(sys.argv) != 3:
    print """Estimate the ages of the nonsyntenic tandem duplications using annotated
NEXUS trees from BEAST/TreeAnnotator. Be sure that the raw output from
TreeANnotator has the ampersands (&) removed, because they cause the Biopython
parser to choke. Assume that the distances between tandem duplicate genes
represents their age, and that substitutions are clock-like. Takes two
arguments:
    1) Nonsyntenic tandem duplicate assignment file
    2) Directory of fixed annotated trees"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2])
