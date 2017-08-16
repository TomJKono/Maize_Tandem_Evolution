#!/usr/bin/env python
"""Estimate the ages of the syntenic tandem duplications using annotated NEXUS
trees from BEAST/TreeAnnotator. Be sure that the raw output from TreeAnnotator
has had the ampersands (&) removed, because they cause the Biopython parser to
choke. We will assume that the distances beteen genes in the same tandem
duplicate cluster represents their age, and that substitutions accumulate in a
clock fashion. Takes two arguments:
    1) Synteny tandem duplicate assignment file
    2) Directory of annotated trees
"""

import sys
import os
try:
    from Bio import Phylo
except ImportError:
    print 'This script requires Biopython.'
    exit(1)


def parse_synteny(synteny):
    """Parse the synteny file and return it as a list of lists. We want to keep
    track of the syntenic duplication ID, the genes in each of the subgenomes
    of both genotypes, and the ancestral gene."""
    syn = []
    with open(synteny, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                syn.append(line.strip().split())
    return syn


def tree_distances(tree, genes):
    """Given a parsed tree object and a list of gene IDs, return the distance
    between them. If one gene ID is given, return NA. Calculate nearest gene
    distances, and return as a list"""
    gene_ids = genes.split(',')
    # Get the tip names
    tree_names = [n.name for n in tree.find_clades(terminal=True)]
    if len(gene_ids) == 1:
        return ['NA']
    # This assumes that the gene duplication order goes unidirectionally across
    # the chromosome. This may be a bad assumption, but it's the only way to be
    # sure we don't over-count the duplication history by treating each
    # pairwise combination of genes as an indepdendent duplication event.
    gene_ids = sorted(gene_ids)
    distances = []
    for i in range(0, len(gene_ids)-1):
        j = i+1
        g1 = [x for x in tree_names if x.startswith(gene_ids[i])][0]
        g2 = [x for x in tree_names if x.startswith(gene_ids[j])][0]
        distances.append(str(tree.distance(g1, g2)/2))
    return distances


def main(synteny, tree_dir):
    """Main function."""
    # Get the full path to the tree directory
    t_dir = os.path.abspath(os.path.expanduser(tree_dir))
    syntenic_rel = parse_synteny(synteny)
    # Print a header
    print 'Duplicate_ID\tAncestral_Gene\tB1_Genes\tB2_Genes\tP1_Genes\tP2_Genes\tB1_Ages\tB2_Ages\tP1_Ages\tP2_Ages'
    for r in syntenic_rel:
        dup_id = r[0]
        anc = r[1]
        b1_genes = r[6]
        b2_genes = r[7]
        p1_genes = r[8]
        p2_genes = r[9]
        # Parse the tree
        try:
            tfile = os.path.join(t_dir, dup_id + '_Ann.tree')
            tree = Phylo.read(tfile, 'nexus')
        except:
            continue
        # Then, get the distances
        b1_distances = tree_distances(tree, b1_genes)
        b2_distances = tree_distances(tree, b2_genes)
        p1_distances = tree_distances(tree, p1_genes)
        p2_distances = tree_distances(tree, p2_genes)
        # Format for printing
        b1_distances = ','.join(b1_distances)
        b2_distances = ','.join(b2_distances)
        p1_distances = ','.join(p1_distances)
        p2_distances = ','.join(p2_distances)
        # Print it
        print '\t'.join([
            dup_id,
            anc,
            ','.join(sorted(b1_genes.split(','))),
            ','.join(sorted(b2_genes.split(','))),
            ','.join(sorted(p1_genes.split(','))),
            ','.join(sorted(p2_genes.split(','))),
            b1_distances,
            b2_distances,
            p1_distances,
            p2_distances])
    return


if len(sys.argv) != 3:
    print """Estimate the ages of the syntenic tandem duplications using annotated NEXUS
trees from BEAST/TreeAnnotator. Be sure that the raw output from TreeAnnotator
has had the ampersands (&) removed, because they cause the Biopython parser to
choke. We will assume that the distances beteen genes in the same tandem
duplicate cluster represents their age, and that substitutions accumulate in a
clock fashion. Takes two arguments:
    1) Synteny tandem duplicate assignment file
    2) Directory of annotated trees"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2])
