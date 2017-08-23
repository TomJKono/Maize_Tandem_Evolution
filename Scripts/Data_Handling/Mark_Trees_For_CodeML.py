#!/home/tomkono/anaconda-ete/bin/python
"""Mark a Newick tree for testing with CodeML. This script is for estimating
functional outcomes of maize tandem duplicates. The tree will be marked with
three partitions:
    #1: Maize tandem branches
    #2: Maize genes generally
    (No mark): Background branches
Takes two arguments:
    1) Directory of trees
    2) Orthogroup-Tandem ID table
"""

import sys
import os
try:
    from ete3 import EvolTree
except ImportError:
    print 'This script requires the ete3 toolkit to be installed.'
    exit(1)


def list_trees(d):
    """List a directory and return paths to each tree file."""
    fd = os.path.abspath(os.path.expanduser(d))
    tree_paths = [t
        for t
        in os.listdir(fd)
        if t.startswith('RAxML_bestTree')]
    return (fd, tree_paths)


def parse_assignment(tandems):
    """Parse the OG-tandem gene ID file and return it as a dictionary."""
    og_t = {}
    with open(tandems, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            og_t[tmp[0]] = tmp[1].split(',')
    return og_t


def parse_tree(t):
    """Simple function to read a tree file off the disk and return it as a
    Tree object. Also calls the render() method on it so that it gets node IDs
    that we can use to apply labels/marks."""
    treeobj = EvolTree(t)
    treeobj.render('')
    return treeobj


def find_tandems(tree, tandem_names):
    """Use the gene IDs to get the node names that lead to tandem duplicates.
    We need this because orthofinder adds external data to the names of the
    genes in the trees."""
    tandem_nodes = []
    for node in tree:
        # the format is SPECIES|GENEID_TXID_COUNTER
        if node.name.split('|')[1].split('_')[0] in tandem_names:
            tandem_nodes.append(node.name)
    return tandem_nodes


def label_branches(tree, tandem_ids):
    """Use the tandem gene IDs to apply labels to the branches that go to maize
    tandem duplicates. We also want to apply labels to the internal branches
    that go to only tandem duplicates. Also use the naming scheme of the maize
    genes to label the maize branches. We leave every other branch alone to
    serve as the background. We do not distinguish between B73 and PH207,
    because the data are not appropriate for that distinction."""
    # Keep lists of nodes for marking
    tandem = []
    maize = []
    # First, label the branches as tandem or not
    for leaf in tree:
        if leaf.name in tandem_ids:
            leaf.add_features(cat='Tandem')
        elif leaf.name.startswith('B73') or leaf.name.startswith('PH207'):
            leaf.add_features(cat='Maize')
        else:
            leaf.add_features(cat='Grass')
    # Then find monophyletic chunks within those marked as tandem
    for i in tree.get_monophyletic(values=['Tandem'], target_attr='cat'):
        # For each monophyletic sub-tree, we want to label the branches
        for node in i.traverse():
            tandem.append(node.node_id)
    # Same for maize
    for i in tree.get_monophyletic(values=['Maize'], target_attr='cat'):
        for node in i.traverse():
            maize.append(node.node_id)
    # Then apply the marks!
    tree.mark_tree(tandem + maize, marks=['#1']*len(tandem) + ['#2']*len(maize))
    return tree


def write_marked(tree, og):
    """Write the marked tree into a standarized filename."""
    fname = og + '_Marked.tree'
    tree.write(features=[], outfile=fname)
    return


def write_ctl(og):
    """Write a CodeML control file to go along with the tree file."""
    fname = og + '.ctl'
    handle = open(fname, 'w')
    # Set the input and output paths
    handle.write('seqfile = ' + og + '_BKT.fa\n')
    handle.write('outfile = ' + og + '_codeml\n')
    handle.write('treefile = ' + og + '_Marked.tree\n')
    # Set the run verbosity parameters
    handle.write('noisy = 9\n')
    handle.write('verbose = 1\n')
    handle.write('runmode = 0\n')
    # Set broad parameters
    #   seqtype=1: codons; 2: AA; 3: codons->AAs
    handle.write('seqtype = 1\n')
    #   CodonFreq=1: 1/61 each; 2: Based on nuc freq; 3: nuc freq at each pos; 4: free parameters
    handle.write('CodonFreq = 2\n')
    #   The number of datasets
    handle.write('ndata = 1\n')
    #   clock=0: no clock all branches free; 1: strict global clock; 2: local clock rates
    handle.write('clock = 0\n')
    #   aaDist=0: equal AA distances for al subs; 1: use Grantham matrix
    handle.write('aaDist = 0\n')
    #   We use model=3 and NSites=2 for clade model C
    handle.write('model = 3\n')
    handle.write('NSsites = 2\n')
    #   icode=0: Standard genetic code
    handle.write('icode = 0\n')
    #   Mgene=0: No partitioning of sites
    handle.write('Mgene = 0\n')
    #   Do not fix kappa, estimate it, start at 2
    handle.write('fix_kappa = 0\n')
    handle.write('kappa = 2\n')
    #   Same with omega (dN/dS), start at 0.5
    handle.write('fix_omega = 0\n')
    handle.write('omega = 0.5\n')
    #   Fix alpha and set to 0, because we are using NSsites
    handle.write('fix_alpha = 1\n')
    handle.write('alpha = 0\n')
    handle.write('Malpha = 0\n')
    handle.write('ncatG = 3\n')
    #   Fix rho (rates at adjacent sites)
    handle.write('fix_rho = 1\n')
    handle.write('rho = 0\n')
    #   We want standard errors
    handle.write('getSE = 1\n')
    #   We do not want ancestral state reconstruction nor ancestral rates
    handle.write('RateAncestor = 0\n')
    #   Set the convergence tolerance
    handle.write('Small_Diff = 0.5e-6\n')
    #   Treat gaps as missing
    handle.write('cleandata = 1\n')
    #   Do not fix branch lengths
    handle.write('fix_blength = 0\n')
    #   Estimate all at once
    handle.write('method = 0\n')
    handle.close()
    return


def main(tree_dir, tandems):
    """Main function."""
    fd, trees = list_trees(tree_dir)
    og_tandem = parse_assignment(tandems)
    for og in sorted(og_tandem):
        t_file = os.path.join(fd, 'RAxML_bestTree.' + og)
        tree = parse_tree(t_file)
        tandem_names = find_tandems(tree, og_tandem[og])
        lab_tree = label_branches(tree, tandem_names)
        # Write the marked tree
        write_marked(lab_tree, og)
        # And a control file
        write_ctl(og)
    return


if len(sys.argv) != 3:
    print """Mark a Newick tree for testing with CodeML. This script is for estimating
functional outcomes of maize tandem duplicates. The tree will be marked with
three partitions:
    #1: Maize tandem branches
    #2: Maize genes generally
    (No mark): Background branches
Takes two arguments:
    1) Directory of trees
    2) Orthogroup-Tandem ID table"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2])
