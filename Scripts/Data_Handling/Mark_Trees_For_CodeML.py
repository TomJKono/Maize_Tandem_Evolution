#!/home/tomkono/anaconda-ete/bin/python
"""Mark a Newick tree for testing with CodeML. This script is for estimating
functional outcomes of maize tandem duplicates. The tree will be marked with
several schemes to test relative rates of duplicates. See the markdown file for
details on how they will be marked. Takes three arguments:
    1) Directory of trees
    2) Orthogroup-Tandem ID table
    3) Output directory
"""

import sys
import os
import copy
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


def label_branches(tree, tandem_ids, hyp):
    """Use the tandem gene IDs to apply labels to the branches that go to maize
    tandem duplicates. We also want to apply labels to the internal branches
    that go to only tandem duplicates. Also use the naming scheme of the maize
    genes to label the maize branches. We leave every other branch alone to
    serve as the background. We do not distinguish between B73 and PH207,
    because the data are not appropriate for that distinction."""
    # Copy the tree because we need to modify it multiple times, and don't want
    # to clobber the marks
    new_tree = copy.deepcopy(tree)
    # First, label the branches as tandem or not
    for leaf in new_tree:
        if leaf.name in tandem_ids:
            leaf.add_features(cat='Tandem')
        elif leaf.name.startswith('B73') or leaf.name.startswith('PH207'):
            leaf.add_features(cat='Maize')
        else:
            leaf.add_features(cat='Grass')
    # Then, depending on which hypothesis we are testing, we mark the trees
    # differently.
    if hyp == 'H0':
        # Null hypothesis: no variation among branches. No marks.
        return new_tree
    elif hyp == 'Ha1':
        # Keep lists of nodes for marking
        maize = []
        # This hypothesis marks all maize genes (both tandem and non-tandem) as
        # different.
        for i in new_tree.get_monophyletic(values=['Tandem'], target_attr='cat'):
            # For each monophyletic sub-tree, we want to label the branches
            for node in i.traverse():
                maize.append(node.node_id)
        for i in new_tree.get_monophyletic(values=['Maize'], target_attr='cat'):
            for node in i.traverse():
                maize.append(node.node_id)
        # Then apply the marks!
        new_tree.mark_tree(maize, marks=['#1']*len(maize))
        return new_tree
    elif hyp == 'Ha2':
        # Keep lists of nodes for marking
        tandem = []
        # This hypothesis marks maize tandems as different from all other genes
        for i in new_tree.get_monophyletic(values=['Tandem'], target_attr='cat'):
            # For each monophyletic sub-tree, we want to label the branches
            for node in i.traverse():
                tandem.append(node.node_id)
        # Then apply the marks!
        new_tree.mark_tree(tandem, marks=['#1']*len(tandem))
        return new_tree
    elif hyp == 'Ha3':
        # Keep lists of nodes for marking
        tandem = []
        maize = []
        # This hypothesis marks tandems, maize non-tandems, and other grass
        # genes with different labels.
        for i in new_tree.get_monophyletic(values=['Maize'], target_attr='cat'):
            for node in i.traverse():
                maize.append(node.node_id)
        for i in new_tree.get_monophyletic(values=['Tandem'], target_attr='cat'):
            # For each monophyletic sub-tree, we want to label the branches
            for node in i.traverse():
                tandem.append(node.node_id)
        # Then apply the marks!
        new_tree.mark_tree(tandem + maize, marks=['#1']*len(tandem) + ['#2']*len(maize))
        return new_tree
    else:
        return None


def make_dirs(og, out_dir):
    """Make the directory structure for the PAML tests. PAML writes a bunch of
    files into the directory it runs, so we segregate all the analyses."""
    fullpath = os.path.abspath(os.path.expanduser(out_dir))
    ogpath = os.path.join(fullpath, og)
    # Make one for each model
    os.makedirs(os.path.join(ogpath, 'Null'))
    os.makedirs(os.path.join(ogpath, 'Ha1'))
    os.makedirs(os.path.join(ogpath, 'Ha2'))
    os.makedirs(os.path.join(ogpath, 'Ha3'))
    return ogpath


def write_marked(tree, og, hyp, outdir):
    """Write the marked tree into a standarized filename."""
    hyp_path = os.path.join(outdir, hyp)
    fname = os.path.join(hyp_path, og + '_Marked.tree')
    tree.write(features=[], outfile=fname)
    return


def write_ctl(og, hyp, outdir):
    """Write a CodeML control file to go along with the tree file."""
    hyp_path = os.path.join(outdir, hyp)
    fname = os.path.join(hyp_path, og + '.ctl')
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
    #       Null model gets model=0 and NSsites=22
    if hyp == 'Null':
        handle.write('model = 0\n')
        handle.write('NSsites = 22\n')
    else:
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
    #   Don't try to get SEs - they cause problems when estimating sub rates
    handle.write('getSE = 0\n')
    #   We do not want ancestral state reconstruction nor ancestral rates
    handle.write('RateAncestor = 0\n')
    #   Set the convergence tolerance
    handle.write('Small_Diff = 1e-6\n')
    #   Keep gapping structure - do not remove sites with gaps
    handle.write('cleandata = 0\n')
    #   Do not fix branch lengths
    handle.write('fix_blength = 0\n')
    #   Estimate all at once
    handle.write('method = 0\n')
    handle.close()
    return


def main(tree_dir, tandems, outdir):
    """Main function."""
    fd, trees = list_trees(tree_dir)
    og_tandem = parse_assignment(tandems)
    for og in sorted(og_tandem):
        t_file = os.path.join(fd, 'RAxML_bestTree.' + og)
        if not os.path.isfile(t_file):
            continue
        else:
            tree = parse_tree(t_file)
            tandem_names = find_tandems(tree, og_tandem[og])
            # Mark the tree for the four models
            lab_tree_null = label_branches(tree, tandem_names, 'H0')
            lab_tree_ha1 = label_branches(tree, tandem_names, 'Ha1')
            lab_tree_ha2 = label_branches(tree, tandem_names, 'Ha2')
            lab_tree_ha3 = label_branches(tree, tandem_names, 'Ha3')
            # Then make the directory structure
            og_out = make_dirs(og, outdir)
            # Write the trees
            write_marked(lab_tree_null, og, 'Null', og_out)
            write_marked(lab_tree_ha1, og, 'Ha1', og_out)
            write_marked(lab_tree_ha2, og, 'Ha2', og_out)
            write_marked(lab_tree_ha3, og, 'Ha3', og_out)
            # Write the control files
            write_ctl(og, 'Null', og_out)
            write_ctl(og, 'Ha1', og_out)
            write_ctl(og, 'Ha2', og_out)
            write_ctl(og, 'Ha3', og_out)
    return


if len(sys.argv) != 4:
    print """Mark a Newick tree for testing with CodeML. This script is for estimating
functional outcomes of maize tandem duplicates. The tree will be marked with
several schemes to test relative rates of duplicates. See the markdown file for
details on how they will be marked. Takes three arguments:
    1) Directory of trees
    2) Orthogroup-Tandem ID table
    3) Output directory"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2], sys.argv[3])
