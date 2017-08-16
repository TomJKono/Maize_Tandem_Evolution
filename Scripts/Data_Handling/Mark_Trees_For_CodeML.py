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
    from ete3 import Tree
except ImportError:
    print 'This script requires the ete3 toolkit to be installed.'
    exit(1)


def list_trees(d):
    """List a directory and return paths to each tree file."""
    fd = os.path.abspath(os.path.expanduser(d))
    tree_paths = [t
        for t
        in os.listdir(fd)
        if t.endswith('.tree')]
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
    treeobj = Tree(t)
    treeobj.render('')
    return treeobj


def label_branches(tree, tandem_ids):
    """Use the tandem gene IDs to apply labels to the branches that go to maize
    tandem duplicates. We also want to apply labels to the internal branches
    that go to only tandem duplicates. Also use the naming scheme of the maize
    genes to label the maize branches. We leave every other branch alone to
    serve as the background. We do not distinguish between B73 and PH207,
    because the data are not appropriate for that distinction."""
    pass


def main(tree_dir, tandems):
    """Main function."""
    fd, trees = list_trees(tree_dir)
    og_tandem = parse_assignment(tandems)
    for og in sorted(og_tandem):
        t_file = os.path.join(fd, og + '.tree')
        tree = parse_tree(t_file)
        lab_tree = label_branches(tree, og_tandem[og])
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
