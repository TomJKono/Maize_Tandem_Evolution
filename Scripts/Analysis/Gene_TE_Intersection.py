#!/usr/bin/env python
"""Generate a table that counts genes in TEs, and TEs in genes, for each major
class of TE. Will separate tandem from non-tandem, and not double-count
genes for nested TEs. Takes two arguments:
    1) B73 tandem duplicate clusters
    2) Directory of TE-Gene intersection results (See workflow)
"""

import sys
import os


def parse_clusters(clust):
    """Return a list of tandem duplicate gene IDs."""
    c = []
    with open(clust, 'r') as f:
        for line in f:
            genes = line.strip().split('\t')[1]
            for g in genes.split(','):
                c.append(g)
    return c


def find_files(tedir):
    """Return the paths to the various files that are required. There should be
    10 files."""
    file_dict = {}
    # Generate the full path
    fpath = os.path.abspath(os.path.expanduser(tedir))
    for fname in os.listdir(fpath):
        # Split the filename on underscores - we are interested in the first
        # and third parts
        parts = fname.split('_')
        outer = parts[0]
        inner =  parts[2].split('.')[0]
        if parts[0] not in file_dict:
            file_dict[outer] = {inner: os.path.join(fpath, fname)}
        else:
            file_dict[outer].update({inner: os.path.join(fpath, fname)})
    return file_dict


def separate_tandems(overlap_file, tandem_ids):
    """Parse the overlap file (which looks like a GFF) and return a count of the
    number of genes that are tandem and that are not tandem."""
    tandem = []
    non_tandem = []
    with open(overlap_file, 'r') as f:
        for line in f:
            tmp = line.strip().split('\t')
            # Get the gene ID out of the 9th field. It is labeled as 'ID=...'
            # so remove the first three characters
            gene_id = tmp[8].split(';')[0][3:]
            if gene_id in tandem_ids:
                if gene_id in tandem:
                    continue
                else:
                    tandem.append(gene_id)
            else:
                if gene_id in non_tandem:
                    continue
                else:
                    non_tandem.append(gene_id)
    return (tandem, non_tandem)


def main(clust, tedir):
    """Main function."""
    tandems = parse_clusters(clust)
    files = find_files(tedir)
    for inner in sorted(files):
        for outer in sorted(files[inner]):
            key = outer + '_in_' + inner
            t, nt = separate_tandems(files[inner][outer], tandems)
            print key + '\t' + str(len(t)) + '\t' + str(len(nt))
    return


if len(sys.argv) != 3:
    print """Generate a table that counts genes in TEs, and TEs in genes, for each major
class of TE. Will separate tandem from non-tandem, and not double-count
genes for nested TEs. Takes two arguments:
    1) B73 tandem duplicate clusters
    2) Directory of TE-Gene intersection results"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2])
