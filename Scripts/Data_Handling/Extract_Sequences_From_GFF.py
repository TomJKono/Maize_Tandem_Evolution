#!/usr/bin/env python
"""Generates a FASTA file with sequences for each feature listed in a GFF. The
GFF that contains the TE annotations has strands that are not part of the
original GFF specification. In particular, some are asterisks (*). These will
be treated as + features. Takes two arguments:
    1) GFF file with features that are to have sequences
    2) Reference sequence associated with the GFF
"""

import sys

try:
    from Bio import SeqIO
    from Bio import SeqFeature
    from Bio.SeqFeature import FeatureLocation
except ImportError:
    print "This script requires the Biopython library to be installed."""
    exit(1)


def read_ref(fasta):
    """Read the supplied FASTA file and return it as a dictionary, keyed on
    chromosome or contig name."""
    f = SeqIO.parse(fasta, 'fasta')
    return SeqIO.to_dict(f)


def read_gff(gff):
    """Read the supplied GFF file and return it as a dictionary with the
    following structure:
    {
        feature_name: (chromosome, SeqFeature),
        feature_name: (chromosome, SeqFeature),
        ...
    }
    """
    seqfeat = {}
    with open(gff, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                tmp = line.strip().split('\t')
                chromosome = tmp[0]
                start = int(tmp[3])
                stop = int(tmp[4])
                if tmp[6] == '*' or tmp[6] == '+':
                    strand = 1
                elif tmp[6] == '-':
                    strand = -1
                metadata = tmp[8]
                #   Split the metadata field on semicolons and iterate through
                #   the fields. We want to save the ID field, to use as the
                #   sequece name
                for m in metadata.split(';'):
                    if m.startswith('ID='):
                        featname = m[3:]
                        break
                #   Build the SeqFeature object
                feat = SeqFeature.SeqFeature(
                    FeatureLocation(start, stop, strand=strand),
                    id=featname)
                #   And put it in the dictionary
                seqfeat[featname] = (chromosome, feat)
    return seqfeat


def main(gff, fasta):
    """Main function."""
    refseq = read_ref(fasta)
    features = read_gff(gff)
    for feat in features.iteritems():
        c = feat[1][0]
        f = feat[1][1]
        if c not in refseq:
            continue
        teseq = f.extract(refseq[c])
        teseq.name = feat[0]
        teseq.id = feat[0]
        teseq.description = ''
        SeqIO.write(teseq, sys.stdout, 'fasta')
    return


main(sys.argv[1], sys.argv[2])
