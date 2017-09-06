#!/usr/bin/env python
"""Filter sites in a nucleotide alignment that have missing data greater than
a supplied proportion. Sites failing the filter are dropped from the alignment.
Takes two arguments:
    1) Alignment
    2) Maxiumum proportion of missing data per site
"""

import sys
try:
    from Bio import AlignIO
    from Bio import SeqIO
    from Bio import Seq
    from Bio import SeqRecord
except ImportError:
    print 'This script requires Biopython.'
    exit(1)


def missing_cols(aln, prop):
    """Iterate over columns of the alignment and return indices for sites that
    have a proportion of missing data that is greater than the supplied
    proportion."""
    miss = []
    # How many sequences, and how many alignment columns?
    nseq = float(len(aln))
    ncol = len(aln[0])
    # For each column, get the proportion of missing data
    for column in range(0, ncol):
        nmiss = aln[:,column].count('-')
        propmiss = nmiss/nseq
        if propmiss > prop:
            miss.append(column)
    return miss


def filter_columns(aln, cols):
    """Remove columns from the alignment that fail the missing data filter."""
    flt_seqs = []
    for s in aln:
        t_s = Seq.Seq('')
        for index, base in enumerate(s):
            if index in cols:
                continue
            else:
                t_s += base
        flt_seqs.append(SeqRecord.SeqRecord(
            t_s,
            id=s.id,
            description=''))
    return flt_seqs


def main(aln, prop):
    """Main function."""
    try:
        fprop = float(prop)
    except ValueError:
        print 'Please supply a numeric proportion'
        exit(1)
    # Parse the FASTA alignment
    alignment = AlignIO.read(aln, 'fasta')
    mcols = missing_cols(alignment, fprop)
    flt_aln = filter_columns(alignment, mcols)
    # Write the filtered alignment to stdout, omitting sequences that are all
    # gapped sites.
    nogap = []
    for s in flt_aln:
        is_gap = [
            True
            if base == '-'
            else False
            for base in s.seq]
        if all(is_gap) or len(s.seq) == 0:
            continue
        else:
            nogap.append(s)
    SeqIO.write(nogap, sys.stdout, 'fasta')
    return


if len(sys.argv) != 3:
    print """Filter sites in a nucleotide alignmetn that have missing data greater than
a supplied proportion. Sites failing the filter are dropped from the alignment.
Takes two arguments:
    1) Alignment
    2) Maximum proportion of missing data per site"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2])
