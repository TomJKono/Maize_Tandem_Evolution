#!/usr/bin/env python
"""Subset a GFF based on a list of suppled transcript IDs. We do it this way
instead of with grep because grep is very slow when comparing files of this
size. Takes two arguments:
    1) Tx ID file
    2) GFF to subset
"""

import sys


def parse_tx(transcript_ids):
    """Parse the transcript IDs and return a list of gene IDs and transcript IDs
    to search in the GFF."""
    genes = []
    txs = []
    with open(transcript_ids, 'r') as f:
        for line in f:
            tmp = line.strip()
            genes.append(tmp.split('_')[0])
            txs.append(tmp)
    return (genes, txs)


def trim_gff(gff, genes, transcripts):
    """Read through the GFF and print out the annotations for the longest
    transcripts in each gene."""
    with open(gff, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                tmp = line.strip().split('\t')
                # Print out the gene lines if they have transcripts we are
                # interested in
                if tmp[2] == 'gene':
                    geneid = tmp[8].split(';')[0].split('=')[1]
                    if geneid in genes:
                        print line.strip()
                    else:
                        continue
                elif tmp[2] in ['CDS', 'mRNA']:
                    # For other gene parts, we want to print them out, only if
                    # they are in the list of longest transcript IDs.
                    feat_id = tmp[8].split(';')[0].split('=')[1]
                    if feat_id in transcripts:
                        print line.strip()
                    else:
                        continue
                elif tmp[2] == 'exon':
                    parent = tmp[8].split(';')[1].split('=')[1]
                    if parent in transcripts:
                        print line.strip()
                    else:
                        continue
    return


def main(tx, gff):
    """Main function."""
    g_ids, tx_ids = parse_tx(tx)
    trim_gff(gff, g_ids, tx_ids)
    return


if len(sys.argv) != 3:
    print """Subset a GFF based on a list of suppled transcript IDs. We do it this way
instead of with grep because grep is very slow when comparing files of this
size. Takes two arguments:
    1) Tx ID file
    2) GFF to subset"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2])
