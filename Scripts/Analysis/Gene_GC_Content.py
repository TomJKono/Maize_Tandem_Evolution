#!/usr/bin/env python
"""Calculate GC fraction of all genes in a FASTA file. Requires Biopython.
Takes one argument:
    1) FASTA file of nucleotide seqeucnes.
"""

import sys
try:
    from Bio import SeqIO
except ImportError:
    print 'This script requires Biopython.'
    exit(1)

if len(sys.argv) != 2:
    print """Calculate GC fraction of all genes in a FASTA file. Requires Biopython.
Takes one argument:
    1) FASTA file of nucleotide seqeucnes."""
    exit(1)
else:
    seqs = SeqIO.parse(sys.argv[1], 'fasta')
    for s in seqs:
        sname = s.id.split('_')[0]
        slen = float(len(s))
        sgc = s.seq.count('G') + s.seq.count('C')
        print sname, sgc/slen
