#!/usr/bin/env python
"""Convert from FASTA to NEXUS format. This is designed for the maize tandem
duplicate evolution project, where we have a directory of FASTA alignments of
potential tandem duplicates. We convert them to NEXUS files for divergence
dating with BEAST. Takes two arguments:
    1) Directory of FASTA files
    2) Output directory
"""

import sys
import os
try:
    from Bio import SeqIO
    from Bio import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import IUPAC
except ImportError:
    print 'You need the Biopython library installed to use this script.'
    exit(1)


def list_dir(fasta_dir):
    """Return a list of the FASTA files in the input directory."""
    fdir = os.path.abspath(os.path.expanduser(fasta_dir))
    contents = os.listdir(fasta_dir)
    full_contents = [
        os.path.join(fdir, f)
        for f
        in contents
        if f.endswith('fasta')]
    return full_contents


def conv_nexus(fasta_file, out_dir):
    """Convert the FASTA file to a NEXUS file."""
    fname = os.path.basename(fasta_file)
    fname = fname.replace('fasta', 'nex')
    odir = os.path.abspath(os.path.expanduser(out_dir))
    oname = os.path.join(odir, fname)
    # Then read the sequnce data to convert it to a NEXUS file
    fastaseq = list(SeqIO.parse(fasta_file, 'fasta'))
    # We have to re-define the SeqRecord because there is no alphabet associated
    # with them right now.
    newseq = []
    for s in fastaseq:
        n = Seq.Seq(str(s.seq), alphabet=IUPAC.unambiguous_dna)
        newseq.append(
            SeqRecord(
                n,
                id=s.id,
                description=''))
    handle = open(oname, 'w')
    SeqIO.write(newseq, handle, 'nexus')
    handle.close()
    return


def main(fasta_dir, output_dir):
    """Main function."""
    fasta_seqs = list_dir(fasta_dir)
    for s in fasta_seqs:
        conv_nexus(s, output_dir)
    return


if len(sys.argv) != 3:
    print """
Convert from FASTA to NEXUS format. This is designed for the maize tandem
duplicate evolution project, where we have a directory of FASTA alignments of
potential tandem duplicates. We convert them to NEXUS files for divergence
dating with BEAST. Takes two arguments:
    1) Directory of FASTA files
    2) Output directory"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2])
