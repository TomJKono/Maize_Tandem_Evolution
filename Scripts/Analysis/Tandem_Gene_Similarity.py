#!/usr/bin/env python
"""Script to generate alignments for tandem duplicates. Translates the sequence
to amino acids, aligns them with clustal-omega, then back-translates them. Takes
three arguments:
    1) Tandem duplicates CSV
    2) Representative transcripts FASTA
    3) Output directory
"""

import sys
import tempfile
import os
import itertools
from Bio import SeqIO
from Bio import SeqRecord
from Bio import Seq
from Bio.Align.Applications import ClustalOmegaCommandline


def parse_transcripts(trans):
    """Return a dictionary of sequences."""
    s = SeqIO.parse(trans, 'fasta')
    seq_dict = SeqIO.to_dict(s)
    # Remove the _whatever at the end
    seq_dict_nosuff = {}
    for seqid in seq_dict:
        seq_dict_nosuff[seqid.split('_')[0]] = seq_dict[seqid]
    return seq_dict_nosuff


def get_cds(geneid, seqdict):
    """Return the amino acid sequence of a gene with a given ID."""
    nuc_seq = seqdict[geneid]
    # Translate it
    aa_seq = nuc_seq.seq.translate()
    # Decorate it like you would a full SeqRecord object
    aa_seq_rec = SeqRecord.SeqRecord(
        aa_seq,
        id=geneid,
        description='')
    return aa_seq_rec


def align_genes(gene1, gene2):
    """Align the two genes with clustal-omega"""
    # Make temp files for clustal in and out
    clust_in = tempfile.NamedTemporaryFile(
        prefix='CO_in_',
        suffix='.fasta',
        mode='w+t')
    clust_out = tempfile.NamedTemporaryFile(
        prefix='CO_out_',
        suffix='.fasta',
        mode='w+t')
    # Write the sequences into the temp file
    SeqIO.write([gene1, gene2], clust_in, 'fasta')
    # Seek to the beginning else the file will appear empty
    clust_in.seek(0)
    # Run the command
    cline = ClustalOmegaCommandline(
        infile=clust_in.name,
        outfile=clust_out.name,
        seqtype='protein',
        force=True,
        iterations=10,
        distmat_full=True,
        distmat_full_iter=True)
    cline()
    clust_in.close()
    # Return the handle to the output file
    return clust_out


def back_translate(aln_file, seqdict):
    """Back-translate the aligned amino acids, using the original CDS
    nucleotide sequences."""
    aln = SeqIO.parse(aln_file.name, 'fasta')
    bt_seq = []
    for prot_seq in aln:
        codon = 0
        bt = ''
        nuc = seqdict[prot_seq.id]
        for aa in prot_seq:
            if aa == '-':
                bt += '---'
            else:
                bt += nuc[codon*3:(codon*3)+3]
                codon += 1
        bt_seq.append(bt)
    return bt_seq


def write_alignment(nuc_aln, outdir):
    """Write the alignment in FASTA format to the alignment directory."""
    # Get the gene IDs from the alignment data
    gids = [s.id.split('_')[0] for s in nuc_aln]
    # Join them to make a filename
    fname = '-'.join(gids)
    # Generate an output filename
    abs_outdir = os.path.abspath(os.path.expanduser(outdir))
    outname = os.path.join(abs_outdir, fname + '.fasta')
    # Then write the data into the file
    handle = open(outname, 'w')
    SeqIO.write(nuc_aln, handle, 'fasta')
    handle.close()
    # Write a lil message to stderr that says we finished
    sys.stderr.write('Wrote ' + fname + '\n')
    return


def main(tandem, transcripts, outdir):
    """Main function."""
    tx_seqs = parse_transcripts(transcripts)
    with open(tandem, 'r') as f:
        for line in f:
            # We want to take two genes at a time.
            genes = line.strip().split(',')
            gene_pairs = itertools.combinations(genes, 2)
            for gp in gene_pairs:
                g1 = get_cds(gp[0], tx_seqs)
                g2 = get_cds(gp[1], tx_seqs)
                aln = align_genes(g1, g2)
                nuc_aln = back_translate(aln, tx_seqs)
                write_alignment(nuc_aln, outdir)
                # Be good, and clean up our open handles
                aln.close()
    return


main(sys.argv[1], sys.argv[2], sys.argv[3])
