#!/usr/bin/env python
"""Script to generate alignments for adjacent genes. Translates the sequence to
amino acids, aligns them with clustal-omega, then back-translates them. This is
to find the mean pairwise similarity of adjacent genes. Takes three arguments:
    1) Sorted genes GFF
    2) Representative transcripts FASTA
    3) Output directory
"""

import sys
import tempfile
import os
import multiprocessing
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


def parse_gff(gff):
    """Parse the GFF and return a dictionary of lists of the structure
    {
        chr1: [gene1, gene2, gene3,...].
        ...
    }
    Assumes that the gnees are sorted."""
    gff_dat = {}
    with open(gff, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            gene_id = tmp[8].split(';')[0][3:]
            chrom = tmp[0]
            if chrom not in gff_dat:
                gff_dat[chrom] = [gene_id]
            else:
                gff_dat[chrom].append(gene_id)
    # Filter to those that only have at least two genes, since that will mess
    # up our calculations.
    gff_dat_flt = {}
    for chrom in gff_dat:
        if len(gff_dat[chrom]) < 2:
            continue
        else:
            gff_dat_flt[chrom] = gff_dat[chrom]
    return gff_dat_flt


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


def work_func(args):
    """Function that does the work by calling the functions above. We do it this
    way because the multiprocessing module requires it."""
    # Unpack the arguments
    gene1, gene2, tx, output = args
    # Extract the gene sequences
    g1 = get_cds(gene1, tx)
    g2 = get_cds(gene2, tx)
    # Align them
    p_aln = align_genes(g1, g2)
    # Backtranslate
    n_aln = back_translate(p_aln, tx)
    # Write them into the output directory
    write_alignment(n_aln, output)
    # Close those handles
    p_aln.close()
    return


if __name__ == '__main__':
    gff = sys.argv[1]
    transcripts = sys.argv[2]
    outdir = sys.argv[3]
    tx_seqs = parse_transcripts(transcripts)
    sorted_genes = parse_gff(gff)
    # Build up the list of arguments to send
    worker_args = []
    for chrom in sorted(sorted_genes):
        for index, gene in enumerate(sorted_genes[chrom][:-1]):
            worker_args.append((gene, sorted_genes[chrom][index+1], tx_seqs, outdir))
    # Then send them to the map function
    p = multiprocessing.Pool()
    p.map(work_func, worker_args)
