#!/usr/bin/env python
"""Generate alignments of nonsyntenic duplications in B73 and PH207. Takes
one argument:
    1) Output directory.
"""

import sys
import os
import subprocess
import tempfile
try:
    from Bio import SeqIO
    from Bio import SeqRecord
    from Bio import Seq
    from Bio.Align.Applications import ClustalOmegaCommandline
    from Bio.Seq import IUPAC
except ImportError:
    print 'This script requires Biopython.'
    exit(1)

# Define paths to files that house the data as constants
B_CLUST = '/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Filtering/B73_True_Tandem_Clusters.txt'
P_CLUST = '/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Filtering/PH207_True_Tandem_Clusters.txt'
GENE_KEY = '/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Data/ABB_Synteny/B73-PH207_gene-key.txt'
B_CDS = '/Volumes/LaCie/Maize_Tandem_Evolution/References/B73-rep-cDNA.fa'
P_CDS = '/Volumes/LaCie/Maize_Tandem_Evolution/References/ZmaysPH207_443_v1.1.cds.fa'


def read_seqs(b, p):
    """Parse the sequence IDs of the CDS files and return a dictionary that
    relates gene_id => transcript_id."""
    gdict = {}
    b_seqs = [s.id for s in SeqIO.parse(b, 'fasta')]
    p_seqs = [s.id for s in SeqIO.parse(p, 'fasta')]
    for t in b_seqs:
        g = t.split('_')[0]
        gdict[g] = t
    for t in p_seqs:
        g = t.split('_')[0]
        gdict[g] = t
    return gdict


def parse_clusters(clst):
    """Parse the cluster assignment, and store as a gene=>cluster dictionary."""
    clust = {}
    with open(clst, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            c_id = tmp[0]
            genes = tmp[1].split(',')
            for g in genes:
                clust[g] = c_id
    return clust


def parse_nonsyn(nonsyn):
    """Parse the gene key, and save only the nonsyntenic relationships. The
    syntenic relationships will have the value 'curated-syntenic' in them. Every
    other line is nonsyntenic."""
    ns = []
    with open(nonsyn, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            rel_class = tmp[8]
            if rel_class == 'curated-syntenic':
                continue
            else:
                b_gene = tmp[3]
                p_gene = tmp[7]
                ns.append([b_gene, p_gene])
    # Sort it on the B73 gene
    ns = sorted(ns, key=lambda x: x[0])
    return ns


def identify_homologous_dups(bdup, pdup, gkey):
    """Using the nonsyntenic gene key, identify clusters of genes that are
    homologous outside of syntenic regions."""
    # Print a header
    print 'Cluster_ID\tB_Cluster\tP_Cluster\tB_Genes\tP_Genes'
    hom = []
    # Iterate over relationships
    c_num = 0
    for rel in gkey:
        b = rel[0]
        p = rel[1]
        # If the B73 and PH207 genes are duplicated, then we process the
        # clusters
        b_clust = bdup.get(b, None)
        p_clust = pdup.get(p, None)
        # If there are no clusters in B or P, skip it
        if not any([b_clust, p_clust]):
            continue
        else:
            # Then, if the genes are duplicated, pull those in
            b_genes = []
            p_genes = []
            if b_clust:
                for g, c in bdup.iteritems():
                    if b_clust == c:
                        b_genes.append(g)
            else:
                b_clust = 'NA'
            if p_clust:
                for g, c in pdup.iteritems():
                    if p_clust == c:
                        p_genes.append(g)
            else:
                p_clust = 'NA'
            # Add the syntenic genes to the cluster, but only if they are
            # annotated gene models, and not fills,nor NA
            # Build the syntenic cluster genes
            if b.startswith('Zm') and b not in b_genes:
                b_genes.append(b)
            if p.startswith('Zm') and p not in p_genes:
                p_genes.append(p)
            nonsyn_clust = b_genes + p_genes
            # If the nonsyn cluster is already in the list of clusters that
            # we are keeping track of, we skip it
            if nonsyn_clust in hom:
                continue
            else:
                # Print out the ancestral and cluster assignments, for reference!
                # If b_genes or p_genes is empty, drop in NA
                if not b_genes:
                    b_genes = ['NA']
                if not p_genes:
                    p_genes = ['NA']
                c_name = 'Nonsyntenic_Dup_' + str(c_num).zfill(4)
                print '\t'.join([
                    c_name,
                    b_clust,
                    p_clust,
                    ','.join(b_genes),
                    ','.join(p_genes)])
                # Tack them onto the list of syntenic clusters
                hom.append(nonsyn_clust)
                c_num += 1
    return hom


def extract_sequences(gene_ids, trans):
    """Extract gene sequences given their IDs. Will perform the necessary
    transformations to make sure that we don't have weird sequence-not-found
    issues."""
    # get the CDS IDs from the gene IDs
    c_ids = [trans[x] for x in gene_ids if x in trans]
    # Figure out which  genes are from which FASTA
    b_ids = [c for c in c_ids if c.startswith('Zm00001d')]
    p_ids = [c for c in c_ids if c.startswith('Zm00008a')]
    # Open a temporary file to save the sequences
    unaligned = tempfile.NamedTemporaryFile(
        prefix='Nonsyntenic_Tandem_',
        suffix='.fasta',
        mode='w+t')
    # And build the commands
    b_cmd = ['samtools', 'faidx', B_CDS] + b_ids
    p_cmd = ['samtools', 'faidx', P_CDS] + p_ids
    # Execute the commands! Save the output in the temp file
    subprocess.call(b_cmd, shell=False, stdout=unaligned)
    subprocess.call(p_cmd, shell=False, stdout=unaligned)
    # Seek to the beginning so we can read the sequence data out of it later
    unaligned.seek(0)
    # Return the handle
    return unaligned


def align(sequences):
    """Translate, then align, then back-translate the sequences."""
    # First, start a new tempfile
    translated = tempfile.NamedTemporaryFile(
        prefix='Translated_',
        suffix='.fasta',
        mode='w+t')
    # And write the translated sequences into it
    nuc_seqs = {}
    for seq in SeqIO.parse(sequences, 'fasta'):
        s = SeqRecord.SeqRecord(
            seq.seq.translate(),
            id=seq.id,
            description='')
        SeqIO.write(s, translated, 'fasta')
        nuc_seqs[s.id] = str(seq.seq)
    # Seek to the beginning to read them with Clustal Omega
    translated.seek(0)
    # Open a temp file for clustal output
    aligned = tempfile.NamedTemporaryFile(
        prefix='Aligned_',
        suffix='.fasta',
        mode='w+t')
    # And align them
    co_cmd = ClustalOmegaCommandline(
        infile=translated.name,
        outfile=aligned.name,
        seqtype='protein',
        force=True,
        iterations=10,
        distmat_full=True,
        distmat_full_iter=True)
    co_cmd()
    # Close the translated unaligned handle. We are done with it
    translated.close()
    # Then, we want to back-translate the sequences
    backtrans = tempfile.NamedTemporaryFile(
        prefix='Backtranslated_',
        suffix='.fasta',
        mode='w+t')
    aligned.seek(0)
    aln = SeqIO.parse(aligned.name, 'fasta')
    for prot_seq in aln:
        bt = ''
        codon = 0
        nuc = nuc_seqs[prot_seq.id]
        for aa in prot_seq:
            if aa == '-':
                bt += '---'
            else:
                bt += nuc[codon*3:(codon*3)+3]
                codon += 1
        # Make it a SeqRecord to write into the output file
        b = SeqRecord.SeqRecord(Seq.Seq(bt), id=prot_seq.id, description='')
        SeqIO.write(b, backtrans, 'fasta')
    # Again, seek to the beginning so we can read it later
    backtrans.seek(0)
    # Close the aligned and unaligned handles; we are done with them
    aligned.close()
    sequences.close()
    return backtrans


def write_aln(alignment, outdir, outname):
    """Write the alignment as FASTA format into the output directory."""
    outpath = os.path.abspath(os.path.expanduser(outdir))
    outname = os.path.join(outpath, outname)
    # Parse the sequence data again, convert to NEXUS
    newseq = []
    for s in SeqIO.parse(alignment, 'fasta'):
        n = Seq.Seq(str(s.seq), alphabet=IUPAC.unambiguous_dna)
        newseq.append(
            SeqRecord.SeqRecord(
                n,
                id=s.id,
                description=''))
    handle = open(outname, 'w')
    SeqIO.write(newseq, handle, 'fasta')
    handle.close()
    # We are also done with the alignment handle
    alignment.close()
    return


def main(outdir):
    """Main function."""
    gene_name_trans = read_seqs(B_CDS, P_CDS)
    b73_clusters = parse_clusters(B_CLUST)
    ph207_clusters = parse_clusters(P_CLUST)
    ns_gene_key = parse_nonsyn(GENE_KEY)
    hom_clusters = identify_homologous_dups(b73_clusters, ph207_clusters, ns_gene_key)
    for index, cluster in enumerate(hom_clusters):
        # Build an output name based on the index
        outname = 'Nonsyntenic_Dup_' + str(index).zfill(4) + '.fasta'
        un = extract_sequences(cluster, gene_name_trans)
        if not un:
            continue
        bt_aln = align(un)
        write_aln(bt_aln, outdir, outname)
    return

if len(sys.argv) != 2:
    print """Generate alignments of nonsyntenic duplications in B73 and PH207. Takes
one argument:
    1) Output directory."""
    exit(1)
else:
    main(sys.argv[1])
