#!/usr/bin/env python
"""Generate alignments from syntenic tandem duplications in B73 and PH207, but
keep track of the cases where there are syntenic single-copy genes. Takes one
argument:
    1) Output directory
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
    print 'This script requires Biopython'
    exit(1)

# Define constants. These are paths to files that house the data.
B_CLUST = '/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Filtering/B73_True_Tandem_Clusters.txt'
P_CLUST = '/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Filtering/PH207_True_Tandem_Clusters.txt'
SYNTENY = '/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Data/ABB_Synteny/b73-ph207-synteny-bt2corrected-nah-coord.txt'
B_CDS = '/Volumes/LaCie/Maize_Tandem_Evolution/References/B73-rep-cDNA.fa'
P_CDS = '/Volumes/LaCie/Maize_Tandem_Evolution/References/ZmaysPH207_443_v1.1.cds.fa'
SB_CDS = '/Volumes/LaCie/Maize_Tandem_Evolution/References/Sbicolor_313_v3.1.cds.fa'
OS_CDS = '/Volumes/LaCie/Maize_Tandem_Evolution/References/Osativa_323_v7.0.cds.fa'


def read_seqs(b, p, sorg, oryz):
    """Parse the sqeuence IDs of the CDS files and return a dictionary that
    gives gene_id => transcript ID."""
    gdict = {}
    b_seqs = [s.id for s in SeqIO.parse(b, 'fasta')]
    p_seqs = [s.id for s in SeqIO.parse(p, 'fasta')]
    sb_seqs = [s.id for s in SeqIO.parse(sorg, 'fasta')]
    os_seqs = [s.id for s in SeqIO.parse(oryz, 'fasta')]
    # Then from each sequecne ID (which should be CDS IDs), derive the
    # gene ID and put it in the dictionary
    # B73 CDS are named like GENE_T00X
    for t in b_seqs:
        g = t.split('_')[0]
        gdict[g] = t
    # PH207 CDS are named like GENE_C0X
    for t in p_seqs:
        g = t.split('_')[0]
        gdict[g] = t
    # Sorghum bicolor CDS are named like Sobic.Gene.X
    # and there are multiple transcripts listed per gene. We just take the first
    # one.
    for t in sb_seqs:
        if t in gdict:
            continue
        else:
            g = t.rsplit('.', 1)[0]
            gdict[g] = t
    # Oryza sativa CDS are named like LOC_GENE.X
    # and there are multiple transcripts listed per gene. Take just the first.
    for t in os_seqs:
        if t in gdict:
            continue
        else:
            g = t.rsplit('.')[0]
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


def parse_synteny(syn):
    """Parse the synteny file and return it as a big ugly list. We do it this
    way because we will have to look up by either of the B73 or PH207 columns,
    and a dictionary will not help us here."""
    s = []
    with open(syn, 'r') as f:
        for line in f:
            s.append(line.strip().split()[0:5])
    # sort on the ancestral gene
    s = sorted(s, key=lambda x: x[1])
    return s


def identify_homologous_dups(bdup, pdup, syn):
    """Using the synteny table, identify clusters of genes that are homologous.
    Return them as a list of gene IDs."""
    # Print a header
    print 'Cluster_ID\tAnc_Gene\tB1_Cluster\tB2_Cluster\tP1_Cluster\tP2_Cluster\tB1_Genes\tB2_Genes\tP1_Genes\tP2_Genes'
    hom = []
    # Iterate over syntenic relationships
    c_num = 0
    for rel in syn:
        anc = rel[0]
        b1 = rel[1]
        b2 = rel[2]
        p1 = rel[3]
        p2 = rel[4]
        # If the B73 and PH207 genes are duplicated, then we process the
        # clusters
        b1_clust = bdup.get(b1, None)
        b2_clust = bdup.get(b2, None)
        p1_clust = pdup.get(p1, None)
        p2_clust = pdup.get(p2, None)
        # If there are no clusters in B1, B2, P1, or P2, skip it
        if not any([b1_clust, b2_clust, p1_clust, p2_clust]):
            continue
        else:
            # Then, if the syntenic genes are duplicated, pull those in
            b1_genes = []
            b2_genes = []
            p1_genes = []
            p2_genes = []
            if b1_clust:
                for g, c in bdup.iteritems():
                    if b1_clust == c:
                        b1_genes.append(g)
            else:
                b1_clust = 'NA'
            if b2_clust:
                for g, c in bdup.iteritems():
                    if b2_clust == c:
                        b2_genes.append(g)
            else:
                b2_clust = 'NA'
            if p1_clust:
                for g, c in pdup.iteritems():
                    if p1_clust == c:
                        p1_genes.append(g)
            else:
                p1_clust = 'NA'
            if p2_clust:
                for g, c in pdup.iteritems():
                    if p2_clust == c:
                        p2_genes.append(g)
            else:
                p2_clust = 'NA'
            # Add the syntenic genes to the cluster, but only if they are
            # annotated gene models, and not fills,nor NA
            # Build the syntenic cluster genes
            if b1.startswith('Zm') and b1 not in b1_genes:
                b1_genes.append(b1)
            if b2.startswith('Zm') and b2 not in b2_genes:
                b2_genes.append(b2)
            if p1.startswith('Zm') and p1 not in p1_genes:
                p1_genes.append(p1)
            if p2.startswith('Zm') and p2 not in p2_genes:
                p2_genes.append(p2)
            syn_clust = b1_genes + b2_genes + p1_genes + p2_genes + [anc]
            # Print out the ancestral and cluster assignments, for reference!
            # If b_genes or p_genes is empty, drop in NA
            if not b1_genes:
                b1_genes = ['NA']
            if not b2_genes:
                b2_genes = ['NA']
            if not p1_genes:
                p1_genes = ['NA']
            if not p2_genes:
                p2_genes = ['NA']
            c_name = 'Syntenic_Dup_' + str(c_num).zfill(4)
            print '\t'.join([
                c_name,
                anc,
                b1_clust,
                b2_clust,
                p1_clust,
                p2_clust,
                ','.join(b1_genes),
                ','.join(b2_genes),
                ','.join(p1_genes),
                ','.join(p2_genes)])
            # Tack them onto the list of syntenic clusters
            hom.append(syn_clust)
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
    sb_ids = [c for c in c_ids if c.startswith('Sobic')]
    os_ids = [c for c in c_ids if c.startswith('LOC_Os')]
    # Open a temporary file to save the sequences
    unaligned = tempfile.NamedTemporaryFile(
        prefix='Syntenic_Tandem_',
        suffix='.fasta',
        mode='w+t')
    # And build the commands
    b_cmd = ['samtools', 'faidx', B_CDS] + b_ids
    p_cmd = ['samtools', 'faidx', P_CDS] + p_ids
    if sb_ids:
        anc_cmd = ['samtools', 'faidx', SB_CDS] + sb_ids
    elif os_ids:
        anc_cmd = ['samtools', 'faidx', OS_CDS] + os_ids
    else:
        anc_cmd = []
    # Execute the commands! Save the output in the temp file
    subprocess.call(b_cmd, shell=False, stdout=unaligned)
    subprocess.call(p_cmd, shell=False, stdout=unaligned)
    if not anc_cmd:
        return None
    else:
        subprocess.call(anc_cmd, shell=False, stdout=unaligned)
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
    gene_name_trans = read_seqs(B_CDS, P_CDS, SB_CDS, OS_CDS)
    b73_clusters = parse_clusters(B_CLUST)
    ph207_clusters = parse_clusters(P_CLUST)
    synteny = parse_synteny(SYNTENY)
    hom_clusters = identify_homologous_dups(b73_clusters, ph207_clusters, synteny)
    for index, cluster in enumerate(hom_clusters):
        # Build an output name based on the index
        outname = 'Syntenic_Dup_' + str(index).zfill(4) + '.fasta'
        un = extract_sequences(cluster, gene_name_trans)
        if not un:
            continue
        bt_aln = align(un)
        write_aln(bt_aln, outdir, outname)
    return


if len(sys.argv) != 2:
    print """Generate alignments from syntenic tandem duplications in B73 and PH207, but
keep track of the cases where there are syntenic single-copy genes. Takes one
argument:
    1) Output directory"""
    exit(1)
else:
    main(sys.argv[1])
