#!/usr/bin/env python
"""Generate alignments for syntenic and nonsyntenic duplications in B73 and
PH207. Keep track of cases where there are single-copy genes. This script
differs from the syntenic and nonsyntenic specific ones in that it tries to
do a better job at distinguishing them: if there is a syentenic gene in the
tandem duplicate cluster, it is called syntenic. Else it is nonsyntenic. Takes
four arguments:
    1) Syntenic assignment output file
    2) Nonsyntenic assignment output file
    3) Syntenic alignment output directory
    4) Nonsyntenic alignment output directory.
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

# Define constants that are paths to the data files. We do these as constants
# rather than arguments because there are a lot of them.
B_CLUST = '/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Filtering/B73_True_Tandem_Clusters.txt'
P_CLUST = '/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Results/Filtering/PH207_True_Tandem_Clusters.txt'
SYNTENY = '/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Data/ABB_Synteny/b73-ph207-synteny-bt2corrected-nah-coord.txt'
GENE_KEY = '/Users/tomkono/Dropbox/GitHub/Maize_Tandem_Evolution/Data/ABB_Synteny/B73-PH207_gene-key.txt'
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
    s = sorted(s, key=lambda x: x[0])
    return s


def parse_key(genekey):
    """Parse the gene key and identify genes that are in syntenic relationships.
    We will return a dictionary of lists, with the key being the genotype, and
    the value bein gthe gene IDs that are syntenic."""
    syntenic = {'B73': [], 'PH207': []}
    with open(genekey, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            if tmp[8] == 'curated-syntenic' or tmp[8] == 'curated-syntenic-fusion':
                syntenic['B73'].append(tmp[3])
                syntenic['PH207'].append(tmp[7])
    return syntenic


def find_syntenic(dups, genekey, genotype):
    """Separate the clusters into syntenic and nonsyntenic, based on whether or
    not one of the genes appears in a syntenic relationship."""
    syntenic = []
    nonsyntenic = []
    with open(dups, 'r') as f:
        for line in f:
            tmp = line.strip().split('\t')
            c_id = tmp[0]
            genes = tmp[1].split(',')
            # Find out if any of the genes in the cluster show up in syntenic
            # relationships
            in_syn = [True if g in genekey[genotype] else False for g in genes]
            # If any of them do, then it's syntenic
            if any(in_syn):
                syntenic.append(c_id)
            else:
                nonsyntenic.append(c_id)
    return (syntenic, nonsyntenic)


def link_syntenic_dups(bclust, pclust, synteny, syn_output):
    """Read through the syntenic B73 duplicates and the PH207 duplicates and
    return lists of gene IDs that should be aligned. Also write the synteny
    assignment file."""
    handle = open(syn_output, 'w')
    handle.write('Duplicate_ID\tAnc_Gene\tB1_Cluster\tB2_Cluster\tP1_Cluster\tP2_Cluster\tB1_Genes\tB2_Genes\tP1_Genes\tP2_Genes\n')
    syn_clust = 0
    clusters = []
    for rel in synteny:
        # Unpack the syntenic relationship
        anc = rel[0]
        b1 = rel[1]
        b2 = rel[2]
        p1 = rel[3]
        p2 = rel[4]
        # Get clusters for the syntenic genes
        b1_clust = bclust.get(b1, None)
        b2_clust = bclust.get(b2, None)
        p1_clust = pclust.get(p1, None)
        p2_clust = pclust.get(p2, None)
        # Get the genes that are part of the clusters
        b1_genes = []
        b2_genes = []
        p1_genes = []
        p2_genes = []
        if b1_clust:
            for g in sorted(bclust):
                if bclust[g] == b1_clust:
                    b1_genes.append(g)
        if b2_clust:
            for g in sorted(bclust):
                if bclust[g] == b2_clust:
                    b2_genes.append(g)
        if p1_clust:
            for g in sorted(pclust):
                if pclust[g] == p1_clust:
                    p1_genes.append(g)
        if p2_clust:
            for g in sorted(pclust):
                if pclust[g] == p2_clust:
                    p2_genes.append(g)
        # If there are no tandem dups in either M1 or M2, then we continue
        if not any([b1_genes, b2_genes, p1_genes, p2_genes]):
            continue
        else:
            # Next, we need at add in the non-duplicated syntenic genes
            if b1 not in b1_genes and b1 != 'NA' and b1.startswith('Zm'):
                b1_genes.append(b1)
            if b2 not in b2_genes and b2 != 'NA' and b2.startswith('Zm'):
                b2_genes.append(b2)
            if p1 not in p1_genes and p1 != 'NA' and p1.startswith('Zm'):
                p1_genes.append(p1)
            if p2 not in p2_genes and p2 != 'NA' and p2.startswith('Zm'):
                p2_genes.append(p2)
            # Sort them
            b1_genes = sorted(b1_genes)
            b2_genes = sorted(b2_genes)
            p1_genes = sorted(p1_genes)
            p2_genes = sorted(p2_genes)
            # put all the genes together
            clusters.append([anc] + b1_genes + b2_genes + p1_genes + p2_genes)
            # Turn empty lists into NA
            if not b1_clust:
                b1_clust = 'NA'
            if not b2_clust:
                b2_clust = 'NA'
            if not p1_clust:
                p1_clust = 'NA'
            if not p2_clust:
                p2_clust = 'NA'
            if not b1_genes:
                b1_genes = ['NA']
            if not b2_genes:
                b2_genes = ['NA']
            if not p1_genes:
                p1_genes = ['NA']
            if not p2_genes:
                p2_genes = ['NA']
            # Write the table data
            dup_id = 'Syntenic_Dup_' + str(syn_clust).zfill(4)
            handle.write('\t'.join([
                dup_id,
                anc,
                b1_clust,
                b2_clust,
                p1_clust,
                p2_clust,
                ','.join(b1_genes),
                ','.join(b2_genes),
                ','.join(p1_genes),
                ','.join(p2_genes)]) + '\n')
            syn_clust += 1
    handle.flush()
    handle.close()
    return clusters


def link_nonsyntenic_dups(bclust, pclust, b_s, p_s, ns_out):
    """Assocate nonsyntenic duplicates like we would the syntenic dups. For each
    nonsyntenic cluster in B73, find its homologous gene and/or cluster in
    PH207. Unfortunately, we have to re-parse the gene key file to associate
    non-syntenic genes with each other."""
    handle = open(ns_out, 'w')
    # Identify nonsyntenic relationships
    rels = []
    with open(GENE_KEY, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            if 'curated' in tmp[8]:
                continue
            else:
                rels.append((tmp[3], tmp[7]))
    # For each nonsyntenic relationship, ask if B or P are duplicated. Because
    # this goes on a gene-by-gene basis, we can encounter double-counted entires
    # (if B and P are both duplicated)
    ns_dups = []
    counted = []
    ns_clust = 0
    handle.write('Duplicate_ID\tB_Cluster\tP_Cluster\tB_Genes\tP_Genes\n')
    for ns in rels:
        b = ns[0]
        p = ns[1]
        b_dup = bclust.get(b, None)
        p_dup = pclust.get(p, None)
        # If they are not duplicated, we continue
        if not b_dup and not p_dup:
            continue
        elif b_dup in b_s or p_dup in p_s:
            # Similarly, if they are in the list of syntenic duplicates, we
            # continue
            continue
        elif (b_dup, p_dup) in counted:
            # If we've seen this combination of duplicate IDs before, we
            # continue
            continue
        else:
            counted.append((b_dup, p_dup))
            # Get the genes in each duplicate
            b_genes = []
            p_genes = []
            if b_dup:
                for g in sorted(bclust):
                    if bclust[g] == b_dup:
                        b_genes.append(g)
            if p_dup:
                for g in sorted(pclust):
                    if pclust[g] == p_dup:
                        p_genes.append(g)
            # Append the genes under consideration if they are not a part of
            # the lists
            if b != 'NA' and b not in b_genes and b.startswith('Zm'):
                b_genes.append(b)
            if p != 'NA' and p not in p_genes and p.startswith('Zm'):
                p_genes.append(p)
            # Put them in the list of nonsyntenic dups
            ns_dups.append(b_genes + p_genes)
            # Write the table
            dup_id = 'Nonsyntenic_Dup_' + str(ns_clust).zfill(4)
            if not b_dup:
                b_dup = 'NA'
            if not p_dup:
                p_dup = 'NA'
            if not b_genes:
                b_genes = ['NA']
            if not p_genes:
                p_genes = ['NA']
            handle.write('\t'.join([
                dup_id,
                b_dup,
                p_dup,
                ','.join(b_genes),
                ','.join(p_genes)]) + '\n')
            ns_clust += 1
    handle.flush()
    handle.close()
    return ns_dups


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
    # If sb and os IDs are empty, then we are looking at a nonsyntenic cluster
    if not sb_ids and not os_ids:
        unaligned = tempfile.NamedTemporaryFile(
            prefix='Syntenic_Tandem_',
            suffix='.fasta',
            mode='w+t')
        b_cmd = ['samtools', 'faidx', B_CDS] + b_ids
        p_cmd = ['samtools', 'faidx', P_CDS] + p_ids
        subprocess.call(b_cmd, shell=False, stdout=unaligned)
        subprocess.call(p_cmd, shell=False, stdout=unaligned)
        unaligned.seek(0)
        return unaligned
    else:
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
            # Seek to the beginning so we can read the sequence data out of it
            # later
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


def main(syn_out, nonsyn_out, syn_aln, nonsyn_aln):
    """Main function."""
    gene_name_trans = read_seqs(B_CDS, P_CDS, SB_CDS, OS_CDS)
    b73_clusters = parse_clusters(B_CLUST)
    ph207_clusters = parse_clusters(P_CLUST)
    synteny = parse_synteny(SYNTENY)
    gkey = parse_key(GENE_KEY)
    b_s_dups, b_ns_dups = find_syntenic(B_CLUST, gkey, 'B73')
    p_s_dups, p_ns_dups = find_syntenic(P_CLUST, gkey, 'PH207')
    syn_dups = link_syntenic_dups(b73_clusters, ph207_clusters, synteny, syn_out)
    ns_dups = link_nonsyntenic_dups(b73_clusters, ph207_clusters, b_s_dups, p_s_dups, nonsyn_out)
    # Iterate through the syntenic dups and align them
    for index, cluster in enumerate(syn_dups):
        # build an output name based on the index. The numbers should sync up
        # with the numbers written into the table
        outname = 'Syntenic_Dup_' + str(index).zfill(4) + '.fa'
        un = extract_sequences(cluster, gene_name_trans)
        if not un:
            continue
        bt_aln = align(un)
        write_aln(bt_aln, syn_aln, outname)
    # Do the same for the nonsyntenic
    for index, cluster in enumerate(ns_dups):
        outname = 'Nonsyntenic_Dup_' + str(index).zfill(4) + '.fa'
        un = extract_sequences(cluster, gene_name_trans)
        if not un:
            continue
        bt_aln = align(un)
        write_aln(bt_aln, nonsyn_aln, outname)
    return


if len(sys.argv) != 5:
    print """Generate alignments for syntenic and nonsyntenic duplications in B73 and
PH207. Keep track of cases where there are single-copy genes. This script
differs from the syntenic and nonsyntenic specific ones in that it tries to
do a better job at distinguishing them: if there is a syentenic gene in the
tandem duplicate cluster, it is called syntenic. Else it is nonsyntenic. Takes
four arguments:
    1) Syntenic assignment output file
    2) Nonsyntenic assignment output file
    3) Syntenic alignment output directory
    4) Nonsyntenic alignment output directory."""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
