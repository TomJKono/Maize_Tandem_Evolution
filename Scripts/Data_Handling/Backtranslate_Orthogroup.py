#!/usr/bin/env python
"""Back-translate the protein sequences used to generate the orthologous groups
to nucleotides. This is useful for calculating dN/dS for individual genes in the
orthologue alignment. Requires Biopython and samtools. Takes two arguments:
    1) Directory of FASTA files with nucleotide CDS sequences
    2) Amino acid alignment

There are a couple weird things we have to keep track of for making sure we can
get the sequence from the database. These species have '.p' at the end of their
amino acid sequence names that need to be removed for fetching the CDS:
    - Panicum virgatum
    - Brachypodium distachyon
    - Sorghum bicolor
    - Setaria viridis
    - Setaria_italica
Triticum urartu needs to have '-P1' replaced with '-T1'
B73 and PH207 names end in '_1' and the CDS does not.
"""

import sys
import subprocess
import os

try:
    from Bio import SeqIO
except ImportError:
    print 'This script requires Biopython.'
    exit(1)


def list_files(directory):
    """Get the FASTA files that are present in the supplied directory."""
    try:
        abpath = os.path.abspath(directory)
        filepaths = os.listdir(abpath)
        filepaths = [
            os.path.join(abpath, f)
            for f
            in filepaths
            if f.endswith('fa')
            ]
    except OSError, IOError:
        print 'The directory supplied is not readable, or does not exist!'
        exit(1)
    return filepaths


def fix_seqname(sname):
    """Perform the 'fixes' for sequence fetching from the indexed CDS files.
    This just makes sure the name matches the ones listed in the FASTA."""
    #   The species name is the first part of the pipe
    species, protid = sname.split('|')
    #   Apply our corrections
    remove_dotp = [
        'Panicum_virgatum',
        'Brachypodium_distachyon',
        'Sorghum_bicolor',
        'Setaria_viridis',
        'Setaria_italica'
        ]
    if species in remove_dotp:
        protid = protid.rstrip('.p')
    elif species == 'Triticum_urartu':
        protid = protid.replace('-P1', '-T1')
    elif species == 'B73' or species == 'PH207':
        protid = protid[:-2]
    return (species, protid)


def extract_cds(file_list, species, prot_id):
    """Given the protein ID from the alignment and the file list, get the
    correct FASTA to extract from, then use samtools to fetch the CDS sequence.
    """
    #   Get the filename that we want to query
    for f in file_list:
        if species in f:
            break
    #   Then, buld the command line
    cmd = ['samtools', 'faidx', f, prot_id]
    proc = subprocess.Popen(
        cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
    out, err = proc.communicate()
    #   Then, process the output. We will remove the first line, since it is
    #   just the sequence ID. We will then put all the nucleotides into a single
    #   long string.
    lines = out.split('\n')
    cds = ''.join(lines[1:])
    return cds


def backtranslate(p_seq, n_seq):
    """Iterate through the aligned protein sequence, and replace the amino acids
    with codon triplets from the CDS file."""
    #   Keep track of the new sequence. Also keep track of which codon we are
    #   actually processing (gaps don't count)
    newseq = ''
    codon = 0
    for aa in p_seq:
        if aa == '-':
            newseq += '---'
        else:
            newseq += n_seq[codon*3:(codon*3) + 3]
            codon += 1
    return newseq


def main(db_dir, msa):
    """Main function."""
    #   Get the paths of the CDS files
    paths = list_files(db_dir)
    #   Parse the alignment
    aln = SeqIO.parse(msa, 'fasta')
    for sequence in aln:
        s, p = fix_seqname(sequence.id)
        cdsseq = extract_cds(paths, s, p)
        bt_seq = backtranslate(sequence.seq, cdsseq)
        #   Print out the sequence in FASTA format
        print '>' + sequence.id + '\n' + bt_seq
    return


main(sys.argv[1], sys.argv[2])
