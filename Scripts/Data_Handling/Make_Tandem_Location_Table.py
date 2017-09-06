#!/usr/bin/env python
"""Make a summary of the tandem duplicates genomic locations. For 1-Mb segments
of the genome, report the number of tandem duplicate genes, the number of
RNA TEs, the number of DNA TEs, and which subgenome the majority of the window
is assigned to. This script only works for B73v4. Takes four arguments:
    1) Tandem duplicates table
    2) Sorted genes GFF
    3) TE GFF
    4) Syntenic blocks
"""

import sys

# Define the chromosome sizes as a constant
CHR_SIZES = {
    '1': 307041717,
    '2': 244442276,
    '3': 235667834,
    '4': 246994605,
    '5': 223902240,
    '6': 174033170,
    '7': 182381542,
    '8': 181122637,
    '9': 159769782,
    '10': 150982314}
# Define RNA and DNA TE types as constants, too
RNA_TES = ['LINE_element', 'LTR_retrotransposon', 'solo_LTR', 'SINE_element']
DNA_TES = ['helitron', 'terminal_inverted_repeat_element']


def parse_tandem(tand):
    """Return a list of tandem duplicates."""
    tandems = []
    with open(tand, 'r') as f:
        for line in f:
            t = line.strip().split()[1]
            genes = t.split(',')
            tandems += genes
    return tandems


def parse_genes(gff, tand):
    """Parse the GFF of sorted genes and return two dictionaries: one for all
    genes and one for tandems. Be sure they are sorted."""
    genomewide = {}
    tandems = {}
    with open(gff, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            chrom = tmp[0]
            start = int(tmp[3])
            gene_id = tmp[8].split(';')[0].split('=')[1]
            if chrom in genomewide:
                genomewide[chrom].append(start)
            else:
                genomewide[chrom] = [start]
            if gene_id in tand:
                if chrom in tandems:
                    tandems[chrom].append(start)
                else:
                    tandems[chrom] = [start]
    return (genomewide, tandems)


def parse_tes(te):
    """Return dictionaries of TEs, sorted on start positions. Treat RNA TEs and
    DNA TEs separately. The TE GFF should be sorted on start position already,
    so we assume that it's in the proper order."""
    rna_te = {}
    dna_te = {}
    with open(te, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                tmp = line.strip().split()
                chrom = tmp[0]
                start = int(tmp[3])
                te_type = tmp[2]
                if te_type in RNA_TES:
                    if chrom in rna_te:
                        rna_te[chrom].append(start)
                    else:
                        rna_te[chrom] = [start]
                elif te_type in DNA_TES:
                    if chrom in dna_te:
                        dna_te[chrom].append(start)
                    else:
                        dna_te[chrom] = [start]
                else:
                    continue
    return (rna_te, dna_te)


def parse_blocks(blk):
    """Return a dictionary of subgenome assignments."""
    blocks = {}
    with open(blk, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            chrom = tmp[0]
            start = int(tmp[1])
            end = int(tmp[2])
            subg = tmp[4]
            if chrom in blocks:
                blocks[chrom].append((start, end, subg))
            else:
                blocks[chrom] = [(start, end, subg)]
    return blocks


def build_intervals(step):
    """Return a list of intervals that cover the genome in the specified size.
    The intervals will not overlap."""
    # This is super ugly, but Python doesn't have a built-in natural sort func
    # and I don't want to write one right now.
    chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10']
    intervals = []
    for c in chroms:
        i_start = 0
        end = CHR_SIZES[c]
        while i_start < end:
            i_end = i_start + step
            if i_end > end:
                i_end = end
            intervals.append((c, i_start, i_end))
            i_start += step
    return intervals


def count_features(interval, feats):
    """Count the number of features in the supplied interval."""
    chrom = interval[0]
    i_start = interval[1]
    i_end = interval[2]
    chr_feats = feats[chrom]
    int_feats = [
        f
        for f
        in chr_feats
        if f >= i_start
        and f <= i_end]
    return len(int_feats)


def get_block_assign(interval, blocks):
    """Return an identifier for the major subgenome that the supplied interval
    occupies."""
    chrom = interval[0]
    i_start = interval[1]
    i_end = interval[2]
    chr_blocks = blocks[chrom]
    # Do the easy test: if the interval is fully within a block, then return
    # the subgenome it's in
    in_block = [
        b
        for b
        in chr_blocks
        if b[0] <= i_start
        and b[1] >= i_end]
    if in_block:
        return in_block[0][2]
    else:
        # Next, do the harder tests: if the interval runs off the right hand
        # side of the block
        over_right = [
            b
            for b
            in chr_blocks
            if b[0] <= i_start
            and b[1] <= i_end]
        if over_right:
            overlap = (over_right[0][1] - i_start) / float(i_end - i_start)
            if overlap > 0.5:
                return over_right[0][2]
            else:
                return 'Nonsyntenic'
        else:
            # If the interval runs of the left side of the block
            over_left = [
                b
                for b
                in chr_blocks
                if b[0] >= i_start
                and b[1] >= i_end]
            if over_left:
                overlap = (over_left[0][1] - i_end) / float(i_end - i_start)
                if overlap > 0.5:
                    return over_left[0][2]
                else:
                    return 'Nonsyntenic'
            else:
                return'Nonsyntenic'

def main(tandem, genes, tes, syn):
    """Main function."""
    tandem_gene_ids = parse_tandem(tandem)
    gw, tandem_pos = parse_genes(genes, tandem_gene_ids)
    rna, dna = parse_tes(tes)
    blocks = parse_blocks(syn)
    intervals = build_intervals(1000000)
    # Print a header
    print 'Interval\tN_Genes\tN_Tandem\tN_RNATE\tN_DNATE\tSubgenome'
    # For each interval:
    for i in intervals:
        # We want to get the number of features in each interval
        n_genes = count_features(i, gw)
        n_tand = count_features(i, tandem_pos)
        n_rna = count_features(i, rna)
        n_dna = count_features(i, dna)
        b = get_block_assign(i, blocks)
        # make it print friendly
        c_i = i[0] + ':' + str(i[1]) + '-' + str(i[2])
        print '\t'.join([
            c_i,
            str(n_genes),
            str(n_tand),
            str(n_rna),
            str(n_dna),
            b])
    return


if len(sys.argv) != 5:
    print """Make a summary of the tandem duplicates genomic locations. For 1-Mb segments
of the genome, report the number of tandem duplicate genes, the number of
RNA TEs, the number of DNA TEs, and which subgenome the majority of the window
is assigned to. This script only works for B73v4. Takes four arguments:
    1) Tandem duplicates table
    2) Sorted genes GFF
    3) TE GFF
    4) Syntenic blocks"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
