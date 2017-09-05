#!/usr/bin/env python
"""Try to identify tandemly duplicated gene cassettes by identifying clusters
that match certain criteria:
    - Clusters are overlapping, but not nested
    - Clusters have similar "width"

Takes two arguments:
    1) Tandem duplicate clusters file
    2) Sorted GFF
"""

import sys
import pprint


def usage():
    """Print a nice usage message."""
    print """Usage:
Tandem_Cassettes.py Tandem_Clusters.csv genes.gff

will identify tandem duplicate clusters that look like potential gene cassette
duplications. It will do this by identifying overlapping tandem duplicate
clusters that have similar intervals between duplicates, in terms of the number
of intervening genes. These can be filtered by synteny with Sorghum to obtain
a more refined list.

The clusters.csv file has one line per cluster, with member IDs separated by
commas. The GFF must be sorted by chromosome, then base pair start position.
For simplicity, it is assumed to have gene annotations only."""
    exit(1)
    return


def parse_gff(gff):
    """Read the GFF, and return a dictionary of the following form:
    {
        chr: [gene1, gene2, gene3, ...],
        chr: [ ... ],
        ...
    }
    The genes are assumed to be ordered in the GFF, by chromosome start pos."""
    gff_data = {}
    with open(gff, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                tmp = line.strip().split('\t')
                chrom = tmp[0]
                #   Assume that gene ID is the first piece of metadata, and that
                #   the field has this form:
                #   ID=Zm00001d027230;biotype=protein_coding;de...
                geneid = tmp[-1].split(';')[0][3:]
                if chrom not in gff_data:
                    gff_data[chrom] = [geneid]
                else:
                    gff_data[chrom].append(geneid)
    return gff_data


def lookup_cluster(clst, gff_dat):
    """Takes a list of gene IDs and looks up their gene coordinate in the
    genome, and sorts them by this coordinate. Returns a list of tuples:
    [(geneA, chr, idx), (geneB, chr, idx), ... ]
    """
    #   First, isolate the chromosome. We assume that all tandem duplicates are
    #   present on the same chromosome
    for c in gff_dat:
        if clst[0] in gff_dat[c]:
            chromgenes = gff_dat[c]
            break
    indices = [chromgenes.index(g) for g in clst]
    #   We have to make a long list of the chromosome, else zip() won't handle
    #   it properly
    return zip(clst, [c]*len(indices), indices)


def order_clusters(clst_list):
    """Takes a list of clusters and orders them by the occurence of the first
    member of the cluster. The input is a list of list of tuples:
    [
        [(geneA, chr, idx), (geneB, chr, idx), ... ],
        [(geneC, chr, idx), (geneD, chr, idx), ...],
        ...
    ]
    and the output is a list of tuples, sorted by idx of the first element."""
    tmp = []
    #   First order the inner lists - the individual clusters. We want to put
    #   the genes in single clusters in order from first to last, over the
    #   chromosome.
    for c in clst_list:
        #   First sort by chromosome, then by chromosome index
        tmp.append(sorted(c, key=lambda x: (x[1], x[2])))
    #   Then sort the outer list - the whole cluster. We just want to sort by
    #   the first gene in the cluster
    return sorted(tmp, key=lambda x: (x[0][1], x[0][2]))


def find_overlapping(clst_list):
    """Takes a sorted list of lists of tuples and identifies elements that are
    overlapping, but not nested. Since they are sorted by start coordinate of
    the genomic position, this should be pretty simple - just check if the
    previous element encapsulates the current one in a loop. Returns a list
    of overlapping groups of clusters:
    [
        [
            [(geneA, chr, idx), (geneB, chr, idx), ...],
            [(geneC, chr, idx), ...],
            ...
        ], 
        [
            [(geneX, chr, idx), (geneY, chr, idx), ...],
            [(geneW, chr, idx), (geneZ, chr, idx), ...],
            ...
        ]
    ]
    """
    overlap = []
    #   Define a sub-function to test whether one interval is a subset of
    #   another
    def is_subset(i1, i2):
        """Returns True if i1 is completely contained within i2, False
        otherwise."""
        if i1[0] > i2[0] and i1[1] < i2[1]:
            return True
        else:
            return False

    #   Check each cluster for overlapping, but not nesting. We have to respect
    #   the chromosome assignment, too. We use range() so that we start the
    #   count from the second cluster.
    for i in range(1, len(clst_list)):
        j = i - 1
        #   A temporary list to hold overlapping clusters
        o = []
        #   Check if we are on a new chromosome or not. If we are, contiunue
        if clst_list[i][0][1] != clst_list[j][0][1]:
            continue
        else:
            #   Start comparing clusters. Calculate the intervals for the two
            #   adjacent (and potentially overlapping) clusters. i and j index
            #   the clusters, and x indexes the genes within the cluster
            curr_int = [
                (clst_list[i][x-1][2], clst_list[i][x][2])
                for x
                in range(1, len(clst_list[i]))]
            prev_int = [
                (clst_list[j][x-1][2], clst_list[j][x][2])
                for x
                in range(1, len(clst_list[j]))]
            #   First, ask if the intervals overlap. If not, then we continue.
            #   If the start of the current cluster is beyond the end of the
            #   previous cluster, then we do not overlap
            if curr_int[0][0] > prev_int[-1][1]:
                continue
            else:
                #   Now, check if any interval from the current can completely
                #   fit inside an interval from the previous.
                subsets = []
                for current_interval in curr_int:
                    for previous_interval in prev_int:
                        subsets.append(
                            is_subset(current_interval, previous_interval)
                            )
                #   If any of the intervals are subsets, reject them
                if any(subsets):
                    continue
                else:
                    o.append(clst_list[j])
                    o.append(clst_list[i])
        overlap.append(o)
    #   Now, we have a list of lists, with each sub-list being two elements
    #   long. We want to merge these if they have common clusters. Start with
    #   the first cluster in overlaps
    merged = [overlap[0]]
    #   And iterate over the rest
    for i in overlap[1:]:
        #   If the rightmost cluster of the previous chunk is the same as the
        #   leftmost cluster of the current chunk:
        if merged[-1][1] == i[0]:
            #   Append the rightmost cluster from the current chunk to the
            #   merged cluster list
            merged[-1].append(i[1])
        else:
            #   Else, keep the overlapping clusters as discrete chunks
            merged.append(i)
    return merged


def find_similar_size(clst_list, thresh=2):
    """Takes a list of overlapping tandem duplicate clusters and filters out
    those that have large differences in cluster size. For now, we will allow
    +/- 2 genes in size."""
    passing = []
    #   For each group of overlapping clusters:
    for oc in clst_list:
        #   Generate the interval sizes
        sizes = [len(c) for c in oc]
        #    Then, calculate all pairwise differences in sizes
        sizediff = []
        for i, x in enumerate(sizes):
            for y in sizes[i+1:]:
                sizediff.append(abs(y - x))
        passthresh = [True if a <= thresh else False for a in sizediff]
        #   If any interval passes the wobble filter, keep it
        if any(passthresh):
            passing.append(oc)
    return passing

    #   Generate the firt cluster
def report_clusters(clst_list):
    """Prints the report. Since it is a somewhat complicated data structure, we
    are opting to do all the printing and formatting in a seprate function, to
    keep main() from getting difficult to follow."""
    print 'Clusters\tGeneNumber\tChromosome\tSortedClusterIDs'
    for c in clst_list:
        clusters = []
        coords = []
        for clust in c:
            clusters.append(','.join([s[0] for s in clust]))
            coords.append(','.join([str(s[2]) for s in clust]))
        sorted_clusters = sorted(
            [
                g
                for ct
                in clusters
                for g
                in ct.split(',')
            ])
        toprint = '\t'.join([
            ';'.join(clusters),
            ';'.join(coords),
            c[0][0][1],
            ','.join(sorted_clusters)
            ])
        print toprint
    return


def main(tandems, gff):
    """Main function."""
    g_data = parse_gff(gff)
    clusters = []
    with open(tandems, 'r') as f:
        for line in f:
            clust = line.strip().split('\t')[1].split(',')
            clusters.append(lookup_cluster(clust, g_data))
    #   Then, sort all of them
    ord_clusters = order_clusters(clusters)
    overlapping = find_overlapping(ord_clusters)
    sim_size = find_similar_size(overlapping)
    report_clusters(sim_size)


if len(sys.argv) != 3:
    usage()
else:
    main(sys.argv[1], sys.argv[2])
