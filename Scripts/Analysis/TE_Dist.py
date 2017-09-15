#!/usr/bin/env python
"""Calculate the upstream and downstream distance to the nearest DNA TE for a
specified list of tandem duplicate clusters. Takes three arguments:
    1) Tandem duplicate clusters
    2) TE GFF
    3) Gene GFF
"""

import sys
import pprint


def parse_gff(gff, te=False):
    """Read the gff and return a dictionary of the following form:
    {
        chr: {
            geneID: (start, end),
            geneID: (start, end),
            ...
        chr: ...,
    }
    """
    gff_data = {}
    with open(gff, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                tmp = line.strip().split('\t')
                chrom = tmp[0]
                startpos = int(tmp[3])
                endpos = int(tmp[4])
                if te:
                    geneid = tmp[-1].split(';')[1][6:]
                else:
                    geneid = tmp[-1].split(';')[0][3:]
                if chrom not in gff_data:
                    gff_data[chrom] = {geneid: (startpos, endpos)}
                else:
                    gff_data[chrom][geneid] = (startpos, endpos)
    return gff_data


def calculate_distance(geneid, genes, tes):
    """Calculate the upstream distance and downstream distance to the nearest
    DNA TE from the specified gene."""
    #   Get which chromosome
    for c in genes:
        if geneid in genes[c]:
            chromosome = c
            break
    #   Get the gene position
    genestart, geneend = genes[chromosome][geneid]
    #   if the gene chromosome does not have any TEs, return NA
    if chromosome not in tes:
        return ('NA', 'NA', 'NA', 'NA')
    #   Then get the TE that is closest downstream. We do this by iterating
    #   forwards, then breaking when we find a TE that starts after the end of
    #   the gene.
    for downstream in sorted(list(tes[chromosome].iteritems()), key=lambda x: x[1][0]):
        if downstream[1][0] > geneend:
            break
    #   And get the TE that is closest upstream. We do this by doing the same
    #   strategy as above, but in reverse
    for upstream in reversed(sorted(list(tes[chromosome].iteritems()), key=lambda x: x[1][0])):
        if upstream[1][1] < genestart:
            break
    #   Then calculate the distances
    #   Gene start - TE end and TE start - gene end
    dist_upstream = genestart - upstream[1][1]
    dist_downstream = downstream[1][0] - geneend
    #   then return it all
    return (upstream[0], downstream[0], dist_upstream, dist_downstream)


def main(tandems, tegff, genegff):
    """Main function."""
    genes = parse_gff(genegff)
    tes = parse_gff(tegff, te=True)
    #   Print a header
    print 'Gene\tUpstreamTE\tDownstreamTE\tSameFam\tUpstreamBP\tDownstreamBP'
    with open(tandems, 'r') as f:
        for line in f:
            tandem_dups = line.strip().split()[1].split(',')
            for t in tandem_dups:
                u_te, d_te, dist_u, dist_d = calculate_distance(t, genes, tes)
                if dist_u < 0:
                    dist_u = 'NA'
                    u_te = 'NA'
                if dist_d < 0:
                    dist_d = 'NA'
                    d_te = 'NA'
                #   Figure out if the upstream and downstream TEs are the same
                #   family or not. Split on 'v' and take the first part. Then
                #   remove the last three characters ("B73") to get the family
                d_te_fam = d_te.split('v')[0][:-3]
                u_te_fam = u_te.split('v')[0][:-3]
                if d_te_fam == u_te_fam and (d_te != 'NA' or u_te != 'NA'):
                    samefam = "1"
                elif d_te == 'NA' or u_te == 'NA':
                    samefam = 'NA'
                else:
                    samefam = "0"
                print '\t'.join([
                    t,
                    u_te,
                    d_te,
                    samefam,
                    str(dist_u),
                    str(dist_d)
                    ])
    return


main(sys.argv[1], sys.argv[2], sys.argv[3])
