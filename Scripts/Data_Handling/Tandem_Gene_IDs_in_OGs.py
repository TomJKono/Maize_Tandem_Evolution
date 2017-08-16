#!/usr/bin/env python
"""Take the output from Tandem_Orthogroups.py and generate a table that gives
a two-column table. The first column is the orthogroup ID, and the second
column is the maize tandem gene IDs that are present in it. Takes three
arguments:
    1) Tandem_Orthogroups output table
    2) B73 tandem clusters file
    3) PH207 tandem clusters file.
"""

import sys
import pprint


def parse_tandems(tand):
    """Parse the tandem duplicate clusters file and return it as a dictionary.
    The key is the cluster ID and the value is a list of gene IDs."""
    t = {}
    with open(tand, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            t[tmp[0]] = tmp[1].split(',')
    return t


def parse_table(tab):
    """Parse the output table from Tandem_Orthogroups.py and store it as a
    dictionary. It will be a somewhat messy dictionary:
    {
        OG_ID: {'B73': [cluster, cluster, ...],
                'PH207': [cluster, cluster, ...]
                },
        ...
    }"""
    og_dict = {}
    with open(tab, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split('\t')
                genotype = tmp[1]
                c_id = tmp[3]
                og_id = tmp[5]
                if c_id == 'NA' or og_id == 'NA':
                    continue
                if og_id in og_dict:
                    if genotype in og_dict[og_id]:
                        if c_id in og_dict[og_id][genotype]:
                            continue
                        else:
                            og_dict[og_id][genotype].append(c_id)
                    else:
                        og_dict[og_id][genotype] = [c_id]
                else:
                    og_dict[og_id] = {}
                    og_dict[og_id][genotype] = [c_id]
    return og_dict


def make_gene_list(table, tandems):
    """Build the table of tandem gene IDs from the parsed table and the dicts
    of tandem duplicate IDs."""
    og_t = {}
    for og in sorted(table):
        # Unpack the B73 and PH207 clusters
        b_clusters = table[og].get('B73', ['NA'])
        p_clusters = table[og].get('PH207', ['NA'])
        tand_genes = []
        # Tack the B73 and PH207 tandem genes onto the list
        if b_clusters != ['NA']:
            for cluster in b_clusters:
                tand_genes += tandems['B73'][cluster]
        if p_clusters != ['NA']:
            for cluster in p_clusters:
                tand_genes += tandems['PH207'][cluster]
        og_t[og] = ','.join(sorted(list(set(tand_genes))))
    return og_t


def main(tandem_ogs, b_tand, p_tand):
    """Main function."""
    b_clust = parse_tandems(b_tand)
    p_clust = parse_tandems(p_tand)
    og_table = parse_table(tandem_ogs)
    # Put the B73 and PH207 clusters together into a big dictionary
    t_dict = {'B73': b_clust, 'PH207': p_clust}
    og_tandem = make_gene_list(og_table, t_dict)
    # Print it out
    for o in sorted(og_tandem):
        print o + '\t' + og_tandem[o]
    return


if len(sys.argv) != 4:
    print """Take the output from Tandem_Orthogroups.py and generate a table that gives
a two-column table. The first column is the orthogroup ID, and the second
column is the maize tandem gene IDs that are present in it. Takes three
arguments:
    1) Tandem_Orthogroups output table
    2) B73 tandem clusters file
    3) PH207 tandem clusters file."""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2], sys.argv[3])
