#!/usr/bin/env python
"""Parse the cassette file and the tandem duplication dates files to report the
dates of cassette duplications. Takes three arguments:
    1) Cassettes file
    2) Syntenic ages
    3) Nonsyntenic ages
"""

import sys


def parse_syn(syntenic):
    """Parse the syntenic table and store as a dictionary that associates a
    tandem duplicate with its age."""
    syn_ages = {}
    with open(syntenic, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split('\t')
                # keep only the tandems
                b_genes = []
                p_genes = []
                if ',' in tmp[2]:
                    b_genes += tmp[2].split(',')
                if ',' in tmp[3]:
                    b_genes += tmp[3].split(',')
                if ',' in tmp[4]:
                    p_genes += tmp[4].split(',')
                if ',' in tmp[5]:
                    p_genes += tmp[5].split(',')
                b_ages = []
                p_ages = []
                if tmp[6] != 'NA':
                    b_ages += tmp[6].split(',')
                if tmp[7] != 'NA':
                    b_ages += tmp[7].split(',')
                if tmp[8] != 'NA':
                    p_ages += tmp[8].split(',')
                if tmp[9] !='NA':
                    p_ages += tmp[9].split(',')
                # If there are B73 tandem duplicates, then we iterate through
                # and build the dictionary
                if b_ages:
                    ct = 1
                    for age in b_ages:
                        b1 = b_genes[ct-1]
                        b2 = b_genes[ct]
                        syn_ages[(b1, b2)] = age
                        ct += 1
                if p_ages:
                    ct = 1
                    for age in p_ages:
                        p1 = p_genes[ct-1]
                        p2 = p_genes[ct]
                        syn_ages[(p1, p2)] = age
                        ct += 1
    return syn_ages


def parse_nonsyn(nonsyn):
    """Parse the nonsyntenic table and store as a dictionary in the same format
    as the syntenic tandem duplicates."""
    nonsyn_ages = {}
    with open(nonsyn, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split('\t')
                b_genes = []
                p_genes = []
                if ',' in tmp[1]:
                    b_genes += tmp[1].split(',')
                if ',' in tmp[2]:
                    p_genes += tmp[2].split(',')
                b_ages = []
                p_ages = []
                if tmp[3] != 'NA':
                    b_ages += tmp[3].split(',')
                if tmp[4] != 'NA':
                    p_ages += tmp[4].split(',')
                if b_ages:
                    ct = 1
                    for age in b_ages:
                        b1 = b_genes[ct-1]
                        b2 = b_genes[ct]
                        nonsyn_ages[(b1, b2)] = age
                        ct += 1
                if p_ages:
                    ct = 1
                    for age in p_ages:
                        p1 = p_genes[ct-1]
                        p2 = p_genes[ct]
                        nonsyn_ages[(p1, p2)] = age
                        ct += 1
    return nonsyn_ages


def main(cass, syn, nonsyn):
    """Main function."""
    # Parse the ages
    s_ages = parse_syn(syn)
    ns_ages = parse_nonsyn(nonsyn)
    # And then iterate through the cassettes
    with open(cass, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                print line.strip() + '\t' + 'Ages'
            else:
                tmp = line.strip().split('\t')
                c_dups = tmp[0].split(';')
                ages = []
                for c in c_dups:
                    # Get the tuples with the genes for lookup
                    genes = sorted(c.split(','))
                    idx = 1
                    for i in genes[1:]:
                        lookup = (genes[idx-1], genes[idx])
                        if lookup in s_ages:
                            age = s_ages[lookup]
                        elif lookup in ns_ages:
                            age = ns_ages[lookup]
                        else:
                            age = 'NA'
                        ages.append(age)
                        idx += 1
                print line.strip() + '\t' + ','.join(ages)


if len(sys.argv) != 4:
    print """Parse the cassette file and the tandem duplication dates files to report the
dates of cassette duplications. Takes three arguments:
    1) Cassettes file
    2) Syntenic ages
    3) Nonsyntenic ages"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2], sys.argv[3])
