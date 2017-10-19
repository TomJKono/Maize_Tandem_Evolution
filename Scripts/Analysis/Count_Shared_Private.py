#!/usr/bin/env python
"""Count the number of duplications that are shared/private between B73 and
PH207. Takes four arguments:
    1) Syntenic assignment table
    2) Nonsyntenic assignment table
    3) B73 clusters
    4) PH207 clusters
"""

import sys


def parse_clust(clust):
    """Return a list of cluster IDs."""
    cids = []
    with open(clust, 'r') as f:
        for line in f:
            cids.append(line.strip().split()[0])
    return cids


def main(syn, nonsyn, b, p):
    """Main function."""
    # Get the "master" list of B73 and PH207 clusters
    b_clusters = parse_clust(b)
    p_clusters = parse_clust(p)
    # Keep track of clusters that appear in duplications
    b_dups = []
    p_dups = []
    # Keep counts of private and shared
    syn_b1_priv = 0
    syn_p1_priv = 0
    syn_b2_priv = 0
    syn_p2_priv = 0
    syn_m1_shared = 0
    syn_m2_shared = 0
    nonsyn_b_priv = 0
    nonsyn_p_priv = 0
    nonsyn_shared = 0
    with open(syn, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split()
                b1 = tmp[2]
                b2 = tmp[3]
                p1 = tmp[4]
                p2 = tmp[5]
                # This is going to be ugly. We want to test the combinations of
                # NA and not-NA
                dup_state = []
                if b1 == 'NA':
                    dup_state.append(False)
                else:
                    dup_state.append(True)
                    b_dups.append(b1)
                if b2 == 'NA':
                    dup_state.append(False)
                else:
                    dup_state.append(True)
                    b_dups.append(b2)
                if p1 == 'NA':
                    dup_state.append(False)
                else:
                    dup_state.append(True)
                    p_dups.append(p1)
                if p2 == 'NA':
                    dup_state.append(False)
                else:
                    dup_state.append(True)
                    p_dups.append(p2)
                # Then, test the combinations.
                if dup_state == [True, True, True, True]:
                    syn_m1_shared += 1
                    syn_m2_shared += 1
                elif dup_state == [True, True, True, False]:
                    syn_m1_shared += 1
                    syn_b2_priv += 1
                elif dup_state == [True, True, False, True]:
                    syn_m2_shared += 1
                    syn_b1_priv += 1
                elif dup_state == [True, False, True, True]:
                    syn_m1_shared += 1
                    syn_p2_priv += 1
                elif dup_state == [False, True, True, True]:
                    syn_m2_shared += 1
                    syn_p1_priv += 1
                elif dup_state == [True, True, False, False]:
                    syn_b1_priv += 1
                    syn_b2_priv += 1
                elif dup_state == [True, False, True, False]:
                    syn_m1_shared += 1
                elif dup_state == [True, False, False, True]:
                    syn_b1_priv += 1
                    syn_p2_priv += 1
                elif dup_state == [False, True, True, False]:
                    syn_p1_priv += 1
                    syn_b2_priv += 1
                elif dup_state == [False, True, False, True]:
                    syn_m2_shared += 1
                elif dup_state == [False, False, True, True]:
                    syn_p1_priv += 1
                    syn_p2_priv += 1
                elif dup_state == [True, False, False, False]:
                    syn_b1_priv += 1
                elif dup_state == [False, True, False, False]:
                    syn_b2_priv += 1
                elif dup_state == [False, False, True, False]:
                    syn_p1_priv += 1
                elif dup_state == [False, False, False, True]:
                    syn_p2_priv += 1
                else:
                    continue
    with open(nonsyn, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split()
                b = tmp[1]
                p = tmp[2]
                if b in b_dups or p in p_dups:
                    continue
                if b == 'NA' and p != 'NA':
                    nonsyn_p_priv += 1
                if p == 'NA' and b != 'NA':
                    nonsyn_b_priv += 1
                if b != 'NA' and p != 'NA':
                    nonsyn_shared += 1
                if b != 'NA':
                    b_dups.append(b)
                if p != 'NA':
                    p_dups.append(p)
    # Then count up the clusters that do not have homologues
    b_no_hom = set(b_clusters) - set(b_dups)
    p_no_hom = set(p_clusters) - set(p_dups)
    print 'M1 Shared:', syn_m1_shared
    print 'M2 Shared:', syn_m2_shared
    print 'B1 Private:', syn_b1_priv
    print 'B2 Private:', syn_b2_priv
    print 'P1 Private:', syn_p1_priv
    print 'P2 Private:', syn_p2_priv
    print 'Nonsyntenic Shared:', nonsyn_shared
    print 'B Nonsyn Private:', nonsyn_b_priv
    print 'P Nonsyn Private:', nonsyn_p_priv
    print 'B No Homologues:', len(b_no_hom)
    print 'P No Homologues:', len(p_no_hom)
    return


if len(sys.argv) != 5:
    print """Count the number of duplications that are shared/private between B73 and
PH207. Takes four arguments:
    1) Syntenic assignment table
    2) Nonsyntenic assignment table
    3) B73 clusters
    4) PH207 clusters"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
