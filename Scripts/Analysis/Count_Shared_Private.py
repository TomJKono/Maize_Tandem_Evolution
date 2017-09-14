#!/usr/bin/env python
"""Count the number of duplications that are shared/private between B73 and
PH207. Takes two arguments:
    1) Syntenic assignment table
    2) Nonsyntenic assignment table
"""

import sys


def main(syn, nonsyn):
    """Main function."""
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
                if b2 == 'NA':
                    dup_state.append(False)
                else:
                    dup_state.append(True)
                if p1 == 'NA':
                    dup_state.append(False)
                else:
                    dup_state.append(True)
                if p2 == 'NA':
                    dup_state.append(False)
                else:
                    dup_state.append(True)
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
                if b == 'NA' and p != 'NA':
                    nonsyn_p_priv += 1
                if p == 'NA' and b != 'NA':
                    nonsyn_b_priv += 1
                if b != 'NA' and p != 'NA':
                    nonsyn_shared += 1
    print syn_m1_shared, syn_m2_shared, syn_b1_priv, syn_b2_priv, syn_p1_priv, syn_p2_priv
    print nonsyn_b_priv, nonsyn_p_priv, nonsyn_shared
    return


if len(sys.argv) != 3:
    print """Count the number of duplications that are shared/private between B73 and
PH207. Takes two arguments:
    1) Syntenic assignment table
    2) Nonsyntenic assignment table"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2])
