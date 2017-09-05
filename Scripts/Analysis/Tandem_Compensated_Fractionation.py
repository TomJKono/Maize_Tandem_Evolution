#!/usr/bin/env python
"""Calculate the proportion of "compensated" fractionation states. I.e., if B73
maize1 is single copy, and PH207 maize1 is duplicated, do the maize2 homeologs
show retention in B2 and loss in P2? Takes one argument:
    1) Syntenic cluster assignment table
"""

import sys

# Define lists of vectors that key out various retention/duplication scenarios
subg_comp = [
    (2, 0, 1, 1),
    (0, 2, 1, 1),
    (1, 1, 0, 2),
    (1, 1, 2, 0)
    ]
shared_dup_frac = [
    (2, 0, 2, 0),
    (0, 2, 0, 2)
    ]
shared_dup_ret = [
    (2, 1, 2, 1),
    (1, 2, 1, 2)
    ]
no_comp = [
    (2, 0, 1, 0),
    (0, 2, 0, 1),
    (1, 0, 2, 0),
    (0, 1, 0, 2),
    (2, 1, 1, 1),
    (1, 2, 1, 1),
    (1, 1, 2, 1),
    (1, 1, 1, 2)
    ]
n_comp = 0
n_shared_frac = 0
n_shared_ret = 0
n_no_comp = 0
n_tot = 0
with open(sys.argv[1], 'r') as f:
    for index, line in enumerate(f):
        if index == 0:
            continue
        else:
            tmp = line.strip().split()
            # Set some flags:
            #   0: gene lost
            #   1: single copy
            #   2: duplicated
            if tmp[2] == 'NA' and tmp[6] == 'NA':
                b1 = 0
            elif tmp[2] == 'NA' and tmp[6].startswith('Z'):
                b1 = 1
            elif tmp[2] != 'NA':
                b1 = 2
            if tmp[3] == 'NA' and tmp[7] == 'NA':
                b2 = 0
            elif tmp[3] == 'NA' and tmp[7].startswith('Z'):
                b2 = 1
            elif tmp[3] != 'NA':
                b2 = 2
            if tmp[4] == 'NA' and tmp[8] == 'NA':
                p1 = 0
            elif tmp[4] == 'NA' and tmp[8].startswith('Z'):
                p1 = 1
            elif tmp[4] != 'NA':
                p1 = 2
            if tmp[5] == 'NA' and tmp[9] == 'NA':
                p2 = 0
            elif tmp[5] == 'NA' and tmp[9].startswith('Z'):
                p2 = 1
            elif tmp[5] != 'NA':
                p2 = 2
            # Then, count up the cases of compensation:
            #   1 1 2 0
            #   1 1 0 2
            #   2 0 1 1
            #   0 2 1 1
            # and not compensation: Everything else
            if (b1, b2, p1, p2) in subg_comp:
                n_comp += 1
            if (b1, b2, p1, p2) in shared_dup_frac:
                n_shared_frac += 1
            if (b1, b2, p1, p2) in shared_dup_ret:
                n_shared_ret += 1
            if (b1, b2, p1, p2) in no_comp:
                n_no_comp += 1
            n_tot += 1

print 'Private dup, compensating:', n_comp
print 'Private dup, non-compensating:', n_no_comp
print 'Shared dup, fractionated:', n_shared_frac
print 'Shared dup, retained:', n_shared_ret
print 'Total:', n_tot
