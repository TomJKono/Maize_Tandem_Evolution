#!/usr/bin/env python
"""Count the number of duplications that are shared/private to B73 and PH207.
Parses the output from the Homologous_Cluster_Sizes.py script. Takes one
argument:
    1) Homologous size table
"""

import sys


def main(homologous):
    """Main function."""
    # Keep counts of private and shared
    b_priv = 0
    p_priv = 0
    shared = 0
    with open(homologous, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split()
                b = int(tmp[0])
                p = int(tmp[1])
                if b>1 and p>1:
                    shared += 1
                elif b<=1 and p>1:
                    p_priv += 1
                elif p<=1 and b>1:
                    b_priv += 1
                else:
                    continue
    print b_priv, p_priv, shared
    return


if len(sys.argv) != 2:
    print """Count the number of duplications that are shared/private to B73 and PH207.
Parses the output from the Homologous_Cluster_Sizes.py script. Takes one
argument:
    1) Homologous size table"""
    exit(1)
else:
    main(sys.argv[1])
