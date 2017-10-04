#!/usr/bin/env python
"""Repport whether the maize tandems are evolving under constraint or relaxed
constraint with respect to maize non-tandem duplicate genes. The reporting is
based on the compiled output from codeml tests. Takes one argument:
    1) Model summary from PAML
"""

import sys


def main(paml_model_summary):
    """Main function."""
    with open(paml_model_summary, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                # Print a header
                print 'Orthogroup\tTandem_Evolution'
            else:
                tmp = line.strip().split('\t')
                og_id = tmp[0]
                best_model = tmp[1]
                mz_omega = tmp[3]
                td_omega = tmp[4]
                # if the best model is Null or Ha1, then there is no evidence
                # that the tandems are evolving at different rates
                if best_model == 'Null' or best_model == 'Ha1':
                    rep = [og_id, 'No Difference']
                elif best_model == 'Ha2' or best_model == 'Ha3':
                    # If we get here, then there is evidence that tandems are
                    # evolving at a different rate from non-tandem duplicates.
                    # We have to compare the omega values to figure out the
                    # selective regime
                    m2 = float(mz_omega.split(',')[2])
                    t2 = float(td_omega.split(',')[2])
                    if m2 > 10 or t2 > 10:
                        # If we have a really high omega for tandems or maize,
                        # then there are probably no synonymous substitutions
                        # among the genes, and we can't say anything meaningful
                        rep = [og_id, 'No dS']
                    elif m2 > t2:
                        # If maize has higher omega than tandems, then there is
                        # more constraint on tandems
                        rep = [og_id, 'More constraint']
                    elif t2 > m2:
                        # If tandems have higher omega, then there is divergent
                        # evolution of the genes
                        rep = [og_id, 'Divergent']
                print '\t'.join(rep)
    return


if len(sys.argv) != 2:
    print """Repport whether the maize tandems are evolving under constraint or relaxed
constraint with respect to maize non-tandem duplicate genes. The reporting is
based on the compiled output from codeml tests. Takes one argument:
    1) Model summary from PAML"""
    exit(1)
else:
    main(sys.argv[1])
