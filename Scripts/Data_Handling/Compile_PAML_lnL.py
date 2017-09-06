#!/usr/bin/env python
"""Compile the results from the PAML tests. This script will gather the log-
likelihood values from the 'rst' files in the PAML results directories, and
compare them to the alternate log-likelihood values with a chi-squared test.
Takes one argument:
    1) PAML results directory
"""

import sys
import os
try:
    from scipy.stats import chi2
except ImportError:
    print 'This script requires Scipy.'
    exit(1)


def find_rst(pamldir):
    """Return absolute paths to the rst files."""
    fullpath = os.path.abspath(os.path.expanduser(pamldir))
    # Then, use the os.walk() function to recursively list all contents in this
    # directory
    rsts = []
    for path, subdirs, files in os.walk(fullpath):
        if 'rst' in files:
            rsts.append(os.path.join(path, 'rst'))
    return rsts


def extract_meta(rst_path):
    """Simple function to extract the OG and hypothesis from the path, because
    they are standarized."""
    parts = rst_path.split('/')
    og = parts[-3]
    hyp = parts[-2]
    return (og, hyp)


def extract_lnl(rst):
    """Parse the rst file and extract the log-likelihood value."""
    # Define a lnl value up here so that we don't get an undefined value error
    lnl = None
    with open(rst, 'r') as f:
        for line in f:
            if 'lnL' in line:
                lnl = float(line.strip().split()[2])
                break
    return lnl


def print_table(res):
    """Calculate the p-values for the model comparisons and the log-
    likelihood values for the models."""
    print 'Orthogroup\tNull_lnL\tHa1_lnL\tHa2_lnL\tHa3_lnL\tHa1_P\tHa2_P\tHa3_P'
    for o in sorted(res):
        ha1 = res[o].get('Ha1', None)
        ha2 = res[o].get('Ha2', None)
        ha3 = res[o].get('Ha3', None)
        null = res[o].get('Null', None)
        # Do the comparisons
        if not null or not ha1:
            ha1_v_null = 'NA'
        else:
            ha1_v_null = 1 - chi2.cdf(2 * (ha1 - null), 1)
        if not null or not ha2:
            ha2_v_null = 'NA'
        else:
            ha2_v_null = 1 - chi2.cdf(2 * (ha2 - null), 1)
        if not null or not ha3:
            ha3_v_null = 'NA'
        else:
            ha3_v_null = 1 - chi2.cdf(2 * (ha3 - null), 2)
        # Print the lines, but turn NoneType into NA
        if not ha1:
            ha1 = 'NA'
        if not ha2:
            ha2 = 'NA'
        if not ha3:
            ha3 = 'NA'
        if not null:
            null = 'NA'
        print '\t'.join([
            o,
            str(null),
            str(ha1),
            str(ha2),
            str(ha3),
            str(ha1_v_null),
            str(ha2_v_null),
            str(ha3_v_null)])
    return


def main(pamldir):
    """Main function."""
    outfiles = find_rst(pamldir)
    results = {}
    for r in outfiles:
        og, h = extract_meta(r)
        lnl = extract_lnl(r)
        if og in results:
            results[og].update({h: lnl})
        else:
            results[og] = {h: lnl}
    print_table(results)
    return


if len(sys.argv) != 2:
    print """Compile the results from the PAML tests. This script will gather the log-
likelihood values from the 'rst' files in the PAML results directories, and
compare them to the alternate log-likelihood values with a chi-squared test.
Takes one argument:
    1) PAML results directory"""
    exit(1)
else:
    main(sys.argv[1])
