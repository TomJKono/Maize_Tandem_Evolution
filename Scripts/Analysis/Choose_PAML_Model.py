#!/usr/bin/env python
"""Parse the compiled PAML table and choose which model fits each orthogroup
the best. First, compare each model to the null, then compare the nested
alternative models. For each orthogroup, report the best fitting model, and the
omega values and proportions at each omega. Takes one argument:
    1) PAML table.
"""

import sys


def choose_model(paml_res_line, alpha):
    """Compare the p-values of the model comparisons to determine which model
    fits the data best. Do a two-step comparison:
        Compare Ha1 v Null
        Compare Ha2 v Null
        Compare Ha3 v Null
    If Ha1 or Ha2 is signifcant and Ha3 is significant, then compare Ha3 to
    the simpler models.
    If Ha3 v Ha1,Ha2 is significant, report Ha3
    else report the simpler model"""
    try:
        ha1_v_null = float(paml_res_line[5])
        ha2_v_null = float(paml_res_line[6])
        ha3_v_null = float(paml_res_line[7])
        ha3_v_ha1 = float(paml_res_line[8])
        ha3_v_ha2 = float(paml_res_line[9])
    except ValueError:
        # If any of the p-values are NA, then we cannot choose a best fitting
        # model. Return NA.
        return 'NA'
    # Ask if the model comparisons are significant or not
    if ha1_v_null < alpha:
        ha1n = True
    else:
        ha1n = False
    if ha2_v_null < alpha:
        ha2n = True
    else:
        ha2n = False
    if ha3_v_null < alpha:
        ha3n = True
    else:
        ha3n = False
    if ha3_v_ha1 < alpha:
        ha3ha1 = True
    else:
        ha3ha1 = False
    if ha3_v_ha2 < alpha:
        ha3ha2 = True
    else:
        ha3ha2 = False
    # Do the comparisons
    if not ha1n and not ha2n and not ha3n:
        return 'Null'
    else:
        if ha1n:
            if ha3n:
                if ha3ha1:
                    return 'Ha3'
                else:
                    return 'Ha1'
            else:
                return 'Ha1'
        elif ha2n:
            if ha3n:
                if ha3ha2:
                    return 'Ha3'
                else:
                    return 'Ha2'
            else:
                return 'Ha2'
        elif ha3n:
            return 'Ha3'
        else:
            return 'Null'


def get_omegas(paml_res_line, model):
    """Print the omegas for the various branch types, depending on the model
    that is chosen to be best-fitting. If the null model fits the best, then
    report the same omega for the grass/maize/tandem paritions. If Ha1 fits
    the best, then report the same omega for the maize/tandem partitions. If
    Ha2 fits the best, then report the same omega for the grass/maize
    partitions."""
    # Define a dictionary so we can decide which columns to get for the
    # appropriate models
    # The nine columes are for the following, in order:
    #   Grass Omega 0
    #   Grass Omega 1
    #   Grass Omega 2
    #   Maize Omega 0
    #   Maize Omega 1
    #   Maize Omega 2
    #   Tandem Omega 0
    #   Tandem Omega 1
    #   Tandem Omega 2
    # For the null model, all species omegas are the same
    # For Ha1, maize=tandem
    # For Ha2, maize=grass
    # For Ha3, all distinct
    omega_cols = {
        'Null': [[11, 13, 15], [11, 13, 15], [11, 13, 15]],
        'Ha1': [[18, 21, 24], [19, 22, 25], [19, 22, 25]],
        'Ha2': [[28, 31, 34], [28, 31, 34], [29, 32, 35]],
        'Ha3': [[38, 42, 46], [39, 43, 47], [40, 44, 48]]
        }
    # Same for omega proportions
    prop_cols = {
        'Null': [12, 14, 16],
        'Ha1': [20, 23, 26],
        'Ha2': [30, 33, 36],
        'Ha3': [41, 45, 49]
        }
    if model == 'NA':
        om = ['NA,NA,NA']*3
        pr = ['NA'] * 3
    else:
        om = [
            ','.join([
                paml_res_line[o]
                for o
                in x])
            for x
            in omega_cols[model]]
        pr = [paml_res_line[p] for p in prop_cols[model]]
    return (om, pr)


def main(paml_res):
    """Main function."""
    # Print a header
    print '\t'.join([
            'Orthogroup',
            'BestModel',
            'GrassOmegas',
            'MaizeOmegas',
            'TandemOmegas',
            'PropOmega0',
            'PropOmega1',
            'PropOmega2'])
    with open(paml_res, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split()
                best_model = choose_model(tmp, 0.01)
                omegas, props = get_omegas(tmp, best_model)
                toprint = [tmp[0], best_model] + omegas + props
                print '\t'.join(toprint)
    return


if len(sys.argv) != 2:
    print """Parse the compiled PAML table and choose which model fits each orthogroup
the best. First, compare each model to the null, then compare the nested
alternative models. For each orthogroup, report the best fitting model, and the
omega values and proportions at each omega. Takes one argument:
    1) PAML table."""
    exit(1)
else:
    main(sys.argv[1])
