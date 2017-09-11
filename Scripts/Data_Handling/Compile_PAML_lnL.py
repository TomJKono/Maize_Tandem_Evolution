#!/usr/bin/env python
"""Compile the results from the PAML tests. This script will pull out some
interesting information from the '_codeml' files are produced as primary output
from the codeml program. Models will be compared with likelihood ratio tests,
and estimated omega values will be printed, if they are available. This script
requires Biopython 1.58 or greater and scipy. Takes one argument:
    1) PAML results directory
"""

import sys
import os
import pprint
try:
    from Bio.Phylo.PAML import codeml
    from scipy.stats import chi2
except ImportError:
    print 'This script requires Scipy.'

    exit(1)


def find_rst(pamldir):
    """Return absolute paths to the rst files."""
    fullpath = os.path.abspath(os.path.expanduser(pamldir))
    # Then, use the os.walk() function to recursively list all contents in this
    # directory
    outs = []
    for path, subdirs, files in os.walk(fullpath):
        for fi in files:
            if fi.endswith('_codeml'):
                outs.append(os.path.join(path, fi))
    return outs


def extract_meta(out_path):
    """Simple function to extract the OG and hypothesis from the path, because
    they are standarized."""
    parts = out_path.split('/')
    og = parts[-3]
    hyp = parts[-2]
    return (og, hyp)


def read_results(out_file, hyp):
    """Use the Biopython codeml output parser to read the output file and
    extract all the relevant information."""
    res = codeml.read(out_file)
    # extract the values that we are interested in, depending on which model
    # we have fit.
    if hyp == 'Null':
        lnl = res['NSsites'][22].get('lnL', 'NA')
        if 'parameters' not in res['NSsites'][22]:
            dat = {
                'null_lnl': lnl,
                'null_kappa': 'NA',
                'null_omega0dnds': 'NA',
                'null_omega1dnds': 'NA',
                'null_omega2dnds': 'NA',
                'null_omega0prop': 'NA',
                'null_omega1prop': 'NA',
                'null_omega2prop': 'NA'}
        else:
            kappa = res['NSsites'][22]['parameters'].get('kappa', 'NA')
            if 'site classes' not in res['NSsites'][22]['parameters']:
                omega_0_dnds = 'NA'
                omega_1_dnds = 'NA'
                omega_2_dnds = 'NA'
                omega_0_prop = 'NA'
                omega_1_prop = 'NA'
                omega_2_prop = 'NA'
            else:
                omega_0_dnds = res['NSsites'][22]['parameters']['site classes'][0].get('omega', 'NA')
                omega_1_dnds = res['NSsites'][22]['parameters']['site classes'][1].get('omega', 'NA')
                omega_2_dnds = res['NSsites'][22]['parameters']['site classes'][2].get('omega', 'NA')
                omega_0_prop = res['NSsites'][22]['parameters']['site classes'][0].get('proportion', 'NA')
                omega_1_prop = res['NSsites'][22]['parameters']['site classes'][1].get('proportion', 'NA')
                omega_2_prop = res['NSsites'][22]['parameters']['site classes'][2].get('proportion', 'NA')
            dat = {
                'null_lnl': lnl,
                'null_kappa': str(kappa),
                'null_omega0dnds': str(omega_0_dnds),
                'null_omega1dnds': str(omega_1_dnds),
                'null_omega2dnds': str(omega_2_dnds),
                'null_omega0prop': str(omega_0_prop),
                'null_omega1prop': str(omega_1_prop),
                'null_omega2prop': str(omega_2_prop)}
    elif hyp == 'Ha1':
        # For this model, #1 is maize/tandem, #0 is grass
        lnl = res['NSsites'][2].get('lnL', 'NA')
        if 'parameters' not in res['NSsites'][2]:
            dat = {
                'ha1_lnl': lnl,
                'ha1_kappa': 'NA',
                'ha1_omegagrass0dnds': 'NA',
                'ha1_omegamaize0dnds': 'NA',
                'ha1_omegagrass1dnds': 'NA',
                'ha1_omegamaize1dnds': 'NA',
                'ha1_omegagrass2dnds': 'NA',
                'ha1_omegamaize2dnds': 'NA',
                'ha1_omega0prop': 'NA',
                'ha1_omega1prop': 'NA',
                'ha1_omega2prop': 'NA'}
        else:
            kappa = res['NSsites'][2]['parameters'].get('kappa', 'NA')
            if 'site classes' not in res['NSsites'][2]['parameters']:
                omega_grass_0_dnds = 'NA'
                omega_maize_0_dnds = 'NA'
                omega_grass_1_dnds = 'NA'
                omega_maize_1_dnds = 'NA'
                omega_grass_2_dnds = 'NA'
                omega_maize_2_dnds = 'NA'
                omega_0_prop = 'NA'
                omega_1_prop = 'NA'
                omega_2_prop = 'NA'
            else:
                omega_grass_0_dnds = res['NSsites'][2]['parameters']['site classes'][0]['branch types'].get(0, 'NA')
                omega_maize_0_dnds = res['NSsites'][2]['parameters']['site classes'][0]['branch types'].get(1, 'NA')
                omega_grass_1_dnds = res['NSsites'][2]['parameters']['site classes'][1]['branch types'].get(0, 'NA')
                omega_maize_1_dnds = res['NSsites'][2]['parameters']['site classes'][1]['branch types'].get(1, 'NA')
                omega_grass_2_dnds = res['NSsites'][2]['parameters']['site classes'][2]['branch types'].get(0, 'NA')
                omega_maize_2_dnds = res['NSsites'][2]['parameters']['site classes'][2]['branch types'].get(1, 'NA')
                omega_0_prop = res['NSsites'][2]['parameters']['site classes'][0].get('proportion', 'NA')
                omega_1_prop = res['NSsites'][2]['parameters']['site classes'][1].get('proportion', 'NA')
                omega_2_prop = res['NSsites'][2]['parameters']['site classes'][2].get('proportion', 'NA')
            dat = {
                'ha1_lnl': lnl,
                'ha1_kappa': str(kappa),
                'ha1_omegagrass0dnds': str(omega_grass_0_dnds),
                'ha1_omegamaize0dnds': str(omega_maize_0_dnds),
                'ha1_omegagrass1dnds': str(omega_grass_1_dnds),
                'ha1_omegamaize1dnds': str(omega_maize_1_dnds),
                'ha1_omegagrass2dnds': str(omega_grass_2_dnds),
                'ha1_omegamaize2dnds': str(omega_maize_2_dnds),
                'ha1_omega0prop': str(omega_0_prop),
                'ha1_omega1prop': str(omega_1_prop),
                'ha1_omega2prop': str(omega_2_prop)
            }
    elif hyp == 'Ha2':
        # For this model, #1 is tandem, #0 is grass/maize
        lnl = res['NSsites'][2].get('lnL', 'NA')
        if 'parameters' not in res['NSsites'][2]:
            dat = {
                'ha2_lnl': lnl,
                'ha2_kappa': 'NA',
                'ha2_omegagrass0dnds': 'NA',
                'ha2_omegatandem0dnds': 'NA',
                'ha2_omegagrass1dnds': 'NA',
                'ha2_omegatandem1dnds': 'NA',
                'ha2_omegagrass2dnds': 'NA',
                'ha2_omegatandem2dnds': 'NA',
                'ha2_omega0prop': 'NA',
                'ha2_omega1prop': 'NA',
                'ha2_omega2prop': 'NA'
            }
        else:
            kappa = res['NSsites'][2]['parameters'].get('kappa', 'NA')
            if 'site classes' not in res['NSsites'][2]['parameters']:
                omega_grass_0_dnds = 'NA'
                omega_tandem_0_dnds = 'NA'
                omega_grass_1_dnds = 'NA'
                omega_tandem_1_dnds = 'NA'
                omega_grass_2_dnds = 'NA'
                omega_tandem_2_dnds = 'NA'
                omega_0_prop = 'NA'
                omega_1_prop = 'NA'
                omega_2_prop = 'NA'
            else:
                omega_grass_0_dnds = res['NSsites'][2]['parameters']['site classes'][0]['branch types'].get(0, 'NA')
                omega_tandem_0_dnds = res['NSsites'][2]['parameters']['site classes'][0]['branch types'].get(1, 'NA')
                omega_grass_1_dnds = res['NSsites'][2]['parameters']['site classes'][1]['branch types'].get(0, 'NA')
                omega_tandem_1_dnds = res['NSsites'][2]['parameters']['site classes'][1]['branch types'].get(1, 'NA')
                omega_grass_2_dnds = res['NSsites'][2]['parameters']['site classes'][2]['branch types'].get(0, 'NA')
                omega_tandem_2_dnds = res['NSsites'][2]['parameters']['site classes'][2]['branch types'].get(1, 'NA')
                omega_0_prop = res['NSsites'][2]['parameters']['site classes'][0].get('proportion', 'NA')
                omega_1_prop = res['NSsites'][2]['parameters']['site classes'][1].get('proportion', 'NA')
                omega_2_prop = res['NSsites'][2]['parameters']['site classes'][2].get('proportion', 'NA')
            dat = {
                'ha2_lnl': lnl,
                'ha2_kappa': str(kappa),
                'ha2_omegagrass0dnds': str(omega_grass_0_dnds),
                'ha2_omegatandem0dnds': str(omega_tandem_0_dnds),
                'ha2_omegagrass1dnds': str(omega_grass_1_dnds),
                'ha2_omegatandem1dnds': str(omega_tandem_1_dnds),
                'ha2_omegagrass2dnds': str(omega_grass_2_dnds),
                'ha2_omegatandem2dnds': str(omega_tandem_2_dnds),
                'ha2_omega0prop': str(omega_0_prop),
                'ha2_omega1prop': str(omega_1_prop),
                'ha2_omega2prop': str(omega_2_prop)
            }
    elif hyp == 'Ha3':
        # Recall that for this model, #1 is tandem and #2 is maize
        lnl = res['NSsites'][2].get('lnL', 'NA')
        if 'parameters' not in res['NSsites'][2]:
            dat = {
                'ha3_lnl': lnl,
                'ha3_kappa': 'NA',
                'ha3_omegagrass0dnds': 'NA',
                'ha3_omegatandem0dnds': 'NA',
                'ha3_omegamaize0dnds': 'NA',
                'ha3_omegagrass1dnds': 'NA',
                'ha3_omegatandem1dnds': 'NA',
                'ha3_omegamaize1dnds': 'NA',
                'ha3_omegagrass2dnds': 'NA',
                'ha3_omegatandem2dnds': 'NA',
                'ha3_omegamaize2dnds': 'NA',
                'ha3_omega0prop': 'NA',
                'ha3_omega1prop': 'NA',
                'ha3_omega2prop': 'NA'
            }
        else:
            kappa = res['NSsites'][2]['parameters'].get('kappa', 'NA')
            if 'site classes' not in res['NSsites'][2]['parameters']:
                omega_grass_0_dnds = 'NA'
                omega_tandem_0_dnds = 'NA'
                omega_maize_0_dnds = 'NA'
                omega_grass_1_dnds = 'NA'
                omega_tandem_1_dnds = 'NA'
                omega_maize_1_dnds = 'NA'
                omega_grass_2_dnds = 'NA'
                omega_tandem_2_dnds = 'NA'
                omega_maize_2_dnds = 'NA'
                omega_0_prop = 'NA'
                omega_1_prop = 'NA'
                omega_2_prop = 'NA'
            else:
                omega_grass_0_dnds = res['NSsites'][2]['parameters']['site classes'][0]['branch types'].get(0, 'NA')
                omega_tandem_0_dnds = res['NSsites'][2]['parameters']['site classes'][0]['branch types'].get(1, 'NA')
                omega_maize_0_dnds = res['NSsites'][2]['parameters']['site classes'][0]['branch types'].get(2, 'NA')
                omega_grass_1_dnds = res['NSsites'][2]['parameters']['site classes'][1]['branch types'].get(0, 'NA')
                omega_tandem_1_dnds = res['NSsites'][2]['parameters']['site classes'][1]['branch types'].get(1, 'NA')
                omega_maize_1_dnds = res['NSsites'][2]['parameters']['site classes'][1]['branch types'].get(2, 'NA')
                omega_grass_2_dnds = res['NSsites'][2]['parameters']['site classes'][2]['branch types'].get(0, 'NA')
                omega_tandem_2_dnds = res['NSsites'][2]['parameters']['site classes'][2]['branch types'].get(1, 'NA')
                omega_maize_2_dnds = res['NSsites'][2]['parameters']['site classes'][2]['branch types'].get(2, 'NA')
                omega_0_prop = res['NSsites'][2]['parameters']['site classes'][0].get('proportion', 'NA')
                omega_1_prop = res['NSsites'][2]['parameters']['site classes'][1].get('proportion', 'NA')
                omega_2_prop = res['NSsites'][2]['parameters']['site classes'][2].get('proportion', 'NA')
            dat = {
                'ha3_lnl': lnl,
                'ha3_kappa': str(kappa),
                'ha3_omegagrass0dnds': str(omega_grass_0_dnds),
                'ha3_omegatandem0dnds': str(omega_tandem_0_dnds),
                'ha3_omegamaize0dnds': str(omega_maize_0_dnds),
                'ha3_omegagrass1dnds': str(omega_grass_1_dnds),
                'ha3_omegatandem1dnds': str(omega_tandem_1_dnds),
                'ha3_omegamaize1dnds': str(omega_maize_1_dnds),
                'ha3_omegagrass2dnds': str(omega_grass_2_dnds),
                'ha3_omegatandem2dnds': str(omega_tandem_2_dnds),
                'ha3_omegamaize2dnds': str(omega_maize_2_dnds),
                'ha3_omega0prop': str(omega_0_prop),
                'ha3_omega1prop': str(omega_1_prop),
                'ha3_omega2prop': str(omega_2_prop)
            }
    return dat


def print_table(res):
    """Calculate the p-values from the LRTs, and print the table. The test that
    we will do are Ha1vNull, Ha2vNull, Ha3vNull, Ha3vHa1, and Ha3vHa2. This
    will be a huge table."""
    # Print a header
    header = ['OrthogroupID', 'Null_lnL', 'Ha1_lnL', 'Ha2_lnL', 'Ha3_lnL',
              'Ha1_v_Null', 'Ha2_v_Null', 'Ha3_v_Null', 'Ha3_v_Ha1', 'Ha3_v_Ha2',
              'Null_Kappa',
              'Null_Omega0', 'Null_PropOmega0',
              'Null_Omega1', 'Null_PropOmega1',
              'Null_Omega2', 'Null_PropOmega2',
              'Ha1_Kappa',
              'Ha1_Omega0_Grass', 'Ha1_Omega0_Maize', 'Ha1_PropOmega0',
              'Ha1_Omega1_Grass', 'Ha1_Omega1_Maize', 'Ha1_PropOmega1',
              'Ha1_Omega2_Grass', 'Ha1_Omega2_Maize', 'Ha1_PropOmega2',
              'Ha2_Kappa',
              'Ha2_Omega0_Grass', 'Ha2_Omega0_Tandem', 'Ha2_PropOmega0',
              'Ha2_Omega1_Grass', 'Ha2_Omega1_Tandem', 'Ha2_PropOmega1',
              'Ha2_Omega2_Grass', 'Ha2_Omega2_Tandem', 'Ha2_PropOmega2',
              'Ha3_Kappa',
              'Ha3_Omega0_Grass', 'Ha3_Omega0_Tandem','Ha3_Omega0_Maize', 'Ha3_PropOmega0',
              'Ha3_Omega1_Grass', 'Ha3_Omega1_Tandem','Ha3_Omega1_Maize', 'Ha3_PropOmega1',
              'Ha3_Omega2_Grass', 'Ha3_Omega2_Tandem','Ha3_Omega2_Maize', 'Ha3_PropOmega2']
    print '\t'.join(header)
    for og in sorted(res):
        # Calculate the p-values. The LRT statistic is chisq-distributed, with
        # df equal to the difference in the number of parameters of the two
        # models.
        if res[og]['ha1_lnl'] == 'NA' or res[og]['null_lnl'] == 'NA':
            ha1_v_null = 'NA'
        else:
            ha1_v_null = str(1 - chi2.cdf(2 * (res[og]['ha1_lnl'] - res[og]['null_lnl']), 1))
        if res[og]['ha2_lnl'] == 'NA' or res[og]['null_lnl'] == 'NA':
            ha2_v_null = 'NA'
        else:
            ha2_v_null = str(1 - chi2.cdf(2 * (res[og]['ha2_lnl'] - res[og]['null_lnl']), 1))
        if res[og]['ha3_lnl'] == 'NA' or res[og]['null_lnl'] == 'NA':
            ha3_v_null = 'NA'
        else:
            ha3_v_null = str(1 - chi2.cdf(2 * (res[og]['ha3_lnl'] - res[og]['null_lnl']), 2))
        if res[og]['ha3_lnl'] == 'NA' or res[og]['ha1_lnl'] == 'NA':
            ha3_v_ha1 = 'NA'
        else:
            ha3_v_ha1 = str(1 - chi2.cdf(2 * (res[og]['ha3_lnl'] - res[og]['ha1_lnl']), 1))
        if res[og]['ha3_lnl'] == 'NA' or res[og]['ha2_lnl'] == 'NA':
            ha3_v_ha2 = 'NA'
        else:
            ha3_v_ha2 = str(1 - chi2.cdf(2 * (res[og]['ha3_lnl'] - res[og]['ha2_lnl']), 1))
        # Now, print out this big gross table
        toprint = [og, str(res[og]['null_lnl']), str(res[og]['ha1_lnl']),
                   str(res[og]['ha2_lnl']), str(res[og]['ha3_lnl']), ha1_v_null,
                   ha2_v_null, ha3_v_null, ha3_v_ha1, ha3_v_ha2,
                   res[og]['null_kappa'],
                   res[og]['null_omega0dnds'], res[og]['null_omega0prop'],
                   res[og]['null_omega1dnds'], res[og]['null_omega1prop'],
                   res[og]['null_omega2dnds'], res[og]['null_omega2prop'],
                   res[og]['ha1_kappa'],
                   res[og]['ha1_omegagrass0dnds'], res[og]['ha1_omegamaize0dnds'], res[og]['ha1_omega0prop'],
                   res[og]['ha1_omegagrass1dnds'], res[og]['ha1_omegamaize1dnds'], res[og]['ha1_omega1prop'],
                   res[og]['ha1_omegagrass2dnds'], res[og]['ha1_omegamaize2dnds'], res[og]['ha1_omega2prop'],
                   res[og]['ha2_kappa'],
                   res[og]['ha2_omegagrass0dnds'], res[og]['ha2_omegatandem0dnds'], res[og]['ha2_omega0prop'],
                   res[og]['ha2_omegagrass1dnds'], res[og]['ha2_omegatandem1dnds'], res[og]['ha2_omega1prop'],
                   res[og]['ha2_omegagrass2dnds'], res[og]['ha2_omegatandem2dnds'], res[og]['ha2_omega2prop'],
                   res[og]['ha3_kappa'],
                   res[og]['ha3_omegagrass0dnds'], res[og]['ha3_omegatandem0dnds'], res[og]['ha3_omegamaize0dnds'], res[og]['ha3_omega0prop'],
                   res[og]['ha3_omegagrass1dnds'], res[og]['ha3_omegatandem1dnds'], res[og]['ha3_omegamaize1dnds'], res[og]['ha3_omega1prop'],
                   res[og]['ha3_omegagrass2dnds'], res[og]['ha3_omegatandem2dnds'], res[og]['ha3_omegamaize2dnds'], res[og]['ha3_omega2prop']]
        print '\t'.join(toprint)
    return


def main(pamldir):
    """Main function."""
    outfiles = find_rst(pamldir)
    results = {}
    for r in outfiles:
        og, h = extract_meta(r)
        model_res = read_results(r, h)
        if og in results:
            results[og].update(model_res)
        else:
            results[og] = model_res
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
