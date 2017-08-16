#!/usr/bin/env python
"""Generate two CodeML control files that specify a null hypothesis (all genes
evolve at the same rate) and an alternate Clade Model C hypothesis (tandems !=
maize != grass genes), run the two models, and return their likelihoods. The
CMC models get parameters model = 3 and Nssites = 2.
