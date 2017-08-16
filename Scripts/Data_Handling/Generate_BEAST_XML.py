#!/usr/bin/env python
"""This will be a big ugly script. It will generate a BEAST XML file for running
a tandem duplication dating analysis. Takes one argument:
    1) NEXUS alignment of the tandem duplicates
"""

import sys
import os
import xml.etree.ElementTree as ET
import xml.dom.minidom
try:
    from Bio import SeqIO
except ImportError:
    print 'This script requires Biopython.'
    exit(1)



class BEASTXML(object):
    """A class that stores the BEAST XML control file that will be used in
    the run. It looks like BEAST has trouble with multiple data sections, so
    one XML file per gene will be the way we go here."""

    def __init__(self):
        """Initialize the XML file with the headers and namespace etc."""
        self.beast = ET.Element(
            'beast',
            attrib={
                'beautitemplate': 'Standard',
                'beautistatus': '',
                'namespace': 'beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood',
                'version': '2.4',
                'required': ''})
        self.locname = ''
        self.maize = []
        return

    def add_data(self, aln):
        """Add the data section to the XML file."""
        # we want to set the name of the locus/alignment
        fname = os.path.basename(aln)
        self.locname = fname.rsplit('.', 1)[0]
        dat = ET.Element('data', id=self.locname, name='alignment')
        aln_seqs = SeqIO.parse(aln, 'fasta')
        for s in aln_seqs:
            attribs = {
                'id': 'seq_' + s.id,
                'taxon': s.id,
                'totalcount': '4',
                'value': str(s.seq)}
            dat.append(ET.Element('sequence', attrib=attribs))
            # Also keep track of the maize IDs for our date prior.
            if s.id.startswith('Zm'):
                self.maize.append(s.id)
        self.beast.append(dat)
        return

    def add_map(self):
        """Add the distributions in the <map> tags. This will be pretty easy
        since it has no external data to depend on."""
        unif = ET.Element('map', attrib={'name': 'Uniform'})
        unif.text = 'beast.math.distributions.Uniform'
        exp = ET.Element('map', attrib={'name': 'Exponential'})
        exp.text = 'beast.math.distributions.Exponential'
        lognorm = ET.Element('map', attrib={'name': 'LogNormal'})
        lognorm.text = 'beast.math.distributions.LogNormalDistributionModel'
        norm = ET.Element('map', attrib={'name': 'Normal'})
        norm.text = 'beast.math.distributions.Normal'
        beta = ET.Element('map', attrib={'name': 'Beta'})
        beta.text = 'beast.math.distributions.Beta'
        gamma = ET.Element('map', attrib={'name': 'Gamma'})
        gamma.text = 'beast.math.distributions.Gamma'
        laplace = ET.Element('map', attrib={'name': 'LaplaceDistribution'})
        laplace.text = 'beast.math.distributions.LaplaceDistribution'
        prior = ET.Element('map', attrib={'name': 'prior'})
        prior.text = 'beast.math.distributions.Prior'
        invgamma = ET.Element('map', attrib={'name': 'InverseGamma'})
        invgamma.text = 'beast.math.distributions.InverseGamma'
        oneonx = ET.Element('map', attrib={'name': 'OneOnX'})
        oneonx.text = 'beast.math.distributions.OneOnX'
        self.beast.append(unif)
        self.beast.append(exp)
        self.beast.append(lognorm)
        self.beast.append(norm)
        self.beast.append(beta)
        self.beast.append(gamma)
        self.beast.append(laplace)
        self.beast.append(prior)
        self.beast.append(invgamma)
        self.beast.append(oneonx)
        return

    def subst_model(self, state_elem):
        """Puts the parameters of the nucleotide substitution model into the 
        <state> element. We will hard-code in a GTR model."""
        freq_param = ET.Element(
            'parameter',
            attrib={
                'id': 'freqParameter.s:' + self.locname,
                'dimension': '4',
                'lower': '0.0',
                'name': 'stateNode',
                'upper': '1.0'}
                )
        freq_param.text = '0.25'
        ac_rate = ET.Element(
            'parameter',
            attrib={
                'id': 'rateAC.s:' + self.locname,
                'lower': '0.0',
                'name': 'stateNode'}
                )
        ac_rate.text = '1.0'
        ag_rate = ET.Element(
            'parameter',
            attrib={
                'id': 'rateAG.s:' + self.locname,
                'lower': '0.0',
                'name': 'stateNode'}
                )
        ag_rate.text = '1.0'
        at_rate = ET.Element(
            'parameter',
            attrib={
                'id': 'rateAT.s:' + self.locname,
                'lower': '0.0',
                'name': 'stateNode'}
                )
        at_rate.text = '1.0'
        cg_rate = ET.Element(
            'parameter',
            attrib={
                'id': 'rateCG.s:' + self.locname,
                'lower': '0.0',
                'name': 'stateNode'}
                )
        cg_rate.text = '1.0'
        gt_rate = ET.Element(
            'parameter',
            attrib={
                'id': 'rateGT.s:' + self.locname,
                'lower': '0.0',
                'name': 'stateNode'}
                )
        gt_rate.text = '1.0'
        ct_rate = ET.Element(
            'parameter',
            attrib={
                'id': 'rateCT.s:' + self.locname,
                'lower': '0.0',
                'name': 'stateNode'}
                )
        ct_rate.text = '1.0'
        mut_rate = ET.Element(
            'parameter',
            attrib={
                'id': 'mutationRate.s:' + self.locname,
                'name': 'stateNode'}
                )
        mut_rate.text = '1.0'
        clock_rate = ET.Element(
            'parameter',
            attrib={
                'id': 'clockrates.c:' + self.locname,
                'dimension': '8',
                'name': 'stateNode'}
                )
        clock_rate.text = '1.0'
        birth_rate = ET.Element(
            'parameter',
            attrib={
                'id': 'birthRateY.t:' + self.locname,
                'name': 'stateNode'}
                )
        birth_rate.text = '1.0'
        mean_clock = ET.Element(
            'parameter',
            attrib={
                'id': 'meanClockRate.c:' + self.locname,
                'name': 'stateNode'}
                )
        mean_clock.text = '1.0'
        statenode = ET.Element(
            'stateNode',
            attrib={
                'id': 'Indicators.c:' + self.locname,
                'spec': 'parameter.BooleanParameter',
                'dimension': '8'}
                )
        statenode.text = 'false'
        # Append all this stuff to the <state> elem
        state = state_elem
        state.append(freq_param)
        state.append(ac_rate)
        state.append(ag_rate)
        state.append(at_rate)
        state.append(cg_rate)
        state.append(gt_rate)
        state.append(ct_rate)
        state.append(mut_rate)
        state.append(statenode)
        state.append(clock_rate)
        state.append(birth_rate)
        state.append(mean_clock)
        return state

    def add_init(self, state):
        """Add the <init> element. This is a somewhat small section for our use
        case. We use a random starting tree."""
        i = ET.Element(
            'init',
            attrib={
                'id': 'RandomTree.t:' + self.locname,
                'spec': 'beast.evolution.tree.RandomTree',
                'estimate': 'false',
                'initial': '@Tree.t:' + self.locname,
                'taxa': '@' + self.locname}
                )
        pop = ET.Element(
            'populationModel',
            attrib={
                'id': 'ConstantPopulation0.t:' + self.locname,
                'spec': 'ConstantPopulation'}
                )
        param = ET.Element(
            'parameter',
            attrib={
                'id': 'randomPopSize.t:' + self.locname,
                'name': 'popSize'}
                )
        param.text = '1.0'
        pop.append(param)
        i.append(pop)
        state_el = state
        state_el.append(i)
        return state_el

    def add_distr(self, state):
        """Add the big ugly <distribution> elements to the <run> element."""
        # The first one is the posterior
        posterior = ET.Element(
            'distribution',
            attrib={
                'id': 'posterior',
                'spec': 'util.CompoundDistribution'}
                )
        # A prior is nested within the posterior element
        prior = ET.Element(
            'distribution',
            attrib={
                'id': 'prior',
                'spec': 'util.CompoundDistribution'}
                )
        # And our priors are listed in this element. We are using the Calibrated
        # Yule model.
        cal_yule = ET.Element(
            'distribution',
            attrib={
                'id': 'CalibratedYuleModel.t:' + self.locname,
                'spec': 'beast.evolution.speciation.CalibratedYuleModel',
                'birthRate': '@birthRateY.t:' + self.locname,
                'tree': '@Tree.t:' + self.locname}
                )
        # Put the calibrated Yule model into the priors
        prior.append(cal_yule)
        # Then add a prior on rate changes
        rrate_changes = ET.Element(
            'prior',
            attrib={
                'id': 'RRateChangesPrior.c:' + self.locname,
                'name': 'distribution'}
                )
        x = ET.Element(
            'x',
            attrib={
                'id': 'RRateChanges.c:' + self.locname,
                'spec': 'util.Sum'}
                )
        x.append(ET.Element('arg', attrib={'idref': 'Indicators.c:' + self.locname}))
        rrate_changes.append(x)
        # Then, we make a section for a poisson distribution. Not sure where the
        # value comes from, but it's in the XML.
        distr = ET.Element(
            'distr',
            attrib={
                'id': 'Poisson.0',
                'spec': 'beast.math.distributions.Poisson'}
                )
        param = ET.Element(
            'parameter',
            attrib={
                'id': 'RealParameter.13',
                'estimate': 'false',
                'name': 'lambda'}
                )
        param.text = '0.6931471805599453'
        distr.append(param)
        rrate_changes.append(distr)
        prior.append(rrate_changes)
        # After the RRateChanges prior, we have one on the Cal. Yule Birth Rate
        cy_br = ET.Element(
            'prior',
            attrib={
                'id': 'CalibratedYuleBirthRatePrior.t:' + self.locname,
                'name': 'distribution',
                'x': '@birthRateY.t:' + self.locname}
                )
        cy_br.append(ET.Element(
            'Uniform',
            attrib={
                'id': 'Uniform.4',
                'name': 'distr',
                'upper': '1000.0'}
                ))
        prior.append(cy_br)
        # Next, we have an RRatesPrior. These will have mutational rates priors
        rrates_prior = ET.Element(
            'prior',
            attrib={
                'id': 'RRatesPrior.c:s' + self.locname,
                'name': 'distribution',
                'x': '@clockrates.c:' + self.locname}
                )
        gamma = ET.Element(
            'Gamma',
            attrib={'id': 'Gamma.0', 'name': 'distr'})
        alpha = ET.Element(
            'parameter',
            attrib={
                'id': 'RealParameter.1',
                'estimate': 'false',
                'name': 'alpha'}
                )
        alpha.text = '0.05'
        beta = ET.Element(
            'parameter',
            attrib={
                'id': 'RealParameter.2',
                'estimate': 'false',
                'name': 'beta'}
                )
        beta.text = '10.0'
        gamma.append(alpha)
        gamma.append(beta)
        rrates_prior.append(gamma)
        prior.append(rrates_prior)
        prior.append(ET.Element(
            'prior',
            attrib={
                'id': 'RateACPrior.s:' + self.locname,
                'distr': '@Gamma.0',
                'name': 'distribution',
                'x': '@rateAC.s:' + self.locname}))
        # And we have to repeat that for each of AC, AT, GC, CT, GT
        ag_rate = ET.Element(
            'prior',
            attrib={
                'id': 'RateAGPrior.s:' + self.locname,
                'name': 'distribution',
                'x': '@rateAG.s:' + self.locname}
                )
        gamma = ET.Element(
            'Gamma',
            attrib={'id': 'Gamma.1', 'name': 'distr'})
        alpha = ET.Element(
            'parameter',
            attrib={
                'id': 'RealParameter.3',
                'estimate': 'false',
                'name': 'alpha'}
                )
        alpha.text = '0.05'
        beta = ET.Element(
            'parameter',
            attrib={
                'id': 'RealParameter.4',
                'estimate': 'false',
                'name': 'beta'}
                )
        beta.text = '20.0'
        gamma.append(alpha)
        gamma.append(beta)
        ag_rate.append(gamma)
        prior.append(ag_rate)
        at_rate = ET.Element(
            'prior',
            attrib={
                'id': 'RateATPrior.s:' + self.locname,
                'name': 'distribution',
                'x': '@rateAT.s:' + self.locname}
                )
        gamma = ET.Element(
            'Gamma',
            attrib={'id': 'Gamma.2', 'name': 'distr'})
        alpha = ET.Element(
            'parameter',
            attrib={
                'id': 'RealParameter.5',
                'estimate': 'false',
                'name': 'alpha'}
                )
        alpha.text = '0.05'
        beta = ET.Element(
            'parameter',
            attrib={
                'id': 'RealParameter.6',
                'estimate': 'false',
                'name': 'beta'}
                )
        beta.text = '10.0'
        gamma.append(alpha)
        gamma.append(beta)
        at_rate.append(gamma)
        prior.append(at_rate)
        cg_rate = ET.Element(
            'prior',
            attrib={
                'id': 'RateCGPrior.s:' + self.locname,
                'name': 'distribution',
                'x': '@rateCG.s:' + self.locname}
                )
        gamma = ET.Element(
            'Gamma',
            attrib={'id': 'Gamma.3', 'name': 'distr'})
        alpha = ET.Element(
            'parameter',
            attrib={
                'id': 'RealParameter.7',
                'estimate': 'false',
                'name': 'alpha'}
                )
        alpha.text = '0.05'
        beta = ET.Element(
            'parameter',
            attrib={
                'id': 'RealParameter.8',
                'estimate': 'false',
                'name': 'beta'}
                )
        beta.text = '10.0'
        gamma.append(alpha)
        gamma.append(beta)
        cg_rate.append(gamma)
        prior.append(cg_rate)
        ct_rate = ET.Element(
            'prior',
            attrib={
                'id': 'RateCTPrior.s:' + self.locname,
                'name': 'distribution',
                'x': '@rateCT.s:' + self.locname}
                )
        gamma = ET.Element(
            'Gamma',
            attrib={'id': 'Gamma.4', 'name': 'distr'})
        alpha = ET.Element(
            'parameter',
            attrib={
                'id': 'RealParameter.9',
                'estimate': 'false',
                'name': 'alpha'}
                )
        alpha.text = '0.05'
        beta = ET.Element(
            'parameter',
            attrib={
                'id': 'RealParameter.10',
                'estimate': 'false',
                'name': 'beta'}
                )
        beta.text = '20.0'
        gamma.append(alpha)
        gamma.append(beta)
        ct_rate.append(gamma)
        prior.append(ct_rate)
        gt_rate = ET.Element(
            'prior',
            attrib={
                'id': 'RateGTPrior.s:' + self.locname,
                'name': 'distribution',
                'x': '@rateGT.s:' + self.locname}
                )
        gamma = ET.Element(
            'Gamma',
            attrib={'id': 'Gamma.5', 'name': 'distr'})
        alpha = ET.Element(
            'parameter',
            attrib={
                'id': 'RealParameter.11',
                'estimate': 'false',
                'name': 'alpha'}
                )
        alpha.text = '0.05'
        beta = ET.Element(
            'parameter',
            attrib={
                'id': 'RealParameter.12',
                'estimate': 'false',
                'name': 'beta'}
                )
        beta.text = '10.0'
        gamma.append(alpha)
        gamma.append(beta)
        gt_rate.append(gamma)
        prior.append(gt_rate)
        # Then, we add the Maize1-Maize2 TMRCA prior
        tmrca_prior = ET.Element(
            'distribution',
            attrib={
                'id': 'M1-M2.prior',
                'spec': 'beast.math.distributions.MRCAPrior',
                'monophyletic': 'true',
                'tree': '@Tree.t:' + self.locname}
                )
        # Within this TMRCA prior, we add the taxonset
        tx_set = ET.Element('taxonset', attrib={'id': 'M1-M2', 'spec': 'TaxonSet'})
        # All the maize genes are within this taxonset
        for m_gene in self.maize:
            tx_set.append(ET.Element('taxon', attrib={'id': m_gene, 'spec': 'Taxon'}))
        # Tack the taxonset onto the TMRCA
        tmrca_prior.append(tx_set)
        # Then, add a normal distribution to the TMRCA prior. We use a mean of
        # 11.9 MYA and a sd of 1 MY for the prior
        norm = ET.Element(
            'Normal',
            attrib={
                'id': 'Normal.0',
                'name': 'distr'}
                )
        n_mean = ET.Element(
            'parameter',
            attrib={
                'id': 'RealParameter.14',
                'estimate': 'false',
                'name': 'mean'}
                )
        n_mean.text = '11.9'
        n_sd = ET.Element(
            'parameter',
            attrib={
                'id': 'RealParameter.15',
                'estimate': 'false',
                'name': 'sigma'}
                )
        n_sd.text = '1.0'
        # Add the parameters to the distribution
        norm.append(n_mean)
        norm.append(n_sd)
        # Add the distribution to the TMRCA prior
        tmrca_prior.append(norm)
        # And add the TMRCA prior to the big prior
        prior.append(tmrca_prior)
        posterior.append(prior)
        # Next, we have to make a new distribution for the likelihood.
        likelihood = ET.Element(
            'distribution',
            attrib={
                'id': 'likelihood',
                'spec': 'util.CompoundDistribution',
                'useThreads': 'true'}
                )
        # Add a treeLikelihood to this one
        tree_likelihood = ET.Element(
            'distribution',
            attrib={
                'id': 'treeLikelihood.' + self.locname,
                'spec': 'ThreadedTreeLikelihood',
                'data': '@' + self.locname,
                'tree': '@Tree.t:' + self.locname}
                )
        #  Within the tree likelihood, we have a site model
        site_model = ET.Element(
            'siteModel',
            attrib={
                'id': 'SiteModel.s:' + self.locname,
                'spec': 'SiteModel',
                'gammaCategoryCount': '1',
                'mutationRate': '@mutationRate.s:' + self.locname}
                )
        # Add parameters to the site model
        shape = ET.Element(
            'parameter',
            attrib={
                'id': 'gammaShape.s:' + self.locname,
                'estimate': 'false',
                'name': 'shape'}
                )
        shape.text = '1.0'
        prop_inv = ET.Element(
            'parameter',
            attrib={
                'id': 'proportionInvariant.s:' + self.locname,
                'estimate': 'false',
                'lower': '0.0',
                'name': 'proportionInvariant',
                'upper': '1.0'}
                )
        prop_inv.text = '0.0'
        # Add these to the gamma element
        site_model.append(shape)
        site_model.append(prop_inv)
        # Next, a bit about the substitution model
        subst = ET.Element(
            'substModel',
            attrib={
                'id': 'gtr.s:' + self.locname,
                'spec': 'GTR',
                'rateAC': '@rateAC.s:' + self.locname,
                'rateAG': '@rateAG.s:' + self.locname,
                'rateAT': '@rateAT.s:' + self.locname,
                'rateCG': '@rateCG.s:' + self.locname,
                'rateCT': '@rateCT.s:' + self.locname,
                'rateGT': '@rateGT.s:' + self.locname}
                )
        subst.append(ET.Element(
            'frequencies',
            attrib={
                'id': 'estimatedFreqs.s:' + self.locname,
                'spec': 'Frequencies',
                'frequencies': '@freqParameter.s:' + self.locname}
                ))
        # Add the substitution model to the site model
        site_model.append(subst)
        # And add the site model to the tree likelihood
        tree_likelihood.append(site_model)
        # Add a branchratemodel to the tree likelihood. This is where the
        # local clock is defined.
        tree_likelihood.append(ET.Element(
            'branchRateModel',
            attrib={
                'id': 'RandomLocalClock.c:' + self.locname,
                'spec': 'beast.evolution.branchratemodel.RandomLocalClockModel',
                'clock.rate': '@meanClockRate.c:' + self.locname,
                'indicators': '@Indicators.c:' + self.locname,
                'rates': '@clockrates.c:' + self.locname,
                'tree': '@Tree.t:' + self.locname}
                ))
        # Add the tree likelihood to the posterior
        likelihood.append(tree_likelihood)
        # And put the likelihood into the posterior
        posterior.append(likelihood)
        state_el = state
        state_el.append(posterior)
        return state_el

    def add_operators(self, run_elem):
        """Add the <operator> elements to the run element.""" 
        freq_exchanger = ET.Element(
            'operator',
            attrib={
                'id': 'FrequenciesExchanger.s:' + self.locname,
                'spec': 'DeltaExchangeOperator',
                'delta': '0.01',
                'weight': '0.1'}
                )
        freq_exchanger.append(ET.Element(
            'parameter',
            attrib={
                'idref': 'freqParameter.s:' + self.locname}))
        ac_scaler = ET.Element(
            'operator',
            attrib={
                'id': 'RateACScaler.s:' + self.locname,
                'spec': 'ScaleOperator',
                'parameter': '@rateAC.s:' + self.locname,
                'scaleFactor': '0.5',
                'weight': '0.1'})
        ag_scaler = ET.Element(
            'operator',
            attrib={
                'id': 'RateAGScaler.s:' + self.locname,
                'spec': 'ScaleOperator',
                'parameter': '@rateAG.s:' + self.locname,
                'scaleFactor': '0.5',
                'weight': '0.1'})
        at_scaler = ET.Element(
            'operator',
            attrib={
                'id': 'RateATScaler.s:' + self.locname,
                'spec': 'ScaleOperator',
                'parameter': '@rateAT.s:' + self.locname,
                'scaleFactor': '0.5',
                'weight': '0.1'})
        cg_scaler = ET.Element(
            'operator',
            attrib={
                'id': 'RateCGScaler.s:' + self.locname,
                'spec': 'ScaleOperator',
                'parameter': '@rateCG.s:' + self.locname,
                'scaleFactor': '0.5',
                'weight': '0.1'})
        gt_scaler = ET.Element(
            'operator',
            attrib={
                'id': 'RateGTScaler.s:' + self.locname,
                'spec': 'ScaleOperator',
                'parameter': '@rateGT.s:' + self.locname,
                'scaleFactor': '0.5',
                'weight': '0.1'})
        ct_scaler = ET.Element(
            'operator',
            attrib={
                'id': 'RateCTScaler.s:' + self.locname,
                'spec': 'ScaleOperator',
                'parameter': '@rateCT.s:' + self.locname,
                'scaleFactor': '0.5',
                'weight': '0.1'})
        mean_mut_op = ET.Element(
            'operator',
            attrib={
                'id': 'FixMeanMutationRatesOperator',
                'spec': 'DeltaExchangeOperator',
                'delta': '0.75',
                'weight': '2.0'})
        mean_mut_op.append(ET.Element('parameter', attrib={'idref': 'mutationRate.s:' + self.locname}))
        wv = ET.Element(
            'weightvector',
            attrib={
                'id': 'weightparameter',
                'spec': 'parameter.IntegerParameter',
                'estimate': 'false',
                'lower': '0',
                'upper': '0'})
        wv.text = '3264'
        mean_mut_op.append(wv)
        bitflip = ET.Element(
            'operator',
            attrib={
                'id': 'IndicatorsBitFlip.c:' + self.locname,
                'spec': 'BitFlipOperator',
                'parameter': '@Indicators.c:' + self.locname,
                'weight': '15.0'})
        clock_scaler = ET.Element(
            'operator',
            attrib={
                'id': 'ClockRateScaler.c:' + self.locname,
                'spec': 'ScaleOperator',
                'parameter': '@clockrates.c:' + self.locname,
                'scaleFactor': '0.5',
                'weight': '15.0'})
        tree_scaler = ET.Element(
            'operator',
            attrib={
                'id': 'CalibratedYuleModelTreeScaler.t:' + self.locname,
                'spec': 'ScaleOperator',
                'scaleFactor': '0.5',
                'tree': '@Tree.t:' + self.locname,
                'weight': '3.0'})
        tree_root_scaler = ET.Element(
            'operator',
            attrib={
                'id': 'CalibratedYuleModelTreeRootScaler.t:' + self.locname,
                'spec': 'ScaleOperator',
                'rootOnly': 'true',
                'scaleFactor': '0.5',
                'tree': '@Tree.t:' + self.locname,
                'weight': '3.0'})
        unif = ET.Element(
            'operator',
            attrib={
                'id': 'CalibratedYuleModelUniformOperator.t:' + self.locname,
                'spec': 'Uniform',
                'tree': '@Tree.t:' + self.locname,
                'weight': '30.0'})
        subtree_slide = ET.Element(
            'operator',
            attrib={
                'id': 'CalibratedYuleModelSubtreeSlide.t:' + self.locname,
                'spec': 'SubtreeSlide',
                'tree': '@Tree.t:' + self.locname,
                'weight': '15.0'})
        narrow = ET.Element(
                'operator',
                attrib={
                    'id': 'CalibratedYuleModelNarrow.t:' + self.locname,
                    'spec': 'Exchange',
                    'tree': '@Tree.t:' + self.locname,
                    'weight': '15.0'})
        wide = ET.Element(
            'operator',
            attrib={
                'id': 'CalibratedYuleModelWide.t:' + self.locname,
                'spec': 'Exchange',
                'isNarrow': 'false',
                'tree': '@Tree.t:' + self.locname,
                'weight': '3.0'})
        wilsonbalding = ET.Element(
            'operator',
            attrib={
                'id': 'CalibratedYuleModelWilsonBalding.t:' + self.locname,
                'spec': 'WilsonBalding',
                'tree': '@Tree.t:' + self.locname,
                'weight': '3.0'})
        br_scaler = ET.Element(
            'operator',
            attrib={
                'id': 'CalibratedYuleBirthRateScaler.t:' + self.locname,
                'spec': 'ScaleOperator',
                'parameter': '@birthRateY.t:' + self.locname,
                'scaleFactor': '0.75',
                'weight': '3.0'})
        rand_clock_scaler = ET.Element(
            'operator',
            attrib={
                'id': 'randomClockScaler.c:' + self.locname,
                'spec': 'ScaleOperator',
                'parameter': '@meanClockRate.c:' + self.locname,
                'scaleFactor': '0.5',
                'weight': '1.0'})
        clock_up_down = ET.Element(
            'operator',
            attrib={
                'id': 'randomClockUpDownOperator.c:' + self.locname,
                'spec': 'UpDownOperator',
                'scaleFactor': '0.75',
                'weight': '3.0'})
        clock_up_down.append(ET.Element('up', attrib={'idref': 'meanClockRate.c:' + self.locname}))
        clock_up_down.append(ET.Element('down', attrib={'idref': 'Tree.t:' + self.locname}))
        # Put all the operators into the run element
        r_elem = run_elem
        r_elem.append(freq_exchanger)
        r_elem.append(ac_scaler)
        r_elem.append(ag_scaler)
        r_elem.append(at_scaler)
        r_elem.append(cg_scaler)
        r_elem.append(gt_scaler)
        r_elem.append(ct_scaler)
        r_elem.append(mean_mut_op)
        r_elem.append(bitflip)
        r_elem.append(clock_scaler)
        r_elem.append(tree_scaler)
        r_elem.append(tree_root_scaler)
        r_elem.append(unif)
        r_elem.append(subtree_slide)
        r_elem.append(narrow)
        r_elem.append(wide)
        r_elem.append(wilsonbalding)
        r_elem.append(br_scaler)
        r_elem.append(rand_clock_scaler)
        r_elem.append(clock_up_down)
        return r_elem

    def add_loggers(self, run_elem):
        """Finally, add the logger elements to the run element."""
        # Start with the trace log
        tracelog = ET.Element(
            'logger',
            attrib={
                'id': 'tracelog',
                'fileName': 'Beast_GTR_' + self.locname + '.log',
                'logEvery': '1000',
                'model': '@posterior',
                'sanitiseHeaders': 'true',
                'sort': 'smart'})
        tracelog.append(ET.Element('log', idref='posterior'))
        tracelog.append(ET.Element('log', idref='likelihood'))
        tracelog.append(ET.Element('log', idref='prior'))
        tracelog.append(ET.Element('log', idref='treeLikelihood.' + self.locname))
        tracelog.append(ET.Element('log', id='TreeHeight.t:' + self.locname, spec='beast.evolution.tree.TreeHeightLogger', tree='@Tree.t:' + self.locname))
        tracelog.append(ET.Element('log', idref='rateAC.s:' + self.locname))
        tracelog.append(ET.Element('log', idref='rateAG.s:' + self.locname))
        tracelog.append(ET.Element('log', idref='rateAT.s:' + self.locname))
        tracelog.append(ET.Element('log', idref='rateCG.s:' + self.locname))
        tracelog.append(ET.Element('log', idref='rateGT.s:' + self.locname))
        tracelog.append(ET.Element('log', idref='rateCT.s:' + self.locname))
        tracelog.append(ET.Element('log', idref='mutationRate.s:' + self.locname))
        tracelog.append(ET.Element('log', idref='clockrates.c:' + self.locname))
        tracelog.append(ET.Element('log', idref='RRateChanges.c:' + self.locname))
        tracelog.append(ET.Element('log', idref='CalibratedYuleModel.t:' + self.locname))
        tracelog.append(ET.Element('log', idref='birthRateY.t:' + self.locname))
        tracelog.append(ET.Element('log', idref='M1-M2.prior'))
        tracelog.append(ET.Element('log', idref='meanClockRate.c:' + self.locname))
        tracelog.append(ET.Element('log', idref='freqParameter.s:' + self.locname))
        tracelog.append(ET.Element('log', idref='Indicators.c:' + self.locname))
        # Then add the screen log
        screenlog = ET.Element(
            'logger',
            attrib={
                'id': 'screenlog',
                'logEvery': '1000'})
        screenlog.append(ET.Element('log', idref='posterior'))
        screenlog.append(ET.Element('log', id='ESS.0', spec='util.ESS', arg='@posterior'))
        screenlog.append(ET.Element('log', idref='likelihood'))
        screenlog.append(ET.Element('log', idref='prior'))
        # And lastly, the tree log
        treelog = ET.Element(
            'logger',
            attrib={
                'id': 'treelog.t:' + self.locname,
                'fileName': '$(tree).trees',
                'logEvery': '1000',
                'mode': 'tree'})
        treelog.append(ET.Element(
            'log',
            attrib={
                'id': 'TreeWithMetaDataLogger.t:' + self.locname,
                'spec': 'beast.evolution.tree.TreeWithMetaDataLogger',
                'branchratemodel': '@RandomLocalClock.c:' + self.locname,
                'tree': '@Tree.t:' + self.locname}))
        r_elem = run_elem
        r_elem.append(tracelog)
        r_elem.append(screenlog)
        r_elem.append(treelog)
        return r_elem

    def add_run(self):
        """Add the run element, which defines the MCMC parameters, in addition
        to a bunch of distributions and priors. This is a somewhat complicated
        section of the BEAST XML."""
        # Start the big <run> element
        r = ET.Element(
            'run',
            id='mcmc',
            spec='MCMC',
            chainLength='10000000')
        # There is a <state> element at the top of the <run>
        s = ET.Element(
            'state',
            attrib={'id': 'state', 'storeEvery': '5000'})
        # Then, there is a <tree> under the state
        t = ET.Element(
            'tree',
            attrib={'id': 'Tree.t:' + self.locname, 'name': 'stateNode'})
        # THEN, under the <tree>, there is a <taxonset>
        tx = ET.Element(
            'taxonset',
            attrib={'id': 'TaxonSet.' + self.locname, 'spec': 'TaxonSet'})
        # And finally, the <alignment> tag under that
        tx.append(ET.Element('alignment', attrib={'idref': self.locname}))
        # Put the taxonset into the tree, and the tree into the state
        t.append(tx)
        s.append(t)
        # Then, put the <parameter> tags
        s_param = self.subst_model(s)
        r.append(s_param)
        # Put the initial tree into the run
        r_init = self.add_init(r)
        # Add the prior and posterior
        r_dist = self.add_distr(r_init)
        # Then, we have to add various <operator> elements
        r_op = self.add_operators(r_dist)
        # And finally, the logger
        r_log = self.add_loggers(r_op)
        self.beast.append(r_log)
        return


def main(aln):
    """Main function."""
    b = BEASTXML()
    b.add_data(aln)
    b.add_map()
    b.add_run()
    x = xml.dom.minidom.parseString(ET.tostring(b.beast))
    print x.toprettyxml()
    return


if len(sys.argv) != 2:
    print """This will be a big ugly script. It will generate a BEAST XML file for running
a tandem duplication dating analysis. Takes one argument:
    1) NEXUS alignment of the tandem duplicates"""
    exit(1)
else:
    main(sys.argv[1])
