import xml.etree.ElementTree as ET
import functions as f


def define_operators(jointPrior,joint, MCMC, fasta, priors_table, partition_file, rootHeight_mean, rootHeight_stdev, rootHeight_offset):
    subs_HKY(jointPrior, partition_file)
    clock_rate(jointPrior, partition_file)
    treeModel(jointPrior, rootHeight_mean, rootHeight_stdev, rootHeight_offset)
    tip(jointPrior, fasta, priors_table)
    skygrid(jointPrior)
    strict_clock(jointPrior, partition_file)
    likelihood(joint, MCMC, partition_file)


def subs_HKY(jointPrior, partition_file):
    # Substitution Model Priors
    for p, exclude, every3, susb in f.get_partition_info(partition_file):
        p = f.Partition(p, exclude, every3, susb)

        if 'HKY' in p.subs_model:
            if p.every3:
                for i in range(1, 4):
                    jointLogNormal = ET.SubElement(jointPrior, 'logNormalPrior',
                                                   attrib={'mu': '1.0', 'sigma': '1.25', 'offset': '0.0'})
                    ET.SubElement(jointLogNormal, 'parameter', attrib={'idref': '{}.CP{}.kappa'.format(p.p, i)})
                for i in range(1, 4):
                    jointDirichlet = ET.SubElement(jointPrior, 'dirichletPrior',
                                                   attrib={'alpha': '1.0', 'sumsTo': '1.0'})
                    ET.SubElement(jointDirichlet, 'parameter', attrib={'idref': '{}.CP{}.frequencies'.format(p.p, i)})
            else:
                jointLogNormal = ET.SubElement(jointPrior, 'logNormalPrior',
                                               attrib={'mu': '1.0', 'sigma': '1.25', 'offset': '0.0'})
                ET.SubElement(jointLogNormal, 'parameter', attrib={'idref': '{}.kappa'.format(p.p)})
                jointDirichlet = ET.SubElement(jointPrior, 'dirichletPrior', attrib={'alpha': '1.0', 'sumsTo': '1.0'})
                ET.SubElement(jointDirichlet, 'parameter', attrib={'idref': '{}.frequencies'.format(p.p)})
        if 'G' in p.subs_model:
            if p.every3:
                for i in range(1, 4):
                    jointExponential = ET.SubElement(jointPrior, 'exponentialPrior',
                                                     attrib={'mean': '0.5', 'offset': '0.0'})
                    ET.SubElement(jointExponential, 'parameter', attrib={'idref': '{}.CP{}.alpha'.format(p.p, i)})
            else:
                jointExponential = ET.SubElement(jointPrior, 'exponentialPrior',
                                                 attrib={'mean': '0.5', 'offset': '0.0'})
                ET.SubElement(jointExponential, 'parameter', attrib={'idref': '{}.alpha'.format(p.p)})
        if 'I' in p.subs_model:
            if p.every3:
                for i in range(1, 4):
                    jointUniform = ET.SubElement(jointPrior, 'uniformPrior', attrib={'lower': '0.0', 'upper': '1.0'})
                    ET.SubElement(jointUniform, 'parameter', attrib={'idref': '{}.CP{}.pInv'.format(p.p, i)})
            else:
                jointUniform = ET.SubElement(jointPrior, 'uniformPrior', attrib={'lower': '0.0', 'upper': '1.0'})
                ET.SubElement(jointUniform, 'parameter', attrib={'idref': '{}.pInv'.format(p.p)})


def clock_rate(jointPrior, partition_file):
    # Clock Model Priors
    for p, exclude, every3, susb in f.get_partition_info(partition_file):
        p = f.Partition(p, exclude, every3, susb)

        if p.every3:
            dirichletPrior = ET.SubElement(jointPrior, 'dirichletPrior', attrib={'alpha': '1.0', 'sumsTo': '1.0'})
            ET.SubElement(dirichletPrior, 'parameter', attrib={'idref': '{}.allNus'.format(p.p)})
        ctmcScalePrior = ET.SubElement(jointPrior, 'ctmcScalePrior')
        ctmcScale = ET.SubElement(ctmcScalePrior, 'ctmcScale')
        ET.SubElement(ctmcScale, 'parameter', attrib={'idref': '{}.clock.rate'.format(p.p)})
        ET.SubElement(ctmcScalePrior, 'treeModel', attrib={'idref': 'treeModel'})


def treeModel(jointPrior, rootHeight_mean, rootHeight_stdev, rootHeight_offset):
    # Module 11.3 --- Tree model priors
    rootLogNormal = ET.SubElement(jointPrior, 'logNormalPrior', attrib={'mean': str(rootHeight_mean),
                                                                        'stdev': str(rootHeight_stdev),
                                                                        'offset': str(rootHeight_offset)})
    ET.SubElement(rootLogNormal, 'parameter', attrib={'idref': 'treeModel.rootHeight'})
    rootGamma = ET.SubElement(jointPrior, 'gammaPrior', attrib={'shape': '0.001', 'scale': '1000.0', 'offset': '0.0'})
    ET.SubElement(rootGamma, 'parameter', attrib={'idref': 'skygrid.precision'})


def tip(jointPrior, fasta, priors_table):
    # Tip samples priors
    for sample in f.get_ND_taxa(fasta):
        if f.get_tip_priors(priors_table)[sample][0] == 'logNormal':
            tipPriors = ET.SubElement(jointPrior, 'logNormalPrior',
                                      attrib={'mu': str(f.get_tip_priors(priors_table)[sample][1]),
                                              'sigma': str(f.get_tip_priors(priors_table)[sample][2]),
                                              'offset': str(f.get_tip_priors(priors_table)[sample][3])})
        if f.get_tip_priors(priors_table)[sample][0] == 'uniform':
            tipPriors = ET.SubElement(jointPrior, 'uniformPrior',
                                      attrib={'lower': str(f.get_tip_priors(priors_table)[sample][1]),
                                              'upper': str(f.get_tip_priors(priors_table)[sample][2]),
                                              'offset': str(f.get_tip_priors(priors_table)[sample][3])})
        ET.SubElement(tipPriors, 'parameter', attrib={'idref': 'age({})'.format(sample)})


def skygrid(jointPrior):
    # Population model priors
    ET.SubElement(jointPrior, 'gmrfSkyGridLikelihood', attrib={'idref': 'skygrid'})


def strict_clock(jointPrior, partition_file):
    # Clock rate priors
    for p in f.get_partition_name(partition_file):
        ET.SubElement(jointPrior, 'strictClockBranchRates', attrib={'idref': '{}.branchRates'.format(p)})


def likelihood(joint, MCMC, partition_file):
    # Likelihoods
    likelihood = ET.SubElement(joint, 'likelihood', attrib={'id': 'likelihood'})
    for p in f.get_partition_name(partition_file):
        ET.SubElement(likelihood, 'treeDataLikelihood', attrib={'idref': '{}.treeLikelihood'.format(p)})
    ET.SubElement(MCMC, 'operators', attrib={'idref': 'operators'})