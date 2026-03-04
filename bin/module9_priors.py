import xml.etree.ElementTree as ET
import functions as f


def define_operators(jointPrior,joint, MCMC, fasta, priors_table, partition_file, split_partition,
                     rootHeight_mean, rootHeight_stdev, rootHeight_offset):
    subs_HKY(jointPrior, partition_file, split_partition)
    clock_rate(jointPrior, partition_file, split_partition)
    treeModel(jointPrior, rootHeight_mean, rootHeight_stdev, rootHeight_offset)
    tip(jointPrior, fasta, priors_table)
    skygrid(jointPrior)
    strict_clock(jointPrior, partition_file, split_partition)
    likelihood(joint, MCMC, partition_file, split_partition)


def subs_HKY(jointPrior, partition_file, split_partition):
    # Substitution Model Priors
    if split_partition:
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
    else:
        for p, exclude, every3, susb in f.get_partition_info(partition_file):
            p = f.Partition(p, exclude, every3, susb)

        if 'HKY' in p.subs_model:
            jointLogNormal = ET.SubElement(jointPrior, 'logNormalPrior',
                                           attrib={'mu': '1.0', 'sigma': '1.25', 'offset': '0.0'})
            ET.SubElement(jointLogNormal, 'parameter', attrib={'idref': 'kappa'})
            jointDirichlet = ET.SubElement(jointPrior, 'dirichletPrior', attrib={'alpha': '1.0', 'sumsTo': '1.0'})
            ET.SubElement(jointDirichlet, 'parameter', attrib={'idref': 'frequencies'})
        if 'G' in p.subs_model:
            jointExponential = ET.SubElement(jointPrior, 'exponentialPrior',
                                             attrib={'mean': '0.5', 'offset': '0.0'})
            ET.SubElement(jointExponential, 'parameter', attrib={'idref': 'alpha'})
        if 'I' in p.subs_model:
            jointUniform = ET.SubElement(jointPrior, 'uniformPrior', attrib={'lower': '0.0', 'upper': '1.0'})
            ET.SubElement(jointUniform, 'parameter', attrib={'idref': 'pInv'})





def clock_rate(jointPrior, partition_file, split_partition):
    # Clock Model Priors
    if split_partition:
        for p, exclude, every3, susb in f.get_partition_info(partition_file):
            p = f.Partition(p, exclude, every3, susb)

            if p.every3:
                dirichletPrior = ET.SubElement(jointPrior, 'dirichletPrior', attrib={'alpha': '1.0', 'sumsTo': '1.0'})
                ET.SubElement(dirichletPrior, 'parameter', attrib={'idref': '{}.allNus'.format(p.p)})
            ctmcScalePrior = ET.SubElement(jointPrior, 'ctmcScalePrior')
            ctmcScale = ET.SubElement(ctmcScalePrior, 'ctmcScale')
            ET.SubElement(ctmcScale, 'parameter', attrib={'idref': '{}.clock.rate'.format(p.p)})
            ET.SubElement(ctmcScalePrior, 'treeModel', attrib={'idref': 'treeModel'})
    else:
        ctmcScalePrior = ET.SubElement(jointPrior, 'ctmcScalePrior')
        ctmcScale = ET.SubElement(ctmcScalePrior, 'ctmcScale')
        ET.SubElement(ctmcScale, 'parameter', attrib={'idref': 'clock.rate'})
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
        if f.get_tip_priors(priors_table, fasta)[sample][0] == 'logNormal':
            tipPriors = ET.SubElement(jointPrior, 'logNormalPrior',
                                      attrib={'mu': str(f.get_tip_priors(priors_table, fasta)[sample][1]),
                                              'sigma': str(f.get_tip_priors(priors_table, fasta)[sample][2]),
                                              'offset': str(f.get_tip_priors(priors_table, fasta)[sample][3])})
        if f.get_tip_priors(priors_table, fasta)[sample][0] == 'uniform':
            tipPriors = ET.SubElement(jointPrior, 'uniformPrior',
                                      attrib={'lower': str(f.get_tip_priors(priors_table, fasta)[sample][1]),
                                              'upper': str(f.get_tip_priors(priors_table, fasta)[sample][2]),
                                              'offset': str(f.get_tip_priors(priors_table, fasta)[sample][3])})
        if f.get_tip_priors(priors_table, fasta)[sample][0] == 'normal':
            tipPriors = ET.SubElement(jointPrior, 'normalPrior',
                                      attrib={'mean': str(f.get_tip_priors(priors_table, fasta)[sample][1]),
                                              'stdev': str(f.get_tip_priors(priors_table, fasta)[sample][2]),
                                              'offset': str(f.get_tip_priors(priors_table, fasta)[sample][3])})
        ET.SubElement(tipPriors, 'parameter', attrib={'idref': 'age({})'.format(sample)})


def skygrid(jointPrior):
    # Population model priors
    ET.SubElement(jointPrior, 'gmrfSkyGridLikelihood', attrib={'idref': 'skygrid'})


def strict_clock(jointPrior, partition_file, split_partition):
    # Clock rate priors
    if split_partition:
        for p in f.get_partition_name(partition_file):
            ET.SubElement(jointPrior, 'strictClockBranchRates', attrib={'idref': '{}.branchRates'.format(p)})
    else:
        ET.SubElement(jointPrior, 'strictClockBranchRates', attrib={'idref': 'branchRates'})



def likelihood(joint, MCMC, partition_file, split_partition):
    # Likelihoods
    likelihood = ET.SubElement(joint, 'likelihood', attrib={'id': 'likelihood'})
    if split_partition:
        for p in f.get_partition_name(partition_file):
            ET.SubElement(likelihood, 'treeDataLikelihood', attrib={'idref': '{}.treeLikelihood'.format(p)})
    else:
        ET.SubElement(likelihood, 'treeDataLikelihood', attrib={'idref': 'treeLikelihood'})
    ET.SubElement(MCMC, 'operators', attrib={'idref': 'operators'})