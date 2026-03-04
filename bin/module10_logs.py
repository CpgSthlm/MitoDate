import xml.etree.ElementTree as ET
import functions as f


def to_screen(MCMC, partition_file, log_every, split_partition):
    log = ET.SubElement(MCMC, 'log', attrib={'id': 'screenLog', 'logEvery': str(log_every)})
    jointCol = ET.SubElement(log, 'column', attrib={'label': 'Joint', 'dp': '4', 'width': '12'})
    ET.SubElement(jointCol, 'joint', attrib={'idref': 'joint'})
    priorCol = ET.SubElement(log, 'column', attrib={'label': 'Prior', 'dp': '4', 'width': '12'})
    ET.SubElement(priorCol, 'prior', attrib={'idref': 'prior'})
    likelihoodCol = ET.SubElement(log, 'column', attrib={'label': 'Likelihood', 'dp': '4', 'width': '12'})
    ET.SubElement(likelihoodCol, 'likelihood', attrib={'idref': 'likelihood'})
    ageCol = ET.SubElement(log, 'column', attrib={'label': 'age(root)', 'sf': '6', 'width': '12'})
    ET.SubElement(ageCol, 'tmrcaStatistic', attrib={'idref': 'age(root)'})

    if split_partition:
        for p in f.get_partition_name(partition_file):
            pCol = ET.SubElement(log, 'column', attrib={'label': '{}.clock.rate'.format(p), 'sf': '6', 'width': '12'})
            ET.SubElement(pCol, 'parameter', attrib={'idref': '{}.clock.rate'.format(p)})
    else:
        pCol = ET.SubElement(log, 'column', attrib={'label': 'clock.rate', 'sf': '6', 'width': '12'})
        ET.SubElement(pCol, 'parameter', attrib={'idref': 'clock.rate'})

    f.make_comment(MCMC, ' Write log to screen ')


def to_file(MCMC, log_every, log_name, taxon_set, partition_file, fasta, split_partition):
    logFile = ET.SubElement(MCMC, 'log',
                            attrib={'id': 'fileLog', 'logEvery': str(log_every), 'fileName': '{}.log'.format(log_name),
                                    'overwrite': 'false'})
    ET.SubElement(logFile, 'joint', attrib={'idref': 'joint'})
    ET.SubElement(logFile, 'prior', attrib={'idref': 'prior'})
    ET.SubElement(logFile, 'likelihood', attrib={'idref': 'likelihood'})
    ET.SubElement(logFile, 'parameter', attrib={'idref': 'treeModel.rootHeight'})
    ET.SubElement(logFile, 'tmrcaStatistic', attrib={'idref': 'age(root)'})
    ET.SubElement(logFile, 'treeLengthStatistic', attrib={'idref': 'treeLength'})
    for set in taxon_set:
        set_name = [x.strip('.') for x in set]
        ET.SubElement(logFile, 'tmrcaStatistic', attrib={'idref': 'tmrca({})'.format('_and_'.join(set_name))})
    for set in taxon_set:
        set_name = [x.strip('.') for x in set]
        ET.SubElement(logFile, 'tmrcaStatistic', attrib={'idref': 'age({})'.format('_and_'.join(set_name))})
    ET.SubElement(logFile, 'parameter', attrib={'idref': 'skygrid.precision'})
    ET.SubElement(logFile, 'parameter', attrib={'idref': 'skygrid.logPopSize'})
    ET.SubElement(logFile, 'parameter', attrib={'idref': 'skygrid.cutOff'})

    if split_partition:
        for p, exclude, every3, susb in f.get_partition_info(partition_file):
            p = f.Partition(p, exclude, every3, susb)
            if 'HKY' in p.subs_model:
                if p.every3:
                    for i in range(1, 4):
                        ET.SubElement(logFile, 'parameter', attrib={'idref': '{}.CP{}.kappa'.format(p.p, i)})
                        ET.SubElement(logFile, 'parameter', attrib={'idref': '{}.CP{}.frequencies'.format(p.p, i)})
                else:
                    ET.SubElement(logFile, 'parameter', attrib={'idref': '{}.kappa'.format(p.p)})
                    ET.SubElement(logFile, 'parameter', attrib={'idref': '{}.frequencies'.format(p.p)})
            if 'G' in p.subs_model:
                if p.every3:
                    for i in range(1, 4):
                        ET.SubElement(logFile, 'parameter', attrib={'idref': '{}.CP{}.alpha'.format(p.p, i)})
                else:
                    ET.SubElement(logFile, 'parameter', attrib={'idref': '{}.alpha'.format(p.p)})
            if 'I' in p.subs_model:
                if p.every3:
                    for i in range(1, 4):
                        ET.SubElement(logFile, 'parameter', attrib={'idref': '{}.CP{}.pInv'.format(p.p, i)})
                else:
                    ET.SubElement(logFile, 'parameter', attrib={'idref': '{}.pInv'.format(p.p)})
    else:
        for p, exclude, every3, susb in f.get_partition_info(partition_file):
            p = f.Partition(p, exclude, every3, susb)
        if 'HKY' in p.subs_model:
            ET.SubElement(logFile, 'parameter', attrib={'idref': 'kappa'})
            ET.SubElement(logFile, 'parameter', attrib={'idref': 'frequencies'})
        if 'G' in p.subs_model:
            ET.SubElement(logFile, 'parameter', attrib={'idref': 'alpha'})
        if 'I' in p.subs_model:
            ET.SubElement(logFile, 'parameter', attrib={'idref': 'pInv'})

    if split_partition:
        for p, exclude, every3, susb in f.get_partition_info(partition_file):
            p = f.Partition(p, exclude, every3, susb)
            if p.every3:
                ET.SubElement(logFile, 'compoundParameter', attrib={'idref': 'CDS.allNus'})
                for i in range(1, 4):
                    ET.SubElement(logFile, 'statistic', attrib={'idref': '{}.CP{}.mu'.format(p.p, i)})
            ET.SubElement(logFile, 'parameter', attrib={'idref': '{}.clock.rate'.format(p.p)})
    else:
        ET.SubElement(logFile, 'parameter', attrib={'idref': 'clock.rate'})

    if split_partition:
        for p in f.get_partition_name(partition_file):
            ET.SubElement(logFile, 'parameter', attrib={'idref': '{}.meanRate'.format(p)})
        for sample in f.get_ND_taxa(fasta):
            ET.SubElement(logFile, 'parameter', attrib={'idref': 'age({})'.format(sample)})
        for p in f.get_partition_name(partition_file):
            ET.SubElement(logFile, 'treeDataLikelihood', attrib={'idref': '{}.treeLikelihood'.format(p)})
        for p in f.get_partition_name(partition_file):
            ET.SubElement(logFile, 'strictClockBranchRates', attrib={'idref': '{}.branchRates'.format(p)})
    else:
        ET.SubElement(logFile, 'parameter', attrib={'idref': 'meanRate'})
        for sample in f.get_ND_taxa(fasta):
            ET.SubElement(logFile, 'parameter', attrib={'idref': 'age({})'.format(sample)})
        ET.SubElement(logFile, 'treeDataLikelihood', attrib={'idref': 'treeLikelihood'})
        ET.SubElement(logFile, 'strictClockBranchRates', attrib={'idref': 'branchRates'})




    ET.SubElement(logFile, 'gmrfSkyGridLikelihood', attrib={'idref': 'skygrid'})

    f.make_comment(MCMC, ' Write log to file ')


def trees_to_file(MCMC, log_every, log_name, partition_file, split_partition):
    logTree = ET.SubElement(MCMC, 'logTree', attrib={'id': 'treeFileLog', 'logEvery': str(log_every),
                                                     'nexusFormat': 'true',
                                                     'fileName': '{}.trees'.format(log_name),
                                                     'sortTranslationTable': 'true'})
    ET.SubElement(logTree, 'treeModel', attrib={'idref': 'treeModel'})
    if split_partition:
        for p in f.get_partition_name(partition_file):
            trait = ET.SubElement(logTree, 'trait', attrib={'name': 'rate', 'tag': '{}.rate'.format(p)})
            ET.SubElement(trait, 'strictClockBranchRates', attrib={'idref': '{}.branchRates'.format(p)})
        ET.SubElement(logTree, 'joint', attrib={'idref': 'joint'})
    else:
        trait = ET.SubElement(logTree, 'trait', attrib={'name': 'rate', 'tag': 'rate'})
        ET.SubElement(trait, 'strictClockBranchRates', attrib={'idref': 'branchRates'})
        ET.SubElement(logTree, 'joint', attrib={'idref': 'joint'})


    f.make_comment(MCMC, ' Write tree log to file ')


def report(beast):
    report = ET.SubElement(beast, 'report')
    property = ET.SubElement(report, 'property', attrib={'name': 'timer'})
    ET.SubElement(property, 'mcmc', attrib={'idref': 'mcmc'})
