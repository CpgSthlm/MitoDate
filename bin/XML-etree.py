import xml.etree.ElementTree as ET
import pandas as pd
from Bio import SeqIO
import math


def pretty_xml(element, level=0):
    i = '\n' + '\t' * (level + 1)
    if element:
        if (element.text is None) or element.text.isspace():
            element.text = i
        else:
            element.text = i + element.text.strip() + i
    temp = list(element)
    for subelement in temp:
        if temp.index(subelement) < (len(temp) - 1):
            subelement.tail = i
        else:
            subelement.tail = '\n' + '\t' * level
        pretty_xml(subelement, level=level + 1)


class Taxa:  # should contain date, sequence, taxon set, tip priors
    def __init__(self, name, date, direction='backwards', units='years'):
        self.name = name
        self.date = date
        self.direction = direction
        self.units = units


class Partition:  # should contain excluding, every three, clock model, substitution model
    def __init__(self, p, exclude, every3, subs_model):
        self.p = p
        self.exclude = exclude
        self.every3 = every3
        self.subs_model = subs_model


class Xml:  # should contain all the tags (the structure), population model, root height ...
    def __init__(self, population_model, root_height_lower, root_height_mean, root_height_stdev, root_height_offset='0.0', units='years', ):
        self.units = units
        self.population_model = population_model
        self.root_height_lower = root_height_lower
        self.root_height_mean = root_height_mean
        self.root_height_stdev = root_height_stdev
        self.root_height_offset = root_height_offset


def get_taxa_name(fasta):
    taxa_name_list = []
    with open(fasta, 'r') as fasta:
        for line in fasta:
            if line.startswith('>'):
                taxa_name_list.append(str(line).strip('>').strip('\n'))
    return taxa_name_list


def get_taxa_date(fasta, priors_table):
    taxa_date_list = []
    with open(fasta, 'r') as fasta:
        for line in fasta:
            if line.startswith('>'):
                taxa_name = str(line).strip('>').strip('\n')
                if taxa_name.split('_')[-1] == 'ND':
                    taxa_date_list.append(tip_date_table(priors_table)[taxa_name])
                else:
                    taxa_date_list.append(taxa_name.split('_')[-1])
    return taxa_date_list


def get_ND_taxa(fasta):
    taxa_ND_list = []
    with open(fasta, 'r') as fasta:
        for line in fasta:
            if line.startswith('>'):
                taxa_name = str(line).strip('>').strip('\n')
                if taxa_name.split('_')[-1] == 'ND':
                    taxa_ND_list.append(taxa_name)
    return taxa_ND_list


def get_taxa_sequence(fasta):
    fa_dict = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
    return fa_dict


def get_tip_priors(priors_table):
    with open(priors_table, 'r') as priors:
        priors_table = pd.read_csv(priors, header=0)
        tip_taxa = priors_table['taxa'].values.tolist()
        tip_prior = priors_table['prior'].values.tolist()
        tip_mu_lower = priors_table['mu/lower'].values.tolist()
        tip_sigma_upper = priors_table['sigma/upper'].values.tolist()
        tip_offset = priors_table['offset'].values.tolist()

        tip_priors_list = []
        for i in range(0, len(tip_taxa)):
            priors = [tip_prior[i], tip_mu_lower[i], tip_sigma_upper[i], tip_offset[i]]
            tip_priors_list.append(priors)
    tip_priors_dict = dict(zip(tip_taxa, tip_priors_list))
    return tip_priors_dict


def tip_date_table(priors_table):
    with open(priors_table, 'r') as priors:
        priors_table = pd.read_csv(priors, header=0)
        tip_taxa = priors_table['taxa'].values.tolist()
        tip_date = priors_table['date'].values.tolist()
        tip_date_dict = dict(zip(tip_taxa, tip_date))
    return tip_date_dict


def read_partition(gff):
    global seq_type, start, end, df3
    with open(gff, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                gene_future = pd.read_table(f, header=None)
                gene_future.columns = ['genome_id', 'source',
                                       'type', 'start', 'end',
                                       'uk', 'uk1', 'uk1', 'more_info']
                df1 = gene_future.loc[:, ['type', 'start']]
                df2 = gene_future.loc[:, ['type', 'end']]
                df3 = gene_future.loc[:, ['type', 'genome_id']]

                seq_type = df3['type'].values.tolist()
                start = df1['start'].values.tolist()
                end = df2['end'].values.tolist()
    return seq_type, start, end


def get_partition(taxa, partition, fasta, gff):
    read_partition(gff)
    p_start = []
    p_end = []
    for a in range(0, len(seq_type)):
        if seq_type[a] == partition:
            p_start.append(start[a])
            p_end.append(end[a])
    if len(p_end) == 1:
        s = [p_start[0] - 1]
        e = [p_end[0]]
    else:
        s = [p_start[0] - 1]
        e = [p_end[0]]
        for i in range(1, len(p_end)):
            if p_start[i] > p_end[i - 1]:
                s.append(p_start[i] - 1)
                e.append(p_end[i])
            else:
                s.pop()
                e.pop()
                s.append(p_start[i - 1] - 1)
                e.append(p_end[i])
    p_list = []
    for i in range(0, len(s)):
        s_num = s[i]
        e_num = e[i]
        p = str(get_taxa_sequence(fasta)[taxa][s_num:e_num].seq).upper()
        p_list.append(p)
    return ''.join(p_list)


def partition_exclude(taxa, partition_name, partition_exclude, fasta, gff):
    read_partition(gff)
    p_start = []
    p_end = []
    for a in range(0, len(seq_type)):
        if seq_type[a] == partition_name:
            p_start = int(start[a] - 1)
            p_end = int(end[a])
        if seq_type[a] == partition_exclude:
            s_VNTR = int(start[a])
            e_VNTR = int(end[a] - 1)
    s = [p_start, e_VNTR]
    e = [s_VNTR, p_end]
    p_list = []
    for i in range(0, len(s)):
        s_num = s[i]
        e_num = e[i]
        partition = str(get_taxa_sequence(fasta)[taxa][s_num:e_num].seq).upper()
        p_list.append(partition)
    return ''.join(p_list)


def get_partition_info(partition_file):
    global p_name, p_exclude, p_every3, susb_model
    with open(partition_file, 'rb') as p:
        partition_table = pd.read_excel(p, header=0)
        p_name = partition_table['Partition'].values.tolist()
        p_exclude = partition_table['Exclude'].values.tolist()
        p_every3 = partition_table['Every3'].values.tolist()
        susb_model = partition_table['SubstitutionModel'].values.tolist()
    return


def make_comment(tag, text, line=-1):
    comments = ET.Comment(text)
    tag.insert(line, comments)


## Input
fasta = '/Users/vanssyli/Master/SU/CPG/BEASTdata/dating_comparisons/DPVetal2021/input_fasta/one_sample/DVPetal-mitogenomes_NDexcluded_date1.fasta'
priors_table = '/Users/vanssyli/Master/SU/CPG/Scripts/beast1XMLgenerater/data/Priors_DPVetal-mitogenomes_5oldest_5youngest.csv'
gff = '/Users/vanssyli/Master/SU/CPG/Scripts/beast1XMLgenerater/data/NC_007596.2.liftoff.gff3'
partition_file = '/Users/vanssyli/Master/SU/CPG/Scripts/beast1XMLgenerater/data/partitions.xlsx'
log_name = 'DVPetal-mitogenomes_NDexcluded_date1-etree_run5'

taxon_set = [['L.'],
             ['M.'],
             ['M.', 'E.'],
             ['ND']]
ND_list = ['ND', 'P.']

partition_list = ['tRNA', 'rRNA', 'CDS', 'D_loop']

beast_setting = Xml(units='years',
                    population_model='skygrid',
                    root_height_lower=1000000.0,
                    root_height_mean=5300000.0,
                    root_height_stdev=500000.0,
                    root_height_offset=0.0)

skygrid_cutoff = (int(list(str(beast_setting.root_height_mean))[0]) + 1) * \
                 (math.pow(10, len(str(round(beast_setting.root_height_mean))) - 1))

chainLength = 10000000
log_every = 1000



## Module 0 --- Build the basic structure of the XML file
beast = ET.Element('beast')
tree = ET.ElementTree(beast)
# listA = ['taxa', 'alignment', 'patterns', 'constantSize', 'coalescentSimulator', 'treeModel', 'treeLengthStatistic',
#          'tmrcaStatistic', 'gmrfSkyGridLikelihood', 'strictClockBranchRates', 'rateStatistic', 'HKYModel', 'siteModel',
#          'statistic', 'compoundParameter', 'treeDataLikelihood', 'operators', 'mcmc', 'report']
# for item in listA:
#     item = ET.SubElement(beast, item)

make_comment(beast, ' The list of taxa to be analysed (can also include dates/ages). ')


## Module 1 ---- TAXA
taxa = ET.SubElement(beast, 'taxa')
taxa.set('id', 'taxa')
for name, date in zip(get_taxa_name(fasta), get_taxa_date(fasta, priors_table)):
    sample = Taxa(name, str(float(date)))
    taxon = ET.SubElement(taxa, 'taxon')
    taxon.set('id', sample.name)
    date = ET.SubElement(taxon, 'date')
    date.set('value', sample.date)
    date.set('direction', sample.direction)
    date.set('units', sample.units)


# Module 1.1 --- Taxon Set
for set in taxon_set:
    set_name = [x.strip('.') for x in set]
    taxa_set = ET.SubElement(beast, 'taxa')
    taxa_set.set('id', '_and_'.join(set_name))
    for sample in get_taxa_name(fasta):
        for item in set:
            if item in sample:
                taxon = ET.SubElement(taxa_set, 'taxon')
                taxon.set('idref', sample)


## Module 2 --- Sequence
get_partition_info(partition_file)
for p, exclude, every3, susb in zip(p_name, p_exclude, p_every3, susb_model):
    p = Partition(p, exclude, every3, susb)
    print('alignment{}: '.format(str(int(p_name.index(p.p)) + 1)) + str(p.p) + ', excluded: ' + str(p.exclude) + ', every3: ' + str(p.every3),
          ', substitution model: ' + str(p.subs_model))

    alignment = ET.SubElement(beast, 'alignment', attrib={'id': 'alignment{}'.format(str(int(p_name.index(p.p)) + 1)), 'dataType': 'nucleotide'})
    make_comment(beast, ' The sequence alignment (each sequence refers to a taxon above). ')
    make_comment(beast, ' ntax={} nchar={} '.format(len(get_taxa_name(fasta)), len(get_partition(get_taxa_name(fasta)[0], p.p, fasta, gff))))

    for sample in get_taxa_name(fasta):
        sequence = ET.SubElement(alignment, 'sequence')
        ET.SubElement(sequence, 'taxon').set('idref', sample)
        if not p.exclude:
            sequence.text = get_partition(sample, p.p, fasta, gff)
        else:
            sequence.text = partition_exclude(sample, p.p, p.exclude, fasta, gff)



# Module 2.1 --- Patterns
for p, exclude, every3, susb in zip(p_name, p_exclude, p_every3, susb_model):
    p = Partition(p, exclude, every3, susb)
    if p.every3:
        for i in range(1, 4):
            patterns = ET.SubElement(beast, 'patterns', attrib={'id': '{0}.CP{1}.patterns'.format(p.p, i),
                                                                'from': str(i), 'every': '3',
                                                                'strip': 'false'})
            ET.SubElement(patterns, 'alignment').set('idref', 'alignment{}'.format(str(int(p_name.index(p.p)) + 1)))
            make_comment(beast, ' The unique patterns for codon position {} '.format(i))
    else:
        patterns = ET.SubElement(beast, 'patterns', attrib={'id': '{}.patterns'.format(p.p), 'from': '1', 'strip': 'false'})
        ET.SubElement(patterns, 'alignment').set('idref', 'alignment{}'.format(str(int(p_name.index(p.p)) + 1)))
        make_comment(beast, ' The unique patterns from 1 to end ')


## Module 3 --- Initial the tree model
constantSize = ET.SubElement(beast, 'constantSize', attrib={'id': 'initialDemo', 'units': beast_setting.units})
populationSize = ET.SubElement(constantSize, 'populationSize')
popSize_parameter = ET.SubElement(populationSize, 'parameter', attrib={'id': 'initialDemo.popSize', 'value': '100.0'})
make_comment(beast, ' This is a simple constant population size coalescent model ')
make_comment(beast, ' to generate an initial tree for the chain ')

coalescentSimulator = ET.SubElement(beast, 'coalescentSimulator', attrib={'id': 'startingTree'})
ET.SubElement(coalescentSimulator, 'taxa').set('idref', 'taxa')
ET.SubElement(coalescentSimulator, 'constantSize').set('idref', 'initialDemo')
make_comment(beast, ' Generate a random starting tree under the coalescent process ')


## Module 4 --- Generate the tree model
treeModel = ET.SubElement(beast, 'treeModel', attrib={'id': 'treeModel'})
ET.SubElement(treeModel, 'coalescentTree').set('idref', 'startingTree')
make_comment(beast, ' Generate a tree model ')

rootHeight = ET.SubElement(treeModel, 'rootHeight')
ET.SubElement(rootHeight, 'parameter').set('id', 'treeModel.rootHeight')

nodeHeights = ET.SubElement(treeModel, 'nodeHeights')
nodeHeights.set('internalNodes', 'true')
ET.SubElement(nodeHeights, 'parameter').set('id', 'treeModel.internalNodeHeights')

nodeHeights = ET.SubElement(treeModel, 'nodeHeights')
nodeHeights.set('internalNodes', 'true')
nodeHeights.set('rootNode', 'true')
ET.SubElement(nodeHeights, 'parameter').set('id', 'treeModel.allInternalNodeHeights')

# Module 4.1 --- Tip date sampling
for sample in get_ND_taxa(fasta):
    leafHeight = ET.SubElement(treeModel, 'leafHeight')
    leafHeight.set('taxon', sample)
    ET.SubElement(leafHeight, 'parameter').set('id', 'age({})'.format(sample))

# Module 4.2 --- Statistic sum of branch lengths of the tree
treeLengthStatistic = ET.SubElement(beast, 'treeLengthStatistic', attrib={'id': 'treeLength'})
ET.SubElement(treeLengthStatistic, 'treeModel').set('idref', 'treeModel')
make_comment(beast, ' Statistic for sum of the branch lengths of the tree (tree length) ')

tmrcaStatistic = ET.SubElement(beast, 'tmrcaStatistic', attrib={'id': 'age(root)', 'absolute': 'true'})
ET.SubElement(tmrcaStatistic, 'treeModel').set('idref', 'treeModel')
make_comment(beast, ' Statistic for time of most recent common ancestor of tree ')

# Module 4.3 --- Statistic sum of branch lengths for taxon sets
for set in taxon_set:
    set_name = [x.strip('.') for x in set]
    tmrcaStatistic_taxon = ET.SubElement(beast, 'tmrcaStatistic',
                                         attrib={'id': 'tmrca({})'.format('_and_'.join(set_name)), 'absolute': 'false', 'includeStem': 'false'})
    mrca = ET.SubElement(tmrcaStatistic_taxon, 'mrca')
    ET.SubElement(mrca, 'taxa').set('idref', '_and_'.join(set_name))
    ET.SubElement(tmrcaStatistic_taxon, 'treeModel').set('idref', 'treeModel')
    tmrcaStatistic_taxon = ET.SubElement(beast, 'tmrcaStatistic',
                                         attrib={'id': 'age({})'.format('_and_'.join(set_name)), 'absolute': 'true', 'includeStem': 'false'})
    mrca = ET.SubElement(tmrcaStatistic_taxon, 'mrca')
    ET.SubElement(mrca, 'taxa').set('idref', '_and_'.join(set_name))
    ET.SubElement(tmrcaStatistic_taxon, 'treeModel').set('idref', 'treeModel')
make_comment(beast, ' Taxon Sets', line=-(len(taxon_set)*2))



## Module 5 --- Population Model (skygrid for now)
gmrfSkyGridLikelihood = ET.SubElement(beast, 'gmrfSkyGridLikelihood')
gmrfSkyGridLikelihood.set('id', 'skygrid')
make_comment(beast, ' Generate a gmrfSkyGridLikelihood for the Bayesian SkyGrid process ')

skygrid_popSizes = ET.SubElement(gmrfSkyGridLikelihood, 'populationSizes')
ET.SubElement(skygrid_popSizes, 'parameter', attrib={'id': 'skygrid.logPopSize', 'dimension': '50', 'value': '1.0'})

precisionParameter = ET.SubElement(gmrfSkyGridLikelihood, 'precisionParameter')
ET.SubElement(precisionParameter, 'parameter', attrib={'id': 'skygrid.precision', 'value': '0.1', 'lower': '0.0'})

numGridPoints = ET.SubElement(gmrfSkyGridLikelihood, 'numGridPoints')
ET.SubElement(numGridPoints, 'parameter', attrib={'id': 'skygrid.numGridPoints', 'value': '49.0'})

cutOff = ET.SubElement(gmrfSkyGridLikelihood, 'cutOff')
ET.SubElement(cutOff, 'parameter', attrib={'id': 'skygrid.cutOff', 'value': str(skygrid_cutoff)})

populationTree = ET.SubElement(gmrfSkyGridLikelihood, 'populationTree')
ET.SubElement(populationTree, 'treeModel', attrib={'idref': 'treeModel'})


## Module 6 --- Clock Model (Strict clock for now)
for p, exclude, every3, susb in zip(p_name, p_exclude, p_every3, susb_model):
    p = Partition(p, exclude, every3, susb)

    strictClockBranchRates = ET.SubElement(beast, 'strictClockBranchRates', attrib={'id': '{}.branchRates'.format(p.p)})
    rate = ET.SubElement(strictClockBranchRates, 'rate')
    ET.SubElement(rate, 'parameter', attrib={'id': '{}.clock.rate'.format(p.p), 'value': '1.0', 'lower': '0.0'})

    rateStatistic = ET.SubElement(beast, 'rateStatistic', attrib={'id': '{}.meanRate'.format(p.p), 'name': '{}.meanRate'.format(p.p), 'mode': 'mean', 'internal': 'true', 'external': 'true'})
    ET.SubElement(rateStatistic, 'treeModel', attrib={'idref': 'treeModel'})
    ET.SubElement(rateStatistic, 'strictClockBranchRates', attrib={'idref': '{}.branchRates'.format(p.p)})

    make_comment(beast, ' The strict clock (Uniform rates across branches) ', line=(-2))



## Module 7 --- Substitution Model (HKY+G+I and HKY+I for now)
for p, exclude, every3, susb in zip(p_name, p_exclude, p_every3, susb_model):
    p = Partition(p, exclude, every3, susb)

    if p.every3:
        for i in range(1, 4):
            HKYModel = ET.SubElement(beast, 'HKYModel', attrib={'id': '{0}.CP{1}.hky'.format(p.p, i)})
            frequencies = ET.SubElement(HKYModel, 'frequencies')
            frequencyModel = ET.SubElement(frequencies, 'frequencyModel', attrib={'dataType': 'nucleotide'})
            frequencies_parameter = ET.SubElement(frequencyModel, 'frequencies')
            ET.SubElement(frequencies_parameter, 'parameter',
                          attrib={'id': '{0}.CP{1}.frequencies'.format(p.p, i), 'value': '0.25 0.25 0.25 0.25'})
            kappa = ET.SubElement(HKYModel, 'kappa')
            ET.SubElement(kappa, 'parameter', attrib={'id': '{0}.CP{1}.kappa'.format(p.p, i), 'value': '2.0', 'lower': '0.0'})
        for i in range(1, 4):
            siteModel = ET.SubElement(beast, 'siteModel', attrib={'id': '{0}.CP{1}.siteModel'.format(p.p, i)})
            substitutionModel = ET.SubElement(siteModel, 'substitutionModel')
            ET.SubElement(substitutionModel, 'HKYModel', attrib={'idref': '{0}.CP{1}.hky'.format(p.p, i)})
            relativeRate = ET.SubElement(siteModel, 'relativeRate', attrib={'weight': '3.0'})
            ET.SubElement(relativeRate, 'parameter', attrib={'id': '{}.CP{}.nu'.format(p.p, i), 'value': '0.3333333333333333', 'lower': '0.0', 'upper': '1.0'})
            if 'G' in p.subs_model:
                gammaShape = ET.SubElement(siteModel, 'gammaShape', attrib={'gammaCategories': '4'})
                ET.SubElement(gammaShape, 'parameter',
                              attrib={'id': '{0}.CP{1}.alpha'.format(p.p, i), 'value': '0.5', 'lower': '0.0'})
            if 'I' in p.subs_model:
                proportionInvariant = ET.SubElement(siteModel, 'proportionInvariant')
                ET.SubElement(proportionInvariant, 'parameter',
                              attrib={'id': '{}.CP{}.pInv'.format(p.p, i), 'value': '0.5', 'lower': '0.0', 'upper': '1.0'})
            subs_statistic = ET.SubElement(beast, 'statistic', attrib={'id': '{}.CP{}.mu'.format(p.p, i), 'name': 'mu'})
            ET.SubElement(subs_statistic, 'siteModel', attrib={'idref': '{}.CP{}.siteModel'.format(p.p, i)})
        make_comment(beast, ' The HKY substitution model (Hasegawa, Kishino & Yano, 1985) ', -9)

    else:
        HKYModel = ET.SubElement(beast, 'HKYModel', attrib={'id': '{}.hky'.format(p.p)})
        frequencies = ET.SubElement(HKYModel, 'frequencies')
        frequencyModel = ET.SubElement(frequencies, 'frequencyModel', attrib={'dataType': 'nucleotide'})
        frequencies_parameter = ET.SubElement(frequencyModel, 'frequencies')
        ET.SubElement(frequencies_parameter, 'parameter', attrib={'id': '{}.frequencies'.format(p.p), 'value': '0.25 0.25 0.25 0.25'})
        kappa = ET.SubElement(HKYModel, 'kappa')
        ET.SubElement(kappa, 'parameter', attrib={'id': '{}.kappa'.format(p.p), 'value': '2.0', 'lower': '0.0'})

        siteModel = ET.SubElement(beast, 'siteModel', attrib={'id': '{}.siteModel'.format(p.p)})
        substitutionModel = ET.SubElement(siteModel, 'substitutionModel')
        ET.SubElement(substitutionModel, 'HKYModel', attrib={'idref': '{}.hky'.format(p.p)})
        if 'G' in p.subs_model:
            gammaShape = ET.SubElement(siteModel, 'gammaShape', attrib={'gammaCategories': '4'})
            ET.SubElement(gammaShape, 'parameter', attrib={'id': '{}.alpha'.format(p.p), 'value': '0.5', 'lower': '0.0'})
        if 'I' in p.subs_model:
            proportionInvariant = ET.SubElement(siteModel, 'proportionInvariant')
            ET.SubElement(proportionInvariant, 'parameter', attrib={'id': '{}.pInv'.format(p.p), 'value': '0.5', 'lower': '0.0', 'upper': '1.0'})

        subs_statistic = ET.SubElement(beast, 'statistic', attrib={'id': '{}.mu'.format(p.p), 'name': 'mu'})
        ET.SubElement(subs_statistic, 'siteModel', attrib={'idref': '{}.siteModel'.format(p.p)})

        make_comment(beast, ' The HKY substitution model (Hasegawa, Kishino & Yano, 1985) ', -3)

for p, exclude, every3, susb in zip(p_name, p_exclude, p_every3, susb_model):
    p = Partition(p, exclude, every3, susb)
    if p.every3:
        compoundParameter = ET.SubElement(beast, 'compoundParameter', attrib={'id': '{}.allNus'.format(p.p)})
        for i in range(1, 4):
            ET.SubElement(compoundParameter, 'parameter', attrib={'idref': '{}.CP{}.nu'.format(p.p, i)})


## Module 8 --- Likelihood for tree given sequence data
for p, exclude, every3, susb in zip(p_name, p_exclude, p_every3, susb_model):
    p = Partition(p, exclude, every3, susb)

    treeDataLikelihood = ET.SubElement(beast, 'treeDataLikelihood', attrib={'id': '{}.treeLikelihood'.format(p.p), 'useAmbiguities': 'false'})
    make_comment(beast, ' Likelihood for tree given sequence data ')

    if p.every3:
        for i in range(1, 4):
            treeLikelihood_partition = ET.SubElement(treeDataLikelihood, 'partition')
            ET.SubElement(treeLikelihood_partition, 'patterns', attrib={'idref': '{}.CP{}.patterns'.format(p.p, i)})
            ET.SubElement(treeLikelihood_partition, 'siteModel', attrib={'idref': '{}.CP{}.siteModel'.format(p.p, i)})
    else:
        treeLikelihood_partition = ET.SubElement(treeDataLikelihood, 'partition')
        ET.SubElement(treeLikelihood_partition, 'patterns', attrib={'idref': '{}.patterns'.format(p.p)})
        ET.SubElement(treeLikelihood_partition, 'siteModel', attrib={'idref': '{}.siteModel'.format(p.p)})
    ET.SubElement(treeDataLikelihood, 'treeModel', attrib={'idref': 'treeModel'})
    ET.SubElement(treeDataLikelihood, 'strictClockBranchRates', attrib={'idref': '{}.branchRates'.format(p.p)})


## Module 9 --- Define operators
operators = ET.SubElement(beast, 'operators', attrib={'id': 'operators', 'optimizationSchedule': 'log'})
make_comment(beast, ' Define operators ')

# Modeule 9.1 --- Substitution Model Operators
for p, exclude, every3, susb in zip(p_name, p_exclude, p_every3, susb_model):
    p = Partition(p, exclude, every3, susb)

    if 'HKY' in p.subs_model:
        if p.every3:
            for i in range(1, 4):
                scaleOperator = ET.SubElement(operators, 'scaleOperator', attrib={'scaleFactor': '0.75', 'weight': '1'})
                ET.SubElement(scaleOperator, 'parameter', attrib={'idref': '{}.CP{}.kappa'.format(p.p, i)})
            for i in range(1, 4):
                deltaExchange = ET.SubElement(operators, 'deltaExchange', attrib={'delta': '0.01', 'weight': '1'})
                ET.SubElement(deltaExchange, 'parameter', attrib={'idref': '{}.CP{}.frequencies'.format(p.p, i)})
        else:
            scaleOperator = ET.SubElement(operators, 'scaleOperator', attrib={'scaleFactor': '0.75', 'weight': '1'})
            ET.SubElement(scaleOperator, 'parameter', attrib={'idref': '{}.kappa'.format(p.p)})
            deltaExchange = ET.SubElement(operators, 'deltaExchange', attrib={'delta': '0.01', 'weight': '1'})
            ET.SubElement(deltaExchange, 'parameter', attrib={'idref': '{}.frequencies'.format(p.p)})

    if 'G' in p.subs_model:
        if p.every3:
            for i in range(1, 4):
                gammaScaleOperator = ET.SubElement(operators, 'scaleOperator', attrib={'scaleFactor': '0.75', 'weight': '1'})
                ET.SubElement(gammaScaleOperator, 'parameter', attrib={'idref': '{}.CP{}.alpha'.format(p.p, i)})
        else:
            gammaScaleOperator = ET.SubElement(operators, 'scaleOperator', attrib={'scaleFactor': '0.75', 'weight': '1'})
            ET.SubElement(gammaScaleOperator, 'parameter', attrib={'idref': '{}.alpha'.format(p.p)})

    if 'I' in p.subs_model:
        if p.every3:
            for i in range(1, 4):
                randomWalkOperator = ET.SubElement(operators, 'randomWalkOperator', attrib={'windowSize': '0.75', 'weight': '1', 'boundaryCondition': 'logit'})
                ET.SubElement(randomWalkOperator, 'parameter', attrib={'idref': '{}.CP{}.pInv'.format(p.p, i)})
        else:
            randomWalkOperator = ET.SubElement(operators, 'randomWalkOperator', attrib={'windowSize': '0.75', 'weight': '1', 'boundaryCondition': 'logit'})
            ET.SubElement(randomWalkOperator, 'parameter', attrib={'idref': '{}.pInv'.format(p.p)})

# Module 10.2 --- Define clock model operators
for p, exclude, every3, susb in zip(p_name, p_exclude, p_every3, susb_model):
    p = Partition(p, exclude, every3, susb)

    clockScaleOperator = ET.SubElement(operators, 'scaleOperator', attrib={'scaleFactor': '0.75', 'weight': '3'})
    ET.SubElement(clockScaleOperator, 'parameter', attrib={'idref': '{}.clock.rate'.format(p.p)})
    upDownOperator = ET.SubElement(operators, 'upDownOperator', attrib={'scaleFactor': '0.75', 'weight': '3'})
    up = ET.SubElement(upDownOperator, 'up')
    ET.SubElement(up, 'parameter', attrib={'idref': 'treeModel.allInternalNodeHeights'})
    down = ET.SubElement(upDownOperator, 'down')
    ET.SubElement(down, 'parameter', attrib={'idref': '{}.clock.rate'.format(p.p)})
    if p.every3:
        clockDeltaExchange = ET.SubElement(operators, 'deltaExchange', attrib={'delta': '0.01', 'weight': '3'})
        ET.SubElement(clockDeltaExchange, 'parameter', attrib={'idref': '{}.allNus'.format(p.p)})

# Module 10.3 --- Define tree model operators
subtreeSlide = ET.SubElement(operators, 'subtreeSlide', attrib={'size': '1.0', 'gaussian': 'true', 'weight': '30'})
ET.SubElement(subtreeSlide, 'treeModel', attrib={'idref': 'treeModel'})
narrowExchange = ET.SubElement(operators, 'narrowExchange', attrib={'weight': '30'})
ET.SubElement(narrowExchange, 'treeModel', attrib={'idref': 'treeModel'})
wideExchange = ET.SubElement(operators, 'wideExchange', attrib={'weight': '3'})
ET.SubElement(wideExchange, 'treeModel', attrib={'idref': 'treeModel'})
wilsonBalding = ET.SubElement(operators, 'wilsonBalding', attrib={'weight': '3'})
ET.SubElement(wilsonBalding, 'treeModel', attrib={'idref': 'treeModel'})
treeScaleOperator = ET.SubElement(operators, 'scaleOperator', attrib={'scaleFactor': '0.75', 'weight': '3'})
ET.SubElement(treeScaleOperator, 'parameter', attrib={'idref': 'treeModel.rootHeight'})
treeUniformOperator = ET.SubElement(operators, 'uniformOperator', attrib={'weight': '30'})
ET.SubElement(treeUniformOperator, 'parameter', attrib={'idref': 'treeModel.internalNodeHeights'})

# Module 10.4 --- Define population model operators (skygrid here)
gmrfGridBlockUpdateOperator = ET.SubElement(operators, 'gmrfGridBlockUpdateOperator', attrib={'scaleFactor': '1.0', 'weight': '2'})
ET.SubElement(gmrfGridBlockUpdateOperator, 'gmrfSkyrideLikelihood', attrib={'idref': 'skygrid'})
popModelScaleOperator = ET.SubElement(operators, 'scaleOperator', attrib={'scaleFactor': '0.75', 'weight': '1'})
ET.SubElement(popModelScaleOperator, 'parameter', attrib={'idref': 'skygrid.precision'})

# Module 10.5 ---- Define tipping operators
for sample in get_ND_taxa(fasta):
    tipOperator = ET.SubElement(operators, 'uniformOperator', attrib={'weight': '2'})
    ET.SubElement(tipOperator, 'parameter', attrib={'idref': 'age({})'.format(sample)})


## Module 11 --- Define MCMC
MCMC = ET.SubElement(beast, 'mcmc', attrib={'id': 'mcmc', 'chainLength': str(chainLength), 'autoOptimize': 'true'})
joint = ET.SubElement(MCMC, 'joint', attrib={'id': 'joint'})
jointPrior = ET.SubElement(joint, 'prior', attrib={'id': 'prior'})

make_comment(beast,  ' Define MCMC ')

# Module 11.1 --- Substitution Model Priors
for p, exclude, every3, susb in zip(p_name, p_exclude, p_every3, susb_model):
    p = Partition(p, exclude, every3, susb)

    if 'HKY' in p.subs_model:
        if p.every3:
            for i in range(1, 4):
                jointLogNormal = ET.SubElement(jointPrior, 'logNormalPrior', attrib={'mu': '1.0', 'sigma': '1.25', 'offset': '0.0'})
                ET.SubElement(jointLogNormal, 'parameter', attrib={'idref': '{}.CP{}.kappa'.format(p.p, i)})
            for i in range(1, 4):
                jointDirichlet = ET.SubElement(jointPrior, 'dirichletPrior', attrib={'alpha': '1.0', 'sumsTo': '1.0'})
                ET.SubElement(jointDirichlet, 'parameter', attrib={'idref': '{}.CP{}.frequencies'.format(p.p, i)})
        else:
            jointLogNormal = ET.SubElement(jointPrior, 'logNormalPrior', attrib={'mu': '1.0', 'sigma': '1.25', 'offset': '0.0'})
            ET.SubElement(jointLogNormal, 'parameter', attrib={'idref': '{}.kappa'.format(p.p)})
            jointDirichlet = ET.SubElement(jointPrior, 'dirichletPrior', attrib={'alpha': '1.0', 'sumsTo': '1.0'})
            ET.SubElement(jointDirichlet, 'parameter', attrib={'idref': '{}.frequencies'.format(p.p)})
    if 'G' in p.subs_model:
        if p.every3:
            for i in range(1, 4):
                jointExponential = ET.SubElement(jointPrior, 'exponentialPrior', attrib={'mean': '0.5', 'offset': '0.0'})
                ET.SubElement(jointExponential, 'parameter', attrib={'idref': '{}.CP{}.alpha'.format(p.p, i)})
        else:
            jointExponential = ET.SubElement(jointPrior, 'exponentialPrior', attrib={'mean': '0.5', 'offset': '0.0'})
            ET.SubElement(jointExponential, 'parameter', attrib={'idref': '{}.alpha'.format(p.p)})
    if 'I' in p.subs_model:
        if p.every3:
            for i in range(1, 4):
                jointUniform = ET.SubElement(jointPrior, 'uniformPrior', attrib={'lower': '0.0', 'upper': '1.0'})
                ET.SubElement(jointUniform, 'parameter', attrib={'idref': '{}.CP{}.pInv'.format(p.p, i)})
        else:
            jointUniform = ET.SubElement(jointPrior, 'uniformPrior', attrib={'lower': '0.0', 'upper': '1.0'})
            ET.SubElement(jointUniform, 'parameter', attrib={'idref': '{}.pInv'.format(p.p)})

# Module 11.2 --- Clock Model Priors
for p, exclude, every3, susb in zip(p_name, p_exclude, p_every3, susb_model):
    p = Partition(p, exclude, every3, susb)

    if p.every3:
        dirichletPrior = ET.SubElement(jointPrior, 'dirichletPrior', attrib={'alpha': '1.0', 'sumsTo': '1.0'})
        ET.SubElement(dirichletPrior, 'parameter', attrib={'idref': '{}.allNus'.format(p.p)})
    ctmcScalePrior = ET.SubElement(jointPrior, 'ctmcScalePrior')
    ctmcScale = ET.SubElement(ctmcScalePrior, 'ctmcScale')
    ET.SubElement(ctmcScale, 'parameter', attrib={'idref': '{}.clock.rate'.format(p.p)})
    ET.SubElement(ctmcScalePrior, 'treeModel', attrib={'idref': 'treeModel'})

# Module 11.3 --- Tree model priors
rootLogNormal = ET.SubElement(jointPrior, 'logNormalPrior', attrib={'mean': str(beast_setting.root_height_mean),
                                                                   'stdev': str(beast_setting.root_height_stdev),
                                                                   'offset': str(beast_setting.root_height_offset)})
ET.SubElement(rootLogNormal, 'parameter', attrib={'idref': 'treeModel.rootHeight'})
rootGamma = ET.SubElement(jointPrior, 'gammaPrior', attrib={'shape': '0.001', 'scale': '1000.0', 'offset': '0.0'})
ET.SubElement(rootGamma, 'parameter', attrib={'idref': 'skygrid.precision'})

# Module 11.4 --- Tip samples priors
for sample in get_ND_taxa(fasta):
    if get_tip_priors(priors_table)[sample][0] == 'logNormal':
        tipPriors = ET.SubElement(jointPrior, 'logNormalPrior', attrib={'mu': str(get_tip_priors(priors_table)[sample][1]),
                                                                        'sigma': str(get_tip_priors(priors_table)[sample][2]),
                                                                        'offset': str(get_tip_priors(priors_table)[sample][3])})
    if get_tip_priors(priors_table)[sample][0] == 'uniform':
        tipPriors = ET.SubElement(jointPrior, 'uniformPrior', attrib={'lower': str(get_tip_priors(priors_table)[sample][1]),
                                                                      'upper': str(get_tip_priors(priors_table)[sample][2]),
                                                                      'offset': str(get_tip_priors(priors_table)[sample][3])})
    ET.SubElement(tipPriors, 'parameter', attrib={'idref': 'age({})'.format(sample)})

# Module 11.5 --- Population model priors and clock rate priors
skygridPrior = ET.SubElement(jointPrior, 'gmrfSkyGridLikelihood', attrib={'idref': 'skygrid'})
for p in p_name:
    ET.SubElement(jointPrior, 'strictClockBranchRates', attrib={'idref': '{}.branchRates'.format(p)})

# Module 11.6 --- Likelihoods
likelihood = ET.SubElement(joint, 'likelihood', attrib={'id': 'likelihood'})
for p in p_name:
    ET.SubElement(likelihood, 'treeDataLikelihood', attrib={'idref': '{}.treeLikelihood'.format(p)})

ET.SubElement(MCMC, 'operators', attrib={'idref': 'operators'})


## Module 12 --- Write log
# Module 12.1 --- Write log to screen
log = ET.SubElement(MCMC, 'log', attrib={'id': 'screenLog', 'logEvery': str(log_every)})
jointCol = ET.SubElement(log, 'column', attrib={'label': 'Joint', 'dp': '4', 'width': '12'})
ET.SubElement(jointCol, 'joint', attrib={'idref': 'joint'})
priorCol = ET.SubElement(log, 'column', attrib={'label': 'Prior', 'dp': '4', 'width': '12'})
ET.SubElement(priorCol, 'prior', attrib={'idref': 'prior'})
likelihoodCol = ET.SubElement(log, 'column', attrib={'label': 'Likelihood', 'dp': '4', 'width': '12'})
ET.SubElement(likelihoodCol, 'likelihood', attrib={'idref': 'likelihood'})
ageCol = ET.SubElement(log, 'column', attrib={'label': 'age(root)', 'sf': '6', 'width': '12'})
ET.SubElement(ageCol, 'tmrcaStatistic', attrib={'idref': 'age(root)'})
for p in p_name:
    pCol = ET.SubElement(log, 'column', attrib={'label': '{}.clock.rate'.format(p), 'sf': '6', 'width': '12'})
    ET.SubElement(pCol, 'parameter', attrib={'idref': '{}.clock.rate'.format(p)})

make_comment(MCMC, ' Write log to screen ')

# Module 12.2 --- Write log to file
logFile = ET.SubElement(MCMC, 'log', attrib={'id': 'fileLog', 'logEvery': str(log_every), 'fileName': '{}.log'.format(log_name), 'overwrite': 'false'})
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

for p, exclude, every3, susb in zip(p_name, p_exclude, p_every3, susb_model):
    p = Partition(p, exclude, every3, susb)
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

for p, exclude, every3, susb in zip(p_name, p_exclude, p_every3, susb_model):
    p = Partition(p, exclude, every3, susb)
    if p.every3:
        ET.SubElement(logFile, 'compoundParameter', attrib={'idref': 'CDS.allNus'})
        for i in range(1, 4):
            ET.SubElement(logFile, 'statistic', attrib={'idref': '{}.CP{}.mu'.format(p.p, i)})
    ET.SubElement(logFile, 'parameter', attrib={'idref': '{}.clock.rate'.format(p.p)})

for p in p_name:
    ET.SubElement(logFile, 'parameter', attrib={'idref': '{}.meanRate'.format(p)})

for sample in get_ND_taxa(fasta):
    ET.SubElement(logFile, 'parameter', attrib={'idref': 'age({})'.format(sample)})

for p in p_name:
    ET.SubElement(logFile, 'treeDataLikelihood', attrib={'idref': '{}.treeLikelihood'.format(p)})

for p in p_name:
    ET.SubElement(logFile, 'strictClockBranchRates', attrib={'idref': '{}.branchRates'.format(p)})

ET.SubElement(logFile, 'gmrfSkyGridLikelihood', attrib={'idref': 'skygrid'})

make_comment(MCMC, ' Write log to file ')


# Module 12.3 --- Write tree log to file
logTree = ET.SubElement(MCMC, 'logTree', attrib={'id': 'treeFileLog', 'logEvery': str(log_every),
                                                 'nexusFormat': 'true',
                                                 'fileName': '{}.trees'.format(log_name),
                                                 'sortTranslationTable': 'true'})
ET.SubElement(logTree, 'treeModel', attrib={'idref': 'treeModel'})
for p in p_name:
    trait = ET.SubElement(logTree, 'trait', attrib={'name': 'rate', 'tag': '{}.rate'.format(p)})
    ET.SubElement(trait, 'strictClockBranchRates', attrib={'idref': '{}.branchRates'.format(p)})
ET.SubElement(logTree, 'joint', attrib={'idref': 'joint'})

make_comment(MCMC, ' Write tree log to file ')


## Module 13 --- Report
report = ET.SubElement(beast, 'report')
property = ET.SubElement(report, 'property', attrib={'name': 'timer'})
ET.SubElement(property, 'mcmc', attrib={'idref': 'mcmc'})


## Test
pretty_xml(beast)
# ET.dump(beast)

## Store the file
tree.write('{}.xml'.format(log_name), encoding='utf=8', xml_declaration=True)