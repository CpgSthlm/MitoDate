import xml.etree.ElementTree as ET
import pandas as pd
from Bio import SeqIO
import math


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
    with open(partition_file, 'rb') as p:
        partition_table = pd.read_excel(p, header=0)
        p_name = partition_table['Partition'].values.tolist()
        p_exclude = partition_table['Exclude'].values.tolist()
        p_every3 = partition_table['Every3'].values.tolist()
        susb_model = partition_table['SubstitutionModel'].values.tolist()
    return zip(p_name, p_exclude, p_every3, susb_model)


def get_partition_name(partition_file):
    with open(partition_file, 'rb') as p:
        partition_table = pd.read_excel(p, header=0)
        p_name = partition_table['Partition'].values.tolist()
    return p_name


def skygrid_cutoff(a):
    x = (int(list(str(a))[0]) + 1) * (math.pow(10, len(str(round(a))) - 1))
    return x


def make_comment(tag, text, line=-1):
    comments = ET.Comment(text)
    tag.insert(line, comments)


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