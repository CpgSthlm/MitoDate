import xml.etree.ElementTree as ET
import pandas as pd
from Bio import SeqIO
import math


class Taxa:
    """
    Class to represent a Taxon, including its name, date, direction, and units.
    """

    def __init__(self, name, date, direction='backwards', units='years'):
        self.name = name
        self.date = date
        self.direction = direction
        self.units = units


class Partition:
    """
    Class representing a Partition with attributes for the partition, exclusion,
    every three bases, and substitution model.
    """

    def __init__(self, p, exclude, every3, subs_model):
        self.p = p
        self.exclude = exclude
        self.every3 = every3
        self.subs_model = subs_model


class Xml:
    """
    Class for generating XML structure related to population model and root height information.
    """

    def __init__(self, population_model, root_height_lower, root_height_mean, root_height_stdev,
                 root_height_offset='0.0', units='years'):
        self.units = units
        self.population_model = population_model
        self.root_height_lower = root_height_lower
        self.root_height_mean = root_height_mean
        self.root_height_stdev = root_height_stdev
        self.root_height_offset = root_height_offset


def get_taxa_name(fasta):
    """
    Extracts taxa names from a FASTA file.
    """
    taxa_name_list = []
    with open(fasta, 'r') as fasta_file:
        for line in fasta_file:
            if line.startswith('>'):
                taxa_name_list.append(str(line).strip('>').strip('\n'))
    return taxa_name_list


def get_taxa_date(fasta, priors_table):
    """
    Extracts the taxa dates from a FASTA file, using a priors table for ND (No Date) taxa.
    """
    taxa_date_list = []
    with open(fasta, 'r') as fasta_file:
        for line in fasta_file:
            if line.startswith('>'):
                taxa_name = str(line).strip('>').strip('\n')
                if taxa_name.split('_')[-1] == 'ND':  # If ND, use the date from priors_table
                    taxa_date_list.append(tip_date_table(priors_table, fasta)[taxa_name])
                else:
                    taxa_date_list.append(taxa_name.split('_')[-1])  # Extract date from the name
    return taxa_date_list


def get_ND_taxa(fasta):
    """
    Extracts taxa names with 'ND' suffix from a FASTA file.
    """
    taxa_ND_list = []
    with open(fasta, 'r') as fasta_file:
        for line in fasta_file:
            if line.startswith('>'):
                taxa_name = str(line).strip('>').strip('\n')
                if taxa_name.split('_')[-1] == 'ND':
                    taxa_ND_list.append(taxa_name)
    return taxa_ND_list


def get_taxa_sequence(fasta):
    """
    Parses the sequences from a FASTA file and returns them in a dictionary format.
    """
    fa_dict = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
    return fa_dict


def get_tip_priors(priors_table, fasta):
    """
    Reads a priors table and returns a dictionary of tip priors for each taxa.
    """
    with open(priors_table, 'r') as priors_file:
        priors_df = pd.read_csv(priors_file, header=0)
        short_taxa = priors_df['taxa'].values.tolist()
        fasta_taxa = get_taxa_name(fasta)
        
        tip_taxa = []
        matched_indices = []
        
        for i, taxa in enumerate(short_taxa):
            matching_fasta_names = [name for name in fasta_taxa if name.startswith(taxa)]
            if matching_fasta_names:
                tip_taxa.append(matching_fasta_names[0])
                matched_indices.append(i)
    
        tip_prior = priors_df['prior'].values.tolist()
        tip_mu_lower = priors_df['mu/lower'].values.tolist()
        tip_sigma_upper = priors_df['sigma/upper'].values.tolist()
        tip_offset = priors_df['offset'].values.tolist()

        tip_priors_list = [
            [tip_prior[i], tip_mu_lower[i], tip_sigma_upper[i], tip_offset[i]]
            for i in range(len(tip_taxa))
        ]

    return dict(zip(tip_taxa, tip_priors_list))


def tip_date_table(priors_table, fasta):
    """
    Reads a priors table and returns a dictionary of taxa names to their dates.
    """
    with open(priors_table, 'r') as priors_file:
        priors_df = pd.read_csv(priors_file, header=0)
        short_taxa = priors_df['taxa'].values.tolist()
        fasta_taxa = get_taxa_name(fasta)
        
        tip_taxa = []
        matched_indices = []
        
        for i, taxa in enumerate(short_taxa):
            matching_fasta_names = [name for name in fasta_taxa if name.startswith(taxa)]
            if matching_fasta_names:
                tip_taxa.append(matching_fasta_names[0])
                matched_indices.append(i)
        
        tip_date = priors_df['date'].values.tolist()

    return dict(zip(tip_taxa, tip_date))


def get_partition(taxa, partition, fasta, gff):
    """
    Returns the sequence for a specified partition from a FASTA file, excluding others.
    """
    # Read and parse the GFF file into a DataFrame
    gff_df = pd.read_table(gff, comment='#', header=None,
                           names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])

    # Filter the GFF DataFrame to get rows that match the partition type
    partition_df = gff_df[gff_df['type'] == partition]

    # Extract start and end positions for the matching partition
    start_positions = partition_df['start'].values - 1  # Convert to 0-based indexing
    end_positions = partition_df['end'].values

    # Extract the corresponding sequences from the FASTA file
    fa_dict = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
    sequences = [str(fa_dict[taxa][start:end].seq).upper() for start, end in zip(start_positions, end_positions)]

    return ''.join(sequences)


def partition_exclude(taxa, partition_name, partition_exclude, fasta, gff):
    """
    Returns the sequence for a partition excluding certain regions specified by partition_exclude.
    """
    # Read and parse the GFF file into a DataFrame
    gff_df = pd.read_table(gff, comment='#', header=None,
                           names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])

    # Extract the start and end positions for the specified partition
    partition_df = gff_df[gff_df['type'] == partition_name]
    exclude_df = gff_df[gff_df['type'] == partition_exclude]

    # Get the partition start and end positions
    p_start = partition_df['start'].values - 1  # 0-based indexing
    p_end = partition_df['end'].values

    # Get the exclusion region start and end positions
    exclude_start = exclude_df['start'].values
    exclude_end = exclude_df['end'].values - 1  # 0-based indexing

    # Combine the partition and exclusion regions into separate lists
    s = [p_start[0], exclude_end[0]]
    e = [exclude_start[0], p_end[0]]

    # Extract the corresponding sequences from the FASTA file
    fa_dict = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
    sequences = [str(fa_dict[taxa][start:end].seq).upper() for start, end in zip(s, e)]

    return ''.join(sequences)


def get_partition_info(partition_file):
    """
    Reads the partition file and returns the partition information as a tuple.
    """
    with open(partition_file, 'rb') as p:
        partition_df = pd.read_excel(p, header=0)
        p_name = partition_df['Partition'].values.tolist()
        p_exclude = partition_df['Exclude'].values.tolist()
        p_every3 = partition_df['Every3'].values.tolist()
        susb_model = partition_df['SubstitutionModel'].values.tolist()

    return zip(p_name, p_exclude, p_every3, susb_model)


def get_partition_name(partition_file):
    """
    Reads the partition file and returns a list of partition names.
    """
    with open(partition_file, 'rb') as p:
        partition_df = pd.read_excel(p, header=0)
        p_name = partition_df['Partition'].values.tolist()

    return p_name


def skygrid_cutoff(a):
    """
    Returns the cutoff value for skygrid computation based on the input number.
    """
    x = (int(list(str(a))[0]) + 1) * (math.pow(10, len(str(round(a))) - 1))
    return x


def make_comment(tag, text, line=-1):
    """
    Adds a comment to an XML element at the specified line number.
    """
    comments = ET.Comment(text)
    tag.insert(line, comments)


def pretty_xml(element, level=0):
    """
    Formats XML for pretty printing.
    """
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
