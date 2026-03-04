import xml.etree.ElementTree as ET
import functions as f
import pandas as pd
from Bio import SeqIO
import math


class Taxa:  # should contain date, sequence, taxon set, tip priors
    def __init__(self, name, date, direction='backwards', units='years'):
        self.name = name
        self.date = date
        self.direction = direction
        self.units = units


def intro_taxa(beast, fasta, priors_table):
    taxa = ET.SubElement(beast, 'taxa')
    taxa.set('id', 'taxa')
    for name, date in zip(f.get_taxa_name(fasta), f.get_taxa_date(fasta, priors_table)):
        sample = Taxa(name, str(float(date)))
        taxon = ET.SubElement(taxa, 'taxon')
        taxon.set('id', sample.name)
        date = ET.SubElement(taxon, 'date')
        date.set('value', sample.date)
        date.set('direction', sample.direction)
        date.set('units', sample.units)
    return


def intro_taxonSet(beast, taxon_set, fasta):
    for set in taxon_set:
        set_name = [x.strip('.') for x in set]
        taxa_set = ET.SubElement(beast, 'taxa')
        taxa_set.set('id', '_and_'.join(set_name))
        for sample in f.get_taxa_name(fasta):
            for item in set:
                if item in sample:
                    taxon = ET.SubElement(taxa_set, 'taxon')
                    taxon.set('idref', sample)