import xml.etree.ElementTree as ET
import pandas as pd
from Bio import SeqIO
import math
import module1_Taxa as mTaxa
import functions as f




def build_partitionSeq(beast, partition_file, fasta, gff):
    for p, exclude, every3, susb in f.get_partition_info(partition_file):
        p = f.Partition(p, exclude, every3, susb)
        print('alignment{}: '.format(str(int(f.get_partition_name(partition_file).index(p.p))+ 1)) + str(p.p) + ', excluded: ' + str(p.exclude) + ', every3: ' + str(p.every3),
              ', substitution model: ' + str(p.subs_model))

        alignment = ET.SubElement(beast, 'alignment', attrib={'id': 'alignment{}'.format(str(int(f.get_partition_name(partition_file).index(p.p)) + 1)), 'dataType': 'nucleotide'})
        f.make_comment(beast, ' The sequence alignment (each sequence refers to a taxon above). ')
        f.make_comment(beast, ' ntax={} nchar={} '.format(len(f.get_taxa_name(fasta)), len(f.get_partition(f.get_taxa_name(fasta)[0], p.p, fasta, gff))))

        for sample in f.get_taxa_name(fasta):
            sequence = ET.SubElement(alignment, 'sequence')
            ET.SubElement(sequence, 'taxon').set('idref', sample)
            if not p.exclude:
                sequence.text = f.get_partition(sample, p.p, fasta, gff)
            else:
                sequence.text = f.partition_exclude(sample, p.p, p.exclude, fasta, gff)
    return


def build_patterns(beast, partition_file):
    for p, exclude, every3, susb in f.get_partition_info(partition_file):
        p = f.Partition(p, exclude, every3, susb)
        if p.every3:
            for i in range(1, 4):
                patterns = ET.SubElement(beast, 'patterns', attrib={'id': '{0}.CP{1}.patterns'.format(p.p, i),
                                                                    'from': str(i), 'every': '3',
                                                                    'strip': 'false'})
                ET.SubElement(patterns, 'alignment').set('idref', 'alignment{}'.format(str(int(f.get_partition_name(partition_file).index(p.p)) + 1)))
                f.make_comment(beast, ' The unique patterns for codon position {} '.format(i))
        else:
            patterns = ET.SubElement(beast, 'patterns', attrib={'id': '{}.patterns'.format(p.p), 'from': '1', 'strip': 'false'})
            ET.SubElement(patterns, 'alignment').set('idref', 'alignment{}'.format(str(int(f.get_partition_name(partition_file).index(p.p)) + 1)))
            f.make_comment(beast, ' The unique patterns from 1 to end ')
    return