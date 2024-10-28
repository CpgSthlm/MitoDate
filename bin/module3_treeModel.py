import xml.etree.ElementTree as ET
import functions as f
import pandas as pd
from Bio import SeqIO
import math


def build_initialTree(beast, units):
    constantSize = ET.SubElement(beast, 'constantSize', attrib={'id': 'initialDemo', 'units': units})
    populationSize = ET.SubElement(constantSize, 'populationSize')
    popSize_parameter = ET.SubElement(populationSize, 'parameter', attrib={'id': 'initialDemo.popSize', 'value': '100.0'})
    f.make_comment(beast, ' This is a simple constant population size coalescent model ')
    f.make_comment(beast, ' to generate an initial tree for the chain ')

    coalescentSimulator = ET.SubElement(beast, 'coalescentSimulator', attrib={'id': 'startingTree'})
    ET.SubElement(coalescentSimulator, 'taxa').set('idref', 'taxa')
    ET.SubElement(coalescentSimulator, 'constantSize').set('idref', 'initialDemo')
    f.make_comment(beast, ' Generate a random starting tree under the coalescent process ')


def build_treeModel(beast, fasta):
    treeModel = ET.SubElement(beast, 'treeModel', attrib={'id': 'treeModel'})
    ET.SubElement(treeModel, 'coalescentTree').set('idref', 'startingTree')
    f.make_comment(beast, ' Generate a tree model ')

    rootHeight = ET.SubElement(treeModel, 'rootHeight')
    ET.SubElement(rootHeight, 'parameter').set('id', 'treeModel.rootHeight')

    nodeHeights = ET.SubElement(treeModel, 'nodeHeights')
    nodeHeights.set('internalNodes', 'true')
    ET.SubElement(nodeHeights, 'parameter').set('id', 'treeModel.internalNodeHeights')

    nodeHeights = ET.SubElement(treeModel, 'nodeHeights')
    nodeHeights.set('internalNodes', 'true')
    nodeHeights.set('rootNode', 'true')
    ET.SubElement(nodeHeights, 'parameter').set('id', 'treeModel.allInternalNodeHeights')

    for sample in f.get_ND_taxa(fasta):
        leafHeight = ET.SubElement(treeModel, 'leafHeight')
        leafHeight.set('taxon', sample)
        ET.SubElement(leafHeight, 'parameter').set('id', 'age({})'.format(sample))


def calc_tip_branchLen(beast):
    treeLengthStatistic = ET.SubElement(beast, 'treeLengthStatistic', attrib={'id': 'treeLength'})
    ET.SubElement(treeLengthStatistic, 'treeModel').set('idref', 'treeModel')
    f.make_comment(beast, ' Statistic for sum of the branch lengths of the tree (tree length) ')

    tmrcaStatistic = ET.SubElement(beast, 'tmrcaStatistic', attrib={'id': 'age(root)', 'absolute': 'true'})
    ET.SubElement(tmrcaStatistic, 'treeModel').set('idref', 'treeModel')
    f.make_comment(beast, ' Statistic for time of most recent common ancestor of tree ')


def calc_taxonSet_branchLen(beast, taxon_set):
    for set in taxon_set:
        set_name = [x.strip('.') for x in set]
        tmrcaStatistic_taxon = ET.SubElement(beast, 'tmrcaStatistic',
                                             attrib={'id': 'tmrca({})'.format('_and_'.join(set_name)),
                                                     'absolute': 'false', 'includeStem': 'false'})
        mrca = ET.SubElement(tmrcaStatistic_taxon, 'mrca')
        ET.SubElement(mrca, 'taxa').set('idref', '_and_'.join(set_name))
        ET.SubElement(tmrcaStatistic_taxon, 'treeModel').set('idref', 'treeModel')
        tmrcaStatistic_taxon = ET.SubElement(beast, 'tmrcaStatistic',
                                             attrib={'id': 'age({})'.format('_and_'.join(set_name)), 'absolute': 'true',
                                                     'includeStem': 'false'})
        mrca = ET.SubElement(tmrcaStatistic_taxon, 'mrca')
        ET.SubElement(mrca, 'taxa').set('idref', '_and_'.join(set_name))
        ET.SubElement(tmrcaStatistic_taxon, 'treeModel').set('idref', 'treeModel')
    f.make_comment(beast, ' Taxon Sets', line=-(len(taxon_set) * 2))
