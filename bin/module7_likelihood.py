import xml.etree.ElementTree as ET
import functions as f


def calc_likelihood(beast, partition_file, split_partition):
    if split_partition:
        for p, exclude, every3, susb in f.get_partition_info(partition_file):
            p = f.Partition(p, exclude, every3, susb)

            treeDataLikelihood = ET.SubElement(beast, 'treeDataLikelihood',
                                               attrib={'id': '{}.treeLikelihood'.format(p.p), 'useAmbiguities': 'false'})
            f.make_comment(beast, ' Likelihood for tree given sequence data ')
            if p.every3:
                for i in range(1, 4):
                    treeLikelihood_partition = ET.SubElement(treeDataLikelihood, 'partition')
                    ET.SubElement(treeLikelihood_partition, 'patterns', attrib={'idref': '{}.CP{}.patterns'.format(p.p, i)})
                    ET.SubElement(treeLikelihood_partition, 'siteModel',
                                  attrib={'idref': '{}.CP{}.siteModel'.format(p.p, i)})
            else:
                treeLikelihood_partition = ET.SubElement(treeDataLikelihood, 'partition')
                ET.SubElement(treeLikelihood_partition, 'patterns', attrib={'idref': '{}.patterns'.format(p.p)})
                ET.SubElement(treeLikelihood_partition, 'siteModel', attrib={'idref': '{}.siteModel'.format(p.p)})
            ET.SubElement(treeDataLikelihood, 'treeModel', attrib={'idref': 'treeModel'})
            ET.SubElement(treeDataLikelihood, 'strictClockBranchRates', attrib={'idref': '{}.branchRates'.format(p.p)})
    else:
        treeDataLikelihood = ET.SubElement(beast, 'treeDataLikelihood',
                                           attrib={'id': 'treeLikelihood', 'useAmbiguities': 'false'})
        f.make_comment(beast, ' Likelihood for tree given sequence data ')

        treeLikelihood_partition = ET.SubElement(treeDataLikelihood, 'partition')
        ET.SubElement(treeLikelihood_partition, 'patterns', attrib={'idref': 'patterns'})
        ET.SubElement(treeLikelihood_partition, 'siteModel', attrib={'idref': 'siteModel'})
        ET.SubElement(treeDataLikelihood, 'treeModel', attrib={'idref': 'treeModel'})
        ET.SubElement(treeDataLikelihood, 'strictClockBranchRates', attrib={'idref': 'branchRates'})

    return
