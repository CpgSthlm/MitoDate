import xml.etree.ElementTree as ET
import functions as f


def strict_clock(beast, partition_file):
    for p, exclude, every3, susb in f.get_partition_info(partition_file):
        p = f.Partition(p, exclude, every3, susb)

        strictClockBranchRates = ET.SubElement(beast, 'strictClockBranchRates',
                                               attrib={'id': '{}.branchRates'.format(p.p)})
        rate = ET.SubElement(strictClockBranchRates, 'rate')
        ET.SubElement(rate, 'parameter', attrib={'id': '{}.clock.rate'.format(p.p), 'value': '1.0', 'lower': '0.0'})

        rateStatistic = ET.SubElement(beast, 'rateStatistic',
                                      attrib={'id': '{}.meanRate'.format(p.p), 'name': '{}.meanRate'.format(p.p),
                                              'mode': 'mean', 'internal': 'true', 'external': 'true'})
        ET.SubElement(rateStatistic, 'treeModel', attrib={'idref': 'treeModel'})
        ET.SubElement(rateStatistic, 'strictClockBranchRates', attrib={'idref': '{}.branchRates'.format(p.p)})

        f.make_comment(beast, ' The strict clock (Uniform rates across branches) ', line=(-2))


