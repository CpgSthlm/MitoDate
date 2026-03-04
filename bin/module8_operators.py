import xml.etree.ElementTree as ET
import functions as f


def build_operator(beast, partition_file, fasta, split_partition):
    operators = ET.SubElement(beast, 'operators', attrib={'id': 'operators', 'optimizationSchedule': 'log'})
    f.make_comment(beast, ' Define operators ')

    subsModel_HKY(operators, partition_file, split_partition)
    strict_clock(operators, partition_file, split_partition)
    treeModel(operators)
    skygrid(operators)
    tip(operators, fasta)


def subsModel_HKY(operators, partition_file, split_partition):
    if split_partition:
        for p, exclude, every3, susb in f.get_partition_info(partition_file):
            p = f.Partition(p, exclude, every3, susb)

            if 'HKY' in p.subs_model:
                if p.every3:
                    for i in range(1, 4):
                        scaleOperator = ET.SubElement(operators, 'scaleOperator',
                                                      attrib={'scaleFactor': '0.75', 'weight': '1'})
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
                        gammaScaleOperator = ET.SubElement(operators, 'scaleOperator',
                                                           attrib={'scaleFactor': '0.75', 'weight': '1'})
                        ET.SubElement(gammaScaleOperator, 'parameter', attrib={'idref': '{}.CP{}.alpha'.format(p.p, i)})
                else:
                    gammaScaleOperator = ET.SubElement(operators, 'scaleOperator',
                                                       attrib={'scaleFactor': '0.75', 'weight': '1'})
                    ET.SubElement(gammaScaleOperator, 'parameter', attrib={'idref': '{}.alpha'.format(p.p)})

            if 'I' in p.subs_model:
                if p.every3:
                    for i in range(1, 4):
                        randomWalkOperator = ET.SubElement(operators, 'randomWalkOperator',
                                                           attrib={'windowSize': '0.75', 'weight': '1',
                                                                   'boundaryCondition': 'logit'})
                        ET.SubElement(randomWalkOperator, 'parameter', attrib={'idref': '{}.CP{}.pInv'.format(p.p, i)})
                else:
                    randomWalkOperator = ET.SubElement(operators, 'randomWalkOperator',
                                                       attrib={'windowSize': '0.75', 'weight': '1',
                                                               'boundaryCondition': 'logit'})
                    ET.SubElement(randomWalkOperator, 'parameter', attrib={'idref': '{}.pInv'.format(p.p)})


    else:
        for p, exclude, every3, susb in f.get_partition_info(partition_file):
            p = f.Partition(p, exclude, every3, susb)
        if 'HKY' in p.subs_model:
            scaleOperator = ET.SubElement(operators, 'scaleOperator', attrib={'scaleFactor': '0.75', 'weight': '1'})
            ET.SubElement(scaleOperator, 'parameter', attrib={'idref': 'kappa'})
            deltaExchange = ET.SubElement(operators, 'deltaExchange', attrib={'delta': '0.01', 'weight': '1'})
            ET.SubElement(deltaExchange, 'parameter', attrib={'idref': 'frequencies'})
        if 'G' in p.subs_model:
            gammaScaleOperator = ET.SubElement(operators, 'scaleOperator',
                                               attrib={'scaleFactor': '0.75', 'weight': '1'})
            ET.SubElement(gammaScaleOperator, 'parameter', attrib={'idref': 'alpha'})
        if 'I' in p.subs_model:
            randomWalkOperator = ET.SubElement(operators, 'randomWalkOperator',
                                               attrib={'windowSize': '0.75', 'weight': '1',
                                                       'boundaryCondition': 'logit'})
            ET.SubElement(randomWalkOperator, 'parameter', attrib={'idref': 'pInv'})


def strict_clock(operators, partition_file, split_partition=True):
    if split_partition:
        for p, exclude, every3, susb in f.get_partition_info(partition_file):
            p = f.Partition(p, exclude, every3, susb)

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
    else:
        clockScaleOperator = ET.SubElement(operators, 'scaleOperator', attrib={'scaleFactor': '0.75', 'weight': '3'})
        ET.SubElement(clockScaleOperator, 'parameter', attrib={'idref': 'clock.rate'})
        upDownOperator = ET.SubElement(operators, 'upDownOperator', attrib={'scaleFactor': '0.75', 'weight': '3'})
        up = ET.SubElement(upDownOperator, 'up')
        ET.SubElement(up, 'parameter', attrib={'idref': 'treeModel.allInternalNodeHeights'})
        down = ET.SubElement(upDownOperator, 'down')
        ET.SubElement(down, 'parameter', attrib={'idref': 'clock.rate'})


def treeModel(operators):
    # Define treeModel operators
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


def skygrid(operators):
    # Define population model operators (skygrid here)
    gmrfGridBlockUpdateOperator = ET.SubElement(operators, 'gmrfGridBlockUpdateOperator',
                                                attrib={'scaleFactor': '1.0', 'weight': '2'})
    ET.SubElement(gmrfGridBlockUpdateOperator, 'gmrfSkyrideLikelihood', attrib={'idref': 'skygrid'})
    popModelScaleOperator = ET.SubElement(operators, 'scaleOperator', attrib={'scaleFactor': '0.75', 'weight': '1'})
    ET.SubElement(popModelScaleOperator, 'parameter', attrib={'idref': 'skygrid.precision'})


def tip(operators, fasta):
    # Define tipping operators
    for sample in f.get_ND_taxa(fasta):
        tipOperator = ET.SubElement(operators, 'uniformOperator', attrib={'weight': '2'})
        ET.SubElement(tipOperator, 'parameter', attrib={'idref': 'age({})'.format(sample)})










