import xml.etree.ElementTree as ET
import functions as f


def HKY(beast, partition_file, split_partition):
    if split_partition:
        for p, exclude, every3, susb in f.get_partition_info(partition_file):
            p = f.Partition(p, exclude, every3, susb)

            if p.every3:
                for i in range(1, 4):
                    HKYModel = ET.SubElement(beast, 'HKYModel', attrib={'id': '{0}.CP{1}.hky'.format(p.p, i)})
                    frequencies = ET.SubElement(HKYModel, 'frequencies')
                    frequencyModel = ET.SubElement(frequencies, 'frequencyModel', attrib={'dataType': 'nucleotide'})
                    frequencies_parameter = ET.SubElement(frequencyModel, 'frequencies')
                    ET.SubElement(frequencies_parameter, 'parameter',
                                  attrib={'id': '{0}.CP{1}.frequencies'.format(p.p, i), 'value': '0.25 0.25 0.25 0.25'})
                    kappa = ET.SubElement(HKYModel, 'kappa')
                    ET.SubElement(kappa, 'parameter',
                                  attrib={'id': '{0}.CP{1}.kappa'.format(p.p, i), 'value': '2.0', 'lower': '0.0'})
                for i in range(1, 4):
                    siteModel = ET.SubElement(beast, 'siteModel', attrib={'id': '{0}.CP{1}.siteModel'.format(p.p, i)})
                    substitutionModel = ET.SubElement(siteModel, 'substitutionModel')
                    ET.SubElement(substitutionModel, 'HKYModel', attrib={'idref': '{0}.CP{1}.hky'.format(p.p, i)})
                    relativeRate = ET.SubElement(siteModel, 'relativeRate', attrib={'weight': '3.0'})
                    ET.SubElement(relativeRate, 'parameter',
                                  attrib={'id': '{}.CP{}.nu'.format(p.p, i), 'value': '0.3333333333333333', 'lower': '0.0',
                                          'upper': '1.0'})
                    if 'G' in p.subs_model:
                        gammaShape = ET.SubElement(siteModel, 'gammaShape', attrib={'gammaCategories': '4'})
                        ET.SubElement(gammaShape, 'parameter',
                                      attrib={'id': '{0}.CP{1}.alpha'.format(p.p, i), 'value': '0.5', 'lower': '0.0'})
                    if 'I' in p.subs_model:
                        proportionInvariant = ET.SubElement(siteModel, 'proportionInvariant')
                        ET.SubElement(proportionInvariant, 'parameter',
                                      attrib={'id': '{}.CP{}.pInv'.format(p.p, i), 'value': '0.5', 'lower': '0.0',
                                              'upper': '1.0'})
                    subs_statistic = ET.SubElement(beast, 'statistic',
                                                   attrib={'id': '{}.CP{}.mu'.format(p.p, i), 'name': 'mu'})
                    ET.SubElement(subs_statistic, 'siteModel', attrib={'idref': '{}.CP{}.siteModel'.format(p.p, i)})
                f.make_comment(beast, ' The HKY substitution model (Hasegawa, Kishino and Yano, 1985) ', -9)

            else:
                HKYModel = ET.SubElement(beast, 'HKYModel', attrib={'id': '{}.hky'.format(p.p)})
                frequencies = ET.SubElement(HKYModel, 'frequencies')
                frequencyModel = ET.SubElement(frequencies, 'frequencyModel', attrib={'dataType': 'nucleotide'})
                frequencies_parameter = ET.SubElement(frequencyModel, 'frequencies')
                ET.SubElement(frequencies_parameter, 'parameter',
                              attrib={'id': '{}.frequencies'.format(p.p), 'value': '0.25 0.25 0.25 0.25'})
                kappa = ET.SubElement(HKYModel, 'kappa')
                ET.SubElement(kappa, 'parameter', attrib={'id': '{}.kappa'.format(p.p), 'value': '2.0', 'lower': '0.0'})

                siteModel = ET.SubElement(beast, 'siteModel', attrib={'id': '{}.siteModel'.format(p.p)})
                substitutionModel = ET.SubElement(siteModel, 'substitutionModel')
                ET.SubElement(substitutionModel, 'HKYModel', attrib={'idref': '{}.hky'.format(p.p)})
                if 'G' in p.subs_model:
                    gammaShape = ET.SubElement(siteModel, 'gammaShape', attrib={'gammaCategories': '5'})
                    ET.SubElement(gammaShape, 'parameter',
                                  attrib={'id': '{}.alpha'.format(p.p), 'value': '0.5', 'lower': '0.0'})
                if 'I' in p.subs_model:
                    proportionInvariant = ET.SubElement(siteModel, 'proportionInvariant')
                    ET.SubElement(proportionInvariant, 'parameter',
                                  attrib={'id': '{}.pInv'.format(p.p), 'value': '0.5', 'lower': '0.0', 'upper': '1.0'})

                subs_statistic = ET.SubElement(beast, 'statistic', attrib={'id': '{}.mu'.format(p.p), 'name': 'mu'})
                ET.SubElement(subs_statistic, 'siteModel', attrib={'idref': '{}.siteModel'.format(p.p)})

                f.make_comment(beast, ' The HKY substitution model (Hasegawa, Kishino and Yano, 1985) ', -3)

        for p, exclude, every3, susb in f.get_partition_info(partition_file):
            p = f.Partition(p, exclude, every3, susb)
            if p.every3:
                compoundParameter = ET.SubElement(beast, 'compoundParameter', attrib={'id': '{}.allNus'.format(p.p)})
                for i in range(1, 4):
                    ET.SubElement(compoundParameter, 'parameter', attrib={'idref': '{}.CP{}.nu'.format(p.p, i)})

    else:
        for p, exclude, every3, susb in f.get_partition_info(partition_file):
            p = f.Partition(p, exclude, every3, susb)

        HKYModel = ET.SubElement(beast, 'HKYModel', attrib={'id': 'hky'})
        frequencies = ET.SubElement(HKYModel, 'frequencies')
        frequencyModel = ET.SubElement(frequencies, 'frequencyModel', attrib={'dataType': 'nucleotide'})
        frequencies_parameter = ET.SubElement(frequencyModel, 'frequencies')
        ET.SubElement(frequencies_parameter, 'parameter',
                      attrib={'id': 'frequencies', 'value': '0.25 0.25 0.25 0.25'})
        kappa = ET.SubElement(HKYModel, 'kappa')
        ET.SubElement(kappa, 'parameter', attrib={'id': 'kappa', 'value': '2.0', 'lower': '0.0'})

        siteModel = ET.SubElement(beast, 'siteModel', attrib={'id': 'siteModel'})
        substitutionModel = ET.SubElement(siteModel, 'substitutionModel')
        ET.SubElement(substitutionModel, 'HKYModel', attrib={'idref': 'hky'})
        if 'G' in p.subs_model:
            gammaShape = ET.SubElement(siteModel, 'gammaShape', attrib={'gammaCategories': '5'})
            ET.SubElement(gammaShape, 'parameter',
                          attrib={'id': 'alpha', 'value': '0.5', 'lower': '0.0'})
        if 'I' in p.subs_model:
            proportionInvariant = ET.SubElement(siteModel, 'proportionInvariant')
            ET.SubElement(proportionInvariant, 'parameter',
                          attrib={'id': 'pInv', 'value': '0.5', 'lower': '0.0', 'upper': '1.0'})

        subs_statistic = ET.SubElement(beast, 'statistic', attrib={'id': 'mu', 'name': 'mu'})
        ET.SubElement(subs_statistic, 'siteModel', attrib={'idref': 'siteModel'})
        f.make_comment(beast, ' The HKY substitution model (Hasegawa, Kishino and Yano, 1985) ', -3)

    return



def GTR(beast, partition_file, split_partition):
    if split_partition:
        for p, exclude, every3, susb in f.get_partition_info(partition_file):
            p = f.Partition(p, exclude, every3, susb)

            if p.every3:
                for i in range(1, 4):
                    GTRModel = ET.SubElement(beast, 'GTRModel', attrib={'id': '{0}.CP{1}.gtr'.format(p.p, i)})
                    frequencies = ET.SubElement(GTRModel, 'frequencies')
                    frequencyModel = ET.SubElement(frequencies, 'frequencyModel', attrib={'dataType': 'nucleotide'})
                    frequencies_parameter = ET.SubElement(frequencyModel, 'frequencies')
                    ET.SubElement(frequencies_parameter, 'parameter',
                                  attrib={'id': '{0}.CP{1}.frequencies'.format(p.p, i), 'value': '0.25 0.25 0.25 0.25'})
                    rates = ET.SubElement(GTRModel, 'rates')
                    ET.SubElement(rates, 'parameter',
                                  attrib={'id': '{0}.CP{1}.gtr.rates'.format(p.p, i), 'dimension': '6', 'lower': '0.0'})
                for i in range(1, 4):
                    siteModel = ET.SubElement(beast, 'siteModel', attrib={'id': '{0}.CP{1}.siteModel'.format(p.p, i)})
                    substitutionModel = ET.SubElement(siteModel, 'substitutionModel')
                    ET.SubElement(substitutionModel, 'GTRModel', attrib={'idref': '{0}.CP{1}.gtr'.format(p.p, i)})
                    relativeRate = ET.SubElement(siteModel, 'relativeRate', attrib={'weight': '3.0'})
                    ET.SubElement(relativeRate, 'parameter',
                                  attrib={'id': '{}.CP{}.nu'.format(p.p, i), 'value': '0.3333333333333333', 'lower': '0.0',
                                          'upper': '1.0'})
                    if 'G' in p.subs_model:
                        gammaShape = ET.SubElement(siteModel, 'gammaShape', attrib={'gammaCategories': '4'})
                        ET.SubElement(gammaShape, 'parameter',
                                      attrib={'id': '{0}.CP{1}.alpha'.format(p.p, i), 'value': '0.5', 'lower': '0.0'})
                    if 'I' in p.subs_model:
                        proportionInvariant = ET.SubElement(siteModel, 'proportionInvariant')
                        ET.SubElement(proportionInvariant, 'parameter',
                                      attrib={'id': '{}.CP{}.pInv'.format(p.p, i), 'value': '0.5', 'lower': '0.0',
                                              'upper': '1.0'})
                    subs_statistic = ET.SubElement(beast, 'statistic',
                                                   attrib={'id': '{}.CP{}.mu'.format(p.p, i), 'name': 'mu'})
                    ET.SubElement(subs_statistic, 'siteModel', attrib={'idref': '{}.CP{}.siteModel'.format(p.p, i)})
                f.make_comment(beast, ' The GTR substitution model (General Time Reversible) ', -9)

            else:
                GTRModel = ET.SubElement(beast, 'GTRModel', attrib={'id': '{}.gtr'.format(p.p)})
                frequencies = ET.SubElement(GTRModel, 'frequencies')
                frequencyModel = ET.SubElement(frequencies, 'frequencyModel', attrib={'dataType': 'nucleotide'})
                frequencies_parameter = ET.SubElement(frequencyModel, 'frequencies')
                ET.SubElement(frequencies_parameter, 'parameter',
                              attrib={'id': '{}.frequencies'.format(p.p), 'value': '0.25 0.25 0.25 0.25'})
                rates = ET.SubElement(GTRModel, 'rates')
                ET.SubElement(rates, 'parameter', attrib={'id': '{}.gtr.rates'.format(p.p), 'dimension': '6', 'lower': '0.0'})

                siteModel = ET.SubElement(beast, 'siteModel', attrib={'id': '{}.siteModel'.format(p.p)})
                substitutionModel = ET.SubElement(siteModel, 'substitutionModel')
                ET.SubElement(substitutionModel, 'GTRModel', attrib={'idref': '{}.gtr'.format(p.p)})
                if 'G' in p.subs_model:
                    gammaShape = ET.SubElement(siteModel, 'gammaShape', attrib={'gammaCategories': '5'})
                    ET.SubElement(gammaShape, 'parameter',
                                  attrib={'id': '{}.alpha'.format(p.p), 'value': '0.5', 'lower': '0.0'})
                if 'I' in p.subs_model:
                    proportionInvariant = ET.SubElement(siteModel, 'proportionInvariant')
                    ET.SubElement(proportionInvariant, 'parameter',
                                  attrib={'id': '{}.pInv'.format(p.p), 'value': '0.5', 'lower': '0.0', 'upper': '1.0'})

                subs_statistic = ET.SubElement(beast, 'statistic', attrib={'id': '{}.mu'.format(p.p), 'name': 'mu'})
                ET.SubElement(subs_statistic, 'siteModel', attrib={'idref': '{}.siteModel'.format(p.p)})

                f.make_comment(beast, ' The GTR substitution model (General Time Reversible) ', -3)

        for p, exclude, every3, susb in f.get_partition_info(partition_file):
            p = f.Partition(p, exclude, every3, susb)
            if p.every3:
                compoundParameter = ET.SubElement(beast, 'compoundParameter', attrib={'id': '{}.allNus'.format(p.p)})
                for i in range(1, 4):
                    ET.SubElement(compoundParameter, 'parameter', attrib={'idref': '{}.CP{}.nu'.format(p.p, i)})

    else:
        for p, exclude, every3, susb in f.get_partition_info(partition_file):
            p = f.Partition(p, exclude, every3, susb)

        GTRModel = ET.SubElement(beast, 'GTRModel', attrib={'id': 'gtr'})
        frequencies = ET.SubElement(GTRModel, 'frequencies')
        frequencyModel = ET.SubElement(frequencies, 'frequencyModel', attrib={'dataType': 'nucleotide'})
        frequencies_parameter = ET.SubElement(frequencyModel, 'frequencies')
        ET.SubElement(frequencies_parameter, 'parameter',
                      attrib={'id': 'frequencies',  'value': '0.25 0.25 0.25 0.25'})
        rates = ET.SubElement(GTRModel, 'rates')
        ET.SubElement(rates, 'parameter', attrib={'id': 'gtr.rates', 'dimension': '6', 'lower': '0.0'})

        siteModel = ET.SubElement(beast, 'siteModel', attrib={'id': 'siteModel'})
        substitutionModel = ET.SubElement(siteModel, 'substitutionModel')
        ET.SubElement(substitutionModel, 'GTRModel', attrib={'idref': 'gtr'})
        if 'G' in p.subs_model:
            gammaShape = ET.SubElement(siteModel, 'gammaShape', attrib={'gammaCategories': '5'})
            ET.SubElement(gammaShape, 'parameter',
                          attrib={'id': 'alpha', 'value': '0.5', 'lower': '0.0'})
        if 'I' in p.subs_model:
            proportionInvariant = ET.SubElement(siteModel, 'proportionInvariant')
            ET.SubElement(proportionInvariant, 'parameter',
                          attrib={'id': 'pInv', 'value': '0.5', 'lower': '0.0', 'upper': '1.0'})

        subs_statistic = ET.SubElement(beast, 'statistic', attrib={'id': 'mu', 'name': 'mu'})
        ET.SubElement(subs_statistic, 'siteModel', attrib={'idref': 'siteModel'})
        f.make_comment(beast, ' The GTR substitution model (General Time Reversible) ', -3)

    return

