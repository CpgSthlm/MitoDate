import xml.etree.ElementTree as ET
import functions as f


def skygrid(beast, mean):
    gmrfSkyGridLikelihood = ET.SubElement(beast, 'gmrfSkyGridLikelihood')
    gmrfSkyGridLikelihood.set('id', 'skygrid')
    f.make_comment(beast, ' Generate a gmrfSkyGridLikelihood for the Bayesian SkyGrid process ')

    skygrid_popSizes = ET.SubElement(gmrfSkyGridLikelihood, 'populationSizes')
    ET.SubElement(skygrid_popSizes, 'parameter', attrib={'id': 'skygrid.logPopSize', 'dimension': '50', 'value': '1.0'})

    precisionParameter = ET.SubElement(gmrfSkyGridLikelihood, 'precisionParameter')
    ET.SubElement(precisionParameter, 'parameter', attrib={'id': 'skygrid.precision', 'value': '0.1', 'lower': '0.0'})

    numGridPoints = ET.SubElement(gmrfSkyGridLikelihood, 'numGridPoints')
    ET.SubElement(numGridPoints, 'parameter', attrib={'id': 'skygrid.numGridPoints', 'value': '49.0'})

    cutOff = ET.SubElement(gmrfSkyGridLikelihood, 'cutOff')
    ET.SubElement(cutOff, 'parameter', attrib={'id': 'skygrid.cutOff', 'value': str(f.skygrid_cutoff(mean))})

    populationTree = ET.SubElement(gmrfSkyGridLikelihood, 'populationTree')
    ET.SubElement(populationTree, 'treeModel', attrib={'idref': 'treeModel'})