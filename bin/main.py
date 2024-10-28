#!/usr/bin/env python3

import module1_Taxa as mTaxa
import module2_Seq as mSeq
import module3_treeModel as mTreeModel
import module4_popModel as mPop
import module5_clockModel as mClock
import module6_subsModel as mSubs
import module7_likelihood as mLikeli
import module8_operators as mOperator
import module9_priors as mPriors
import module10_logs as mLog
import functions as f
import math
import sys
import xml.etree.ElementTree as ET

## Input files
fasta = sys.argv[1]
priors_table = sys.argv[2]
gff = sys.argv[3]
partition_file = sys.argv[4]

## Parameters

chainLength = int(sys.argv[5])
log_every = int(sys.argv[6])
partition_list = [p.strip() for p in sys.argv[7].split(',')]
ND_list = [nd.strip() for nd in sys.argv[8].split(',')]
taxon_set = [group.strip().strip('()').split(', ') for group in sys.argv[9].split('), (')]
taxon_set = [[item.strip() for item in sublist] for sublist in taxon_set]

# chainLength = 10000
# log_every = 1000
# ND_list = ['ND', 'P.']
# partition_list = ['tRNA', 'rRNA', 'CDS', 'D_loop']
# taxon_set = [['L.', 'P.'], ['M.'], ['M.', 'E.'], ['ND', 'P.']]


# Output file
log_name = fasta.replace('.fasta', '.xml')

beast_setting = f.Xml(units='years',
                      population_model='skygrid',
                      root_height_lower=1000000.0,
                      root_height_mean=5300000.0,
                      root_height_stdev=500000.0,
                      root_height_offset=0.0)

skygrid_cutoff = (int(list(str(beast_setting.root_height_mean))[0]) + 1) * \
                 (math.pow(10, len(str(round(beast_setting.root_height_mean)))))





# Module 1 -- Introduce all the taxa and taxon sets
beast = ET.Element('beast')
tree = ET.ElementTree(beast)

mTaxa.intro_taxa(beast, fasta, priors_table)
mTaxa.intro_taxonSet(beast, taxon_set, fasta)

# Module 2 -- Introduce aligned sequences with/without partitions
mSeq.build_partitionSeq(beast, partition_file, fasta, gff)
mSeq.build_patterns(beast, partition_file)

# Module 3 -- Introduce the tree model
mTreeModel.build_initialTree(beast, units=beast_setting.units)
mTreeModel.build_treeModel(beast, fasta)

mTreeModel.calc_tip_branchLen(beast)
mTreeModel.calc_taxonSet_branchLen(beast, taxon_set)

# Module 4 -- Define the population model
mPop.skygrid(beast, mean=beast_setting.root_height_mean)

# Module 5 -- Define the clock model
mClock.strict_clock(beast, partition_file)

# Module 6 -- Define the substitution model
mSubs.HKY(beast, partition_file)

# Module 7 -- Likelihoods
mLikeli.calc_likelihood(beast, partition_file)

# Module 8 -- Operators
mOperator.build_operator(beast, partition_file, fasta)

# Define MCMC
MCMC = ET.SubElement(beast, 'mcmc', attrib={'id': 'mcmc', 'chainLength': str(chainLength), 'autoOptimize': 'true'})
joint = ET.SubElement(MCMC, 'joint', attrib={'id': 'joint'})
jointPrior = ET.SubElement(joint, 'prior', attrib={'id': 'prior'})
f.make_comment(beast, ' Define MCMC ')

# Module 9 -- Priors
mPriors.define_operators(jointPrior, joint, MCMC, fasta, priors_table, partition_file,
                         rootHeight_mean=beast_setting.root_height_mean,
                         rootHeight_stdev=beast_setting.root_height_stdev,
                         rootHeight_offset=beast_setting.root_height_offset)

# Module 10 -- Write Logs
mLog.to_screen(MCMC, partition_file, log_every)
mLog.to_file(MCMC, log_every, log_name, taxon_set, partition_file, fasta)
mLog.trees_to_file(MCMC, log_every, log_name, partition_file)
mLog.report(beast)

# Test and store the XML file
f.pretty_xml(beast)
# ET.dump(beast)

tree.write(log_name, encoding='utf=8', xml_declaration=True)





