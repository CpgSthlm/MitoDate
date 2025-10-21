#!/usr/bin/env python3

import pandas as pd
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
import xml.etree.ElementTree as ET
import os
import re
import argparse
import sys


def read_metadata(metadata_file, fasta):
    """Extract taxon sets and root height lower bound from metadata."""
    df = pd.read_csv(metadata_file, sep='\t')

    fasta_taxa = set(f.get_taxa_name(fasta))

    # Get unique taxon groups and filter by fasta content
    taxon_set = []
    for group in df['Group-By'].unique():
        has_sequence = any(str(group) in taxon for taxon in fasta_taxa)
        if has_sequence:
            taxon_set.append([group])
        else:
            print(f"Warning: Skipping group '{group}' - no sequences found in fasta file")

    # Get max calibrated age (numeric only)
    root_height_lower = pd.to_numeric(df['Calibrated_yBP'], errors='coerce').max()
    if pd.isna(root_height_lower):
        root_height_lower = 0.0

    return taxon_set, root_height_lower


def create_partition_file_from_args(args, output_dir):
    """Generate partition file based on args parameters."""
    if args.split_partition:
        # Parse multi-partition model: "tRNA:HKY+G+I, rRNA:HKY+G+I, ..."
        partitions = [p.split(':') for p in args.subs_model.split(', ')]
        data = [{'Partition': name, 'Exclude': False, 'Every3': name == 'CDS', 'SubstitutionModel': model}
                for name, model in partitions]
    else:
        # Single partition
        data = [{'Partition': 'all', 'Exclude': False, 'Every3': False, 'SubstitutionModel': args.subs_model}]

    partition_file = os.path.join(output_dir, 'temp_partition.xlsx')
    pd.DataFrame(data).to_excel(partition_file, index=False)
    return partition_file


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate BEAST XML file from configuration",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -f alignment.fasta -m metadata.tsv -p priors.csv --subs_model HKY --root_mean 50000 --root_stdev 5000 --chain_length 10000000 --log_every 1000
  %(prog)s -f alignment.fasta -m metadata.tsv -p priors.csv --subs_model "tRNA:HKY+G+I, rRNA:HKY+G+I, CDS:HKY+G+I" --split_partition --annotation genome.gff --root_mean 50000 --root_stdev 5000 --chain_length 10000000 --log_every 1000 -o output_dir/
        """
    )

    # Required arguments
    parser.add_argument('-f', '--fasta', required=True, help='Path to fasta alignment file')
    parser.add_argument('-m', '--metadata', required=True, help='Path to metadata TSV file')
    parser.add_argument('-p', '--priors', required=True, help='Path to priors table in CSV format')
    parser.add_argument('--subs_model', required=True,
                        help='Substitution model (HKY, GTR, etc.); when --split_partition is enabled, each partition must be specified (e.g. "tRNA:HKY+G+I, rRNA:HKY+G+I, CDS:HKY+G+I, D_loop:HKY+G")')
    parser.add_argument('--root_mean', type=float, required=True, help='Root height mean')
    parser.add_argument('--root_stdev', type=float, required=True, help='Root height stdev')
    parser.add_argument('--chain_length', type=int, required=True, help='Chain length')
    parser.add_argument('--log_every', type=int, required=True, help='Log every N iterations')

    # Optional arguments
    parser.add_argument('--population_model', default='skygrid', help='Population model (default: skygrid)')
    parser.add_argument('--clock_model', default='strict', help='Clock model (default: strict)')
    parser.add_argument('--split_partition', action='store_true', help='Split partition (default: False)')
    parser.add_argument('--annotation', default='', help='Path to annotation file in GFF format, only used when --split_partition is enabled')
    parser.add_argument('--root_offset', type=float, default=0.0, help='Root height offset (default: 0)')
    parser.add_argument('-o', '--output_dir', default=os.getcwd(), help='Output directory for generated files (default: current directory)')

    return parser.parse_args()


def main():
    """Main function."""
    args = parse_arguments()

    print("BEAST XML Generator Ver 2.0")
    print(f"Fasta file: {args.fasta}")
    print()

    # Extract all parameters from args into local variables for module functions
    fasta = args.fasta
    metadata_file = args.metadata
    priors_table = args.priors
    gff = args.annotation
    split_partition = args.split_partition
    population_model = args.population_model
    clock_model = args.clock_model
    root_height_mean = args.root_mean
    root_height_stdev = args.root_stdev
    root_height_offset = args.root_offset
    chainLength = args.chain_length
    log_every = args.log_every

    # Generate log_name: fasta_prefix
    log_name = os.path.basename(fasta).replace('.fasta', '')

    # Validate required files
    required_files = [fasta, priors_table, metadata_file]
    for file_path in required_files:
        if not os.path.exists(file_path):
            print(f"Error: Required file not found: {file_path}")
            sys.exit(1)

    print(f"Processing alignment: {fasta}")
    print(f"Using metadata: {metadata_file}")
    print(f"Using priors: {priors_table}")
    print(f"Output prefix: {log_name}")
    print(f"Population model: {population_model}")
    print(f"Clock model: {clock_model}")
    print(f"Substitution model: {args.subs_model}")
    print(f"Split partition: {split_partition}")
    if split_partition and gff:
        print(f"Annotation file: {gff}")
    print()

    # Read metadata and create partition file
    taxon_set, root_height_lower = read_metadata(metadata_file, fasta)
    partition_file = create_partition_file_from_args(args, args.output_dir)

    # BEAST settings
    beast_setting = f.Xml(
        units='years',
        population_model=population_model,
        root_height_mean=root_height_mean,
        root_height_stdev=root_height_stdev,
        root_height_offset=root_height_offset,
        root_height_lower=root_height_lower
    )

    # Calculate skygrid cutoff
    magnitude = 10 ** (len(str(int(beast_setting.root_height_mean))) - 1)
    skygrid_cutoff = math.ceil(beast_setting.root_height_mean / magnitude) * magnitude

    # Read substitution model from partition file
    with open(partition_file, 'rb') as p:
        subs_model = pd.read_excel(p, header=0)['SubstitutionModel'].iloc[0]

    print("Building XML structure...")

    # Module 1 -- Introduce all the taxa and taxon sets
    beast = ET.Element('beast')
    tree = ET.ElementTree(beast)
    mTaxa.intro_taxa(beast, fasta, priors_table)
    mTaxa.intro_taxonSet(beast, taxon_set, fasta)

    # Module 2 -- Introduce aligned sequences with/without partitions
    mSeq.build_partitionSeq(beast, partition_file, fasta, gff, split_partition)
    mSeq.build_patterns(beast, partition_file, split_partition)

    # Module 3 -- Introduce the tree model
    mTreeModel.build_initialTree(beast, units=beast_setting.units)
    mTreeModel.build_treeModel(beast, fasta)

    mTreeModel.calc_tip_branchLen(beast)
    mTreeModel.calc_taxonSet_branchLen(beast, taxon_set)

    # Module 4 -- Define the population model
    mPop.skygrid(beast, mean=beast_setting.root_height_mean)

    # Module 5 -- Define the clock model
    mClock.strict_clock(beast, partition_file, split_partition)

    # Module 6 -- Define the substitution model
    if 'HKY' in subs_model:
        mSubs.HKY(beast, partition_file, split_partition)
    elif 'GTR' in subs_model:
        mSubs.GTR(beast, partition_file, split_partition)

    # Module 7 -- Likelihoods
    mLikeli.calc_likelihood(beast, partition_file, split_partition)

    # Module 8 -- Operators
    mOperator.build_operator(beast, partition_file, fasta, split_partition)

    # Define MCMC
    MCMC = ET.SubElement(beast, 'mcmc', attrib={'id': 'mcmc', 'chainLength': str(chainLength), 'autoOptimize': 'true'})
    joint = ET.SubElement(MCMC, 'joint', attrib={'id': 'joint'})
    jointPrior = ET.SubElement(joint, 'prior', attrib={'id': 'prior'})
    f.make_comment(beast, ' Define MCMC ')

    # Module 9 -- Priors
    mPriors.define_operators(jointPrior, joint, MCMC, fasta, priors_table, partition_file, split_partition,
                             rootHeight_mean=beast_setting.root_height_mean,
                             rootHeight_stdev=beast_setting.root_height_stdev,
                             rootHeight_offset=beast_setting.root_height_offset)

    # Module 10 -- Write Logs
    mLog.to_screen(MCMC, partition_file, log_every, split_partition)
    mLog.to_file(MCMC, log_every, log_name, taxon_set, partition_file, fasta, split_partition)
    mLog.trees_to_file(MCMC, log_every, log_name, partition_file, split_partition)
    mLog.report(beast)

    # Test and store the XML file
    f.pretty_xml(beast)

    # Generate output XML file path
    xml_output = os.path.join(args.output_dir, f'{log_name}.xml')

    tree.write(xml_output, encoding='utf-8', xml_declaration=True)

    # Clean up fake partition file
    os.remove(partition_file)

    print(f"\nXML file generated successfully: {xml_output}")
    print("Processing complete!")


if __name__ == "__main__":
    main()




