#!/usr/bin/env python3
"""
FASTA Sequence Renaming Tool

Renames FASTA sequences based on metadata from an Excel file.
Matches Sample_IDs and creates new names in format: Sample_ID_Species_Origin_Group-By_Calibrated_yBP
"""

import os
import sys
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pathlib import Path

# ============================================================================
# DIRECT RUN CONFIGURATION
# ============================================================================
# Set your file paths here for direct execution
METADATA_FILE = ""
INPUT_FASTA = ""
OUTPUT_FASTA = '' 

# Set to True for direct run mode, False for command line mode
DIRECT_RUN_MODE = False
# ============================================================================


def read_metadata(metadata_file):
    """Read metadata from TSV file and validate required columns."""
    try:
        df = pd.read_csv(metadata_file, sep='\t')
        print(f"Read {len(df)} records from metadata file")
    
        
        # Check required columns
        required_cols = ['Sample_ID', 'Species', 'Origin', 'Group-By', 'Calibrated_yBP']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            print(f"Error: Missing required columns: {missing_cols}")
            sys.exit(1)
            
        return df
    except Exception as e:
        print(f"Error reading metadata file: {e}")
        sys.exit(1)


def find_matching_sample(seq_name, sample_ids):
    """Find Sample_ID that matches the sequence name (exact or smart matching)."""
    # First try exact matching
    for sample_id in sample_ids:
        sample_id_str = str(sample_id)
        if sample_id_str == seq_name:
            return sample_id
    
    # Then try smart matching - look for complete word boundaries
    for sample_id in sample_ids:
        sample_id_str = str(sample_id)
        # Check if sample_id appears as a complete token in seq_name
        if sample_id_str in seq_name:
            # Make sure it's not a substring of a longer identifier
            # Check boundaries: should be at start/end or surrounded by non-alphanumeric
            start_idx = seq_name.find(sample_id_str)
            end_idx = start_idx + len(sample_id_str)
            
            # Check if it's a complete match (not part of longer string)
            start_ok = (start_idx == 0 or not seq_name[start_idx-1].isalnum())
            end_ok = (end_idx == len(seq_name) or not seq_name[end_idx].isalnum())
            
            if start_ok and end_ok:
                return sample_id
        
        # Also check reverse (seq_name in sample_id) with same logic
        if seq_name in sample_id_str:
            start_idx = sample_id_str.find(seq_name)
            end_idx = start_idx + len(seq_name)
            
            start_ok = (start_idx == 0 or not sample_id_str[start_idx-1].isalnum())
            end_ok = (end_idx == len(sample_id_str) or not sample_id_str[end_idx].isalnum())
            
            if start_ok and end_ok:
                return sample_id
    
    return None


def create_new_name(metadata_row, sample_id):
    """Create new sequence name from metadata row."""
    sample_id_str = str(sample_id)
    species = str(metadata_row.get('Species', ''))
    origin = str(metadata_row.get('Origin', ''))
    group_by = str(metadata_row.get('Group-By', ''))
    calibrated_ybp = str(metadata_row.get('Calibrated_yBP', ''))
    
    # Handle missing values
    if group_by in ['nan', 'NA', 'None', '']:
        group_by = 'NA'
    if calibrated_ybp in ['nan', 'ND', 'None', '']:
        calibrated_ybp = 'ND'
    
    new_name = f"{sample_id_str}_{species}_{origin}_{group_by}_{calibrated_ybp}"
    return new_name.replace(' ', '_')


def process_fasta(input_fasta, output_fasta, metadata_df):
    """Process FASTA file and rename sequences based on metadata."""
    sample_ids = metadata_df['Sample_ID'].tolist()
    metadata_dict = metadata_df.set_index('Sample_ID').to_dict('index')
    
    renamed_sequences = []
    matched_count = 0
    total_sequences = 0
    unmatched_sequences = []
    
    print(f"Processing FASTA file: {input_fasta}")
    
    for record in SeqIO.parse(input_fasta, "fasta"):
        total_sequences += 1
        original_name = record.id
        
        matching_sample = find_matching_sample(original_name, sample_ids)
        
        if matching_sample:
            metadata_row = metadata_dict[matching_sample]
            new_name = create_new_name(metadata_row, matching_sample)
            
            new_record = SeqRecord(record.seq, id=new_name, description="")
            renamed_sequences.append(new_record)
            matched_count += 1
            # print(f"{original_name} -> {new_name}")
        else:
            unmatched_sequences.append(original_name)
    
    # Report results
    if unmatched_sequences:
        print(f"\nWarning: {len(unmatched_sequences)} sequences not found in metadata")
        for seq_name in unmatched_sequences[:5]:  # Show first 5
            print(f"  - {seq_name}")
        if len(unmatched_sequences) > 5:
            print(f"  ... and {len(unmatched_sequences) - 5} more")
    
    if renamed_sequences:
        SeqIO.write(renamed_sequences, output_fasta, "fasta")
        print(f"\nSuccessfully processed {matched_count}/{total_sequences} sequences")
        print(f"Output written to: {output_fasta}")
        
        # Check metadata coverage
        found_samples = {record.id.split('_')[0] for record in renamed_sequences}
        missing_samples = set(sample_ids) - found_samples
        
        if missing_samples:
            print(f"Warning: {len(missing_samples)} metadata samples not found in FASTA")
        else:
            print(f"All {len(sample_ids)} metadata samples were processed")
        
        return True
    else:
        print("Error: No sequences were matched and renamed")
        return False


def validate_files(metadata_file, input_fasta, output_fasta):
    """Validate input files and create output directory if needed."""
    if not os.path.exists(metadata_file):
        print(f"Error: Metadata file not found: {metadata_file}")
        return False
    
    if not os.path.exists(input_fasta):
        print(f"Error: Input FASTA file not found: {input_fasta}")
        return False
    
    output_dir = os.path.dirname(output_fasta)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    return True


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Rename FASTA sequences based on metadata",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -m metadata.xlsx -i input.fasta -o output.fasta
  %(prog)s --metadata data.xlsx --input seqs.fasta --output renamed.fasta

Output format: Sample_ID_Species_Origin_Group-By_Calibrated_yBP
        """
    )
    
    parser.add_argument('-m', '--metadata', required=True,
                       help='Path to metadata Excel file')
    parser.add_argument('-i', '--input', required=True,
                       help='Path to input FASTA file')
    parser.add_argument('-o', '--output', required=True,
                       help='Path to output FASTA file')
    
    return parser.parse_args()


def get_file_paths():
    """Get file paths either from command line args or direct configuration."""
    if DIRECT_RUN_MODE:
        # Direct run mode - use configured paths
        class DirectArgs:
            def __init__(self, metadata, input_file, output):
                self.metadata = metadata
                self.input = input_file
                self.output = output
        
        return DirectArgs(METADATA_FILE, INPUT_FASTA, OUTPUT_FASTA)
    else:
        # Command line mode - parse arguments
        return parse_arguments()


def main():
    """Main function."""
    args = get_file_paths()
    
    print("FASTA Sequence Renaming Tool")
    if DIRECT_RUN_MODE:
        print("Running in DIRECT MODE")
    print(f"Metadata: {args.metadata}")
    print(f"Input: {args.input}")
    print(f"Output: {args.output}")
    print()
    
    # Validate files
    if not validate_files(args.metadata, args.input, args.output):
        sys.exit(1)
    
    # Read metadata
    metadata_df = read_metadata(args.metadata)
    
    # Process FASTA
    success = process_fasta(args.input, args.output, metadata_df)
    
    if success:
        print("\nProcessing complete!")
    else:
        print("\nProcessing failed!")
        sys.exit(1)


if __name__ == "__main__":
    main()
