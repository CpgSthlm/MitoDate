#!/usr/bin/env python3
"""
Filter and rename sequences based on metadata.
Extracts dated sequences (Calibrated_yBP != 'ND') and renames using metadata file.
"""

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import os
import argparse

# Direct run configuration
METADATA_FILE = ""
ALIGNMENT_FILE = ""
DIRECT_RUN_MODE = False


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Filter and rename sequences based on metadata",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example: %(prog)s -m metadata.tsv -i alignment.fasta -o output_dir/"
    )
    parser.add_argument('-m', '--metadata', required=True, help='Metadata TSV file')
    parser.add_argument('-i', '--input', required=True, help='Input FASTA alignment file')
    parser.add_argument('-o', '--output-dir', help='Output directory (default: same as input)')
    return parser.parse_args()


def get_file_paths():
    """Get file paths from command line or direct configuration."""
    if DIRECT_RUN_MODE:
        class DirectArgs:
            def __init__(self):
                self.metadata, self.input, self.output_dir = METADATA_FILE, ALIGNMENT_FILE, None
        return DirectArgs()
    return parse_arguments()


def find_matching_sample(seq_name, sample_ids):
    """Find Sample_ID matching sequence name with exact or boundary matching."""
    # Exact match
    if seq_name in sample_ids:
        return seq_name
    
    # Boundary matching
    for sample_id in sample_ids:
        sid_str = str(sample_id)
        for name, target in [(seq_name, sid_str), (sid_str, seq_name)]:
            if target in name:
                idx = name.find(target)
                start_ok = idx == 0 or not name[idx-1].isalnum()
                end_ok = idx + len(target) == len(name) or not name[idx + len(target)].isalnum()
                if start_ok and end_ok:
                    return sample_id
    return None


def create_new_name(metadata_row, sample_id):
    """Create new sequence name: SampleID_Species_Origin_GroupBy_CalibratedYBP"""
    fields = [
        str(sample_id),
        str(metadata_row.get('Species', '')),
        str(metadata_row.get('Origin', '')),
        str(metadata_row.get('Group-By', 'NA')),
        str(metadata_row.get('Calibrated_yBP', 'ND'))
    ]
    # Clean missing values
    fields = ['NA' if f in ['nan', 'NA', 'None', ''] else f for f in fields]
    return '_'.join(fields).replace(' ', '_')


def process_sequences(input_file, sample_ids, metadata_dict):
    """Process alignment file and rename matching sequences."""
    filtered = []
    for seq in SeqIO.parse(input_file, "fasta"):
        match = find_matching_sample(seq.id, sample_ids)
        if match:
            new_name = create_new_name(metadata_dict[match], match)
            filtered.append(SeqRecord(seq.seq, id=new_name, description=""))
            # print(f"  {seq.id} -> {new_name}")
    return filtered


def main():
    args = get_file_paths()
    
    print(f"Single Dating MSA Processing")
    print(f"Metadata: {args.metadata}")
    print(f"Input: {args.input}\n")
    
    # Setup output directory
    output_dir = args.output_dir if hasattr(args, 'output_dir') and args.output_dir else os.path.dirname(args.input)
    os.makedirs(output_dir, exist_ok=True)
    
    input_basename = os.path.splitext(os.path.basename(args.input))[0]
    output_file = os.path.join(output_dir, f"{input_basename}.DatedSamples.fasta")
    
    # Validate files
    for path, name in [(args.metadata, "Metadata"), (args.input, "Alignment")]:
        if not os.path.exists(path):
            sys.exit(f"Error: {name} file not found: {path}")
    
    # Load metadata
    print("Loading metadata...")
    metadata_df = pd.read_csv(args.metadata, sep='\t')
    print(f"  {len(metadata_df)} rows, columns: {list(metadata_df.columns)}")
    
    # Filter dated sequences
    filtered_meta = metadata_df[metadata_df['Calibrated_yBP'] != 'ND'].copy()
    print(f"  {len(filtered_meta)} dated sequences")
    
    if filtered_meta.empty:
        sys.exit("No dated sequences found")
    
    # Process sequences
    print("\nProcessing sequences...")
    sample_ids = filtered_meta['Sample_ID'].tolist()
    metadata_dict = filtered_meta.set_index('Sample_ID').to_dict('index')
    filtered_seqs = process_sequences(args.input, sample_ids, metadata_dict)
    
    if not filtered_seqs:
        sys.exit("No matching sequences found")
    
    # Write output
    SeqIO.write(filtered_seqs, output_file, "fasta")
    print(f"\nWrote {len(filtered_seqs)} sequences to {output_file}")
    
    # Process TipDating sequences
    tipdating_meta = metadata_df[metadata_df['TipDating'] == 'Y']
    if not tipdating_meta.empty:
        print(f"\nProcessing {len(tipdating_meta)} TipDating sequences...")
        all_seqs = {rec.id: rec for rec in SeqIO.parse(args.input, "fasta")}
        full_meta_dict = metadata_df.set_index('Sample_ID').to_dict('index')
        
        for _, row in tipdating_meta.iterrows():
            sample_id = str(row['Sample_ID'])
            tipdating_file = os.path.join(output_dir, f"{input_basename}.TipDating.{sample_id}.fasta")
            
            # Find matching sequence
            match_id = next((sid for sid in all_seqs if find_matching_sample(sid, [sample_id])), None)
            
            if match_id:
                new_name = create_new_name(full_meta_dict[sample_id], sample_id)
                tipdating_rec = SeqRecord(all_seqs[match_id].seq, id=new_name, description="")
                combined = filtered_seqs + [tipdating_rec]
                SeqIO.write(combined, tipdating_file, "fasta")
                # print(f"  Created {sample_id}: {len(combined)} sequences")
            else:
                print(f"  Warning: {sample_id} not found")
    
    # Summary
    print(f"\nSummary:")
    print(f"  Total metadata rows: {len(metadata_df)}")
    print(f"  Dated sequences: {len(filtered_meta)}")
    print(f"  TipDating sequences: {len(tipdating_meta)}")
    # print(f"  Output: {output_file}")
    print("Complete!")


if __name__ == "__main__":
    main()
