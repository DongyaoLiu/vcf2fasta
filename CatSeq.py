#!/usr/bin/env python3
import os
import re
import argparse
from Bio import SeqIO
from collections import defaultdict
from multiprocessing import Pool, cpu_count

def check_id_format(seq_id):
    """Check if the sequence ID follows the IQ-TREE input requirements."""
    iqtree_id_pattern = re.compile(r'^[A-Za-z0-9_]+$')
    return bool(iqtree_id_pattern.match(seq_id))

def process_fasta_file(file_path, check_alignment, check_id):
    """Process a single FASTA file and return its sequences and lengths."""
    sequences_dict = SeqIO.index(file_path, 'fasta')

    # Check sequence IDs if required
    if check_id:
        for seq_id in sequences_dict:
            if not check_id_format(seq_id):
                raise ValueError(f"Sequence ID {seq_id} in file {file_path} does not match IQ-TREE format.")

    # Check alignment length if required
    if check_alignment:
        seq_lengths = set(len(sequences_dict[seq_id].seq) for seq_id in sequences_dict)
        if len(seq_lengths) > 1:
            raise ValueError(f"Sequences in file {file_path} have varying lengths.")

    # Store sequences and their lengths
    sequences = {seq_id: str(sequences_dict[seq_id].seq) for seq_id in sequences_dict}
    seq_length = len(next(iter(sequences_dict.values())).seq)

    return sequences, seq_length

def process_fasta_files(folder_name, output_fasta, output_partition, check_alignment, check_id):
    # List all FASTA files in the folder
    fasta_files = [f for f in os.listdir(folder_name) if f.endswith('.fasta') or f.endswith('.fa')]

    # Dictionary to store concatenated sequences for each ID
    concatenated_sequences = defaultdict(list)

    # Dictionary to store sequence lengths for each file (for partition file)
    file_lengths = {}

    # Use multiprocessing to process files in parallel
    with Pool(cpu_count()) as pool:
        results = pool.starmap(
            process_fasta_file,
            [(os.path.join(folder_name, fasta_file), check_alignment, check_id) for fasta_file in fasta_files]
        )

    # Combine results from all files
    for (sequences, seq_length), fasta_file in zip(results, fasta_files):
        for seq_id, seq in sequences.items():
            concatenated_sequences[seq_id].append(seq)
        file_lengths[fasta_file] = seq_length

    # Write concatenated sequences to output FASTA file
    with open(output_fasta, 'w') as out_fasta:
        for seq_id, seq_fragments in concatenated_sequences.items():
            concatenated_seq = ''.join(seq_fragments)
            out_fasta.write(f">{seq_id}\n{concatenated_seq}\n")

    # Write partition information to output partition file (if specified)
    if output_partition:
        with open(output_partition, 'w') as out_partition:
            start = 1
            for fasta_file in fasta_files:
                length = file_lengths[fasta_file]
                end = start + length - 1
                out_partition.write(f"DNA, {fasta_file} = {start}-{end}\n")
                start = end + 1

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process FASTA files for IQ-TREE input.")
    parser.add_argument('--folder', required=True, help="Folder containing input FASTA files.")
    parser.add_argument('--output_fasta', required=True, help="Output concatenated FASTA file.")
    parser.add_argument('--output_partition', help="Output partition file (optional).")
    parser.add_argument('--check_alignment', action='store_true', help="Check if sequences in each file are aligned (same length).")
    parser.add_argument('--check_id', action='store_true', help="Check if sequence IDs match IQ-TREE format.")

    # Parse arguments
    args = parser.parse_args()

    # Call the processing function
    process_fasta_files(
        folder_name=args.folder,
        output_fasta=args.output_fasta,
        output_partition=args.output_partition,
        check_alignment=args.check_alignment,
        check_id=args.check_id
    )

if __name__ == "__main__":
    main()
