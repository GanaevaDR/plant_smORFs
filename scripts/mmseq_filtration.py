#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This script filters smORFs from a FASTA file based on MMseqs2 search results.
Only smORFs that do NOT appear in the MMseqs2 output are retained.
"""

from Bio import SeqIO
import pandas as pd
import argparse


def write_fasta(filtered_records, output_path):
    """
    Write a list of SeqRecord objects to a FASTA file.
    """
    try:
        SeqIO.write(filtered_records, output_path, "fasta")
        print(f"Filtered sequences saved to: {output_path}")
    except Exception as e:
        print(f"Error while writing FASTA file: {e}")


if __name__ == "__main__":
    print("*" * 20)
    print("START\n")
    
    # Argument parsing
    parser = argparse.ArgumentParser(description='Filter smORFs based on MMseqs2 search results.')
    parser.add_argument("--mmseq", type=str, required=True, help='Path to MMseqs2 search results.')
    parser.add_argument("--fasta", type=str, required=True, help='Path to the input smORFs FASTA file.')
    parser.add_argument("--out_smorfs", type=str, required=True, help='Path to the output filtered FASTA file.')
    args = parser.parse_args()
    
    # Step 1 - Reading MMseqs2 output
    print("Step 1 - Reading MMseqs2 output")
    try:
        df = pd.read_csv(
            args.mmseq, sep='\t',
            names=['query', 'target', 'identity', 'aln_length', 'mismatches',
                   'number_gaps', 'start_query', 'end_query', 'start_target',
                   'end_target', 'E-value', 'bit_score']
        )
        print(f"MMseqs2 table preview:\n{df.head()}\n")
        smorfs_to_exclude = set(df['query'].unique())  # Using set for faster lookup
        print(f"Number of smORFs to exclude: {len(smorfs_to_exclude)}\n")
        
    except Exception as e:
        print(f"Error reading MMseqs2 output file: {e}")
        exit(1)
    
    # Step 2 - Filtering FASTA file
    print("Step 2 - Filtering FASTA file")
    try:
        filtered_records = [
            record for record in SeqIO.parse(args.fasta, 'fasta')
            if record.id not in smorfs_to_exclude
        ]
        write_fasta(filtered_records, args.out_smorfs)
        print(f"Number of retained smORFs: {len(filtered_records)}\n")
    except Exception as e:
        print(f"Error filtering FASTA file: {e}")
        exit(1)

    print("COMPLETED")
    print("*" * 20)
