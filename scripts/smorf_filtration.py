#!/usr/bin/python
# -*- coding: utf-8 -*-

### Script filters table with MiPepid predictions ###

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import argparse

def write_fasta(df, col_id, col_seq, path):
    """
    Create a FASTA file from a DataFrame.
    """
    try:
        records = []
        for _, row in df.iterrows():
            records.append(SeqRecord(seq=Seq(str(row[col_seq])), id=str(row[col_id]), description=''))
        SeqIO.write(records, path, "fasta")
        print(f"The sequences were saved to {path}")
    except Exception as e:
        print(f"Error while writing FASTA file: {e}")

def translate_dna(sequence):
    """
    Translate a DNA sequence into a protein sequence.
    """
    try:
        seq_obj = Seq(sequence)
        protein = seq_obj.translate(stop_symbol='')  # Remove stop symbols
        return str(protein)
    except Exception as e:
        print(f"Error translating sequence: {sequence} -> {e}")
        return ""

if __name__ == "__main__":
    print("*" * 20)
    print("START")
    
    # Argument parsing
    parser = argparse.ArgumentParser(description='Filter potentially coding smORFs using MiPepid predictions.')
    parser.add_argument("--table", type=str, required=True, help='Path to CSV file with MiPepid predictions.')
    parser.add_argument("--out_bed", type=str, required=True, help='Path to output BED file.')
    parser.add_argument("--out_table", type=str, required=True, help='Path to output CSV table file.')
    parser.add_argument("--out_nucl", type=str, required=True, help='Path to output nucleotide FASTA file.')
    parser.add_argument("--out_prot", type=str, required=True, help='Path to output protein FASTA file.')
    args = parser.parse_args()
    
    # Read arguments
    table_path = args.table
    bed_output = args.out_bed
    table_out = args.out_table
    nucl_out = args.out_nucl
    prot_out = args.out_prot
    LENGTH_LIMIT = 213  # Minimum length for filtering
    
    print("Step 1 - Reading MiPepid output")
    try:
        df = pd.read_csv(table_path)
        df = df[df['classification'] == 'coding']  # Filter only 'coding' classifications
        df['Length,nt'] = df['end_at'] - df['start_at']  # Calculate length
        df = df[df['Length,nt'] > LENGTH_LIMIT]  # Filter by length
        
        print(f"Table preview:\n{df.head()}\n")
        print(f"Statistics:\nThe number of predicted coding smORFs: {df['sORF_ID'].nunique()}\n")
    except Exception as e:
        print(f"Error reading table: {e}")
        exit(1)
    
    print("\nStep 2 - Creating temporary BED file")
    try:
        bed = df[['transcript_DNA_sequence_ID', 'start_at', 'end_at', 'sORF_ID']].copy()
        bed['score'] = 1000
        bed['strand'] = '+'
        bed['start_at'] = bed['start_at'] - 1  # Convert to 0-based indexing
        bed.to_csv(bed_output, index=False, header=False, sep='\t')
        print(f"BED file saved to {bed_output}")
    except Exception as e:
        print(f"Error creating BED file: {e}")
        exit(1)
    
    print("\nStep 3 - Creating FASTA file with nucleotide sequences")
    write_fasta(df, 'sORF_ID', 'sORF_seq', nucl_out)
    print()
    
    print("\nStep 4 - Creating FASTA file with protein sequences")
    try:
        df['Protein'] = df['sORF_seq'].apply(translate_dna)  # Translate DNA to protein
        write_fasta(df, 'sORF_ID', 'Protein', prot_out)
    except Exception as e:
        print(f"Error creating protein FASTA file: {e}")
    print()
    
    print("\nStep 5 - Saving modified MiPepid table")
    try:
        df.to_csv(table_out, index=False)
        print(f"Modified table saved to {table_out}")
    except Exception as e:
        print(f"Error saving modified table: {e}")
    
    print("COMPLETED")
