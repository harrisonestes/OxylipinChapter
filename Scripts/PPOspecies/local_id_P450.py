#!/usr/bin/env python3
# Written by HE, beginning 053025
# This script is used to identify unannotated P450s in faa files
# This will be done by aligning a consensus P450 sequence to the faa files, and keeping any hits with an e-value < 1e-20
# These positive hits will be written in a tsv as genome, number, protein IDs (with these in a single column and comma separated) to enable simple integration with the initial P450 annotation

import subprocess
import sys
import os
from pathlib import Path
import csv

# First, we will find all needed files and extract the genome names to be used in the first column
def find_genome_folders(base_dir):
    """Find all doubledom folders and extract genome names"""
    genome_data = []
    
    for folder_name in os.listdir(base_dir):
        if folder_name.endswith('_doubledom'):
            genome_name = folder_name.replace('_doubledom', '')
            folder_path = os.path.join(base_dir, folder_name)
            faa_file = os.path.join(folder_path, f"{genome_name}_An_peroxidaseonly.faa")
            
            if os.path.exists(faa_file):
                genome_data.append((genome_name, faa_file))
            else:
                print(f"Warning: {faa_file} does not exist. Skipping this genome.")
    return genome_data

def run_blastp(consensus_file, faa_file, evalue_cutoff):
    """Run blastp and return results"""
    cmd = [
        'blastp',
        '-query', consensus_file,
        '-subject', faa_file, 
        '-evalue', str(evalue_cutoff),
        '-outfmt', '6'
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error running blastp: {e}")
        return None

def parse_blast_results(blast_output):
    """Parse blastp output and return hits"""
    if not blast_output or blast_output.strip() == "":
        return []
    
    protein_ids = []
    lines = blast_output.strip().split('\n')
    
    for line in lines:
        if line.strip():
            columns = line.split('\t')
            protein_id = columns[1]
            protein_ids.append(protein_id)
    return protein_ids

def main():
    # Set up file paths
    base_dir = "/Users/hpestes/Programs/Collaborations/Dante/OxylipinChapter/Scripts/PPOspecies/PPO_output"
    consensus_file = "/Users/hpestes/Programs/Collaborations/Dante/OxylipinChapter/Scripts/PPOspecies/sequence_filter/consensus_CD_P450.faa"
    output_file = "unannotated_p450_results.tsv"
    evalue_cutoff = 1e-20
    
    # Find all genome folders
    genome_data = find_genome_folders(base_dir)
    print(f"Found {len(genome_data)} genomes to process")
    
    # Process each genome and collect results
    all_results = []
    
    for genome_name, faa_file in genome_data:
        print(f"Processing {genome_name}...")
        
        # Run blastp
        blast_output = run_blastp(consensus_file, faa_file, evalue_cutoff)
        
        # Parse results
        protein_ids = parse_blast_results(blast_output)
        
        # Format for TSV
        protein_count = len(protein_ids)
        protein_string = ','.join(protein_ids) if protein_ids else ""
        
        all_results.append((genome_name, protein_count, protein_string))
    
    # Write results to TSV using csv module
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        
        # Write header
        writer.writerow(['genome_name', 'protein_count', 'protein_ids'])
        
        # Write each result
        for genome_name, protein_count, protein_string in all_results:
            writer.writerow([genome_name, protein_count, protein_string])
    
    print(f"Results written to {output_file}")
    print(f"Processed {len(all_results)} genomes total")

if __name__ == "__main__":
    main()