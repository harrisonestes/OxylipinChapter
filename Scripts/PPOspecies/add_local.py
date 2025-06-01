#!/usr/bin/env python3
# Written by HE, beginning 053025
# The purpose of this file is to combine the outputs of doubledom.sh and local_id_P450.py
# It will match the genome names from column 1 and add the second column value from local_id_P450.py to the first, thereby accounting for proteins with misidentified p450s as outcases and providing an accurate protein count

import csv
from collections import defaultdict

def merge_p450_files():
    """Merge annotated and unannotated P450 results by genome"""
    
    # Dictionary to store combined data: {genome: [count, protein_list]}
    combined_data = defaultdict(lambda: [0, []])
    
    # Read annotated P450s (complete proteins) - NO HEADER
    print("Reading annotated P450s...")
    with open('all_completeproteins.tsv', 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        # No next(reader) - no header to skip
        
        for row in reader:
            if len(row) >= 3:
                genome = row[0]
                count = int(row[1]) if row[1].isdigit() else 0
                protein_ids = row[2].split(',') if row[2] else []
                
                # Add to combined data
                combined_data[genome][0] += count
                combined_data[genome][1].extend(protein_ids)
    
    # Read unannotated P450s - HAS HEADER
    print("Reading unannotated P450s...")
    with open('unannotated_p450_results.tsv', 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)  # Skip header line
        
        for row in reader:
            if len(row) >= 3:
                genome = row[0]
                count = int(row[1]) if row[1].isdigit() else 0
                protein_ids = row[2].split(',') if row[2] else []
                
                # Add to combined data
                combined_data[genome][0] += count
                combined_data[genome][1].extend(protein_ids)
    
    # Write combined results
    print("Writing combined results...")
    with open('combined_p450_results.tsv', 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        
        # Write header
        writer.writerow(['genome_name', 'total_p450_count', 'all_p450_protein_ids'])
        
        # Write each genome's combined data
        for genome in sorted(combined_data.keys()):
            total_count = combined_data[genome][0]
            all_proteins = combined_data[genome][1]
            
            # Remove empty strings and join
            clean_proteins = [p for p in all_proteins if p.strip()]
            protein_string = ','.join(clean_proteins)
            
            writer.writerow([genome, total_count, protein_string])
    
    print(f"Combined results written to combined_p450_results.tsv")
    print(f"Processed {len(combined_data)} genomes")

if __name__ == "__main__":
    merge_p450_files()