#!/bin/bash
# Written by HE, beginning 053025
# The goal of this script is to identify any proteins containing both an An_peroxidase and p450 domain
# This will be written in a tsv output of genome name, protein count, and protein IDs (in a comma list but single column)
# Additionally, it will find all proteins containing an An_peroxidase without a p450 domain and output those in a separate tsv file containing three columns: genome name, protein count, and protein IDs (in a comma list but single column)
# Finally, it will use the comma list of proteins without p450 domains to reference the .faa genome files and extract the sequences of those proteins, writing them to a .faa file for each genome
# Usage: ./doubledom_edited.sh <genome> <domain1> <domain2>

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <genome> <domain1> <domain2>"
    exit 1
fi

# Set up variables
genome=$1
domain1=$2
domain2=$3

# Prepare the environment and initialize output files
gunzip "${genome}_prot.faa.gz"

hmmer_file="${genome}_hmmer.out" #input file 1
faa_file="${genome}_prot.faa" #input file 2 

# Check if required input files exist
if [ ! -f "$hmmer_file" ]; then
    echo "Error: $hmmer_file not found"
    exit 1
fi

if [ ! -f "$faa_file" ]; then
    echo "Error: $faa_file not found" 
    exit 1
fi

output1="${genome}_completeproteins.tsv" # output file containing proteins with both domains
output2="${genome}_${domain1}only.tsv" # output file containing proteins with only domain1
output3="${genome}_${domain1}only.faa" # output file containing sequences of proteins with only domain1

# Find domain 1
proteins_with_domain1=$(grep "$domain1" "$hmmer_file" | grep -v "^Target_Name" | cut -d$'\t' -f1)

# Find domain 2
proteins_with_domain2=$(grep "$domain2" "$hmmer_file" | grep -v "^Target_Name" | cut -d$'\t' -f1)

# Find proteins with both domains
proteins_both=$(comm -12 <(echo "$proteins_with_domain1" | sort) <(echo "$proteins_with_domain2" | sort))

# Find proteins with only domain1
proteins_only_domain1=$(comm -23 <(echo "$proteins_with_domain1" | sort) <(echo "$proteins_with_domain2" | sort))

# Convert newline-separated to comma-separated and count
if [ -n "$proteins_both" ]; then
    count_both=$(echo "$proteins_both" | wc -l)
    list_both=$(echo "$proteins_both" | tr '\n' ',' | sed 's/,$//')  # Remove trailing comma
else
    count_both=0
    list_both=""
fi

# Convert newline-separated to comma-separated and count for domain1 only
if [ -n "$proteins_only_domain1" ]; then
    count_only_domain1=$(echo "$proteins_only_domain1" | wc -l)
    list_only_domain1=$(echo "$proteins_only_domain1" | tr '\n' ',' | sed 's/,$//')  # Remove trailing comma
else
    count_only_domain1=0
    list_only_domain1=""
fi

# Write proteins with both domains
echo -e "$genome\t$count_both\t$list_both" > "$output1"

# Write proteins with only domain1
echo -e "$genome\t$count_only_domain1\t$list_only_domain1" > "$output2"

# Extract sequences for proteins with only domain1
if [ -n "$proteins_only_domain1" ]; then
    while IFS= read -r protein_id; do
        awk -v id="$protein_id" '
        /^>/ { 
            if ($1 == ">"id) {
                print_seq = 1
            } else {
                print_seq = 0
            }
        }
        print_seq { print }
        ' "$faa_file" >> "$output3"
    done <<< "$proteins_only_domain1"
fi

# Move output files to a folder and compress
mkdir ${genome}_doubledom
mv "$output1" "$output2" "$output3" "${genome}_doubledom/"
tar -czf "${genome}_doubledom.tar.gz" "${genome}_doubledom/"