#!/bin/bash
# Written by HE, beginning 060125
# This script is meant to find lipoxygenase domains in hmmer output files
# It then prints the genome name, the number of lipoxygenases in teh genome, and the protein accessions of the lipoxygenases (with protein names separated by commas) in a tsv
# Usage: ./LPOfind.sh <genome_name> <domain>

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <genome_name> <domain_substring>"
    exit 1
fi

genome=$1
domain=$2
output_file="${genome}_LPOcount.tsv"

# Extract matching lines and process them
awk -v domain="$domain" -v genome="$genome" '
BEGIN {
    FS = "\t";
    count = 0;
    accessions = "";
    domains = "";
}
$0 !~ /^#/ && tolower($2) ~ tolower(domain) {
    count++;
    if (accessions == "") {
        accessions = $1;
        domains = $2;
    } else {
        accessions = accessions "," $1;
        domains = domains "," $2;
    }
}
END {
    print genome "\t" count "\t" accessions "\t" domains;
}' "${genome}_hmmer.out" > "$output_file"
