# Oxylipin Chapter Reproducible Notebook

## Begun by HE, 042325

## Goal

The goal of this project is to create multiple phylogenetic reference points from a species perspective and a protein domain perspective relative to oxylipins. This will require the creation of three trees:

1. A species tree counting the number of PPOs overlaid. ("PPOSpecies" folder)
2. A species tree counting the number of lipoxygenases. ("LPOSpecies" folder)
3. A phylogenetic tree of all PPO proteins from A. nidulans, A. flavus, A. oryzae, A. fumigatus, A. terreus, and A. niger. ("PPO_poteintree" folder)

## Intro

PPOs are defined as a dioxygenase domain fused with a Cytochrome P450. Pfam annotates this dioxygenase as "An_peroxidase" while the P450 is simply "p450".

Lipoxygenases are denoted by Pfam as "Lipoxygenases".

## Used programs

For many of the "mining" scripts, basic bash commands. The anaconda package manager was used in addition to this for phylogenetic tree programs. The environment had the following programs and versions:

- iTOL web server, accessed 062225
- MAFFT v7.526
- trimAl v1.4.rev22
- IQ-TREE multicore version 2.3.6 for MacOS ARM 64-bit

Additionally, HMMER was run against all available genomes, with one representative genome per species being chosen for downstream analyses.

## Species tree and used genomes

The used species tree for overlay was created by Nickles et al 2023.

The genomes used for this study are one representative per species, selected on basis of the highest percentage of BUSCOs present in that species.

## Running HMMER

The following script was used to run HMMER on the UW-Madison CHTC:

Submit script: (/Scripts/hmmer/hmmer.sub)

```bash
# Specify the HTCondor Universe (vanilla is the default)
universe = vanilla
#
log = $(name)_hmmer.log
error = $(name)_hmmer.err
#
# Specify your executable (single binary or a script that runs several
executable = hmmer.sh
arguments = $(name)
#the hmmer outfiles are hugeeeeeeeeeee
#output = $(NewProcess)_hmmer.out
#
#make sure it has gluster
requirements = (Target.HasCHTCStaging == true) && (OpSysMajorVer == 7)
#
# Specify that HTCondor should transfer files to and from the
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files =/home/groups/keller_group/FungalGenomes/Programs/hmmer-3.3.2_compiled.tar.gz,/home/groups/keller_group/FungalGenomes/$(name).tar.gz
#
# Tell HTCondor what amount of compute resources
request_cpus = 1
request_memory = 3GB
request_disk = 5GB
#
queue name from /home/groups/keller_group/FungalGenomes/QueueFile/genomes.txt
```

Script: (/Scripts/hmmer/hmmer.sh)

```bash
#!/bin/sh
#
cp /staging/gnickles/Pfam-A.hmm.gz .
#
tar -xzvf hmmer-3.3.2_compiled.tar.gz
tar -xzvf "$1".tar.gz
gunzip Pfam-A.hmm.gz
gunzip ./"$1"/"$1"_prot.faa.gz
#
export PATH=$(pwd)/hmmer-3.3.2/src:$PATH
#
hmmpress Pfam-A.hmm
hmmsearch -E 1e-4 --cpu 1 --tblout "$1"_hmmer.out Pfam-A.hmm ./"$1"/"$1"_prot.faa
#
#Grant: added the sed statement to remove any lines that have --- in them
cat "$1"_hmmer.out | sed '/^#/d' | awk 'NR>3 {print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' | sort -u -k1,2> temp
echo -e "Target_Name\tDomain\tDomain_Accession\tE-value\tscore\tbias" > topLine
cat topLine temp > "$1"_hmmer.out
rm temp topLine
#
rm *.fa *.gz *.tsv *.txt *h3f *h3i *h3m *h3p *hmm
gzip "$1"_hmmer.out
```

## Counting PPO domains in representative genomes

PPO domains are comprised of a dioxygenase and cytochrome P450, annotated by Pfam as "An_peroxidase" and "p450" respectively. These were counted in 1221 species prior to mapping on a species tree.

### Handling exceptions

Though an exception, the Cytochrome p450 is sometimes missed by HMMER alignment, while the dioxygenase is consistently annotated correctly — including PpoC in Aspergillus fumigatus Af293. This will require an alternative approach than simply counting a HMMER output for co-occurring domains in a single protein.

To accomodate this, a script will harvest the sequences for all proteins containing an annotation with a dioxygenase and without a P450. The NCBI Conserved Domains consensus sequence for a P450, which correctly identified A. fumigatus will be used to align with these proteins. This domain, CYP_LDS-like_C (cd20612), is described as "C-terminal cytochrome P450 domain of linoleate diol synthase and similar cytochrome P450s". This term aligned with the PpoC P450 region with an e-value of 1.67e-127, providing high confidence in this method of identification.

### Harvesting PPO domains

The following script was submitted to the UW-Madison CHTC: (/Scripts/PPOspecies/doubledom_edited.sh)

```bash
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
```

I used the following submit script: (/Scripts/PPOspecies/doubledom_edited.sub)

```bash
# Specify the HTCondor Universe (vanilla is the default)
universe = vanilla
#
# Log Out Error 
log = $(genome)_PPOfind.log
error = $(genome)_PPOfind.err
#output = $(genome)_PPOfind.out
#
# Executable and arguments if any
executable = doubledom_edited.sh
arguments = $(genome) An_peroxidase p450
#
requirements = (Target.HasCHTCStaging == true)
#
# Inputs/Outputs
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = /home/groups/keller_group/FungalGenomes/$(genome)/$(genome)_prot.faa.gz, /home/groups/keller_group/FungalGenomes/$(genome)/$(genome)_hmmer.out
transfer_output_files = $(genome)_doubledom.tar.gz
#
# Tell HTCondor what amount of compute resources
request_cpus = 1
request_memory = 400MB
request_disk = 400MB
# 
# Increasing the number of files that can be run at one time
## Allows it to also go to HTCondor Pools
+WantFlocking = true
## Allows it to also go to OS Pool
+WantGlideIn = true
#
# Query files, use the testing or full file and comment out the other one
#queue genome from /home/hpestes/complete_projects/DanteOxylipin/testlist.txt
queue genome from /home/hpestes/complete_projects/DanteOxylipin/singlegenome_perspecies.txt
```

The queuefile can be found at /Scripts/queuefile_053125.txt.

### Local Identification of ppoC and other atypical P450's

Following completion of the completion of standard domain identification, the files were downloaded locally. Next, the atypical P450 outcase must be identified. To do this, the consensus alignment sequence of the NCBI Conserved Domains CYP_LDS-like_C was saved as /Scripts/PPOspecies/sequence_filter/consensus_CD_P450.faa. This will be blasted against all protein sequences output in the $(genome)_An_peroxidaseonly.faa.

#### Verification of method

To verify this method, the Conserved Domains-identified sequence coordinates were found. The first 692 amino acids, not containing the P450 was blasted against the consensus sequence. This resulted in no similarity being found. When blasted against the entirety of ppoC, it found the p450 successfully, at an e-value of 6e-79 when considering the entire 1121 amino acid protein. Because of this, the cutoff for a blast hit will be 1e-20 for positive identification to account for varying length of proteins.

#### Scripts and output

The following script was written to run BLASTP locally and create an output to combine with the proteins found by the initial method: (/Scripts/PPOspecies/local_id_P450.py)

```python
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
```

### Processing and combining outputs

The output files were combined using bash and a script. The following command was run inside the PPO_output folder:

```bash
cat */GCA_*_completeproteins.tsv > all_completeproteins.tsv
```

Next, a python script was written to combine them: (/Scripts/PPOspecies/add_local.py)

```python
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
```

Due to the handling methods of files in the local method, any genomes only containing non-p450 An_peroxidase domains were removed from the list. To overcome this, the following command was run to find what genomes were not present:

```bash
comm -23 <(sort /Users/hpestes/Programs/Collaborations/Dante/OxylipinChapter/Scripts/queuefile_053125.txt) <(tail -n +2 combined_p450_results.tsv | cut -f1 | sort)
```

Next, the mapfile used for the queuefile was used as a key to replace the genome IDs with their names:

```bash
awk 'NR==FNR{map[$2]=$1; next} {print (map[$1] ? map[$1] : $1), $2}' /Users/hpestes/Programs/Collaborations/Dante/OxylipinChapter/Scripts/fungalmap_oneentry.tsv combined_p450_results.tsv > combined_p450_results_named.tsv
```

These were then overlayed by iTOL to the species tree. The tree was edited for publication in Adobe Illustrator.

## Counting Lipoxygenase domains in representative genomes

Lipoxygenase domains are annotated by Pfam/HMMER correctly, so a basic script will be required to harvest these, saved at LPOharvest.sh and LPOharvest.sub.

```bash
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
}' "${genome}_hmmer_output.txt" > "$output_file"
```

Submit script to HTCondor:

```bash
# Specify the HTCondor Universe (vanilla is the default)
universe = vanilla
#
# Log Out Error 
log = $(genome)_LPOfind.log
error = $(genome)_LPOfind.err
#output = $(genome)_PPOfind.out
#
# Executable and arguments if any
executable = LPOfind.sh
arguments = $(genome) Lipoxygenase
#
requirements = (Target.HasCHTCStaging == true)
#
# Inputs/Outputs
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = /home/groups/keller_group/FungalGenomes/$(genome)/$(genome)_hmmer.out
transfer_output_files = $(genome)_LPOcount.tsv
#
# Tell HTCondor what amount of compute resources
request_cpus = 1
request_memory = 300MB
request_disk = 300MB
# 
# Increasing the number of files that can be run at one time
## Allows it to also go to HTCondor Pools
+WantFlocking = true
## Allows it to also go to OS Pool
+WantGlideIn = true
#
# Query files, use the testing or full file and comment out the other one
queue genome from /home/hpestes/complete_projects/DanteOxylipin/testlist.txt
#queue genome from /home/hpestes/complete_projects/DanteOxylipin/singlegenome_perspecies.txt
```

These outputs were moved to a local destination where the previous command to annotate the genome names was applied. Final overlays were completed via iTOL and Adobe Illustrator.

## Aspergillus PPO tree

Representatives of the Aspergillus genus were selected. These proteins were chosen from the previous harvesting step of PPOs, and trimmed based on InterPro domain identifications to where the protein begins at the beginning of the dioxygenase and ending at the exact end of the p450. These sequences were saved to AspPPOseq.faa.

nidulans    "CBF69709.1,CBF85912.1,CBF76249.1"
flavus  "RAQ42372.1,RAQ46208.1,RAQ53268.1,RAQ45225.1"
oryzae  "OOO07742.1,OOO07799.1,OOO09892.1,OOO11534.1,OOO10188.1"
fumigatus   "EAL84400.2,EAL89712.1,EAL92371.2"
terreus "GES58595.1,GES59806.1,GES60941.1,GES61775.1"
niger   "CAK37756.1,CAK38852.1,CAK41061.1"

Next, mafft was used for alignment of these proteins:

```bash
mafft AspPPOseq.faa > AspPPOseq.aln
```

Then, trimal was used to trim the alignment via the strict gappyout option:

```bash
trimal -gappyout -in AspPPOseq.aln -out AspPPOseq_trimmed.aln
```

Finally, IQTree2 with ModelFinder was run to find the best tree.

```bash
iqtree2 -B 1000 -T 14 -s AspPPOseq_trimmed.aln
```

The output tree was annotated in Adobe Illustrator.
