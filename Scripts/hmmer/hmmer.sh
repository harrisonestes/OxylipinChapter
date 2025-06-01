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