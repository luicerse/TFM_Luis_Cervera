#/bin/bash

cd /home/lab205/Desktop/LLUIS/Analysis/02-miRBase22/precursor_fasta_files

bowtie-build hairpins_cmelo.fa all_sRNA_hairpins

bowtie -f -a -v 0 -l 15 -a -p 2 all_sRNA_hairpins ~/Desktop/LLUIS/Analysis/06-sRNA_alignment/all_sRNA.fasta -S ~/Desktop/LLUIS/Analysis/06-sRNA_alignment/bowtie_sRNAs_19_09_25.sam



