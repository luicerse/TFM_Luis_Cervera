#/bin/bash

cd /home/lab205/Desktop/LLUIS/Analysis/06-sRNA_alignment

bowtie-build fasta_1515_6478.fa bowtie_index/prec_1515_6478

bowtie -f -a -v 0 -l 15 -a -p 2 bowtie_index/prec_1515_6478 ~/Desktop/LLUIS/Analysis/06-sRNA_alignment/all_sRNA.fasta -S ~/Desktop/LLUIS/Analysis/06-sRNA_alignment/bowtie_miR1515_6478.sam
