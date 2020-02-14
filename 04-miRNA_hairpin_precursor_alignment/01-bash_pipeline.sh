#!/bin/bash

#Filtra archivo de hairpins.fa a solo aquellas secuencias que sean de melón

sed 's/>/-\n>/g' ~/Desktop/LLUIS/Analysis/02-miRBase22/precursor_fasta_files/hairpin.fa | sed -n '/>cme/,/^-$/p' | awk '!/^-/' > ~/Desktop/LLUIS/Analysis/02-miRBase22/precursor_fasta_files/hairpins_cmelo.fa

#Now the bowtie alignment: it is end-to-end, performed using version 1 with the following options: -f, -a (all alignments), -vo -l 15 -p2 -s, after building the indeces. First only with cmelo miRBase miRNAs. The fasta files of hairpins have been translated, the "U"s to "T"s,as well as the miRNAs from miRBase--> now in 03-miRNA-hairpin-alignment folder

bowtie-build ~/Desktop/LLUIS/Analysis/03-miRNA-hairpin-alignment/hairpins_cmelo_U_to_T.fa bowtie/hairpins_cmelo

bowtie -f -a -v 0 -l 15 -p 2 -S bowtie/hairpins_cmelo cmelo_U_to_T.fasta 


#Alineamiento de viri y cmelo miRNAs, el archivo viri_miRNAs.fasta se obtuvo filtrado el archivo original de miRBase22 con grep (cmelo y viri) y luego cut.

bowtie-build ~/Desktop/LLUIS/Analysis/03-miRNA-hairpin-alignment/hairpins_cmelo_U_to_T.fa bowtie/hairpins_cmelo

bowtie -f -a -v 0 -l 15 -p 2 -S bowtie/hairpins_cmelo viri_miRNAs.fasta 
