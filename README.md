# TFM_Luis_Cervera
These folder contains all the scripts employed for the Master's Thesis of the MSc in Bioinformatics and Computational Biology at the Universidad Autónoma de Madrid (UAM).
They are distributed as follows:
  
  ## 1. Filtrado_RFAM_y_graficos_cuentas: 
   Contains the initial script that removes any sRNA aligning to the RFAM database using bowtie. This database is composed of sRNAs that cannot perform a function in stress tolerance in plants, such as transfer RNA or ribosomal RNA. It uses as input the cleaned sRNA libraries filtered by quality obtained from the FASTQ files.
   
   ## 2. RPMs: 
   Normalises the read counts associated to the RFAM-filtered sRNA libraries.
   
   ## 3. PCA:
   To correlate biological replicates and different samples we performed a Principal Component Analysis in R, contained in one single script. 
   
   ## 4. miRNA_hairpin_precursor_alignment:
   Starting from a filtered fasta file processed in order to contain only melon miRNAs, mature miRNAs were aligned to miRBase melon hairpins using bowtie aligner. Then, with a Python script, the SAM file was processed so that the regions from which miRNAs were generated could be identified.
   
   ## 5. Heatmaps:
   A simple code in R was devised using the "Morpheus" package from the Broad Institute to represent the direction and the magnitude of the differentially expressed sRNAs extracted at 11 days post-treatment ("Time 4" or "T4"). The master table of sRNAs used as input was obtained with scripts within the directory "7_Differential_expression".
   
   ## 6. sRNA_precursor_alignment:
   With diverse Python routines we have implemented a pipeline to combine read counts and differential expression data to filter relevant sRNA reads, according also to probabilistic filters (FDR).
   
   ## 7. Differential_expression:
   This folder contains one script written in R to compute the differentially expressed sRNAs according to three different methods: EdgeR, DeSeq2 and NOISeq.
   
   ## 8. Conectividad_miRNAs:
   
   ## 9. antagonismo:
