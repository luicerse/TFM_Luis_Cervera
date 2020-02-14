### Extraction and alignment of all sRNA sequences from the libraries

setwd('D://TFM/LLUIS/Analysis/')
library(dplyr)

### 1. IMPORT THE DIRNAMES OF ALL FASTA FILES TO BE READ
#system('find -L ./librerias_con_cuentas_sRNA/0-trimClean-RFAM > tmp_file_names.txt')
#all_files <- read.table('tmp_file_names.txt')
#all_fasta <- as.list(all_files[grepl("fasta$", t(all_files)),])
#system('rm tmp_file_names.txt')

#############################################################################################################################
# 1. CREATE FASTA FILES FOR BOWTIE WITH THE SRNA SEQUENCES

diff_expressed_sRNA <- read.table('04-Tablas_exp_differencial/T4.txt', header = TRUE, stringsAsFactors = FALSE)
subset_sRNA <- diff_expressed_sRNA[,c(1,2,23)]
#miRBase_IDs <- read.table('miRBase22/miRBase22_U_to_T.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t')

# Create a fasta format


#List of unique SNCMel with the unique sRNA seq:
library(dplyr)
unique_sRNAs <- unique(subset_sRNA[,c(1,2)])

indeces_sRNA <- grr::matches((unique_sRNAs[,1]), (subset_sRNA[,1]))

################################################ LOOP FOR THE CREATION OF THE FASTA FILE #########################################
initial <- 1
while (initial <= dim(subset_sRNA)[1]){

  index_match <- indeces_sRNA[initial,1]
  stresses_index <- indeces_sRNA[indeces_sRNA[,1] == index_match,2] ### Devuelve los indices de los estreses para cada sRNA
  stresses <- subset_sRNA[stresses_index,3] ### Devuelve los estreses para cada sRNA
  #rownames(unique_sRNAs) <- c()
  cmelo_id <- unique_sRNAs[index_match, 1] ## Devuelve el identificador de melon correspondiente
  cmelo_id <- gsub(", ", "_", toString(c(cmelo_id,stresses))) ## Junta id cmelo con los estreses en un mismo identificador
  cmelo_id <- gsub(", ", "", toString(c('>',cmelo_id))) ## Añade el ">" para iniciar el ID en fasta

  sRNA_sequence <- unique_sRNAs[index_match,2] ## Selecciona la secuencia del sRNA
  fasta_sequence <- gsub(", ", "\n", toString(c(cmelo_id,sRNA_sequence))) ## Creates the whole fasta sequence (ID+seq)

  ## Para actualizar los indices en el loop
  length_stress <- length(stresses_index)
  initial <- initial + length_stress ## Pasa al siguiente sRNA unico en la lista de subset_sRNA con todos los sRNAs y los estreses
  ##################################################################################################

  print(paste(index_match, initial))

  ### Write to the fasta file
  write(fasta_sequence, '06-sRNA-alignment/all_sRNA_new.fasta', append = TRUE)
}

################################# FINISHED!!!! ################################################################
######################## ALINEAMIENTO CON BOWTIE ###############################################################

### Desde este working directory: /home/lab205/Desktop/LLUIS/Analysis/miRBase22/precursor_fasta_files
#bowtie2-build hairpins_cmelo.fa all_sRNA_hairpins
# bowtie2 --very-sensitive -f ~/Desktop/LLUIS/Analysis/librerias_con_cuentas_sRNA/0-trimClean-RFAM/all_sRNA_libraries/sRNA_fasta.fasta -x all_sRNA_hairpins

