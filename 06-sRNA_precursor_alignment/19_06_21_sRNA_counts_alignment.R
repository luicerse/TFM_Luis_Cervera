############ ALINEAMIENTO DE PRI-MIRNAs con cuentas de SRNA ###################################
#### TABLAS ONLY ############ 
setwd('D://TFM/LLUIS/Analysis/')
new_table <- read.table('~/Desktop/LLUIS/Results/sRNA_hairpin_alignment/19_06_21_new_table_sRNA_precursor_alignment.csv', sep = '\t', header = T)

##############################################################################################################################################
diff_expressed_sRNA <- read.table('Tablas_exp_differencial/T4.txt', header = TRUE, stringsAsFactors = FALSE)
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
  cmelo_id <- unique_sRNAs[index_match, 1] ## Devuelve el identificador de melón correspondiente
  cmelo_id <- gsub(", ", "_", toString(c(cmelo_id,stresses))) ## Junta id cmelo con los estreses en un mismo identificador
  cmelo_id <- gsub(", ", "", toString(c('>',cmelo_id))) ## Añade el ">" para iniciar el ID en fasta
  
  sRNA_sequence <- unique_sRNAs[index_match,2] ## Selecciona la secuencia del sRNA
  fasta_sequence <- gsub(", ", "\n", toString(c(cmelo_id,sRNA_sequence))) ## Creates the whole fasta sequence (ID+seq)
  
  ## Para actualizar los índices en el loop
  length_stress <- length(stresses_index) 
  initial <- initial + length_stress ## Pasa al siguiente sRNA único en la lista de subset_sRNA con todos los sRNAs y los estreses
  ########################## ESTO PROVOCA QUE ESTÉ MAL!!!!!!!!!!!!!!!!!!!!!!!! #####################################################
  
  print(paste(index_match, initial))
  
  ### Write to the fasta file
  write(fasta_sequence, './librerias_con_cuentas_sRNA/0-trimClean-RFAM/all_sRNA_libraries/prova_final_sRNA.fasta', append = TRUE)
}