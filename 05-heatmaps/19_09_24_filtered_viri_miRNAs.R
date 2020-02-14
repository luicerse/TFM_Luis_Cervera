################### VIRIDIPLANTAE MIRNAs##################################################
setwd('D:/TFM/LLUIS/Analysis')
library(tidyr) ## for na removal
library("stringr") ## To access letters within a given substring

diff_expressed_sRNA <- read.table('04-Tablas_exp_differencial/T4.txt', header = TRUE, stringsAsFactors = FALSE)

### Select double and triple stresses by the dot "."
double_triple_stress <- diff_expressed_sRNA[grepl(".", diff_expressed_sRNA$Stress, fixed = T),]
### Decreasing order for the Stress
double_triple_stress <- double_triple_stress[order(double_triple_stress$Stress, decreasing = FALSE),]

### Loading the miRBase collapsed table, to add the miRNAs IDs to the table of doubles and triples and filter by miRNA found
## Awk to process the txt file and substitute Us by Ts to match with the sRNA sequence in this table

#system('awk \'BEGIN {FS = "\t"; OFS = "\t"} {gsub("U","T",$2); print}\' miRBase22/miRBase22-collapsed.txt > miRBase22/miRBase22_U_to_T.txt')
miRBase_IDs <- read.table('02-miRBase22/miRBase22_U_to_T.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t')
# Este archivo miRBase_IDs ha sido modificado a fecha de 10/07/19, añadiendole las secuencias del miR1515 y miR6478 de melón, que son otros candidatos.


## Logical vector representing matches of miRNA ids with SNCMel ids - MATCH
## install package grr
#install.packages('grr')
library(grr)
# List of indexes for x = miRBase_IDs$seq and y = double_triple_stress$sRNA
miRNA_index_matches <- grr::matches(miRBase_IDs$seq, double_triple_stress$sRNA, all.x = FALSE, all.y = FALSE, list = FALSE)

### Join the miRNA ID and the SNC Mel matched in one dataframe.
# Creation of an empty dataframe with proper dimensions
matches_double_triple <- data.frame(matrix(data = NA, nrow = dim(miRNA_index_matches)[1], ncol = length(double_triple_stress) + 2))
index <- 1
for (indexes in t(miRNA_index_matches[2])){ #762 rows
  miRNA_ID <- miRBase_IDs[miRNA_index_matches[index,1],c(1,4)] # Information of the miRNA id (e.g. cme-miR390).
  seq_info <- double_triple_stress[indexes,] ## Information about LOGFC and differential expression data.
  matches_double_triple[index, ] <- cbind(miRNA_ID, seq_info)
  index <- index + 1
}
# Colnames
colnames(matches_double_triple) <- c('miRNA ID', 'miRNA seq', colnames(double_triple_stress))

########################################################################################################
################################### SINGLE STRESSES ####################################################
one_stresses <- diff_expressed_sRNA[!grepl(".", diff_expressed_sRNA$Stress, fixed = T),]

### We will not use the filtered table by highest CPM, just if you want so. All results are required.
### We need to find and compare the same miRNAs

miRNA_single <- grr::matches(miRBase_IDs$seq, one_stresses$sRNA, all.x = FALSE, all.y = FALSE, list = FALSE)
### Join the miRNA ID and the SNC Mel matched in one dataframe.
# Creation of an empty dataframe with proper dimensions
matches_single <- data.frame(matrix(data = NA, nrow = dim(miRNA_index_matches)[1], ncol = length(double_triple_stress) + 2))
index <- 1
for (indexes in t(miRNA_single[2])){ #762 rows
  miRNA_ID <- miRBase_IDs[miRNA_single[index,1],c(1,4)] # Information of the miRNA id (e.g. cme-miR390).
  seq_info <- one_stresses[indexes,] ## Information about LOGFC and differential expression data.
  matches_single[index, ] <- cbind(miRNA_ID, seq_info)
  index <- index + 1
}
# Colnames
colnames(matches_single) <- c('miRNA ID', 'miRNA seq', colnames(one_stresses))

################################# VIRIDIPLANTAE + CMELO ############################################################
######### NOW WE ARE INSERTING DATA FROM SIMPLE STRESSES TO THE DATASET CONTAINING DOUBLES AND THE TRIPLE STRESS.
library(DataCombine)

index <- 1

viri_miRNAs_filt <- matches_double_triple[c('miRNA ID', 'ID', 'sRNA', 'log2FC.E', 'logCPM.E', 'PValue.E', 'FDR.E', 'log2FC.D', 'PValue.D', 'FDR.D', 'log2FC.N', 'prob.N', 'Stress')]

for (miRNA in unique(matches_double_triple$`miRNA ID`)){
  ## store all the matched row (SIMPLE STRESSES)
  aa <- matches_single[matches_single$`miRNA ID` == miRNA & matches_single$Stress %in% c('C', 'D', 'SD', 'SA', 'MON'),
                       c('miRNA ID', 'ID', 'sRNA', 'log2FC.E', 'logCPM.E', 'PValue.E', 'FDR.E', 'log2FC.D', 'PValue.D', 'FDR.D', 'log2FC.N', 'prob.N', 'Stress')]
  print(aa)
  ## Insert the row aa into the dataframe "viri_miRNAs_filt"
  viri_miRNAs_filt <- InsertRow(viri_miRNAs_filt, aa)
  index <- index + 1
}

# Order by miRNA ID
viri_miRNAs_filt <- viri_miRNAs_filt[order(viri_miRNAs_filt$`miRNA ID`),]

################################################################################
######################### FILTERS #############################################
#viri_miRNAs_filt <- viri_miRNAs_filt[viri_miRNAs_filt$logCPM.E >= 5, ] ## logCPM (nº de cuentas)
viri_miRNAs_filt <- viri_miRNAs_filt[viri_miRNAs_filt$FDR.E <= 0.05, ] ## FDR!!!
viri_miRNAs_filt <- viri_miRNAs_filt[viri_miRNAs_filt$FDR.D <= 0.05, ] 
viri_miRNAs_filt <- viri_miRNAs_filt[viri_miRNAs_filt$prob.N >= 0.95, ] 
#labels_rows <- viri_miRNAs_filt[viri_miRNAs_filt$logCPM.E >= log(5,2),]
#viri_miRNAs_filt <- viri_miRNAs_filt[viri_miRNAs_filt$log2FC.E >= 1 | viri_miRNAs_filt$log2FC.E <= -1, ] ## LogFC

## Add a column with length of the sRNAs
list_length <- list()

for (miRNA in viri_miRNAs_filt$sRNA){
  list_length <- append(list_length, nchar(miRNA))
}

viri_miRNAs_filt <- cbind(viri_miRNAs_filt, as.matrix(list_length))
colnames(viri_miRNAs_filt) <- c(colnames(viri_miRNAs_filt)[1:13], 'length')

# last filter by length: only 21- and 22-nt long
viri_miRNAs_filt <- viri_miRNAs_filt[viri_miRNAs_filt$length == 21 | viri_miRNAs_filt$length == 22,]



####################################################################
##### HEATMAPS
#####################################################################

## devtools::install_github('cmap/morpheus.R')
## devtools::install_local(path = '~/Desktop/morpheus.R/') # USAR ESTE COMANDO PARA INSTALAR LA VERSIÓN CORRECTA
library(morpheus)

###################################################################################################################################
###################################################################################################################################
##########################  12/06/19-10/07/19 3. NEW HEATMAP WITH 5P/3P IDs #################################################################
####################################################################################################################################
####################################################################################################################################


#diff_exp_viri_miRNAs <- cbind(rownames(diff_exp_viri_miRNAs), diff_exp_viri_miRNAs)

viri_identifiers <- read.table('D:/TFM/LLUIS/Results/06-miRNA-cmelo-hairpin-alignment/viri_identifiers.txt', header = T, sep = '\t')

unique_stress <- unique(viri_miRNAs_filt$Stress)
diff_exp_viri_miRNAs <- data.frame(matrix(data = NA, nrow = length(unique(viri_miRNAs_filt$`miRNA ID`)), ncol = 1 + length(unique_stress)))
colnames(diff_exp_viri_miRNAs) <- c('5P/3P', unique_stress)

################## CHANGE THE VIRI MIRNA FILT ID, BY THE ONES WITH 5P/3P IN viri_identifiers #####################################
matches_ids_5p_3p_viri <- grr::matches(unique(viri_miRNAs_filt$`miRNA ID`), as.character(viri_identifiers$id.miRNA), all.y = FALSE, indexes = TRUE)
rownames(matches_ids_5p_3p_viri) <- NULL
## all.y = FALSE, since all the NAs in the viri_identifiers, i.e. not matching the diff expressed miRNAs, will be dropped, to have a list of 148 IDs.
### Up to HERE the code is OK
na_indeces <- which(is.na(matches_ids_5p_3p_viri[,2]))
c <- 1
r <- 1
for (viri_id in matches_ids_5p_3p_viri[,2]){
  if (is.na(viri_id)){
    viri_id2 <- matches_ids_5p_3p_viri[na_indeces[c],1]
    matches_ids_5p_3p_viri[na_indeces[c],2] <- as.character(unique(viri_miRNAs_filt$`miRNA ID`)[viri_id2])
    c <- c + 1
    r <- r + 1
  }
  if (!is.na(viri_id)){
    #print(viri_id)
    #print(matches_ids_5p_3p_viri[viri_id,2])
    matches_ids_5p_3p_viri[r,2] <- as.character(viri_identifiers$ID_5P_3P[viri_id])
    #print(matches_ids_5p_3p_viri[viri_id,2])
    r <- r + 1
  }

}

##### Create the original list of rownames (without the new 5P/3P label)
labels_5P_3P <- matches_ids_5p_3p_viri[order(matches_ids_5p_3p_viri[,2]),2]
rownames_no_5P_3P <- list()
for (id in labels_5P_3P){
  if (str_sub(id, -1, -1) == 'P'){  # If the last character is "P", we remove the label "-5P" or "-3P", i.e. the last three characters of the string.
    id <- str_sub(id, 0, -4) # here, up to the fourth to the last (remove the last 3 elements).
  }
  rownames_no_5P_3P <- append(rownames_no_5P_3P, id)
}

###
rownames(diff_exp_viri_miRNAs) <- unique(rownames_no_5P_3P) # Rownames

# There are 3 miRNAs with both labels 5P and 3P due to multiple alignments to precursors at two regions:
labels_no_P <- unique(labels_5P_3P)
counter <- 1
for (id in labels_no_P){
  if (str_sub(id, -1, -1) == 'P'){  # If the last character is "P", we remove the label "-5P" or "-3P", i.e. the last three characters of the string.
    id <- str_sub(id, 0, -4) # here, up to the fourth to the last (remove the last 3 elements).
  }
  labels_no_P[counter] <- id
  counter <- counter + 1
}

duplicated_labels <- labels_no_P[duplicated(labels_no_P)] # length 3

counter <- 1
for (name in labels_5P_3P){
  if (str_sub(name, 0, -4) %in% duplicated_labels){
    labels_5P_3P[counter] <- paste0(str_sub(name, 0, -3), '5P&3P') #str_sub removes the 2-symbol original label (e.g. 5P), and substitutes it by the new one 5P&3P
  }
  counter <- counter + 1
}

diff_exp_viri_miRNAs$`5P/3P` <- unique(labels_5P_3P) # 1st column
######

index <- 1
for (prova in unique(viri_miRNAs_filt$`miRNA ID`)){
  dataf_heatmap <- subset(viri_miRNAs_filt, viri_miRNAs_filt$`miRNA ID` == prova)
  selected_log2FC <- dataf_heatmap$log2FC.E
  selected_stress <- dataf_heatmap$Stress
  for (value in selected_stress){
    index_unique_stress <- match(value, unique_stress)
    index_selected_stress <- match(value, selected_stress)
    diff_exp_viri_miRNAs[index, index_unique_stress + 1] <- selected_log2FC[index_selected_stress]
    ### This is my new value of Log2FC to add to the specific place in the data frame.
  }
  index <- index + 1
}
diff_exp_viri_miRNAs[is.na(diff_exp_viri_miRNAs)] <- 0 ########## Change NA to 0s since morpheus does not accept NAs.

### Representation with morpheus (drop the rownames and substitute by the labels with 5P and 3P)
morpheus_diff_viri_miRNAs <- diff_exp_viri_miRNAs

#FIlter labels_rows with the coincident miRNAs in morpheus_diff_exp_viri_miRNAs
#merge(rownames(morpheus_diff_viri_miRNAs), labels_rows, by=)

rownames(morpheus_diff_viri_miRNAs) <- unique(labels_5P_3P)
morpheus_diff_viri_miRNAs$`5P/3P` <- NULL ## drop the extra column
morpheus::morpheus((morpheus_diff_viri_miRNAs), dendrogram = 'both', 
                   drawValues = TRUE)#rowAnnotations = as.type(labels_rows[1:97,5]))
morpheus_diff_viri_miRNAs <- round(morpheus_diff_viri_miRNAs,2)
# Write diff exp data to a csv file
write.csv(morpheus_diff_viri_miRNAs, file = 'D://TFM/LLUIS/Results/07-heatmaps/morpheus_df.csv')
###########################################################
######## 10/07/19 HEATMAPS POR FAMILIAS ###################

# Extract from miRBase ID the column indicating the miRNA family
colnames(miRBase_IDs)
selection_miRBase <- miRBase_IDs[,c(1,4)]
require(reshape)
require(ggplot2)
# Bucle gordo para sacar un heatmap por familia
for (miR in unique(selection_miRBase[,2])){ ## Unique family names (cmelo and viri in general altogether)
  print(miR)
  list_sRNAs <- subset(selection_miRBase, selection_miRBase$miRNA == miR)[,1]
  ## New subset of diff_exp_viri_miRNAs
  subset_diff_exp_viri <- subset(diff_exp_viri_miRNAs, rownames(diff_exp_viri_miRNAs) %in% list_sRNAs)
  if (dim(subset_diff_exp_viri)[1] == 0){
    next() ### if there is no miRNA from miRBase found in the viri differentially expressed set of miRNAs, the loop passes to the next in the list
  }
  ## Plot a heatmap using only the elements of the subset.
  subset_reshaped <- melt(subset_diff_exp_viri)

  heatmap_ggplot <- ggplot((data.frame(subset_reshaped)), 
                           aes(y = subset_reshaped$`5P/3P`, 
                               x = subset_reshaped$variable,
                               fill = subset_reshaped$value)) + 
    theme(text = element_text(size=15), 
          axis.text.x = element_text(angle=90, hjust=1)) +
    geom_tile() + scale_fill_gradient2(low = 'red', high = 'blue', mid = 'white', midpoint = 0)
  # Save the heatmap in the folder:
  if (!dir.exists('D://TFM/LLUIS/Results/07-heatmaps/heatmaps_miRNA_families/')){
    dir.create('D://TFM/LLUIS/Results/07-heatmaps/heatmaps_miRNA_families/')
  }
  png(filename = paste0('D://TFM/LLUIS/Results/07-heatmaps/heatmaps_miRNA_families/',miR, '_', 'family.png'))
  plot(heatmap_ggplot + labs(x = 'Stress condition', y = 'miRNA', fill = 'Scale of logFC', title = paste(miR,'family')))

  dev.off() # Close and save the plot as a png image

}



########################################  10/07/19 BOXPLOTS OF MIRNA FAMILIES ##########################################
require(car) # Special library for boxplots

# Bucle gordo para sacar un heatmap por familia
for (miR in unique(selection_miRBase[,2])){ ## Unique family names (cmelo and viri in general altogether)
  #print(miR)
  list_sRNAs <- subset(selection_miRBase, selection_miRBase$miRNA == miR)[,1]
  #print(list_sRNAs)
  ## New subset of diff_exp_viri_miRNAs
  subset_diff_exp_viri <- subset(diff_exp_viri_miRNAs, rownames(diff_exp_viri_miRNAs) %in% list_sRNAs)
  if (dim(subset_diff_exp_viri)[1] == 0){
    next() ### if there is no miRNA from miRBase found in the viri differentially expressed set of miRNAs, the loop passes to the next in the list
  }
  print('diff')
  rownames(subset_diff_exp_viri) <- subset_diff_exp_viri$`5P/3P`
  subset_diff_exp_viri$`5P/3P` <- NULL ## drop the extra column

  # Plot a boxplot using only the elements of the subset.

  if (!dir.exists('D://TFM/LLUIS/Results/07-heatmaps/boxplots_miRNA_families/')){
    dir.create('D://TFM/LLUIS/Results/07-heatmaps/boxplots_miRNA_families/')
  }
  png(filename = paste0('D://TFM/LLUIS/Results/07-heatmaps/boxplots_miRNA_families/',miR, '_', 'family.png'), width = 1000, height = 800)
  boxplot_ <- Boxplot(subset_diff_exp_viri, main = paste0(miR, '_', 'family'), 
                      xlab = 'Stress', ylab = 'Log-FC', cex.axis = 1.5,
                      cex.lab = 1.5,
                      id = list(n = 10, location = 'avoid'))
  ###boxplot_ <- boxplot(x = (subset_diff_exp_viri), notch = TRUE, plot = TRUE, border = c('green', 'red'))
  dev.off()
}
##########################################################################################################################
######################  TABLA MIRNAS RESPUESTA (BINARIA -- PRESENTE/NO PRESENTE):
####################### (1) BROAD
####################### (2) INTERMEDIATE
####################### (3) NARROW

################################## CLASIFICAR LOS MIRNA EN LAS TRES CLASES ########

response_range_df <- data.frame(matrix(data = NA, nrow = length(rownames(diff_exp_viri_miRNAs)), ncol = 11)) 

## Change name of columns and rows
colnames(response_range_df) <- colnames(diff_exp_viri_miRNAs)[2:12]
rownames(response_range_df) <- rownames(diff_exp_viri_miRNAs)

#Bucle
for (row in rownames(response_range_df)){
  for (stress in colnames(response_range_df)){
    if (diff_exp_viri_miRNAs[row, stress] != 0){
      response_range_df[row, stress] <- 1
    }
    else{
      response_range_df[row, stress] <- 0
    }
  }
}

## Add sum

sum_column <- apply(response_range_df, 1, function(x) sum(x))
response_range_df$'sum' <- sum_column


######################################################################################
##################### TABLA BINARIA PARA ESTRESES MULTIPLES ########################
## Primero creamos un diccionario con la clase asociada a cada miRNA presente en el objeto diff_exp_viri_miRNAs (102 miembros)
miRNA_list <- rownames(diff_exp_viri_miRNAs)
broad_miRNAs <- c('157', '408', '319', '396', '167', '6478', '156', '858', '393')
intermediate_miRNAs <- c('166', '169', '397', '171', '172', '168', '159', '398')
narrow_miRNAs <- c('165', '394', '1515', '160', '390', '162', '164', '395')

## Create the dictionary

miRNA_split <- str_split(miRNA_list, '-') ## the second substring of each element of the list contains the general miRNA family ID.

## Creation of the dataframe with 1/0
presence_absence_df <- data.frame(matrix(data = NA, nrow = length(miRNA_list), ncol = 5)) # Cols: Class, ID, simple st, double st, triple st.

#miRNA_list_class <- c()

for (i in 1:length(miRNA_split)){
  number_id <- gsub({'\\D'}, '', miRNA_split[[i]][2]) ## Remove all but digits

  if (number_id %in% broad_miRNAs){
    presence_absence_df[i,1] = 'broad'
    presence_absence_df[i,2] = miRNA_list[i]
    #miRNA_list_class <- append(miRNA_list_class, c(miRNA_list[i] = 'broad'))
  }
  if (number_id %in% intermediate_miRNAs){
    presence_absence_df[i,1] = 'intermediate'
    presence_absence_df[i,2] = miRNA_list[i]
    #miRNA_list_class <- append(miRNA_list_class, c(miRNA_list[i] = 'intermediate'))
  }
  if (number_id %in% narrow_miRNAs){
    presence_absence_df[i,1] = 'narrow'
    presence_absence_df[i,2] = miRNA_list[i]
    #miRNA_list_class <- append(miRNA_list_class, c(miRNA_list[i] = 'narrow'))
  }
}

## Find NAs to check which miRNAs are not yet classified

presence_absence_df <- drop_na(presence_absence_df, X1)

## Lists of stresses
simple_stresses <- c('C', 'D', 'MON', 'SA', 'SD')
double_stresses <- c('C.SA', 'C.SD', 'C.D', 'D.MON', 'D.SA')
triple_stress <- c('C.SA.SD')

## Change name of columns
colnames(presence_absence_df) <- c('Class', 'miRNA', 'Simples', 'Doubles', 'Triple')

## Simple stresses
for (i in 1:dim(presence_absence_df)[1]){
  row_values <- diff_exp_viri_miRNAs[i, simple_stresses]
  if (length(row_values[row_values != 0]) == 0){
    presence_absence_df[i,3] <- 0
  }
  else {
    presence_absence_df[i,3] <- 1
  }
}

## Double stresses
for (i in 1:dim(presence_absence_df)[1]){
  row_values <- diff_exp_viri_miRNAs[i, double_stresses]
  if (length(row_values[row_values != 0]) == 0){
    presence_absence_df[i,4] <- 0
  }
  else {
    presence_absence_df[i,4] <- 1
  }
}

## Triple stress
for (i in 1:dim(presence_absence_df)[1]){
  row_values <- diff_exp_viri_miRNAs[i, triple_stress]
  if (length(row_values[row_values != 0]) == 0){
    presence_absence_df[i,5] <- 0
  }
  else {
    presence_absence_df[i,5] <- 1
  }
}

## Add sum

sum_column <- apply(presence_absence_df[,c('Simples', 'Doubles', 'Triple')], 1, function(x) sum(x) - 1)
presence_absence_df$'sum' <- sum_column

## Write the table
write.table(x = presence_absence_df, file = 'D://TFM/LLUIS/Results/07-heatmaps/presence_absence_df.txt', sep = '\t',row.names = FALSE)

############## HCLUSTERING STUDY #######################################

library(stats)
library(dendextend)
library(circlize)

# Change classes to factors
presence_absence_df$Class <- as.factor(presence_absence_df$Class)


#Change Class to colour name:
color = c("green", "orange", "red")[morpheus_diff_viri_miRNAs$Class]

dist_df <- dist(morpheus_diff_viri_miRNAs[,], method = 'binary')
ward_D2 <- as.dendrogram(hclust(dist_df, method = 'ward.D2'))

#Changing colors to branches and width/type of branches
ward_D2 <- ward_D2 %>%
  color_branches(k=3) %>%
  set("branches_lwd", c(3,1,3)) %>%
  set("branches_lty", c(1,2,1)) %>%
  #set("branches_k_color", value = c("pink", "dark green", "light blue"), k =3) %>%
  set("labels_cex", 0.7) %>%
  set("nodes_pch", 19) %>%  # node point type
  set("nodes_cex", 0.5) %>%  # node point size
  set("nodes_col", "blue") # node point color

#Change labels to miRNA ID

list_miRNAs <- c()
list_colors <- c()
for (index in labels(ward_D2)){
  print(index)
  list_miRNAs <- c(list_miRNAs, presence_absence_df$miRNA[as.numeric(index)])
  list_colors <- c(list_colors, color[as.numeric(index)])
}

labels(ward_D2) <- list_miRNAs

# Apply color vector to the dendogram
labels_colors(ward_D2) <- list_colors

#Plot
par(mar = c(10,4,4,0))
plot(ward_D2, main = 'Cluster dendrogram for miRNA classes', 
     ylab = 'Height')
rect.dendrogram(ward_D2, k =3, border = 8, lty = 5, lwd = 2)
abline(h = 1.5, lty = 2, col = 2)

legend(34,1.2, legend = levels(presence_absence_df$Class), 
       fill = c("green", "orange", "red"), cex = 0.5, 
       ncol = 3, text.font = 2)

par(mar = c(0,0,0,0))
circlize_dendrogram(ward_D2,   
                    facing = c("outside", "inside"),
                    labels = TRUE,
                    labels_track_height = 0.4,
                    dend_track_height = 0.3)

legend(34,1.2, legend = levels(presence_absence_df$Class), 
       fill = c("green", "orange", "red"), cex = 0.5, 
       ncol = 3, text.font = 2)
kruskal.test(x=presence_absence_df$sum, g=presence_absence_df$Class)
#SIGNIFICATIVO <0.01
####################################################################################
########################################## CALCULATE MEAN CONNECTIVITY  ######################################################
mean_broad <- subset(presence_absence_df, presence_absence_df$Class == 'broad')
mean(mean_broad$sum)
mean_intermediate <- subset(presence_absence_df, presence_absence_df$Class == 'intermediate')
mean(mean_intermediate$sum)
mean_narrow <- subset(presence_absence_df, presence_absence_df$Class == 'narrow')
mean(mean_narrow$sum)

## Connectivity broad == 1.77
## Connectivity intermediate == 1.46
## Connectivity narrow == 1.10
############################################################################################################################
################################################# ANÁLISIS ESTADÍSTICO #######################################################################################
require(car) # Function of boxplots

#The Mann-Whitney U test for comparing independent data samples: the nonparametric version of the Student t-test.
#The Wilcoxon signed-rank test for comparing paired data samples: the nonparametric version of the paired Student t-test.
#The Kruskal-Wallis H and Friedman tests for comparing more than two data samples: the nonparametric version of the ANOVA and repeated measures ANOVA tests

# Punto 1: Comparar dobles estreses vs simple (valores exp diferencial)
morpheus_diff_viri_miRNAs
# First we order the dataframe by class and id to perform any statistical test
presence_absence_df <- presence_absence_df[order(presence_absence_df$Class, presence_absence_df$miRNA), ]

shapiro.test(t(morpheus_diff_viri_miRNAs)[,1])  # Samples (miRNAs across stresses) are not normal

###### Punto 2

png(filename = paste0('D://TFM/LLUIS/Results/07-heatmaps/boxplot_connectivity.png'), width = 1000, height = 900)
Boxplot(formula = sum ~ Class, data = presence_absence_df, id.method = "y", id = list(n=5, location = 'avoid'))
dev.off()

## test estadístico para distribuciones no normales (aquí, 3 grupos o clases)

kruskal.test(x = presence_absence_df$sum, g = presence_absence_df$Class) # Kruskal Wallis
wilcox.test(x = presence_absence_df$sum, g = presence_absence_df$Class) # Mann-Whitney-Wilcoxon

############### Els dos donen pvalue molt baix -- Resultats significatius #################
# Comparamos con ANOVA (asume normalidad) ##
# Generally, anova is considered to be relatively robust to violations of normality and homogeneity, especially when the sample sizes are equal or nearlhy equal.”
# “ANOVA is generally robust to moderate violations as long as the sample sizes in each group are equal and not unreasonably small (<5)"
# ANOVA se descarta por no ser distribuciones normales.


## Summary of means, sd, members of each group
summary_connectivity <- with(presence_absence_df, by(sum, Class,
              function(x) c(mean = mean(x),
                            sd = sd(x),
                            n = sum(!is.na(x)))
))

## Multiple comparisons of means were performed in the Python script


##################     BARS FOR MEAN CONNECTIVITY  ##########################################3
means_conn <- c(summary_connectivity$broad['mean'],
                summary_connectivity$intermediate['mean'],
                summary_connectivity$narrow['mean'])

sd_conn <- c(summary_connectivity$broad['sd'],
             summary_connectivity$intermediate['sd'],
             summary_connectivity$narrow['sd'])

#Reset par configuration to default

dev.off()
bars <- barplot(height = means_conn,
                main = "Mean connectivity per miRNA class", xlab = 'Class', ylab = "Connectivity",
                ylim = c(0,3),
                col = c('green', 'orange', 'red'), 
                axes = TRUE, axisnames = TRUE, 
                names.arg = c('Broad', 'Intermediate', 'Narrow'))
arrows(bars, means_conn - sd_conn,
       bars, means_conn + sd_conn, angle=90, code=3, col = 'black')
legend("topright", legend = levels(presence_absence_df$Class), 
              fill = c("green", "orange", "red"), cex = 0.7, 
              ncol = 3, text.font = 2)

## For betweenness centrality calculated in Python
means_bet <- c(0.0003236735986084274, 0.000237982944045371, 0.000129297797801827)

sd_bet <- c(0.000126092272573645, 0.00016728369114877174, 0.00016505141586172054)
bars <- barplot(height = means_bet,
                main = "Betweenness centrality per miRNA class", xlab = 'Class', ylab = "Betweenness",
                ylim = c(0,0.0005),
                col = c('green', 'orange', 'red'), 
                axes = TRUE, axisnames = TRUE, 
                names.arg = c('Broad', 'Intermediate', 'Narrow'))
arrows(bars, means_bet - sd_bet,
       bars, means_bet + sd_bet, angle=90, code=3, col = 'black')
legend("topright", legend = levels(presence_absence_df$Class), 
       fill = c("green", "orange", "red"), cex = 0.7, 
       ncol = 3, text.font = 2)
################################ THAT'S ALL FOLKS!!! ################################################################
################################# SEE U VERY SOON (HOPEFULLY) #######################################################


