### 02-RPMs

setwd('D://TFM/LLUIS/Results/01-cuentas_librerias/')

exp <- list.files(path = '.', pattern = 'sorted_reads', recursive = TRUE)

#Create new tables
if (!dir.exists('~/Desktop/LLUIS/Results/02-cuentas_con_RPMs/')){
  dir.create('~/Desktop/LLUIS/Results/02-cuentas_con_RPMs/')
}
if (!dir.exists('~/Desktop/LLUIS/Results/02-cuentas_con_RPMs/T0')){
  dir.create('~/Desktop/LLUIS/Results/02-cuentas_con_RPMs/T0')
}
if (!dir.exists('~/Desktop/LLUIS/Results/02-cuentas_con_RPMs/T1')){
  dir.create('~/Desktop/LLUIS/Results/02-cuentas_con_RPMs/T1')
}
if (!dir.exists('~/Desktop/LLUIS/Results/02-cuentas_con_RPMs/T2')){
  dir.create('~/Desktop/LLUIS/Results/02-cuentas_con_RPMs/T2')
}
if (!dir.exists('~/Desktop/LLUIS/Results/02-cuentas_con_RPMs/T4')){
  dir.create('~/Desktop/LLUIS/Results/02-cuentas_con_RPMs/T4')
}

for (file_ in exp){
  print(file_)
  table <- read.table(file_, stringsAsFactors = F, header = F)
  total <- sum(table[,2])
  colnames(table) <- c('seq', 'counts')
  time <- unlist(strsplit(file_, '/'))[1]
  name <- unlist(strsplit(file_, '/'))[2]
  table$RPM <- apply(data.frame(table[,2]), 2, function(x) x*1000000/total)
  table$percentage <- round(apply(data.frame(table[,2]), 2, function(x) x*100/total), digits = 2)
  print(paste0('Writing table to ', '../02-cuentas_con_RPMs/', time,'/', name))
  write.table(table, file = paste0('../02-cuentas_con_RPMs/', time,'/', name), sep = '\t', row.names = FALSE, quote = FALSE)

}

######################### 10/09/19 GET ONE TABLE PER SAMPLE ######################################
# There will be 8 columns: sRNA seq, Control-2,3,4,5 and S-1, S-2, S-3

exp <- list.files(path = '../02-cuentas_con_RPMs/', pattern = 'sorted_reads', recursive = T)
setwd('../02-cuentas_con_RPMs/')
times <- c('T1', 'T2', 'T4')
stresses <- c('A', 'C', 'D', 'HSVd', 'MON', 'SA', 'SD', 'C.SA', 'C.SD', 'C.D', 'C.SA.SD', 'D.MON', 'D.SA')


### CREATE THE CONTROL TABLE ####

for (time in times){
  controles <- exp[grepl(paste0(time,'/', 'CONTROL', '-'), exp)]
  
  counter <- 1
  for (archivo in controles){
    
    assign(paste0('tabla_', counter), read.table(archivo, header = TRUE) ) # With assign, we can assign a dataframe to a different name object within the for loop.
    counter <- counter + 1 
  }
  
  # Now we merge the replicate dataframes from the same stress condition, according to the sRNA sequence
  if (length(controles) == 3){
    
    # This form of merge is an inner merge, this is, keeping only sequences appearing in the three tables.
    merged_table <- merge(merge(tabla_1, tabla_2, by = 'seq'), tabla_3, by = 'seq') 
    colnames(merged_table) <- c('sRNA', 'counts1', 'CONTROL-1', 'perc1', 'counts2', 'CONTROL-2', 'perc2',
                                'counts3', 'CONTROL-3', 'perc3')
    #New selection of columns
    merged_table <- cbind.data.frame(merged_table$sRNA, merged_table$`CONTROL-1`, merged_table$`CONTROL-2`, merged_table$`CONTROL-3`)
    
  }
  
  if (length(controles) == 4){
    merged_table <- merge(merge(tabla_1, tabla_2, by = 'seq'), tabla_3, by = 'seq')
    colnames(merged_table) <- c('seq', 'counts2', 'CONTROL-2', 'perc2', 'counts3', 'CONTROL-3', 'perc3',
                                'counts4', 'CONTROL-4', 'perc4')
    merged_table <- merge(merged_table, tabla_4, by = 'seq') # Two steps of merge to avoid warning messages about duplicated columns.
    colnames(merged_table) <- c('sRNA', 'counts2', 'CONTROL-2', 'perc2', 'counts3', 'CONTROL-3', 'perc3',
                                'counts4', 'CONTROL-4', 'perc4', 'counts5', 'CONTROL-5', 'perc5')
    #New selection of columns
    merged_table <- cbind.data.frame(merged_table$sRNA, merged_table$`CONTROL-2`, merged_table$`CONTROL-3`, merged_table$`CONTROL-4`, merged_table$`CONTROL-5`)
  }
  
  ## Ahora quiero guardar las tablas de controles en T1, T2 y T4.
  if (time == 'T1'){
    control_T1 <- merged_table
    print('Control T1 table has been created')
  }
  if (time == 'T2'){
    control_T2 <- merged_table
    print('Control T2 table has been created')
  }
  if (time == 'T4'){
    control_T4 <- merged_table
    print('Control T4 table has been created')
  }
}

# Change colnames of control tables
colnames(control_T1) <- c('sRNA', 'CONTROL-1', 'CONTROL-2', 'CONTROL-3')
colnames(control_T2) <- c('sRNA', 'CONTROL-1', 'CONTROL-2', 'CONTROL-3')
colnames(control_T4) <- c('sRNA', 'CONTROL-2', 'CONTROL-3', 'CONTROL-4', 'CONTROL-5')

##############################################################################
### FUSION OF STRESS WITH CONTROL ###

# 1. Create directories

if (!dir.exists('../03-archivos_cuentas_txt/T1')){
  dir.create('../03-archivos_cuentas_txt/T1')
}

if (!dir.exists('../03-archivos_cuentas_txt/T2')){
  dir.create('../03-archivos_cuentas_txt/T2')
}

if (!dir.exists('../03-archivos_cuentas_txt/T4')){
  dir.create('../03-archivos_cuentas_txt/T4')
}

# 2. create stress tables

for (time in times){
  
  # Selection of the correct control library
  
  if (time == 'T1'){
    control_library <- control_T1
  }
  
  if (time == 'T2'){
    control_library <- control_T2
  }
  
  if (time == 'T4'){
    control_library <- control_T4
  }
    
  for (stress in stresses){
    archivos <- exp[grepl(paste0(time,'/', stress, '-'), exp)]
   
    counter <- 1
    for (archivo in archivos){
      assign(paste0('tabla_', counter), read.table(archivo, header = TRUE)) # With assign, we can assign a dataframe to a different name object within the for loop.
      counter <- counter + 1
      #print(counter)
    }
    # Now we merge the replicate dataframes from the same stress condition, according to the sRNA sequence
    if (counter == 1){
      print(paste('There are no', stress, 'data frames', 'at time', time))
      next()
    }
    if (length(archivos == 3)){
      merged_table <- merge(merge(tabla_1, tabla_2, by = 'seq'), tabla_3, by = 'seq')
      merged_table <- cbind.data.frame(merged_table$seq, merged_table$RPM.x, merged_table$RPM.y, merged_table$RPM)
      colnames(merged_table) <- c('sRNA', paste0(stress, '-1'), paste0(stress, '-2'), paste0(stress, '-3'))
      print(paste(time, stress, 'data frames have been joined'))
      
      name_file <- paste0('../03-archivos_cuentas_txt/', time, '/', stress, '.txt')
      final_dataframe <- merge(control_library, merged_table, by = 'sRNA', all = TRUE)
      final_dataframe[is.na(final_dataframe)] <- 0 # Replace NAs by 0s
      write.table(final_dataframe, file = name_file, row.names = FALSE, quote = FALSE, sep = '\t')
      next()
    }
    if (length(archivos == 2)){
      merged_table <- merge.data.frame(tabla_1, tabla_2, by = 'seq')
      merged_table <- cbind.data.frame(merged_table$seq, merged_table$RPM.x, merged_table$RPM.y)
      colnames(merged_table) <- c('sRNA', paste0(stress, '-2'), paste0(stress, '-3'))
      print(paste(time, stress, 'data frames have been joined'))
      
      name_file <- paste0('../03-archivos_cuentas_txt/', time, '/', stress, '.txt')
      final_dataframe <- merge(control_library, merged_table, by = 'sRNA', all = TRUE)
      final_dataframe[is.na(final_dataframe)] <- 0 # Replace NAs by 0s
      write.table(final_dataframe, file = name_file, row.names = FALSE, quote = FALSE, sep = '\t')
                     
                     
    }

  }
}

