##################################### ANÁLISIS DE LAS LIBRERÍAS ##############################################################
##################################### 1er paso: FILTRADO POR RFAM ###########################################################

## En este script se parte de los archivos fasta en los que se han eliminado los adaptadores, y que provienen del filtrado por calidad y procesado de los fastq,
## que a su vez vienen del secuenciador. 
## Vamos a filtrar estas secuencias mediante RFAM, es decir, alinearlas con aquellas secuencias que sean rRNA/tRNA/snRNA/snoRNA y que, por tanto, son irrelevantes
## para nuestro estudio. La secuencias que alineen serán descartadas, obviamente. Representan alrededor de un 3% por librería (+/- 1%).


setwd("D://TFM/LLUIS/Analysis/")

if (!dir.exists("01-Rfam")){  # Creates the directory if it does not exist.
  dir.create("01-Rfam")
}



######################
# LINUX  (YA ESTÁ HECHO)
######################
# # Descargamos las 2686 secuencias de ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/ correspondientes a la versión más actual RFAM
# 
# # Descromprimimos y creamos un solo archivo fasta con todas las secuencias.
# gunzip *
# # Todas empiezan por RF?
# ls | grep -c "^RF"
# cat RF* > rfam.fasta
# # Creamos fichero tabular con cabecera y secuencia
# fasta_formatter -t -i rfam.fasta -o rfam1line.tsv
# # Nos quedamos con aquellas secuencias que contengan rRNA/tRNA/snRNA/snoRNA
# 
# grep -E 'rRNA|tRNA|snRNA|snoRNA|ribosomal|transfer RNA|small nuclear RNA|small nucleolar RNA' rfam1line.tsv | sed 's/#//g' > rfamsel.tsv  

################## Cargamos fichero con todas las secuencias en RFAM, y nos quedamos con aquellos que sean rRNA/tRNA/snRNA/snoRNA 
### * MIRAR PASOS RFAM BOWTIE

rfam <- read.table("01-Rfam/rfamsel.tsv", sep = "\t", header = F, stringsAsFactors = F, quote="\"")
rfam <- rfam[!duplicated(rfam$V1),]
rfam$V1 <- paste0(">", rfam$V1)

if (!dir.exists("01-Rfam/alignment")){  # Creates the directory if it does not exist.
  dir.create("01-Rfam/alignment")
}
write.table(rfam,"01-Rfam/rfamsel.fasta",sep="\n", row.names = F, col.names = F, quote = F)



setwd("D://TFM/LLUIS/Analysis/01-Rfam/alignment/")


exp <- list.files("D://TFM/LLUIS/Analysis/01-Rfam/Rfam_cleaned/", recursive = T)
exp <- exp[!grepl("T4/CONTROL/CONTROL-1.fasta", exp)]
exp <- exp[!grepl("T4/X", exp)]

dtexp <- data.frame(file=exp, 
                    time=unlist(lapply(strsplit(exp ,"/",fixed=TRUE),"[",7)),
                    sample=unlist(lapply(strsplit(exp ,"/",fixed=TRUE),"[",8)), 
                    nsample=gsub(".fasta", "",unlist(lapply(strsplit(exp ,"/",fixed=TRUE),"[",9))),
                    percentage=0, stringsAsFactors = F)

#dtexp <- dtexp[dtexp$time == 'T4',] # Me quedo solo con las de T4

system("bowtie-build ../rfamsel.fasta rfam") # Indexar con bowtie


# usa=1
counts= 0
for(usa in 1:dim(dtexp)[1]){
  e <- dtexp$file[usa]
  # e="/home/lab205/Lab/Cmelo/0-raw_input/outTrim/T4/A/A-2.fasta"
  name <- gsub(".fasta", "", unlist(strsplit(e, "/"))[9])
  time <- gsub(".fasta", "", unlist(strsplit(e, "/"))[7])
  
  #En cada iteración borramos el contenido de la carpeta de alineamientos, ya que al hacer diferentes 
  #tiempos no cambian los nombres y se vuelven a gastar los mismos.
  system("ls -1 | grep -v ^rfam | xargs rm -f") # Aquí borramos
  system(paste0("bowtie rfam -f -v 0 -l 20 --best -k 1 --norc ", e, " -S ", name, ".sam --un ", name, "-rfam.fasta"))
  system(paste0("samtools view -bT ../rfamsel.fasta ", name, ".sam > ", name,".bam"))
  system(paste0("samtools sort ", name, ".bam ", name, ".sorted"))
  system(paste0("samtools index ", name, ".sorted.bam ", name, ".sorted.bai"))
  system(paste0("samtools idxstats ", name, ".sorted.bam > ", name, "-bamstats.tsv")) 
  
  # Sacamos secuencias sin mapear, 
  system(paste0("samtools view -b -f 4 ", name, ".sorted.bam > ", name, ".unmapped.bam"))
  
  # Creamos si no existe carpeta destino
  
  if(!dir.exists(paste0("~/Desktop/LLUIS/Analysis/01-Rfam/Rfam_cleaned/"))){
    dir.create(paste0("~/Desktop/LLUIS/Analysis/01-Rfam/Rfam_cleaned/"))
  }
  
  if(!dir.exists(paste0("~/Desktop/LLUIS/Analysis/01-Rfam/Rfam_cleaned/", unlist(strsplit(e, "/"))[7]))){
    dir.create(paste0("~/Desktop/LLUIS/Analysis/01-Rfam/Rfam_cleaned/", unlist(strsplit(e, "/"))[7]))
  }
  
  if(!dir.exists(paste0("~/Desktop/LLUIS/Analysis/01-Rfam/Rfam_cleaned/", unlist(strsplit(e, "/"))[7], "/", unlist(strsplit(e, "/"))[8]))){
    dir.create(paste0("~/Desktop/LLUIS/Analysis/01-Rfam/Rfam_cleaned/", unlist(strsplit(e, "/"))[7], "/", unlist(strsplit(e, "/"))[8]))
  }
  
  if(!dir.exists(paste0("~/Desktop/LLUIS/Results/01-cuentas_librerias/", time))){
    dir.create(paste0("~/Desktop/LLUIS/Results/01-cuentas_librerias/", time))
  }
  
  dirfasta <- paste0("~/Desktop/LLUIS/Analysis/01-Rfam/Rfam_cleaned/", unlist(strsplit(e, "/"))[7], "/", unlist(strsplit(e, "/"))[8], "/", name, ".fasta")
  
  # sudo apt-get install seqtk    
  system(paste0("samtools bam2fq ", name, ".unmapped.bam | seqtk seq -A - > ", dirfasta))
  
  # Mapeadas samtools view -b -F 4 file.bam > mapped.bam
  # Sacar las secuencias
  # samtools view unmapped.bam | cut -f 10 | sort | uniq -c | sort -nr > .sorted_reads.txt
  system(paste0("samtools view ", name, ".unmapped.bam | cut -f 10 | sort | uniq -c | sort -nr | awk '{print $2\"\\t\"$1}' > ", name, ".sorted_reads.txt "))
  system(paste0("cp *sorted_reads.txt ~/Desktop/LLUIS/Results/01-cuentas_librerias/", time))
  
  # Cargamos tabla 
  tabla <- read.table(paste0(name, ".sorted_reads.txt"), sep="\t", header = F, stringsAsFactors = F)
  # Creamos tabla frecuencia, según estrés, cuentas, sacamos total, y calculamos también en RPM
  dt <- data.frame(nt=20:25, stringsAsFactors = F)
  dt$counts = 0
  dt$total = 0
  dt$RPM = 0
  dt$percentage = 0
  dt$sample = name
  dt$time = unlist(strsplit(e, "/"))[7]
  # Sumamos secuencias por tamaño
  for(ii in 1:dim(dt)[1]){
    # ii=1
    subii <- tabla[nchar(tabla$V1) == dt$nt[ii],] # Filtrar la tabla por aquellas secuencias de longitud 20,21,22,23,24 o 25.
    dt$counts[ii] <- sum(subii$V2)
  }
  sumtotal <- sum(dt$counts)
  dt$total <- sumtotal
  
  # Calculamos RPM
  for(ii in 1:dim(dt)[1]){
    # ii=1
    dt$RPM[ii] <- (dt$counts[ii]*1000000)/sumtotal
    dt$percentage[ii] <- round((dt$counts[ii]/sumtotal)*100, digits = 2)
  }
  
  
  
  bams <-read.table(paste0(name, "-bamstats.tsv"), sep = "\t", header = F, stringsAsFactors = F, quote="\"") 
  colnames(bams) <- c("ID", "length", "mapped","unmapped")
  bams <- bams[!bams$ID == "*",]
  bams <- bams[bams$mapped > 0,]
  bams$description <- "-"
  #Buscamos ID con nombre completo
  for(i in 1:dim(bams)[1]){
    # i=1
    bams$description[i] <- rfam$V1[grep(bams$ID[i], rfam$V1)]
  }
  
  #Cuentas del archivo
  cuentas <- as.numeric(system(paste0("grep -c '^>' ", e), intern = T))
  bams$percentage <- (bams$mapped/cuentas)*100
  print(paste0(dt$time[1], " ", name," ", round(sum(bams$percentage), digits = 3), " %"))
  bams$sample <- name
  bams$time = unlist(strsplit(e, "/"))[7]
  dtexp$percentage[usa] <- round(sum(bams$percentage), digits = 3)
  if(counts == 0){
    bamsfull <- bams
    dtfull <- dt
    counts = counts + 1
  }else{
    bamsfull <- rbind(bamsfull, bams)
    dtfull <- rbind(dtfull, dt)
    
  }
  
  #Limpiamos archivos generados
  system("ls -1 | grep -v ^rfam | xargs rm -f") # Borramos archivos creados
  
}

write.table(dtexp,"../ALL-percentage-sum.tsv",sep="\t", dec = ",", row.names = F, quote = F)

write.table(bamsfull,"../ALL-bams-sum.tsv",sep="\t", dec=",", row.names = F, quote = F)

write.table(dtfull,"../ALL-freq-sum.tsv",sep="\t", dec=",", row.names = F, quote = F)


dfreq <- dtfull

dfreq$condition <- unlist(lapply(strsplit(dfreq$sample ,"-",fixed=TRUE),"[",1))


se <- function(x) sqrt(var(x)/length(x))  # This is the function for calculating the standard deviation



colortable <- data.frame(s=c("T0", "CONTROL", "C", "D", "SA", "SD", "MON", "HSVd", "A", "C.SA", "C.SD", "C.D", "C.SA.SD", "D.MON", "D.SA"),
                         color=c("#E69F00", "#BDBDBD","#0000FF","#CC8877","#66FFFF","#000000","#339900","#FFEE00", "#FF00FF", "#56B4E9","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00"),
                         stringsAsFactors = F) # color table assigning one color to one sample specifically.


ucondi <- unique(dfreq$condition) # These are the unique 14 stress conditions, including T0.
controladd <- 0

for(u in 1:length(ucondi)){
  # u=1
  usubx <- dfreq[dfreq$condition == ucondi[u],]
  uniqtime <- unique(usubx$time)
  for(t in uniqtime){
    # t="T1"
    usub <- usubx[usubx$time == t,]
    
    dt <- data.frame(Sample=rep(ucondi[u], 6), nt=20:25, stringsAsFactors = F) # Here a data frame is created, by replicating 6 times the stress conditions, this is,
                                                                              # one set for each nt length (20,21,22,23,24,25; 6 numbers).
    dt$mean.RPM <- 0
    dt$se.RPM <- 0
    dt$mean.Relative <- 0
    dt$se.Relative <- 0
    dt$Time  <- t
    dt$Color <- colortable$color[colortable$s == ucondi[u]]
    
    for(ui in 1:dim(dt)[1]){
      # ui=1
      #sacamos los valores según tamaños
      dt$mean.RPM[ui] <- round(mean(usub$RPM[usub$nt == dt$nt[ui]]), digits = 2)
      dt$se.RPM[ui] <- round(se(usub$RPM[usub$nt == dt$nt[ui]]), digits = 2)
      dt$mean.Relative[ui] <- round(mean(usub$percentage[usub$nt == dt$nt[ui]]), digits = 2)
      dt$se.Relative[ui] <- round(se(usub$percentage[usub$nt == dt$nt[ui]]), digits = 2)
    }
    
    if(controladd == 0){
      dtfinal <- dt
      controladd = controladd + 1
    }else{
      dtfinal <- rbind(dtfinal, dt)
    }
    
    
    
    
  }
  
}




write.table(dtfinal,"../ALL-freq-sum-FINAL.tsv",sep="\t", dec=",", row.names = F, quote = F)
write.table(dtfinal,"~/Desktop/LLUIS/Results/04-graficos_cuentas_librerias/ALL-freq-sum-FINAL.tsv",sep="\t", dec=",", row.names = F, quote = F)


require(ggplot2)
require(ggExtra)
# Grouped
#dtfinal$Sample[dtfinal$Sample == "T0"] <- "CONTROL"
dtfinal$Sample <- factor(dtfinal$Sample, levels=c("T0", "CONTROL", "C", "D", "SA", "SD", 
                                                  "MON", "HSVd", "A", "C.SA", "C.SD", 
                                                  "C.D", "C.SA.SD", "D.MON", "D.SA"))

#dtfinal$Color <- factor(dtfinal$Color, levels=c("#BDBDBD", "#0000FF", "#CC8877", "#66FFFF", "#000000", "#339900", "#FFEE00", "#FF00FF"))

## Gráficas con RPMs
options(scipen = 10)
ggplot(dtfinal, aes(x=nt, weight=mean.RPM, ymin=mean.RPM-se.RPM, ymax=mean.RPM+se.RPM, fill=Sample)) +
  geom_bar(position=position_dodge(), aes(y=mean.RPM), stat="identity") +
  geom_errorbar (position=position_dodge(width=0.9), colour="black")+
  scale_x_discrete(breaks=c(20,21,22,23,24,25),limits=20:25)+
  scale_fill_manual(values=c("darkgrey", "lightgrey","blue1","darkgoldenrod1","red","antiquewhite1",
                             "salmon4","pink", "green", "cyan","deepskyblue",
                             "cornflowerblue","darkorchid","gold1","darkorange"))+
  theme_bw()+removeGrid( x = T, y = T  )+ xlab("Read size")+ ylab("Reads Per Million (RPM)")+
  geom_segment(aes(x = 19, y = 0, xend = 26, yend = 0), linetype=2, colour = "darkgray", size = 0.4 )+
  facet_wrap( ~ Time, ncol=4)+
  theme(legend.position="bottom", axis.text.x = element_text(angle=65, vjust=0.6), plot.title = element_text(hjust = 0.5)) #panel.grid.major = element_line(colour = "gray",linetype="dashed",size=0.4)
ggsave("~/Desktop/LLUIS/Results/04-graficos_cuentas_librerias/02-RPMs_small.png", height = 7, width = 20, dpi = 300)

# Gráficas con porcentajes
ggplot(dtfinal, aes(x=nt, weight=mean.Relative, ymin=mean.Relative-se.Relative, ymax=mean.Relative+se.Relative, fill=Sample)) +
  geom_bar      (position=position_dodge(), aes(y=mean.Relative), stat="identity") +
  geom_errorbar (position=position_dodge(width=0.9), colour="black")+
  scale_x_discrete(breaks=c(20,21,22,23,24,25),limits=20:25)+
  scale_fill_manual(values=c("darkgrey", "lightgrey","blue1","darkgoldenrod1","red","antiquewhite1",
                             "salmon4","pink", "green", "cyan","deepskyblue",
                             "cornflowerblue","darkorchid","gold1","darkorange"))+
  theme_bw()+removeGrid( x = T, y = T  )+ xlab("Read size")+ ylab("Relative percentage (%)")+
  geom_segment(aes(x = 19, y = 0, xend = 26, yend = 0), linetype=2, colour = "darkgray", size = 0.4 )+
  facet_wrap( ~ Time, ncol=4)+
  theme(legend.position="bottom", axis.text.x = element_text(angle=65, vjust=0.6), plot.title = element_text(hjust = 0.5)) #panel.grid.major = element_line(colour = "gray",linetype="dashed",size=0.4)
ggsave("~/Desktop/LLUIS/Results/04-graficos_cuentas_librerias/04-Percentage_small.png", height = 7, width = 20, dpi = 300)



# Para tabla latex

dfreq
dfreq2 <- dfreq[,c(7,6,1,2,3)]


write.table(dfreq2,"../ALL-freq-for-latex.tsv",sep="\t", dec=".", row.names = F, quote = F)
