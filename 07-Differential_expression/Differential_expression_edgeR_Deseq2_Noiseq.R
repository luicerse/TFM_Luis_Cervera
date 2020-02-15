

# En este documento se describe como se realiza un analisis de expresion diferencial para datos de RNA-Seq.
# 
# En primer lugar se han colapsado las secuencias para cada libreria y se han dado un formato tabular, de modo que se tiene la secuencia y el numero de apariciones para cada una. Se ha renombrado todas las carpetas y ficheros de trabajo segun una nomenclatura.
# 
# El an·lisis de expresion diferencial se va a realizar haciendo estudios pareados entre control y estrÈs.
# 
# 1. Primero creamos la tabla con las cuentas.



###BASH
# sudo apt-get install aptitude
# sudo apt-get install libcurl4-openssl-dev
# sudo apt-get install libxml2-dev
# sudo apt-get install libcairo2-dev
# sudo apt-get install libssl-dev
# 
# rownames(subset(as.data.frame(installed.packages()),Priority=="base"))

#BiocManager::install('limma')
#BiocManager::install('EDASeq')
#BiocManager::install('NOISeq')
#BiocManager::install('ffpe')
#BiocManager::install('DESeq2')
#BiocManager::install('edgeR')
#BiocManager::install('statmod')
#BiocManager::install('genefilter')
#BiocManager::install('calibrate')
#BiocManager::install('RColorBrewer')
#BiocManager::install('BiocParallel')
#BiocManager::install('geneplotter')
#BiocManager::install('gridExtra')
#BiocManager::install('lattice')
#BiocManager::install('fdrtool')
#BiocManager::install('gplots')
#BiocManager::install('ggplot2')
#BiocManager::install('reshape')
#BiocManager::install('venn')


library(limma)
library(EDASeq)
library(NOISeq)
library(ffpe)
library(DESeq2)
library(edgeR)
library(statmod)
library(genefilter)
library(calibrate)
library(RColorBrewer)
library(BiocParallel)
library(geneplotter)
library(gridExtra)
library(lattice)
library(fdrtool)
library(gplots)
library(ggplot2)
library(reshape)
library(venn)

houseDIR="D:/TFM"

setwd(houseDIR)
# setwd("~/Desktop/Lab/01DEN/data20-25groups")


############################################# DISTRIBUCIONES
distr <- read.table("0-distribution/distribution20-25.txt", sep = "\t", stringsAsFactors = F, quote = NULL)

colnames(distr) <- c("Size", "Count", "Time", "Sample")

distr$Size <- as.numeric(as.character(distr$Size))


################ NORMALIZAMOS A RPM
# distr$Count <- distr$Count/1000000


## ponemos porcentajes de secuencias de 0 a 100, la desviacion estandar sera en base a estos numeros

head(distr)

distr$stress <- "-"
distr$replica <- 0


for(j in 1:dim(distr)[1]){
  distr$stress[j] <- strsplit(distr$Sample[j], "-")[[1]][1]
  distr$replica[j] <- strsplit(distr$Sample[j], "-")[[1]][2]
}

########## Normalizamos para cada estrÈs y dentro de cada estrÈs para cada rÈplica

tiemposX <- unique(distr$Time)

#Crearemos una nueva tabla donde iremos metiendo los datos normalizados
#contador par controlar la tabla nueva 
count=0

for(time in tiemposX){
  # time="T1"
  #extraemos aquellos para cada tiempo
  
  casosTiempo <- distr[distr$Time == time,]
  
  #Vemos cuantos estreses tenemos para este tiempo
  estresesX <- unique(casosTiempo$stress)

  for(estr in estresesX){
    # estr="A"
    #extraemos aquellos para cada estrÈs
    
    casosEstres <- casosTiempo[casosTiempo$stress == estr,]
    
    #sacamos replicas para ese estrÈs
    replicasX <- unique(casosEstres$replica)
    for(rep in replicasX){
      # rep="1"
      casosReplica <- casosEstres[casosEstres$replica == rep,]
      
      total=sum(casosReplica$Count)
      
      #Blucle para poner relacion en porcentaje (0-100) 
      #no nos hace falta saber para que tamaÒo es cada cual, simplemente para cada tamaÒo, sacamos la proporcion
      for(nu in 1:length(casosReplica$Size)){
        # nu=1
        casosReplica$Count[nu] <- (casosReplica$Count[nu]/total)*100
      }#Fin bucle proporcion para cada tamaÒo de sRNA
      
      
      #Una vez con los casosReplica con las proporciones vamos a crear una nueva tabla aÒadiendolos
      if(count == 0){
        tablaProp <- casosReplica
        count = count + 1    
      }else if(count >0){
        tablaProp <- rbind(tablaProp, casosReplica)
      }
      
      
    }#Fin bucle para replicas
  }#Fin bucle para estreses
}#Fin bucle para tiempos


head(tablaProp)

#quitamos las dos ultimas columnas y seguimos como antes de hacer las proporciones
tablaProp <- tablaProp[,-c(5,6)]


# Summarize the data :
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


#Reasignamos variable anterior (distr) sin proporciones
distr <- tablaProp

############## grafico pro tiempos

tiempos <- unique(distr$Time)

for(x in tiempos){
  #x="T4"
  ##nos quedamos con aquellos del tiempo en concretro
  distr1.1 <- distr[distr$Time == x,]
  distr2 <- distr1.1[,c(4,1,2)]
  
  #Quitamos numero de muestras
  for(k in 1:dim(distr2)[1]){
    distr2$Sample[k] <- strsplit(distr2$Sample[k], "-")[[1]][1]
  }
  
  head(distr2)
  
  df3 <- data_summary(distr2, varname="Count", 
                      groupnames=c("Sample", "Size"))
  
  ##sacamos por proporciones por tiempos
  # total <- sum(df3$Count)
  # 
  # estresesNormalizar <- unique(df3$Sample)
  
  
  # Convert dose to a factor variable
  df3$Sample =as.factor(df3$Sample)
  head(df3)
  
  
  
  # Standard deviation of the mean as error bar
  p <- ggplot(df3, aes(x=Size, y=Count, fill=Sample, color=Sample)) + 
    geom_line(size=1) + geom_point() +
    geom_errorbar(aes(ymin=Count-sd, ymax=Count+sd), width=.2,
                  position=position_dodge(.1), size = .4, linetype=2, alpha = .5) +
    ggtitle(paste0("Distribution all libraries 20 - 25 NTPs - ", x)) +
    ylab("% of total sRNAs per case") + xlab("Size sRNAs")
  
  # png(filename=paste0("0-distribution/Distribution20-25-", x, ".png"), width = 3000, height = 2000, res = 300)
  p + theme_minimal(base_size = 15) + theme(plot.title=element_text(hjust = .5))
  ggsave(paste0("0-distribution/Distribution20-25-", x, ".png"), width = 10, height = 8, dpi=300)
  
  #dev.off()
} #fin bucle tiempos

######################### Progresion tamaÒo secuencias entre tiempos, para cada estres

dis <- distr
for(l in 1:dim(dis)[1]){
  dis$Sample[l] <- strsplit(dis$Sample[l], "-")[[1]][1]
}


estreses <- unique(dis$Sample)
estreses <- estreses[!grepl("T0", estreses)]

### Hacemos todo igual, solo que la grafica de T0 como estrÈs no servira, ya que solo es un estado, cuando hagamos el resto de graficas, se seleccionara el estres en 
### cuestion y el T0, en todos, y dara los valores de tamaÒos de sRNAs para T0, en todos los estreses.

for(x in estreses){
  #x="A"
  ##nos quedamos con aquellos del tiempo en concreto
  dis1.1 <- dis[dis$Sample == x | dis$Sample == "T0",]
  dis2 <- dis1.1[,c(3,1,2)]
  
  head(dis2)
  
  df3 <- data_summary(dis2, varname="Count", 
                      groupnames=c("Time", "Size"))
  # Convert dose to a factor variable
  df3$Time =as.factor(df3$Time)
  head(df3)
  
  # Standard deviation of the mean as error bar
  ggplot(df3, aes(x=Size, y=Count, fill=Time, color=Time)) + 
    geom_line(size=1) + geom_point() +
    geom_errorbar(aes(ymin=Count-sd, ymax=Count+sd), width=.2,
                  position=position_dodge(.1), size = .4, linetype=2, alpha = .5) +
    ggtitle(paste0("disibution all libraries 20 - 25 NTPs - ", x)) +
    ylab("% of total sRNAs per case") + xlab("Size sRNAs") + theme_minimal(base_size = 15) + theme(plot.title=element_text(hjust = .5))
  
  ggsave(paste0("0-distribution/Distribution20-25-TIME-", x, ".png"), width = 10, height = 8, dpi=300)
  
  
  
  
  # png(filename=paste0("0-distribution/Distribution20-25-TIME", x, ".png"), width = 3000, height = 2000, res = 300)
  # p + theme_minimal(base_size = 15) + theme(plot.title=element_text(hjust = .5))
  # dev.off()
} #fin bucle tiempos


############################################################################################### 
##
##                        INICIO AN¡LISIS DE EXPRESI”N DIFERENCIAL
##
###############################################################################################-

setwd(houseDIR)

#Comprobar si existe la carpeta resultados y borrar
if( dir.exists(paths = "0-results") ){
  unlink("D:/TFM/0-results", recursive = T)
  timeS.DIR <- dir("LLUIS/Results/01-cuentas_librerias/")
  timeS.DIR <- timeS.DIR[!grepl("CUENTAS", timeS.DIR)]
  dir.create("0-results")
}else if( !dir.exists(paths = "0-results")){
  timeS.DIR <- dir("LLUIS/Results/01-cuentas_librerias/")
  timeS.DIR <- timeS.DIR[!grepl("CUENTAS", timeS.DIR)]
  dir.create("0-results")}

#NO hacemos an·lisis de expresion diferencial para T0
timeS.DIR <- timeS.DIR[!timeS.DIR=="T0"]
#timeS.DIR <- c('T1', 'T2', 'T4')

#########################################################
###
###  CREAR TABLAS INPUT EXPRESI”N DIFERENCIAL 
###
#########################################################

setwd(houseDIR)

for (time.DIR in timeS.DIR){
  
  # time.DIR="T1"
  dir.create(paste0("0-results/", time.DIR))
  dir.create(paste0("0-results/", time.DIR,"/data"))
  
  dir.create(paste0("0-results/", time.DIR,"/sRNA"))
  dir.create(paste0("0-results/", time.DIR,"/sRNA/edgeR"))
  dir.create(paste0("0-results/", time.DIR,"/sRNA/DESeq2"))
  dir.create(paste0("0-results/", time.DIR,"/sRNA/NOISeq"))
  dir.create(paste0("0-results/", time.DIR,"/sRNA/D.E.N"))
  
  dir.create(paste0("0-results/", time.DIR,"/images"))
  dir.create(paste0("0-results/", time.DIR,"/images/edgeR"))
  dir.create(paste0("0-results/", time.DIR,"/images/DESeq2"))
  dir.create(paste0("0-results/", time.DIR,"/images/NOISeq"))
  
  dir.create(paste0("0-results/", time.DIR,"/VennDiagrams"))
  dir.create(paste0("0-results/", time.DIR,"/VennDiagrams/DEN"))
  dir.create(paste0("0-results/", time.DIR,"/SAVE"))
  
  # setwd(paste0(houseDIR, "/db/", time.DIR))
  
  DIRECTORIES <- dir(paste0("LLUIS/Results/01-cuentas_librerias/", time.DIR))
  
  #Tenemos una carpeta denominada CONTROL, pero vamos a tratarla por separado
  DIRECTORIES <- DIRECTORIES[!grepl("CONTROL", DIRECTORIES)]
  #DIRECTORIES <- DIRECTORIES[!grepl("txt", DIRECTORIES)]
  
  
  #Queremos un orden concreto para los estreses simples, vamos a comprobar si estan (todos ellos) y los que no estan se pondran a continuacion
  simples <- c("A", "C", "D", "HSVd", "MON", "SA", "SD")
  
  
  if(all(simples %in% DIRECTORIES)){
    #si estan todos, caso tiempo 1,2,4, primero me pones "estos" y luego el resto
    parte1 <- c("C", "D", "SA", "SD", "MON", "HSVd", "A")
    parte2 <- DIRECTORIES[!DIRECTORIES %in% simples]
    ## En este analisis de expresion diferencial no tenemos en cuenta el POLEN
    parte2 <- parte2[!grepl("POLEN", parte2)]
    DIRECTORIES <- c(parte1, parte2)  #en el caso de tiempo 4 se pondran el resto de estreses aqui
  }
  
  #Para hacer la tabla resumen de condiciones, creamos un data frame, que ha medida del experimento se ira llenando
  exp.sumALL <- data.frame(Condition=DIRECTORIES,
                           Lw.D=(1:length(DIRECTORIES)),
                           Up.D=(1:length(DIRECTORIES)),
                           Total.D=(1:length(DIRECTORIES)),
                           Lw.E=(1:length(DIRECTORIES)),
                           Up.E=(1:length(DIRECTORIES)),
                           Total.E=(1:length(DIRECTORIES)),
                           Lw.N=(1:length(DIRECTORIES)),
                           Up.N=(1:length(DIRECTORIES)),
                           Total.N=(1:length(DIRECTORIES)),
                           Lw.D.E.N=(1:length(DIRECTORIES)),
                           Up.D.E.N=(1:length(DIRECTORIES)),
                           Total.D.E.N=(1:length(DIRECTORIES)))
  
  
  exp.BCV.e <- data.frame(Condition=DIRECTORIES,
                          BCV=(1:length(DIRECTORIES)))
  
  ########################################### INI PRUEBA CREAR TABLA
  exp.loop=1
  #bucle para hacerlo en todos los estreses
  DIRECTORIES
  
  for(DIR in DIRECTORIES ){
    #DIR="A"
    #Valores que no queremos borrar, que luego nos serviran para resumen final, o para countinuar con el bucle
    no.rm.1 <- grep("DIR", ls(), value = T)
    no.rm.2 <- grep("^exp.", ls(), value = T)
    no.rm <- c(no.rm.1, no.rm.2)
    
    rm(list=ls()[!(ls()%in%no.rm)])
    corteLgFC = 1
    
    
    #DIR="C"
    workDIR = paste0(houseDIR,"/LLUIS/Results/01-cuentas_librerias/", time.DIR,"/", DIR)
    
    i.edge <- paste0("D://TFM/0-results/", time.DIR, "/images/edgeR/",DIR)
    i.deseq <- paste0("D://TFM/0-results/", time.DIR, "/images/DESeq2/",DIR)
    i.noiseq <- paste0("D://TFM/0-results/", time.DIR, "/images/NOISeq/",DIR)

    #Creamos directorios para el estres
    dir.create(i.edge)
    dir.create(i.deseq)
    dir.create(i.noiseq)
    
    setwd(workDIR)
   
    # directory<-workDIR
    # 
    # 
    files.dir <- dir()
    
    #Queremos asegurarnos orden estres control en la tabla y que sean archivos acabados en ".collapsed20-25.txt" 
    #1. Nos quedamos con solo aquellos que acaben en...
    files.dir.stress <- files.dir[grep(".txt", files.dir)]
    
    #2. Cogemos el control que esta en otra carpeta
    files.dir.control <- dir("../CONTROL")
    
    #Si estamos en el tiempo 4, tendremos 5 controles, eliminamos el primero
    if(time.DIR == "T4"){
      files.dir.control <- files.dir.control[!grepl("CONTROL-1.txt", files.dir.control)]
    }
    
    
    #3. Unimos y generamos orden de archivos para la tabla CONTROL-ESTRES
    files.dir <- c(paste0("../CONTROL/",files.dir.control), files.dir.stress)
    
    
    contador = 0
    for (file in files.dir ) {  
      
      #Como ya sabemos que solo vamos a trabajar con archivos terminados en .collapsed20-25.txt, podemos quitarlo
      #if(grepl(".collapsed20-25.txt", file)){
      # file="C-2.tsv"
      # file="../CONTROL/CONTROL-2.tsv"
      
      #segun el archivo CONTROL o estres (control esta en otra carpeta, tenemos que extraerlo de otra forma)
      if(grepl("CONTROL", file)){
        name <- strsplit(strsplit(file, ".txt")[[1]], "/")[[1]][3]
        
      }else{
        name <- strsplit(file, ".txt")[[1]] 
      }
      
      
  
      
      if (contador == 0){
        opened <-read.table(file,sep="\t", header=F, blank.lines.skip=FALSE,quote = NULL)
        colnames(opened) <- c("sRNA", name)
        contador=contador+1
      }else{
        opened2 <-read.table(file,sep="\t", header=F, blank.lines.skip=FALSE, quote = NULL)
        colnames(opened2) <- c("sRNA", name)
        
        opened <- merge(opened, opened2, by="sRNA", all = TRUE, incomparables = NA)
      }
      #}
    }#Fin para colapsar tablas dentro de un mismo estres
    
    #si son NA ponemos 0
    opened[2:length(opened)][is.na(opened[2:length(opened)])] <- 0
    
    ##################################################################################################
    ##################################################################################################
    #Ya tenemos la tabla con la que haremos el analisis de expresion diferencial, estudio pareado!!!!!
    ##################################################################################################
    ##################################################################################################
    
    #IMPORTANTE, cambiamos a integer nuestra matriz de cuentas
    
    #Vemos que esta como numeric(double), y lo queremos como integer
    #sapply(opened, class)
    #typeof(opened$`C.SA-1`)
    
    #Convertimos a una matriz, eliminando la primera columna o secuencias sRNA que es una cadena
    m <- as.matrix(sapply(opened[-1], as.integer)) 
    #Volvemos a asignar a opened
    opened[-1] <- as.data.frame(m)
    
    head(opened)
    #Vemos que ya esta como integer y no como double ;), deseq ya no necesita hacer el cambio, y puede que tengamos menos problemas en adelante
    #typeof(opened$`Control-2`)
    
    # buscamos los archivos refentes al estres y control de forma separada, para asegurarnos el orden (estres | control)
    # esto es para el diseÒo de DESEQ!! Sin embargo tambien nos sirve para saber cuantos controles tenemos y cuantas
    # condiciones de estres. lo que nos servira para filtrar aquellas secuencias que estan al menos 2 veces en el estres
    # yo 3 veces en los controles.
    
    # Nos sirve lo que hemos hecho anterioremente
    # files.dir.stress
    # files.dir.control
  
    
    #sampleFiles
    ####################### @Filtrado 1
    #### nos quedamos con aquellas que estan en al menos 2 veces en los estreses, o 3 veces en los controles
    sampleFiles.x <- files.dir.control
    sampleFiles.y <- files.dir.stress
    sampleFiles <- c(sampleFiles.x, sampleFiles.y)
    
    
    # A. Control, teniendo en cuenta que la primera fila son los sRNAs
    # primer posicion control --- se restan los estreses para saber hasta que posicion son los controles
    keep1 <- rowSums(opened[ 2:(length(opened)-length(files.dir.stress)) ]>=5) >= length(files.dir.control)-1
    #table(keep1)
    groupControl <- opened[keep1,]
    
    
    # B. Estres, teniendo en cuenta que la primera fila son los sRNAs
    ini <- ((length(opened)-length(files.dir.stress))+1)
    end <- length(opened)
    
    keep2 <- rowSums(opened[ ini:end ]>=5) >= length(files.dir.stress)-1
    #table(keep2)
    groupStress <- opened[keep2,]
    
    
    #juntamos elementos y eliminamos duplicados (aquellos que cumplan tambiÈn la otra condicion)
    g <- rbind(groupControl,groupStress)
    
    g <- g[!duplicated(g), ]
    
    #keep <- rowSums(g[3:7]>=5) >= 3
    #table(keep)
    #g2 <- g[keep,]
    
    ########## ESCRIBIMOS TABLA QUE SERVIR¡ PARA OTROS AN¡LISIS DE EXPRESI”N DIFERENCIAL (si nos lo hacemos de forma interna)
    dataDIR=paste0("0-results/", time.DIR, "/data/", DIR, ".txt")
    
    
    setwd(houseDIR)
    write.table(g, file=dataDIR, sep = "\t", quote = F, row.names = F)
    
    
    #####################################################################################################################
    
    #####################################################################################################################
    
    #####################################################################################################################
    
    #####################################################################################################################
    
    mycounts <- g
    rownames(mycounts) <- g$sRNA
    mycounts <- mycounts[-1]
    
    
    namesCol <- colnames(mycounts)
    namesCol2 <- colnames(mycounts)
    for(k in 1:length(namesCol)){
      namesCol[k] <- strsplit(namesCol[k], "-")[[1]][1]
    }
    
    # namesCol
    
    
    ####################################################################### START NOISEQ
    options(digits=3, width=95)
    
    head(mycounts)
    
    #creamos una forma para poder crear myfactors, frecuencia de controles y estres
    dtCol <- as.data.frame(table(namesCol))
    
    #ordenamos, el control, EN ESTE caso, ser· siempre m·s grande que los estreses
    dtCol <- dtCol[order(dtCol$Freq, decreasing = T),]
    
    rownames(dtCol) <- NULL
    
    cCol <- dtCol$Freq[1]
    sCol <- dtCol$Freq[2]
    
    myfactors = data.frame(Treatment=namesCol, 
                           TreatmentRun=c( paste0(namesCol[1], "_", 1:cCol), paste0(namesCol[length(namesCol)], "_",   1:sCol) ),
                           Run=c( paste0("R", rep(1, cCol)), paste0("R", rep(2, sCol))))
    
    # myfactors
    
    #########################################################################################################
    ###read.data
    mydata <- readData(data=mycounts, factors=myfactors) 
    #mydataU <- readData(data=myUQUA, factors=myfactors)
    #mydataN <- readData(data=myTMM, factors=myfactors)
    #   str(mydata)
    #   head(assayData(mydata)$exprs)
    #   head(pData(mydata))
    
    
    mycountsbio = dat(mydata, factor = NULL, type = "countsbio")
    #mycountsbioN = dat(mydataN, factor = NULL, type = "countsbio")
    #mycountsbioU = dat(mydataU, factor = NULL, type = "countsbio")
    
    i.n <- paste0(i.noiseq, "/", DIR)
    
    i.n.1 <- paste(i.n, "-boxplot-expression.png", sep="")
    png(i.n.1, 1200, 800, pointsize=20)
    explo.plot(mycountsbio, toplot = 1, samples = NULL, plottype = "boxplot")
    #explo.plot(mycountsbioN, toplot = 1, samples = NULL, plottype = "boxplot")
    dev.off()
    
    ## saturation
    mysaturation = dat(mydata, k = 0, ndepth = 6, type = "saturation")
    #mysaturationN = dat(mydataN, k = 0, ndepth = 6, type = "saturation")
    
    i.n.2 <- paste(i.n, "-saturation.png", sep="")
    png(i.n.2, 800, 800, pointsize=20)
    explo.plot(mysaturation, toplot = 1, samples = 1:length(namesCol), yleftlim = NULL, yrightlim = NULL)
    dev.off()
    
    ###################################################
    ### code chunk number 18: fig_boxplot3
    ###################################################
    
    i.n.3 <- paste(i.n, "-fig-boxplot3.png", sep="")
    png(i.n.3, 1200, 800, pointsize=20)
    explo.plot(mycountsbio, toplot = 1, samples = NULL, plottype = "barplot")
    dev.off()
    
    ###################################################
    ### code chunk number 22: fig_countdistr
    ###################################################
    mycd = dat(mydata, type = "cd", norm = FALSE, refColumn = 1)
    
    i.n.4 <- paste(i.n, "-fig-countdistr.png", sep="")
    png(i.n.4, 1600, 800, pointsize=20)
    explo.plot(mycd)
    dev.off()
    
    ###################################################
    ### code chunk number 24: fig_PCA
    ###################################################
    myPCA = dat(mydata, type = "PCA")
    i.n.5 <- paste(i.n, "-PCA.png", sep="")
    png(i.n.5, 800, 800, pointsize=20)
    explo.plot(myPCA, factor = "Treatment")
    dev.off()
    
    ###################################################
    ### code chunk number 26: normalization
    ###################################################
    myTMM = tmm(assayData(mydata)$exprs, long = 1000, lc = 0)
    #myUQUA = uqua(assayData(mydata)$exprs, lc = 0.5, k = 0)
    #head(myTMM[,1:4])
    #head(mycounts[,1:4])
    ##head(myUQUA[,1:4])
    
    ###################################################
    ### code chunk number 27: filtering
    ###################################################
    #myfilt = filtered.data(mycounts, factor = myfactors$Treatment, norm = FALSE, depth = NULL, method = 1, cv.cutoff = 110, cpm = 1, p.adj = "BH")
    
    #myfilt2 = filtered.data(mycounts, factor = myfactors$Treatment, norm = T, depth = NULL, method = 1, cv.cutoff = 100, cpm = 1, p.adj = "fdr")
    
    # head(mynoiseq@results[[1]])
    ###################################################
    ### code chunk number 33: NOISeq.Rnw:875-877
    ###################################################
    mynoiseqbio = noiseqbio(mydata, k = 0.5, norm = "tmm", factor="Treatment", nclust = 10, lc = 0, r = 50, adj = 1.5, plot = FALSE,
                            a0per = 0.9, random.seed = 12345)
    #, filter = 1, cv.cutoff = 95
    #?noiseqbio
    head(mynoiseqbio@results[[1]])
    
    corte = 0.95
    ###################################################
    ### code chunk number 35: NOISeq.Rnw:943-946
    ###################################################
    
    mynoiseq.all = degenes(mynoiseqbio, q = 0, M = NULL)
    
    mynoiseq.deg = degenes(mynoiseqbio, q = corte, M = NULL)
    mynoiseq.deg22 <- mynoiseq.deg[mynoiseq.deg$log2FC >= corteLgFC |  mynoiseq.deg$log2FC <= -corteLgFC,]
    mynoiseq.deg1 = mynoiseq.deg22[mynoiseq.deg22$log2FC > 0,]
    mynoiseq.deg2 = mynoiseq.deg22[mynoiseq.deg22$log2FC < 0,]
    
    #   mynoiseq.deg1 = degenes(mynoiseqbio, q = corte, M = "up")
    #   mynoiseq.deg2 = degenes(mynoiseqbio, q = corte, M = "down")
    
    mynoiseq.all2 <- data.frame(sRNA=rownames(mynoiseq.all))
      
    mynoiseq.all2 <- cbind(mynoiseq.all2, mynoiseq.all)
    
    #Comrprobamos que rownames concuerda con columna sRNA
    #head(mynoiseq.all2)
    
    write.table(mynoiseq.all2, paste0("0-results/",time.DIR, "/sRNA/NOISeq/", DIR,".txt"), sep="\t", row.names = F, quote = F )
    
    
    
    exp.sumALL$Lw.N[exp.loop] <- dim(mynoiseq.deg1)[1]
    exp.sumALL$Up.N[exp.loop] <- dim(mynoiseq.deg2)[1]
    exp.sumALL$Total.N[exp.loop] <- dim(mynoiseq.deg22)[1]
    
    
    
    ###################################################
    ### code chunk number 36: fig_summ_expr
    ###################################################
    
    i.n.6 <- paste(i.n, "-DEplot-expr.png", sep="")
    png(i.n.6, 1000, 1000, pointsize=20)
    DE.plot(mynoiseqbio, q = corte, graphic = "expr", log.scale = TRUE)
    dev.off()
    
    ###################################################
    ### code chunk number 37: fig_summ_MD
    ###################################################
    i.n.7 <- paste(i.n, "-DEplot-MD.png", sep="")
    png(i.n.7, 1000, 1000, pointsize=20)
    DE.plot(mynoiseqbio, q = corte, graphic = "MD")
    dev.off()
    
    #mynoiseq.deg22
    #write.table(resdata2, file="resdata2.tsv", sep="\t", quote = F, dec = ",", row.names = F)
    exp.lw.fil.20 <- subset(mynoiseq.deg22, ( (nchar(rownames(mynoiseq.deg22))) == 20 & mynoiseq.deg22$log2FC <= -corteLgFC) )
    exp.up.fil.20 <- subset(mynoiseq.deg22, ( (nchar(rownames(mynoiseq.deg22))) == 20 & mynoiseq.deg22$log2FC >= corteLgFC) )
    exp.lw.fil.21 <- subset(mynoiseq.deg22, ( (nchar(rownames(mynoiseq.deg22))) == 21 & mynoiseq.deg22$log2FC <= -corteLgFC) )
    exp.up.fil.21 <- subset(mynoiseq.deg22, ( (nchar(rownames(mynoiseq.deg22))) == 21 & mynoiseq.deg22$log2FC >= corteLgFC) )
    exp.lw.fil.22 <- subset(mynoiseq.deg22, ( (nchar(rownames(mynoiseq.deg22))) == 22 & mynoiseq.deg22$log2FC <= -corteLgFC) )
    exp.up.fil.22 <- subset(mynoiseq.deg22, ( (nchar(rownames(mynoiseq.deg22))) == 22 & mynoiseq.deg22$log2FC >= corteLgFC) )
    exp.lw.fil.23 <- subset(mynoiseq.deg22, ( (nchar(rownames(mynoiseq.deg22))) == 23 & mynoiseq.deg22$log2FC <= -corteLgFC) )
    exp.up.fil.23 <- subset(mynoiseq.deg22, ( (nchar(rownames(mynoiseq.deg22))) == 23 & mynoiseq.deg22$log2FC >= corteLgFC) )
    exp.lw.fil.24 <- subset(mynoiseq.deg22, ( (nchar(rownames(mynoiseq.deg22))) == 24 & mynoiseq.deg22$log2FC <= -corteLgFC) )
    exp.up.fil.24 <- subset(mynoiseq.deg22, ( (nchar(rownames(mynoiseq.deg22))) == 25 & mynoiseq.deg22$log2FC >= corteLgFC) )
    exp.lw.fil.25 <- subset(mynoiseq.deg22, ( (nchar(rownames(mynoiseq.deg22))) == 25 & mynoiseq.deg22$log2FC <= -corteLgFC) )
    exp.up.fil.25 <- subset(mynoiseq.deg22, ( (nchar(rownames(mynoiseq.deg22))) == 25 & mynoiseq.deg22$log2FC >= corteLgFC) )
    
    
    d <- data.frame(Size.sRNA=c("20", "20", "21", "21", "22", "22", "23", "23", "24", "24", "25", "25"),
                    State=rep(c("Lw","Up"),6),
                    Count=c(length(exp.lw.fil.20$log2FC), length(exp.up.fil.20$log2FC), length(exp.lw.fil.21$log2FC), length(exp.up.fil.21$log2FC), length(exp.lw.fil.22$log2FC), length(exp.up.fil.22$log2FC), length(exp.lw.fil.23$log2FC), length(exp.up.fil.23$log2FC), length(exp.lw.fil.24$log2FC), length(exp.up.fil.24$log2FC), length(exp.lw.fil.25$log2FC), length(exp.up.fil.25$log2FC)))
    
    name.tab1 <- paste0("0-results/", time.DIR, "/images/NOISeq/",DIR, "/sum20_25-", DIR, ".txt")
    write.table(d, file=name.tab1, quote = F, sep="\t", row.names = F)
    
    #damos nombre al tÌtulo de imagen
    titleplot <- paste0("Sample:", DIR,"\n\nSummary differential expression analysis accounts shown by size\n\nq >= 0.95 and (Log2FC <= -1 or Log2FC >= 1)")
    
    p1<-ggplot(d,aes(x=Size.sRNA,y=Count,fill=State), color=State) +  
      stat_summary(fun.y=mean,position=position_dodge(),geom="bar") +
      ggtitle(titleplot) + 
      ylim(0,8000)+
      theme(plot.title = element_text(lineheight=.8, face="bold"))
    
    i.d.13 <- paste0(i.n,"-counts-diff-exp.png")
    
    #png(i.d.13, width=12, height=8, units='in', res=300)
    par(mar=c(0,0,0,0)) # Remove unnecessary margins
    #png(i.d.13, 1500, 800, pointsize=20)
    p1
    ggsave(i.d.13, width = 10, height = 8, dpi=300)
    #dev.off()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ############################################### START DESEQ2
    
    #Nos quedamos con la parte por delante del "-", el primer elemento de sampleFiles.x (del estrÈs), y el primer elemento de 
    #sampleFiles.y (del control), si es "A-2.collapsed20-25.txt", queremos el "A"
    name.x <- strsplit(sampleFiles.x[1], "-")[[1]]
    name.x <- name.x[1]
    name.y <- strsplit(sampleFiles.y[1], "-")[[1]]
    name.y <- name.y[1]
    
    ### SE PUEDE SALTAR AL EDGE!! 
    
    
    sampleCondition<-c(rep(name.x,length(sampleFiles.x)), rep(name.y,length(sampleFiles.y)))
    
    # The table is created for design
    sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
    
    counts.A <- g[-1]
    rownames(counts.A) <- g$sRNA
    
    colData <- data.frame(row.names=colnames(counts.A),
                          condition=factor( c(rep("untreated",length(sampleFiles.x)), rep("treated",length(sampleFiles.y)))))
    
    #colData
    
    dds <-DESeqDataSetFromMatrix(countData=counts.A, colData = colData, design=~condition)
    #dds
    
    # Small filtering 
    #dds <- dds[ rowSums(counts(dds)) > 1, ]
    
    colData(dds)$condition <- factor(colData(dds)$condition, levels=c("untreated","treated"))
    
    # Estimate size factors
    dds=estimateSizeFactors(dds)
    # ?estimateSizeFactors
    # differential expression test
    dds <- DESeq(dds)#, fitType = 'local')
    
    #plotDispEsts(dds)
    #The variance stabilizing transformation involves integration of the variance as a function of the mean. When the local fit is used, the integration is not a symbolic integration using the parametrized function, but a numeric integration of the local regression curve.
    
    #The Cook's cutoff used here is the .90 quantile 
    p=1 #condition
    m=length(sampleFiles) #en este caso contamos con 8 muestras
    
    cooksCutoff<-qf(0.95, p, m - p)
    
    res <- results(dds, cooksCutoff=cooksCutoff, pAdjustMethod = "BH")
    
    #We can summarize some basic tallies using the summary function.
    summary(res$pvalue <= 0.05)
    
    register(MulticoreParam(10))
    
    
    # We can order our results table by the smallest adjusted p value:
    resOrdered <- res[order(res$pvalue),]
    
    ########################################################  IMPRIMIMOS TABLA COMPLETA
    r.d <- as.data.frame(resOrdered)
    r.d$sRNA <- rownames(r.d)
    
    #Cambiamos el orden de las columnas
    r.d <- r.d[,c(length(r.d),1:(length(r.d)-1))]
    rownames(r.d) <- NULL
    
    i.d.srna <- paste0("0-results/",time.DIR,"/sRNA/DESeq2/", DIR, ".txt")
    write.table(r.d, i.d.srna, row.names = F, quote = F, sep = "\t")
    
    i.d.srna2 <- paste0("0-results/",time.DIR,"/sRNA/DESeq2/", DIR, "-excel.csv")
    write.table(r.d, i.d.srna2, row.names = F, quote = F, sep = ";", dec = ",")
    
    
    # Some plot tests
    
    GeneCounts <- counts(dds)
    idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)})
    sum(idx.nz)
    con <- as.character(colData(dds)$condition)
    
    # plot densities of counts for the different samples
    ### to assess their distributions
    
    i.d <- paste0(i.deseq, "/", DIR)
    
    i.d.1 <- paste(i.d, "-multidensity.png", sep="")
    
    png(i.d.1, 1500, 800, pointsize=20)
    old.par <- par(mfrow=c(1, 2))
    multidensity(counts(dds, normalized = F)[idx.nz,],xlab="mean counts", xlim=c(0, 100))
    multidensity(counts(dds, normalized = T)[idx.nz,],xlab="mean counts", xlim=c(0, 100))
    par(old.par)
    dev.off()
    
    
    i.d.2 <- paste(i.d, "-multiecdf.png", sep="")
    
    png(i.d.2, 1500, 800, pointsize=20)
    old.par <- par(mfrow=c(1, 2))
    multiecdf( counts(dds, normalized = F)[idx.nz ,], xlab="mean counts", xlim=c(0, 100))
    multiecdf( counts(dds, normalized = T)[idx.nz ,], xlab="mean counts", xlim=c(0, 100))
    par(old.par)
    dev.off()
    
    i.d.3 <- paste(i.d, "-DESeq2-BH.png", sep="")
    
    
    maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
      with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
      with(subset(res, (pvalue<=thresh & (log2FoldChange <= -corteLgFC | log2FoldChange >= corteLgFC))), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
    }
    png(i.d.3, 1500, 800, pointsize=20)
    maplot(res, main="plotMA :: DESeq2 - p.value <= 0.05 & (logFC <= -1 | logFC >= 1)")
    dev.off()
    
    
    # we use data.frame(res$baseMean, res$log2FoldChange, res$pvalue <= 0.05),‚Ä¶‚Ä¶ use ‚Äúres‚Äù to plot res$adj values
    #lotMA(data.frame(res$baseMean, res$log2FoldChange, res$pvalue <= 0.05), ylim=c(-10,10),main='plotMA :: DESeq2 - p.value < 0.05')
    #dev.off()
    
    # Plot dispersions
    i.d.4 <- paste(i.d, "-qc-dispersions.png", sep="")
    
    png(i.d.4, 1500, 800, pointsize = 20)
    plotDispEsts(dds, main="Dispersion plot")
    dev.off()
    
    # Regularized log transformation for clustering/heatmaps, etc
    rld <- rlogTransformation(dds)
    #(mycols <- brewer.pal(8, "Dark2")[1:length(unique(dds$condition))])
    mycols <- brewer.pal(8, "Dark2")[1:length(unique(dds$condition))]
    # Sample distance heatmap
    sampleDists <- as.matrix(dist(t(assay(rld))))
    
    i.d.5 <- paste(i.d, "-qc-heatmap-samples.png", sep="")
    png(i.d.5, w=1000, h=1000, pointsize=20)
    heatmap.2(sampleDists, key=F, trace="none",
              srtCol=45, adjCol=c(1,1), offsetCol=-0.1, offsetRow = -0.1,
              col=colorpanel(100, "black", "white"),
              ColSideColors=mycols[dds$condition], RowSideColors=mycols[dds$condition],
              margin=c(10, 10), main="Sample Distance Matrix")
    dev.off()
    
    i.d.6 <- paste(i.d, "-PCA.png", sep="")
    
    #png(i.d.6, w=1000, h=800, pointsize=40)
    DESeq2::plotPCA(rld, intgroup=c("condition"), ntop=500)
    ggsave(i.d.6, width = 10, height = 6, dpi=300)
    #dev.off()
    
    ################################################################ PCA with labels
    # 1 PCA
    rv = rowVars(assay(rld))
    ntop=500
    select = order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
    pca = prcomp(t(assay(rld)[select,]))
    
    # proportion of variance
    variance = pca$sdev^2 / sum(pca$sdev^2)
    variance = round(variance, 3) * 100
    
    d<-plotPCA(rld, returnData=TRUE)
    
    i.d.7 <- paste(i.d, "-PCA-with-names.png", sep="")
    
    #png(i.d.7, w=1000, h=800, pointsize=40)
    
    ggplot(d, aes(x=PC1,y=PC2,col=condition,label=name)) + 
      geom_point() + geom_text(aes(label=name),hjust=0.5, vjust=-0.5) +
      xlab(paste("PC1 (", variance[1], "%)", sep="")) +
      ylab(paste("PC2 (", variance[2], "%)", sep=""))
    
    ggsave(i.d.7, width = 10, height = 8, dpi=300)
    #dev.off()
    
    # 2 PCA
    se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1),
                               colData=colData(dds))
    
    DESeqT.se <- DESeqTransform( se )
    
    i.d.8 <- paste(i.d, "-PCA-DESeqTransform.png", sep="")
    
    #png(i.d.8, w=1000, h=800, pointsize=40)
    plotPCA( DESeqT.se, intgroup = "condition")
    ggsave(i.d.8, width = 10, height = 8, dpi=300)
    #dev.off()
    
    # pca
    rv = rowVars(assay(DESeqT.se))
    ntop=500
    select = order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
    pca = prcomp(t(assay(DESeqT.se)[select,]))
    
    # proportion of variance
    variance = pca$sdev^2 / sum(pca$sdev^2)
    variance = round(variance, 3) * 100
    
    dd<-plotPCA(DESeqT.se, returnData=TRUE)
    
    i.d.9 <- paste(i.d, "-PCA-DESeqTransform-with-names.png", sep="")
    #png(i.d.9, w=1000, h=800, pointsize=40)
    ggplot(dd, aes(x=PC1,y=PC2,col=condition,label=name)) + 
      geom_point() + geom_text(aes(label=name),hjust=0.5, vjust=-0.5) +
      xlab(paste("PC1 (", variance[1], "%)", sep="")) +
      ylab(paste("PC2 (", variance[2], "%)", sep=""))
    ggsave(i.d.9, width = 10, height = 8, dpi=300)
    
    ##dev.off()
    
    #   i.d.99 <- paste(i.d, "-PCA-DESeqTransform-with-names.tiff", sep="")
    #   library(ggplot2)
    #   tiff(i.d.99, width=10, height=5, units="in", res=300)
    #   par(mar=c(0,0,0,0)) # Remove unnecessary margins
    #   ggplot(dd, aes(x=PC1,y=PC2,col=condition,label=name)) + 
    #     geom_point() + geom_text(aes(label=name),hjust=0.6, vjust=1.4) +
    #     xlab(paste("PC1 (", variance[1], "%)", sep="")) +
    #     ylab(paste("PC2 (", variance[2], "%)", sep=""))
    #   dev.off()
    
    ################################################################  Examine plot of p-values
    
    i.d.10 <- paste(i.d, "-hist-nofiltered.png", sep="")
    png(i.d.10, w=1500, h=800, pointsize=20)
    old.par <- par(mfrow=c(1, 2))
    hist(res$pvalue, breaks=50, col="grey", main = "Histogram p.value sRNA", xlab = "p. value")  #, ylim=c(0,1500))
    hist(res$padj, breaks=50, col="grey", main = "Histogram FDR sRNA - p.value < 0.05", xlab = "FDR") #, ylim=c(0,1500))
    par(old.par)
    dev.off()
    
    
    ## SHOW FDR VALUES, WE FILTERED by PVALUE
    filtered.FDR=res$padj[(res$pvalue <= 0.05 & (res$log2FoldChange <= -corteLgFC | res$log2FoldChange >= corteLgFC))]
    filtered.PVALUE=res$pvalue[(res$pvalue <= 0.05 & (res$log2FoldChange <= -corteLgFC | res$log2FoldChange >= corteLgFC))]
    
    i.d.11 <- paste(i.d, "-hist-filtered.png", sep="")
    
    png(i.d.11, w=1500, h=800, pointsize=20)
    old.par <- par(mfrow=c(1, 2))
    hist(filtered.PVALUE, breaks=50, col="grey", main = "Histogram p.value filtered sRNA", xlab = "p. value") #, ylim=c(0,300))
    hist(filtered.FDR, breaks=50, col="grey", main = "Histogram FDR (BH) filtered sRNA - p.value < 0.05", xlab = "FDR") #, ylim=c(0,300))
    par(old.par)
    dev.off()
    
    #mcols(res)$description
    
    rld <- rlog(dds)
    vsd <- varianceStabilizingTransformation(dds)
    
    #head(assay(rld), 3)
    
    DESeq2Res <- res[ !(is.na(res$pvalue)), ]
    #summary(DESeq2Res)
    
    #   i.d.12 <- paste(i.d, "-fdrtool-stat.png", sep="")
    #   
    #   png(i.d.12, w=1000, h=1000, pointsize=20)
    #   FDR.DESeq2Res <- fdrtool(DESeq2Res$stat, statistic= "normal", plot = T)
    #   dev.off()
    
    ################################################################ diff expression results
    # Get differential expression results
    #res <- results(dds, cooksCutoff=cooksCutoff, pAdjustMethod = "BH")
    #table(res$pvalue <= 0.05)
    
    ## Order by adjusted p-value
    res <- res[order(res$pvalue), ]
    
    ## Merge with normalized count data
    resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
    names(resdata)[1] <- "sRNA"
    
    
    #head(resdata)
    
    # restada1 (points in red plotMA), resdata2 wich we will use, to compare tu sRNA miRNAs targets. 
    #resdata1 <- subset(resdata, resdata$pvalue < 0.05 & !is.na(resdata$pvalue))
    resdata2 <- subset(resdata, resdata$pvalue <= 0.05 & (resdata$log2FoldChange >= corteLgFC | resdata$log2FoldChange <= -corteLgFC) & !is.na(resdata$pvalue))
    
    #We keep summary data for this stress
    exp.lw.d <- subset(resdata2, resdata2$log2FoldChange <= -corteLgFC)
    exp.up.d <- subset(resdata2, resdata2$log2FoldChange >= corteLgFC)
    
    #write.table(resdata2, file="resdata2.tsv", sep="\t", quote = F, dec = ",", row.names = F)
    exp.lw.fil.20 <- subset(resdata2, ( (nchar(resdata2$sRNA)) == 20 & resdata2$log2FoldChange <= -corteLgFC) )
    exp.up.fil.20 <- subset(resdata2, ( (nchar(resdata2$sRNA)) == 20 & resdata2$log2FoldChange >= corteLgFC) )
    exp.lw.fil.21 <- subset(resdata2, ( (nchar(resdata2$sRNA)) == 21 & resdata2$log2FoldChange <= -corteLgFC) )
    exp.up.fil.21 <- subset(resdata2, ( (nchar(resdata2$sRNA)) == 21 & resdata2$log2FoldChange >= corteLgFC) )
    exp.lw.fil.22 <- subset(resdata2, ( (nchar(resdata2$sRNA)) == 22 & resdata2$log2FoldChange <= -corteLgFC) )
    exp.up.fil.22 <- subset(resdata2, ( (nchar(resdata2$sRNA)) == 22 & resdata2$log2FoldChange >= corteLgFC) )
    exp.lw.fil.23 <- subset(resdata2, ( (nchar(resdata2$sRNA)) == 23 & resdata2$log2FoldChange <= -corteLgFC) )
    exp.up.fil.23 <- subset(resdata2, ( (nchar(resdata2$sRNA)) == 23 & resdata2$log2FoldChange >= corteLgFC) )
    exp.lw.fil.24 <- subset(resdata2, ( (nchar(resdata2$sRNA)) == 24 & resdata2$log2FoldChange <= -corteLgFC) )
    exp.up.fil.24 <- subset(resdata2, ( (nchar(resdata2$sRNA)) == 25 & resdata2$log2FoldChange >= corteLgFC) )
    exp.lw.fil.25 <- subset(resdata2, ( (nchar(resdata2$sRNA)) == 25 & resdata2$log2FoldChange <= -corteLgFC) )
    exp.up.fil.25 <- subset(resdata2, ( (nchar(resdata2$sRNA)) == 25 & resdata2$log2FoldChange >= corteLgFC) )
    
    
    d <- data.frame(Size.sRNA=c("20", "20", "21", "21", "22", "22", "23", "23", "24", "24", "25", "25"),
                    State=rep(c("Lw","Up"),6),
                    Count=c(length(exp.lw.fil.20$sRNA), length(exp.up.fil.20$sRNA), length(exp.lw.fil.21$sRNA), length(exp.up.fil.21$sRNA), length(exp.lw.fil.22$sRNA), length(exp.up.fil.22$sRNA), length(exp.lw.fil.23$sRNA), length(exp.up.fil.23$sRNA), length(exp.lw.fil.24$sRNA), length(exp.up.fil.24$sRNA), length(exp.lw.fil.25$sRNA), length(exp.up.fil.25$sRNA)))
    
    name.tab1 <- paste0("0-results/",time.DIR,"/images/DESeq2/",DIR, "/sum20_25-", DIR, ".txt")
    write.table(d, file=name.tab1, quote = F, sep="\t", row.names = F)
    
    #damos nombre al tÌtulo de imagen
    titleplot <- paste0("Sample:", DIR, "\n\nSummary differential expression analysis accounts shown by size\n\nP.value <= 0.05 and (Log2FC <= -1 or Log2FC >= 1)")
    
    p1<-ggplot(d,aes(x=Size.sRNA,y=Count,fill=State), color=State) +  
      stat_summary(fun.y=mean,position=position_dodge(),geom="bar") +
      ggtitle(titleplot) +
      ylim(0,8000)+
      theme(plot.title = element_text(lineheight=.8, face="bold"))
    
    i.d.13 <- paste(i.d, "-counts-diff-exp.png", sep="")
    
    #png(i.d.13, width=12, height=8, units='in', res=300)
    par(mar=c(0,0,0,0)) # Remove unnecessary margins
    #png(i.d.13, 1500, 800, pointsize=20)
    p1
    ggsave(i.d.13, width = 10, height = 8, dpi=300)
    #dev.off()
    
    t.d <- subset(resdata, (resdata$pvalue <= 0.05 & (resdata$log2FoldChange >= corteLgFC | resdata$log2FoldChange <= -corteLgFC) & !is.na(resdata$pvalue)))
    
    exp.d <- dim(t.d)[1]
    
    
    i.d.12 <- paste(i.d, "-fdrtool-stat.png", sep="")
    
    png(i.d.12, w=1000, h=1000, pointsize=20)
    FDR.DESeq2Res <- fdrtool(DESeq2Res$stat, statistic= "normal", plot = T)
    dev.off()
    
    
    #llenamos la tabla con la informacion "condition, exp.d, exp.up.d, exp.lw.d, exp.d, exp.lw.d.fil, exp.up.d.fil"
    exp.sumALL$Lw.D[exp.loop] <- dim(exp.lw.d)[1]
    exp.sumALL$Up.D[exp.loop] <- dim(exp.up.d)[1]
    exp.sumALL$Total.D[exp.loop] <- exp.d
    
    
    # options(scipen=10)
    options(scipen=0)
    
    ####
    #####################################################################################################
    ## edgeR
    #####################################################################################################
    ####
    
    #edgeR is an R/Bioconductor package that provides methods for the statistical analysis of count data from comparative experiments on high-throughput sequencing platforms. Particular attention is given to designed multi-factor experiments and to experiments with minimal replication.
    #The package provides methods for assessing differential expression in RNA-Seq, Tag-Seq, SAGE-Seq and other digital gene expression experiments. RNA-Seq is the most common source of gene expression count data (or digital gene expression data), but the methods implemented in edgeR are general and can also be used with ChIP-seq and other genome-scale count data.
    
    #Over-dispersed Poisson count models (primarily the negative binomial) are used to distinguish biological from technical variation. We use generalized linear modelling to handle complex multi-factor experiments. Information sharing techniques ensure rigorous results even for experiments with minimal biological replication.
    
    #Through the use of the negative binomial distribution to model transcript counts we allow the possibility of gene-specific variability, whereby some genes may show more biological variability than others. A measure of this, the biological coefficient of variation (BCV), is inferred from how much the variance of the counts exceeds the variance that would arise from Poisson counts.
    
    #For simple, one-way layout experimental designs (e.g. pairwise comparisons between treatments) conditional likelihood and weighted conditional likelihood are used for estimation of the NB model parameters. This approach permits an exact test to generate p-values for assessing differential expression. Information sharing is used so that genes may take individual values for the BCV, but stabilized towards a common BCV value. This approach profoundly improves inference in small sample experiments.
    
    #For more complex experimental designs we use generalized linear models (GLMs) with the negative binomial distribution to conduct inference on differential expression. A GLM with the same set of explanatory variables but possibly different BCV is fit to the counts for each gene. Cox-Reid adjusted profile likelihood, a well-respected method for adjusting for bias when estimating the variance parameters in non-linear models, is used to estimate the BCV values. Information sharing methods are also applied in the GLM setting. A likelihood ratio test is used to compute p-values for differential expression. Extensive optimization of GLM-fitting routines allows tens of thousands of GLMs to be fit in edgeR in a matter of seconds.
    
    #With integration of the edgeR package with another R/Bioconductor package called GOSeq, differential expression results from edgeR can be related easily to existing annotation databases such as Gene Ontology or the Molecular Signatures Database, while accounting for gene-length bias on differential expression.
    
    #All in all, edgeR represents a very powerful and flexible modular pipeline for the statistical analysis of comparative RNA-Seq experiments and other genome-scale count data.
    
    
    #Leemos los datos tabla guardada, si no van de forma interna
    #g = read.table(dataDIR, sep="\t", header=T)
    
    #preparamos para edgeR
    rownames(g)<-g$sRNA
    gg<-g[-1]
    
    #x <- gg
    
    #reordenamos la columna para edge, queremos control + estrÈs (si mahoma no va a la montaÒa, la montaÒa va a mahoma)
    #x <- gg[,(c(length(sampleFiles.y):length(gg), 1:length(sampleFiles.x)))]
    
    x <- gg
    
    condition.e <- factor(c(rep(name.x, length(sampleFiles.x)), rep(name.y, length(sampleFiles.y))))
    condition.e <- relevel(condition.e, ref="CONTROL")
    design<-model.matrix(~condition.e) 
    
    #group <- c(rep(1,length(sampleFiles.x)),rep(2,length(sampleFiles.y)))
    
    y <- DGEList(counts=x,group=condition.e)
    
    #We filter out lowly expressed genes using the following commands:
    #keep <- rowSums(cpm(y)>1) >= 2
    #y <- y[keep, , keep.lib.sizes=FALSE]
    
    #I assumed that "Healthy" is set as the first level of the factor 'Disease'.  I assumed that because that is how you set it up in the original code you posted.
    
    #You can make this so by typing
    
    #samples$Disease <- relevel(samples$Disease, ref="Healthy")
    
    #The coefficient "Treatmentnpc" gives the effect of npc in healthy patients.
    
    
    ###Calculate Dispersions
    y <- calcNormFactors(y,method=c("TMM"))
    #?calcNormFactors
    design <- model.matrix(~condition.e)
    colnames(design) <- levels(condition.e)
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTrendedDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    fit <- glmFit(y, design)
    
    #Run the likelihood ratio test on all contrasts
    lrt <- glmLRT(fit, coef=ncol(fit$design))
    
    #table(lrt$table$PValue < 0.05)
    
    #?estimateDisp
    #estimateGLMCommonDisp(y,design, verbose=TRUE)
    
    #A likelihood-ratio test might be preferred as being more powerful, and the test statistic might not be a monotone function of the one above.
    #To perform likelihood ratio tests:
    #fit <- glmQLFit(y,design)
    #######lrt <- glmLRT(fit,coef=ncol(fit$design))
    
    #qlf <- glmQLFit(y,design, robust=TRUE)
    
    #Alternatively, one can perform QL F-test to test for differential expression. While the likelihood ratio test is a more obvious choice for inferences with GLMs, the QL F-test is preferred as it reflects the uncertainty in estimating the dispersion for each gene. It provides more robust and reliable error rate control when the number of replicates is small. The QL dispersion estimation and hypothesis testing can be done by using the functions glmQLFit() and glmQLFTest()
    
    #?glmQLFTest
    #coef= integer or character index vector indicating which coefficients of the linear model are to be tested equal to zero. Ignored if contrast is not NULL
    
    
    ##If you want to be conservative and safe, use QL. If you want to be a bit more aggressive and find all the DE genes you can, still with good FDR control, then use glmLRT. You have already said that you prefer the latter, so use glmLRT.
    #qlf <- glmQLFTest(fit,coef=ncol(fit$design))
    
    #topTags(qlf)
    
    #We think that, if variability and uncertainty in estimation needs to be considered, this would be better achieved with the quasi-likelihood framework. Indeed, I use glmQLFit and glmQLFTest for most of my analyses (RNA-seq, ChIP-seq and Hi-C), only using the chi-squared test in glmLRT when I don't have any replicates.
    
    #indeed, the quasi-likelihood framework, i.e., glmQLFit and glmQLFTest, is preferred over the glmLRT in general. In some rare occasions where there are no replicates, one might have to fit GLMs using some prior estimates of NB dispersions and call glmLRT for the DE analysis.
    
    #The glmLRT with test="F" was implemented to account for the uncertainty of the NB dispersion estimates in the LRT framework. Now considering the limited scenarios where it is applicable, we decided to phase out the F-test of the glmLRT in the next release to avoid confusion.
    
    #To perform likelihood ratio tests:
    #fit <- glmFit(y,design)
    #lrt <- glmLRT(fit,coef=2)
    #topTags(lrt)
    
    
    #table(lrt$table$PValue < 0.05)
    
    #y$common.dispersion
    
    #COEFICIENTE VARIACI”N BIOL”GICA
    #sqrt(y$common.dispersion)
    
    exp.BCV.e$BCV[exp.loop]<-sqrt(y$common.dispersion)
    
    #standard_dev_edgeR
    
    #The square root of the common dispersion gives the coefficient of variation of biological variation. Here the common dispersion is found to be 0.6093768, so the coefficient of biological variation is around 0.4.
    
    #The dispersion estimates can be viewed in a BCV plot:
    i.e <- paste0(i.edge, "/", DIR)
    
    i.e.1 <- paste(i.e, "-plotBCV.png", sep="")
    
    png(i.e.1, 1500, 800, pointsize=20)
    plotBCV(y)
    dev.off()
    
    
    #sacamos valores de dentro de qlf y ordenados
    #qlf2 <- as.data.frame(topTags(qlf, n=dim(qlf)[1], adjust.method = "BH", sort.by = "PValue"))
    
    
    lrt2 <- as.data.frame(topTags(lrt, n=dim(lrt)[1], adjust.method = "BH", sort.by = "PValue"))
    #lrt2 <- as.data.frame(topTags(y, n=dim(y)[1], adjust.method = "BH", sort.by = "PValue"))
    
    #o <- order(qlf$table$PValue)
    #cpm(y)[o[1:10],]
    
    #guardamos la tabla ordenada con todos las secuencias y su an·lisis de expresion diferencial
    #table.e <- qlf$table[o,]
    
    
    #We see that all the top genes have consistent tumour vs normal changes for the three patients. The total number of differentially expressed genes at 5% FDR is given by:
    
    #summary(de <- decideTestsDGE(qlf))
    
    #summary(de <- decideTestsDGE(lrt))
    
    
    #Plot log-fold change against log-counts per million, with DE genes highlighted:
    #detags <- rownames(y)[as.logical(de)]
    
    
    
    #etiquetas para pvalor 0.05
    
    lrt.e <- subset(lrt2, lrt2$PValue <= 0.05 & (lrt2$logFC <= -corteLgFC | lrt2$logFC >= corteLgFC))
    detags.e <- rownames(lrt.e)
    
    #lrt.e <- subset(et$table, et$table$PValue <= 0.05 & (et$table$logFC <= -corteLgFC | et$table$logFC >= corteLgFC))
    
    
    #qlf2.e <- subset(qlf2, qlf2$PValue <= 0.05 & (qlf2$logFC <= -corteLgFC | qlf2$logFC >= corteLgFC))
    #detags.e <- rownames(qlf2.e)
    
    
    
    i.e.2 <- paste(i.e, "-plotSmear.png", sep="")
    
    png(i.e.2, 1500, 800, pointsize=20)
    
    
    plotSmear(lrt, de.tags=detags.e)
    abline(h=c(-1, 1), col="blue")
    
    #plotSmear(qlf$table[qlf$table$PValue<=0.05], de.tags=detags)
    #abline(h=c(-1, 1), col="blue")
    dev.off()
    
    
    
    #Data exploration
    #An MDS plot shows the relative similarities of the six samples.
    i.e.3 <- paste(i.e, "-plotMDs.png", sep="")
    
    png(i.e.3, 1500, 800, pointsize=20)
    plotMDS(y)
    dev.off()
    
    
    ###FUTURO. Distances on an MDS plot of a DGEList object correspond to leading log-fold-change between each pair of samples. Leading log-fold-change is the root-mean-square average of the largest log2-fold-changes between each pair of samples. Each pair of samples extracted at each time tend to cluster together, suggesting a batch effect. The hrcc treated samples tend to be below the mock samples for each time, suggesting a treatment effect within each time. The two samples at time 1 are less consistent than at times 2 and 3.
    
    #The square root of dispersion is the coefficient of biological variation (BCV). The common BCV is on the high side, considering that this is a designed experiment using genetically identical plants. The trended dispersion shows a decreasing trend with expression level. At low logCPM, the dispersions are very large indeed.
    
    #Note that only the trended dispersion is used under the quasi-likelihood (QL) pipeline. The tagwise and common estimates are shown here but will not be used further.
    
    #The QL dispersions can be estimated using the glmQLFit function, and then be visualized with the plotQLDisp function.
    
    #i.e.4 <- paste(i.e, "-plotQLDisp.png", sep="")
    
    #png(i.e.4, 1500, 800, pointsize=20)
    #plotQLDisp(fit)
    #dev.off()
    
    
    #Here, a gene is only retained if it is expressed at a count-per-million (CPM) above 0.5 in at least two samples.
    #keep <- rowSums(cpm(y) > 0.5) >= 2
    
    #summary(keep)
    
    #When a negative binomial model is fitted, we need to estimate the BCV(s) before we carry out the analysis. The BCV, as shown in the previous section, is the square root of the dispersion parameter under the negative binomial model. Hence, it is equivalent to estimating the dispersion(s) of the negative binomial model.
    
    #The parallel nature of sequencing data allows some possibilities for borrowing information from the ensemble of genes which can assist in inference about each gene individually. The easiest way to share information between genes is to assume that all genes have the same mean-variance relationship, in other words, the dispersion is the same for all the genes [29]. An extension to this ‚Äúcommon dispersion‚Äù approach is to put a mean-dependent trend on a parameter in the variance function, so that all genes with the same expected count have the same variance.
    
    #However, the truth is that the gene expression levels have non-identical and dependent distribution between genes, which makes the above assumptions too naive. A more general approach that allows genewise variance functions with empirical Bayes shrinkage was intro- duced several years ago [28] and has recently been extended to generalized linear models and thus more complex experimental designs [17]. Only when using tagwise dispersion will genes that are consistent between replicates be ranked more highly than genes that are not. It has been seen in many RNA-Seq datasets that allowing gene-specific dispersion is necessary in order that differential expression is not driven by outliers. Therefore, the tagwise dispersions are strongly recommended in model fitting and testing for differential expression.
    
    #In edgeR, we first estimate a common dispersion for all the tags and then apply an em- pirical Bayes strategy for squeezing the tagwise dispersions towards the common dispersion. The amount of shrinkage is determined by the prior weight given to the common dispersion (or the dispersion trend) and the precision of the tagwise estimates, and can be considered as the prior degrees of freedom. This prior degrees of freedom is estimated by examining the heteroskedasticity of the data [5].
    
    
    ### resumen expr diferencial para cada estrÈs edgeR
    #convertimos rownames de table.e en una columna aparte
    
    ########################################################  datos exp diff
    r.e <- lrt2
    r.e$sRNA <- rownames(r.e)
    
    head(r.e)
    #Cambiamos el orden de las columnas
    r.e <- r.e[,c(length(r.e),1:(length(r.e)-1))]
    rownames(r.e) <- NULL
    
    i.e.srna <- paste0("0-results/",time.DIR,"/sRNA/edgeR/", DIR, ".txt")
    write.table(r.e, i.e.srna, row.names = F, quote = F, sep = "\t")
    
    i.e.srna2 <- paste0("0-results/",time.DIR,"/sRNA/edgeR/", DIR, "-excel.csv")
    
    write.table(r.e, i.e.srna2, row.names = F, quote = F, sep = ";", dec = ",")
    
    
    t.e <- subset(r.e, r.e$PValue <= 0.05 & (r.e$logFC >= corteLgFC | r.e$logFC<= -corteLgFC) & !is.na(r.e$PValue))
    
    #We keep summary data for this stress
    exp.lw.e <- subset(t.e, t.e$logFC <= -corteLgFC)
    exp.up.e <- subset(t.e, t.e$logFC >= corteLgFC)
    
    
    exp.lw.fil.20 <- subset(t.e, ( (nchar(t.e$sRNA)) == 20 & t.e$logFC<= -corteLgFC) )
    exp.up.fil.20 <- subset(t.e, ( (nchar(t.e$sRNA)) == 20 & t.e$logFC>= corteLgFC) )
    exp.lw.fil.21 <- subset(t.e, ( (nchar(t.e$sRNA)) == 21 & t.e$logFC<= -corteLgFC) )
    exp.up.fil.21 <- subset(t.e, ( (nchar(t.e$sRNA)) == 21 & t.e$logFC>= corteLgFC) )
    exp.lw.fil.22 <- subset(t.e, ( (nchar(t.e$sRNA)) == 22 & t.e$logFC<= -corteLgFC) )
    exp.up.fil.22 <- subset(t.e, ( (nchar(t.e$sRNA)) == 22 & t.e$logFC>= corteLgFC) )
    exp.lw.fil.23 <- subset(t.e, ( (nchar(t.e$sRNA)) == 23 & t.e$logFC<= -corteLgFC) )
    exp.up.fil.23 <- subset(t.e, ( (nchar(t.e$sRNA)) == 23 & t.e$logFC>= corteLgFC) )
    exp.lw.fil.24 <- subset(t.e, ( (nchar(t.e$sRNA)) == 24 & t.e$logFC<= -corteLgFC) )
    exp.up.fil.24 <- subset(t.e, ( (nchar(t.e$sRNA)) == 25 & t.e$logFC>= corteLgFC) )
    exp.lw.fil.25 <- subset(t.e, ( (nchar(t.e$sRNA)) == 25 & t.e$logFC<= -corteLgFC) )
    exp.up.fil.25 <- subset(t.e, ( (nchar(t.e$sRNA)) == 25 & t.e$logFC>= corteLgFC) )
    
    
    dd <- data.frame(Size.sRNA=c("20", "20", "21", "21", "22", "22", "23", "23", "24", "24", "25", "25"),
                     State=rep(c("Lw","Up"),6),
                     Count=c(length(exp.lw.fil.20$sRNA), length(exp.up.fil.20$sRNA), length(exp.lw.fil.21$sRNA), length(exp.up.fil.21$sRNA), length(exp.lw.fil.22$sRNA), length(exp.up.fil.22$sRNA), length(exp.lw.fil.23$sRNA), length(exp.up.fil.23$sRNA), length(exp.lw.fil.24$sRNA), length(exp.up.fil.24$sRNA), length(exp.lw.fil.25$sRNA), length(exp.up.fil.25$sRNA)))
    
    name.tab1 <- paste0("0-results/",time.DIR,"/images/edgeR/", DIR,"/sum20_25-", DIR, ".txt")
    write.table(dd, file=name.tab1, quote = F, sep="\t", row.names = F)
    
    #damos nombre al tÌtulo de imagen
    titleplot <- paste0("Sample:", DIR,"\n\nSummary differential expression analysis accounts shown by size\n\nP.value <= 0.05 and (Log2FC <= -1 or Log2FC >= 1)")
    
    p1<-ggplot(dd,aes(x=Size.sRNA,y=Count,fill=State), color=State) +  
      stat_summary(fun.y=mean,position=position_dodge(),geom="bar") +
      ggtitle(titleplot) + 
      ylim(0,8000) + theme(plot.title = element_text(lineheight=.8, face="bold"))
    # p1
    
    i.e.5 <- paste(i.e, "-counts-diff-exp.png", sep="")
    
    
    par(mar=c(0,0,0,0)) # Remove unnecessary margins
    
    p1
    ggsave(i.e.5, width = 10, height = 8, dpi=300)
    #dev.off()
    
    
    
    
    exp.sumALL$Lw.E[exp.loop] <- dim(exp.lw.e)[1]
    exp.sumALL$Up.E[exp.loop] <- dim(exp.up.e)[1]
    exp.sumALL$Total.E[exp.loop] <- dim(t.e)[1]
    
    
    
    
    
    
    #juntamos tablas egde + deseq2 SIN EXPRESIÈìN DIFERENCIAL, TODAS, PARA ALINEAMIENTO MIRBASE
    
    # 
    # t.b <- merge(t.e, t.d, by="sRNA", all = TRUE, incomparables = NA)
    # 
    # t.b.r <- t.b[,c(1,2,3,5,6,7,10:length(t.b))]
    # 
    # colnames(t.b.r[,1:8]) <- c("sRNA", "logFC.E", "logCPM.E", "PValue.E", "baseMean.D", "logFC.D", "PValue.D", "Padj.D")
    # 
    # t.b.r[,c(4,7,8)] <- round(t.b.r[c(4,7,8)],digits=5)
    # t.b.r[,c(2,3,5,6,9:length(t.b.r))] <- round(t.b.r[c(2,3,5,6,9:length(t.b.r))],digits=2)
    # 
    # t.b.r[2:length(t.b.r)][is.na(t.b.r[2:length(t.b.r)])] <- "-"
    # 
    # nameDE.ED <- paste0("0-results/sRNA/DESeq2-edgeR/", DIR, ".tsv")
    # nameDE.ED.2 <- paste0("0-results/sRNA/DESeq2-edgeR/", DIR, "-excel.tsv")
    # 
    # write.table(t.b.r, file=nameDE.ED, quote = F, sep="\t", row.names = F)
    # write.table(t.b.r, file=nameDE.ED.2, quote = F, sep="\t", row.names = F, dec = ",")
    
    
    
    
    ####COMPARARAR LOS TRES como en edge y deseq podemos tener la lista entera... tendremos la de estos, tanto exp diff como no
    ###para todos los del noiseq, que tendr·n menos datos, ya filtrados
    
    
    
    ########################################################################################################
    
    ######## diagramas de ven
    edgeList <- data.frame(sRNA=t.e$sRNA)
    deseqList <- data.frame(sRNA=resdata2$sRNA)
    
    #mynoiseq.deg22 <- mynoiseq.deg[mynoiseq.deg$log2FC >= corteLgFC |  mynoiseq.deg$log2FC <= -corteLgFC,]
    noiseqList <- data.frame(sRNA=rownames(mynoiseq.deg22))
    
    head(list(edgeList))
    
    e <- t.e$sRNA
    d <- resdata2$sRNA
    n <- rownames(mynoiseq.deg22)
    
    x <- list(edgeR=e, DESeq2=d, NOISeq=n)
    
    
    # when x is a list
    png(paste0("0-results/", time.DIR, "/VennDiagrams/", DIR, ".png"), 3200, 2100, pointsize=20, res = 300)
    venn(x, zcolor = "style") 
    dev.off()
    
    
    den <- data.frame(sRNA=d[d %in% e & d %in% n])
    
    #Escribimos unicos tres sitemas expresion diferencial
    write.table(den, paste0("0-results/", time.DIR, "/VennDiagrams/DEN/", DIR, ".txt"), row.names = F, quote = F)
    

    # write.table(edgeList, paste0("0-results/VennDiagrams/", DIR, "/edgeR.txt"), row.names = F, quote = F)
    # write.table(deseqList, paste0("0-results/VennDiagrams/", DIR, "/DESeq2.txt"), row.names = F, quote = F)
    # write.table(noiseqList, paste0("0-results/VennDiagrams/", DIR, "/NOISeq.txt"), row.names = F, quote = F)
    # 
    ########################################################################################################
    
    
    
    r.b <- merge(r.e, resdata, by="sRNA", all = TRUE, incomparables = NA)
    
    r.b.r <- r.b[,c(1,2,3,5,6,7,8,11:length(r.b))]
    
    parte1 <- c("sRNA", "log2FC.E", "logCPM.E", "PValue.E", "FDR.E", "baseMean.D", "log2FC.D", "PValue.D", "FDR.D")
    parte2 <- colnames(r.b.r[,c(10:length(r.b.r))])
    
    dospartes <- c(parte1, parte2)
    colnames(r.b.r) <- dospartes
    
    #redondeamos
    #r.b.r[,c(4,5,8,9)] <- round(r.b.r[c(4,5,8,9)],digits=5)
    
    #r.b.r[,c(2,3,6,7,(10:length(r.b.r)))] <- round(r.b.r[,c(2,3,6,7,10:length(r.b.r))],digits=2)
    
    table(is.na(r.b.r[2:length(r.b.r)]))
    
    #si lo dejamos en NA luego en el reporte HTML no sale
    r.b.r[8:9][is.na(r.b.r[8:9])] <- 1
    
    r.b.r[2:length(r.b.r)][is.na(r.b.r[2:length(r.b.r)])] <- "-"
    
    
    ##Ahora juntamos con NOISeq   mynoiseq.all, sin los de normalizacion, el resto si (ser· limitante para deseq2 y edger)
    mynoiseq.all = degenes(mynoiseqbio, q = 0, M = NULL)
    
    ###Controlamos que el Control estÈ antes que el estrÈs
    #En que posicion est· el control?
    pos <- grep("CONTROL", colnames(mynoiseq.all),value = F)
    
    #Si est· en la segunda posicion, cambiamos el orden y renombramos
    if(pos==2){
      mynoiseq.all <- mynoiseq.all[,c(2,1,3,4,5)]
      colnames(mynoiseq.all) <- c("CONTROL-mean.N", "STRESS-mean.N", "theta.N", "prob.N", "log2FC.N")
    
    #Si est· en la primera posicion tan solo renombramos
    }else if(pos==1){
      colnames(mynoiseq.all) <- c("CONTROL-mean.N", "STRESS-mean.N", "theta.N", "prob.N", "log2FC.N")
    }
    
    mynoiseq.all$sRNA <- rownames(mynoiseq.all)
    rownames(mynoiseq.all) <- NULL
    
    rfinal <- merge(r.b.r, mynoiseq.all, by = "sRNA")
    
    
    nameDE.ED <- paste0("0-results/", time.DIR, "/sRNA/D.E.N/", DIR, ".txt")
    nameDE.ED.2 <- paste0("0-results/", time.DIR, "/sRNA/D.E.N/", DIR, "-excel.csv")
    
    options(scipen=10)
    write.table(rfinal, file=nameDE.ED, quote = F, sep="\t", row.names = F)
    write.table(rfinal, file=nameDE.ED.2, quote = F, sep=";", row.names = F, dec = ",")
    
    options(scipen=0) #RESTORE DEFAULT
    
    #sacamos aquellos que tienen expresion diferencial en los TRES
    edgeList$sRNA <- as.character(edgeList$sRNA)
    deseqList$sRNA <- as.character(deseqList$sRNA)
    noiseqList$sRNA <- as.character(noiseqList$sRNA)
    
    d.e <- rfinal[rfinal$sRNA %in% edgeList$sRNA & rfinal$sRNA %in% deseqList$sRNA & rfinal$sRNA %in% noiseqList$sRNA,]
    
    
    
    #d.e <- subset(r.b.r, (r.b.r$PValue.D <= 0.05 & (r.b.r$logFC.D >= corteLgFC | r.b.r$logFC.D <= -corteLgFC )) & (r.b.r$PValue.E <= 0.05 & (r.b.r$logFC.E >= corteLgFC | r.b.r$logFC.E <= -corteLgFC )) ) 
    
    
    #hacemos grafico secuencias conjuntas
    d.e.lw <- subset(d.e, (d.e$log2FC.D <= -corteLgFC ))#& (d.e$log2FC.E <= -corteLgFC ))
    d.e.up <- subset(d.e , (d.e$log2FC.D >= corteLgFC ))#& (d.e$log2FC.E >= corteLgFC ))
    
    
    exp.sumALL$Total.D.E.N[exp.loop] <- dim(d.e)[1]
    exp.sumALL$Lw.D.E.N[exp.loop] <- dim(d.e.lw)[1]
    exp.sumALL$Up.D.E.N[exp.loop] <- dim(d.e.up)[1]
    
    
    exp.loop=exp.loop+1
    
    setwd(houseDIR)
  }#Fin bucle para cada estrÈs   ### VOLVEMOS A LÈçNEA 352
  
  ##################################################################################
  #fuera del bucle FOR
  
  write.table(exp.sumALL, file=paste0("0-results/",time.DIR, "/sRNA/D.E.N/summary20-25.txt"), quote = F, sep="\t", row.names = F)
  
  #When we have had larger BCV from mouse experiments, we've been able to track the cause down to sex differences or batch effects or RNA 
  #contamination.
  
  write.table(exp.BCV.e, file=paste0("0-results/",time.DIR, "/sRNA/edgeR/BCV20-25.txt"), quote = F, sep="\t", row.names = F)
  
  ##########################################################################################################3
  #DIAGRAMA DE VEN DE 7!!!!
  
  #Solo si tenemos los 7 valores!
  
  files.ven <- dir(paste0("0-results/", time.DIR, "/VennDiagrams/DEN/"))
  namesVEN <- files.ven
  
  #quitamos el txt
  for(j in 1:length(namesVEN)){
    namesVEN[j] <- strsplit(namesVEN[j], ".txt")[[1]][1]
  }
  
  adelaide <- c("C", "D", "SA", "SD", "MON", "HSVd", "A")
  
  if(all(adelaide %in% namesVEN)){
    C=read.table(paste0("0-results/",time.DIR, "/VennDiagrams/DEN/C.txt"), header = T, stringsAsFactors = F)
    D=read.table(paste0("0-results/",time.DIR, "/VennDiagrams/DEN/D.txt"), header = T, stringsAsFactors = F)
    SA=read.table(paste0("0-results/",time.DIR, "/VennDiagrams/DEN/SA.txt"), header = T, stringsAsFactors = F)
    SD=read.table(paste0("0-results/",time.DIR, "/VennDiagrams/DEN/SD.txt"), header = T, stringsAsFactors = F)
    MON=read.table(paste0("0-results/",time.DIR, "/VennDiagrams/DEN/MON.txt"), header = T, stringsAsFactors = F)
    HSVd=read.table(paste0("0-results/",time.DIR, "/VennDiagrams/DEN/HSVd.txt"), header = T, stringsAsFactors = F)
    A=read.table(paste0("0-results/",time.DIR, "/VennDiagrams/DEN/A.txt"), header = T, stringsAsFactors = F)
    
    #Cogemos listas
  
    x <- list(C=C$sRNA, D=D$sRNA, SA=SA$sRNA, SD=SD$sRNA, HSVd=HSVd$sRNA, MON=MON$sRNA, A=A$sRNA)
    
    
    # when x is a list
    
    png(paste0("0-results/",time.DIR, "/VennDiagrams/DEN/Adelaide.png"), 3200, 2100, pointsize=20, res = 300)
    venn(x, zcolor = "style")
    dev.off()
    
  }
  
  
  
  
  save.image(paste0("0-results/", time.DIR,"/SAVE/", time.DIR, "-D-E-N.Rdata"))
  
}##Fin blucle para cada tiempo ###
############################################################# FIN PRUEBA CREAR TABLA


warnings()

################################################################## Juntamos tablas para cada tiempo

# setwd("/media/lab205/hdd/Lab/Cmelo/2-DEN")

#Donde guardaremos los resultados para cada tiempo
dir.create("0-results/ALL")

#Antes de guardar la tabla, aÒadimos ID scmel, segun tabla
SNCMel <- read.table("D://TFM/LLUIS/Analysis/04-Tablas_exp_differencial/T12-4reverse-maestra.tsv", sep = "\t")
colnames(SNCMel) <- c("sRNA", "ID")
SNCMel <- SNCMel[c(2,1)]

for(time.DIR in timeS.DIR){
  # time.DIR = "T4"
  
  controlTablaAdd = 0
  
  files.res <- dir(paste0("0-results/", time.DIR, "/sRNA/D.E.N"))
  
  files.res <- files.res[!grepl("excel", files.res)]
  files.res <- files.res[!grepl("summary", files.res)]
  
  for(file in files.res){
    # file="A.txt"
    reading.table <- read.table(paste0("0-results/", time.DIR, "/sRNA/D.E.N/", file), sep="\t", header = T, stringsAsFactors = F, quote = NULL)
    
    #Arrreglamos nombre columnas para poder juntarlos A.1 A.2 A.3 por STRESS.1 STRESS.2 STRESS.3
    newnames <- colnames(reading.table)
    
    #Como lo unico cierto que tenemos son las 9 primeras y 5 ultimas columnas, para tiempo 1 y 2, siendo rÈplicas de 3 todas no
    #hay problema, se cambia el nombre y ya est·, pero para T4 agro tiene dos muestras, crearemos una nueva en agro pero tendr· 
    #valores de 0, solo es para cuadrar el numero de columnas.
    
    #sabiendo cuantas hay por delante (9) y cuantas por detr·s (5) podemos saber cuantas CONTROLES y estreses tenemos y modificar
    #los nombres
    
    #nombres a cambiar
    chanames <- newnames[c(10:(length(newnames)-5))]
    
    #cuantos son controles
    ncontrol <- length(grep("CONTROL", chanames))
    nNOcontrol <- length(chanames)-ncontrol
   
    #Si son rÈplicas de 3 cambiamos nombre, sino ser· el T4 con 4 controles y 2 en agro
    if((ncontrol == 3 & nNOcontrol == 3) | (ncontrol == 4 & file!="A.txt")){
      chanames <- c(chanames[c(1:ncontrol)], paste0("STRESS.", 1:nNOcontrol))
      
      colnames(reading.table) <- c(newnames[c(1:9)], chanames, newnames[c((length(newnames)-4):(length(newnames)))])
      
      
    }else if (ncontrol == 4 & file=="A.txt"){
      
      parte1 <- reading.table[c(1:(9+ncontrol+nNOcontrol))]
      parte1$x <- 0
      
      #Reescribimos encima
      reading.table <- cbind(parte1, reading.table[c((10+ncontrol+nNOcontrol):length(newnames))] )
      newnames <- colnames(reading.table)
      
      chanames <- c(chanames[c(1:ncontrol)], paste0("STRESS.", 1:(nNOcontrol+1)))
      colnames(reading.table) <- c(newnames[c(1:9)], chanames, newnames[c((length(newnames)-4):(length(newnames)))])
      
    }

    #Separamos estrÈs del nombre fichero y creamos nueva columna
    namestress <- strsplit(file, ".txt")[[1]][1]
    
    reading.table$Stress <- namestress
    reading.table$Time <- time.DIR
    
    #Juntamos tablas
    if( controlTablaAdd == 0){
      
      tabla <- reading.table
      
      controlTablaAdd <- controlTablaAdd + 1
    }else if (controlTablaAdd > 0){
      tabla2 <- reading.table
      
      tabla <- rbind(tabla, tabla2)
    }
    
  }#Bucle para leer, cambiar nombres para cada fichero y juntar tablas para cada Tiempo
  
  
  #Esto lo hacemos por encima del bucle, para no repetirlo en cada tiempo
  # #Antes de guardar la tabla, aÒadimos ID scmel, segun tabla
  # SNCMel <- read.table("../0-ext_data/T12-4reverse-maestra.tsv", sep = "\t")
  # colnames(SNCMel) <- c("sRNA", "ID")
  # SNCMel <- SNCMel[c(2,1)]
  
  t <- merge(SNCMel, tabla, by="sRNA")
  
  t <- t[c(2,1,3:length(t))]
  
  summary(is.na(t$ID))
  
  #Guardamos tabla para cada Tiempo
  write.table(t, paste0("0-results/ALL/", time.DIR, ".txt"), sep="\t", row.names = F, quote = F)
  
}




