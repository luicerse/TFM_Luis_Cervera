################## PCA creation ################################################################
############ Only data from double and triple stress conditions ################################
library(cluster)
library(fpc)
library(rgl)
library(edgeR) #BiocManager::install('edgeR')
library(stats)
## 1. READ THE FILES OF THE sRNA counts

setwd('D://TFM/LLUIS/')

data_C.D <- read.table("Results/03-archivos_cuentas_txt/T4/C.D.txt", sep = "\t", header = T)
data_C.SA <- read.table("Results/03-archivos_cuentas_txt/T4/C.SA.txt", sep = "\t", header = T)
data_C.SA.SD <- read.table("Results/03-archivos_cuentas_txt/T4/C.SA.SD.txt", sep = "\t", header = T)
data_C.SD <- read.table("Results/03-archivos_cuentas_txt/T4/C.SD.txt", sep = "\t", header = T)
data_D.SA <- read.table("Results/03-archivos_cuentas_txt/T4/D.SA.txt", sep = "\t", header = T)
data_D.MON <- read.table("Results/03-archivos_cuentas_txt/T4/D.MON.txt", sep = "\t", header = T)
data_C <- read.table("Results/03-archivos_cuentas_txt/T4/C.txt", sep = "\t", header = T)
data_D <- read.table("Results/03-archivos_cuentas_txt/T4/D.txt", sep = "\t", header = T)
data_MON <- read.table("Results/03-archivos_cuentas_txt/T4/MON.txt", sep = "\t", header = T)
data_SA <- read.table("Results/03-archivos_cuentas_txt/T4/SA.txt", sep = "\t", header = T)
data_SD <- read.table("Results/03-archivos_cuentas_txt/T4/SD.txt", sep = "\t", header = T)
###Primero creamos una columna con todos los datos de controles
head(data_C.D[,c(1:5)])

all_data <- unique(rbind(data_C.D[,c(1:5)], data_C.SA[,c(1:5)], data_C.SA.SD[,c(1:5)], data_C.SD[,c(1:5)], 
                         data_D.MON[,c(1:5)], data_D.SA[,c(1:5)], data_C[,c(1:5)], data_D[,c(1:5)], 
                         data_MON[,c(1:5)], data_SA[,c(1:5)], data_SD[,c(1:5)]))


summary(duplicated(all_data$sRNA))

datax <- merge(all_data, data_C.D[,c(1,6:8)], by = "sRNA", all = T)
datax <- merge(datax, data_C.SA[,c(1,6:8)], by = "sRNA", all = T)
datax <- merge(datax, data_C.SA.SD[,c(1,6:8)], by = "sRNA", all = T)
datax <- merge(datax, data_C.SD[,c(1,6:8)], by = "sRNA", all = T)
datax <- merge(datax, data_D.MON[,c(1,6:8)], by = "sRNA", all = T)
datax <- merge(datax, data_D.SA[,c(1,6:8)], by = "sRNA", all = T)
datax <- merge(datax, data_C[,c(1,6:8)], by = "sRNA", all = T)
datax <- merge(datax, data_D[,c(1,6:8)], by = "sRNA", all = T)
datax <- merge(datax, data_MON[,c(1,6:8)], by = "sRNA", all = T)
datax <- merge(datax, data_SA[,c(1,6:8)], by = "sRNA", all = T)
datax <- merge(datax, data_SD[,c(1,6:8)], by = "sRNA", all = T)

## Change NAs to 0s.
datax[2:38][is.na(datax[2:38])] <- 0

### Write table for PCA

write.table(datax, "Results/05-PCA_cuentas/T4_simples_dobles_triples_stress.txt", row.names = F, quote = F, sep = "\t")

#######################################3

mycounts <- datax
rownames(mycounts) <- datax$sRNA
mycounts <- mycounts[-1]

colnames(mycounts)

condition.e <- factor(rep(c("CONTROL", "C.D", "C.SA", "C.SA.SD", "C.SD", "D.MON", "D.SA", "C", "D", "MON", "SA", "SD"), c(4,3,3,3,3,3,3,3,3,3,3,3)))
condition.e <- relevel(condition.e, ref="CONTROL")
design<-model.matrix(~condition.e) 
design

# library(edgeR)

y <- DGEList(counts=mycounts,group=condition.e)


###Calculate Dispersions
y <- calcNormFactors(y,method=c("TMM"))

normalized.counts <- cpm(y, normalized.lib.sizes=T)

head(normalized.counts)
mycounts2 <- as.data.frame(normalized.counts)
write.table(mycounts2, "Results/05-PCA_cuentas/T4simples_DOB_TRIP-TMM.txt", row.names = T, quote = F, sep = "\t")

# library(stats)
model <- prcomp(t(mycounts2), scale=TRUE)


df <- as.data.frame(model$x)

df$Names <- rownames(df)
df$Names <- gsub("\\.[0-9]", "", df$Names)
df$color <- c(rep("#BDBDBD",4), rep("#009E73",3), rep("#56B4E9",3), 
              rep("#F0E442",3),rep("#FFEE00",3),rep("#0072B2",3), 
              rep("#D55E00",3), rep("#0000FF",3), rep("#CC8877",3), 
              rep("#339900",3), rep("#66FFFF",3), rep("#000000",3))

# Plot PCA

with(df,plot3d(df$PC1, df$PC2, df$PC3, col=df$color,type="s", radius=7, alpha=1,xlab = "", ylab = "", zlab = "", shininess=20))
with(df,aspect3d(1, 1, 1))
with(df,title3d(main = 'PCA of sRNA libraries', sub = NULL, 
                xlab = "PC1", ylab = "PC2", color = "blue",
                zlab = "PC3", font = 2))
with(df,axes3d(col='darkgray'))

# Remove duplicated sample names
gr <- duplicated(df$Names)
empty <- which(gr)
df$Names[empty] = ""
df$Names

# Print names of samples 
with(df,text3d(df$PC1, df$PC2, df$PC3, c(df$Names[1:4], #control
                       rep("",33)),adj = c(1,1.4), cex = 1, font = 2))

with(df,text3d(df$PC1, df$PC2, df$PC3, c(rep("",4), #C.D
               df$Names[5:7], rep("",30)), adj = c(-1,-0.8), cex = 1, font = 2))

with(df,text3d(df$PC1, df$PC2, df$PC3, c(rep("",7), #C.SA
               df$Names[8:10], rep("",27)),adj = c(0.5,-1.4), cex = 1, font = 2))

with(df,text3d(df$PC1, df$PC2, df$PC3, c(rep("",10), #C.SA.SD
               df$Names[11:13], rep("",24)),adj = c(1.5,1.5), cex = 1, font = 2))

with(df,text3d(df$PC1, df$PC2, df$PC3, c(rep("",13), #C.SD
               df$Names[14:16], rep("",21)),adj = c(1.5,2), cex = 1, font = 2))

with(df,text3d(df$PC1, df$PC2, df$PC3, c(rep("",16), #D.MON
               df$Names[17:19], rep("",18)),adj = c(0.9,1.6), cex = 1, font = 2))

with(df,text3d(df$PC1, df$PC2, df$PC3, c(rep("",19), #D.SA
               df$Names[20:22], rep("",15)),adj = c(0.5,-1.4), cex = 1, font = 2))

with(df,text3d(df$PC1, df$PC2, df$PC3, c(rep("",22), #C
               df$Names[23:25], rep("",12)),adj = c(1.7,1), cex = 1, font = 2))

with(df,text3d(df$PC1, df$PC2, df$PC3, c(rep("",25), #D
               df$Names[26:28], rep("",9)),adj = c(0.5,-1.4), cex = 1, font = 2))

with(df,text3d(df$PC1, df$PC2, df$PC3, c(rep("",28), #MON
               df$Names[29:31], rep("",6)),adj = c(0.5,-1.4), cex = 1, font = 2))

with(df,text3d(df$PC1, df$PC2, df$PC3, c(rep("",31), #SA
               df$Names[32:34], rep("",3)),adj = c(0,-1.4), cex = 1, font = 2))

with(df,text3d(df$PC1, df$PC2, df$PC3, c(rep("",34), #SD
               df$Names[35:37]),adj = c(-1.3,1.7), cex = 1, font = 2))

with(df,view3d(theta=45, phi=20, fov=20, zoom=0.7))
with(df,clear3d(type='lights'))
with(df,light3d(-45, 20, ambient='black', diffuse='#dddddd', specular='gray'))
with(df,light3d(60, 30, ambient='#dddddd', diffuse='#dddddd', specular='black'))


# Drawing ellipses of confidence around the group of 3 samples.
# CONTROL
ellipseControl <- ellipse3d(cov(df[1:4,c(1,2,3,4)]), 
          centre = c(mean(df$PC1[1:4]), mean(df$PC2[1:4]), mean(df$PC3[1:4])), 
          level = 0.95)
with(df, shade3d(ellipseControl, alpha = 0.3, col = '#66BBEE', lit = FALSE))

# C.D.
ellipseC.D <- ellipse3d(cov(df[5:7,c(1,2,4,3)]), 
                            centre = c(mean(df$PC1[5:7]), mean(df$PC2[5:7]), mean(df$PC3[5:7])), 
                            level = 0.95)
with(df, shade3d(ellipseC.D, alpha = 0.3, col = '#66BBEE', lit = FALSE))

# C.SA.
ellipseC.SA <- ellipse3d(cov(df[8:10,c(1,2,3)]), 
                         centre = c(mean(df$PC1[8:10]), mean(df$PC2[8:10]), mean(df$PC3[8:10])), 
                         level = 0.95)
with(df, shade3d(ellipseC.SA, alpha = 0.3, col = '#66BBEE', lit = FALSE))

# C.SA.SD.
ellipseC.SA.SD <- ellipse3d(cov(df[11:13,c(4,1,2,3)]), 
                        centre = c(mean(df$PC1[11:13]), mean(df$PC2[11:13]), mean(df$PC3[11:13])), 
                        level = 0.95)
with(df, shade3d(ellipseC.SA.SD, alpha = 0.3, col = '#66BBEE', lit = FALSE))

# C.SD.
ellipseC.SD <- ellipse3d(cov(df[14:16,c(3,2,1,4)]), 
                        centre = c(mean(df$PC1[14:16]), mean(df$PC2[14:16]), mean(df$PC3[14:16])), 
                        level = 0.95)
with(df, shade3d(ellipseC.SD, alpha = 0.3, col = '#66BBEE', lit = FALSE))

# D.MON.
ellipseD.MON <- ellipse3d(cov(df[14:16,c(2,1,3,4)]), 
                        centre = c(mean(df$PC1[17:19]), mean(df$PC2[17:19]), mean(df$PC3[17:19])), 
                        level = 0.95)
with(df, shade3d(ellipseD.MON, alpha = 0.3, col = '#66BBEE', lit = FALSE))

# D.SA.
ellipseD.SA <- ellipse3d(cov(df[20:22,1:3]), 
                        centre = c(mean(df$PC1[20:22]), mean(df$PC2[20:22]), mean(df$PC3[20:22])), 
                        level = 0.95)
with(df, shade3d(ellipseD.SA, alpha = 0.3, col = '#66BBEE', lit = FALSE))

# C
ellipseC <- ellipse3d(cov(df[23:25,c(2,1,3)]), 
                         centre = c(mean(df$PC1[23:25]), mean(df$PC2[23:25]), mean(df$PC3[23:25])), 
                         level = 0.95)
with(df, shade3d(ellipseC, alpha = 0.3, col = '#66BBEE', lit = FALSE))

# D
ellipseD <- ellipse3d(cov(df[26:28, c(3,1,2,4)]), 
                      centre = c(mean(df$PC1[26:28]), mean(df$PC2[26:28]), 
                                 mean(df$PC3[26:28]), mean(df$PC4[26:28])), 
                      level = 0.95)
with(df, shade3d(ellipseD, alpha = 0.3, col = '#66BBEE', lit = FALSE))

# MON
ellipseMON <- ellipse3d(cov(df[29:31,c(2,1,3)]), 
                      centre = c(mean(df$PC1[29:31]), mean(df$PC2[29:31]), mean(df$PC3[29:31])), 
                      level = 0.95)
with(df, shade3d(ellipseMON, alpha = 0.3, col = '#66BBEE', lit = FALSE))

# SA
ellipseSA <- ellipse3d(cov(df[32:34,c(1,3,4,2)]), 
                      centre = c(mean(df$PC1[32:34]), mean(df$PC2[32:34]), mean(df$PC3[32:34])), 
                      level = 0.95)
with(df, shade3d(ellipseSA, alpha = 0.3, col = '#66BBEE', lit = FALSE))

# SD
ellipseSD <- ellipse3d(cov(df[35:37,c(1,2,3,4)]), 
                      centre = c(mean(df$PC1[35:37]), mean(df$PC2[35:37]), mean(df$PC3[35:37])), 
                      level = 0.95)
with(df, shade3d(ellipseSD, alpha = 0.3, col = '#66BBEE', lit = FALSE))

# Reshape the size of the window.
par3d( windowRect=c(20,30, 1000, 1000)  )

browseURL(
  paste("file://", writeWebGL(dir=file.path(paste0("D://TFM/LLUIS/Results/05-PCA_cuentas/"), paste0("PCA3D-simple-double-triple")),
                              width=1000), sep="")
)
