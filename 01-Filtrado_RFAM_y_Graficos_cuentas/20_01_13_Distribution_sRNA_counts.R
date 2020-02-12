library(stringr)

distribution.txt <- read.csv('D://TFM/LLUIS/Analysis/01-Rfam/ALL-freq-for-latex.tsv',
                             sep = '\t', header = T)
#colnames(distribution.txt) <- c('length', 'counts', 'time', 'library')

distribution.txt <- distribution.txt[distribution.txt$time == 'T4', ]
distribution.txt$'library' <- t(data.frame(str_split(distribution.txt$sample, '-')))[,1]


final_df <- aggregate(counts ~ library + nt,
          distribution.txt,
          sum)
final_df <- final_df[order(final_df$library),]
write.csv(final_df, 'D://TFM/LLUIS/Results/04-graficos_cuentas_librerias/df_cuentas.csv')
