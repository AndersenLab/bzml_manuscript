#BZ and ML correlations 
#By: ESS
#last updated 3/6/23, R version 4.1.1

#set the repository and wd
options(repos = "https://cloud.r-project.org/")
setwd("/projects/b1059/projects/Etta/BZML_correlations/")

#install required packages
install.packages('data.table')
install.packages('tidyverse')
install.packages('corrplot')
install.packages('viridis')
install.packages('pheatmap')
install.packages('grDevices')
install.packages('Hmisc')

#load required packages
library(data.table) #fread
library(tidyverse) #reduce function
library(corrplot) #triangle plots
library(viridis) #colors
library(pheatmap) #clustered heatmaps
library(grDevices) #saving pdfs and pngs
library(Hmisc) #correlations and p values

#load BZ dataframes, excluding cv traits
df.albendazole <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files/20221103_Albendazole_traitfile.tsv", drop = c("Albendazole_CV_length"))
df.benomyl <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files/20221103_Benomyl_traitfile.tsv", drop = c("Benomyl_CV_length")) 
df.fenbendazole <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files/20221103_Fenbendazole_traitfile.tsv", drop = c("Fenbendazole_CV_length"))
df.mebendazole <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files/20221103_Mebendazole_traitfile.tsv", drop = c("Mebendazole_CV_length"))
df.thiabendazole <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files/20221103_Thiabendazole_traitfile.tsv", drop = c("Thiabendazole_CV_length"))
df.ricobendazole <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files/20221103_Ricobendazole_traitfile.tsv", drop = c("Ricobendazole_CV_length"))

#load ML dataframes, excluding cv traits
df.abamectin <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files/20221103_Abamectin_traitfile.tsv", drop = c("Abamectin_CV_length"))
df.doramectin <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files/20221103_Doramectin_traitfile.tsv", drop = c("Doramectin_CV_length"))
df.eprinomectin <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files/20221103_Eprinomectin_traitfile.tsv", drop = c("Eprinomectin_CV_length"))
df.ivermectin_0.008 <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files/20221103_Ivermectin_0.008_traitfile.tsv", drop = c("Ivermectin_0.008_CV_length"))
df.ivermectin <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files/20221103_Ivermectin_traitfile.tsv", drop = c("Ivermectin_CV_length"))
df.milbemycin <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files/20221103_Milbemycin_traitfile.tsv", drop = c("Milbemycin_CV_length"))
df.moxidectin <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files/20221103_Moxidectin_traitfile.tsv", drop = c("Moxidectin_CV_length"))
df.selamectin <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files/20221103_Selamectin_traitfile.tsv", drop = c("Selamectin_CV_length"))

#load control dataframes, excluding cv traits
df.dmso <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files/20221103_DMSO_traitfile.tsv", drop = c("DMSO_CV_length"))
df.water <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files/20221103_Water_traitfile.tsv", drop = c("Water_CV_length"))

#makes a data frame from whatever drugs you enter, with strains as row names 
merger <- function (df1, df2, ...) { 
list.int  <-   list(df1, df2, ...)
df.int <- list.int %>% reduce(full_join, by='strain')
column_to_rownames(df.int, 'strain')
}

#merges data frames of interest
df.BZs <- merger(df.albendazole, df.benomyl, df.fenbendazole, df.mebendazole, 
                 df.ricobendazole, df.thiabendazole) #BZs
df.MLs <- merger(df.abamectin, df.doramectin, df.eprinomectin, df.ivermectin, #all MLs
              df.ivermectin_0.008, df.milbemycin, df.moxidectin, df.selamectin) 
df.MLs.alt <- merger(df.abamectin, df.doramectin, df.eprinomectin, #all MLs, omit second ivm dose
                     df.ivermectin, df.milbemycin, df.moxidectin, df.selamectin)
df.conts <- merger(df.dmso, df.water) #controls
df.full <- merger(df.albendazole, df.benomyl, df.fenbendazole, df.mebendazole, 
                  df.ricobendazole, df.thiabendazole, df.abamectin, df.doramectin, df.eprinomectin, df.ivermectin, #all MLs
                  df.ivermectin_0.008, df.milbemycin, df.moxidectin, df.selamectin)
df.full.alt <- merger(df.albendazole, df.benomyl, df.fenbendazole, df.mebendazole, 
                      df.ricobendazole, df.thiabendazole, df.abamectin, df.doramectin, df.eprinomectin, df.ivermectin, #all MLs
                      df.milbemycin, df.moxidectin, df.selamectin) #all drugs, no controls, omit second ivm dose

#renames _length from BZ names, omits na
df.BZs.int  <- rename(df.BZs, 
                      Albendazole = Albendazole_length,
                      Benomyl = Benomyl_length,
                      Fendbendazole = Fenbendazole_length,
                      Mebendazole = Mebendazole_length,
                      Ricobendazole = Ricobendazole_length,
                      Thiabendazole = Thiabendazole_length)
df.BZs <- na.omit(df.BZs.int)
rm(df.BZs.int)

#renames _length from ML names, omits na
df.MLs.int  <- rename(df.MLs, 
                      Abamectin = Abamectin_length,
                      Doramectin = Doramectin_length,
                      Eprinomectin = Eprinomectin_length,
                      Ivermectin = Ivermectin_length,
                      Ivermectin_0.008 = Ivermectin_0.008_length,
                      Milbemycin = Milbemycin_length,
                      Moxidectin = Moxidectin_length,
                      Selamectin = Selamectin_length)
df.MLs <- na.omit(df.MLs.int)
rm(df.MLs.int)

#(no ivm 0.008), renames _length from ML names, omits na
df.MLs.alt.int  <- rename(df.MLs.alt, 
                          Abamectin = Abamectin_length,
                          Doramectin = Doramectin_length,
                          Eprinomectin = Eprinomectin_length,
                          Ivermectin = Ivermectin_length,
                          Milbemycin = Milbemycin_length,
                          Moxidectin = Moxidectin_length,
                          Selamectin = Selamectin_length)
df.MLs.alt <- na.omit(df.MLs.alt.int)
rm(df.MLs.alt.int)

#renames _length from drug names, omits na
df.full.int  <- rename(df.full, 
         Abamectin = Abamectin_length,
         Albendazole = Albendazole_length,
         Benomyl = Benomyl_length,
         Doramectin = Doramectin_length,
         Eprinomectin = Eprinomectin_length,
         Fendbendazole = Fenbendazole_length,
         Ivermectin = Ivermectin_length,
         Ivermectin_0.008 = Ivermectin_0.008_length,
         Mebendazole = Mebendazole_length,
         Milbemycin = Milbemycin_length,
         Moxidectin = Moxidectin_length,
         Ricobendazole = Ricobendazole_length,
         Selamectin = Selamectin_length,
         Thiabendazole = Thiabendazole_length)
df.full <- na.omit(df.full.int)
rm(df.full.int)

#renames _length from drug names, omits na, omits ivm 0.008
df.full.alt.int  <- rename(df.full.alt, 
                       Abamectin = Abamectin_length,
                       Albendazole = Albendazole_length,
                       Benomyl = Benomyl_length,
                       Doramectin = Doramectin_length,
                       Eprinomectin = Eprinomectin_length,
                       Fendbendazole = Fenbendazole_length,
                       Ivermectin = Ivermectin_length,
                       Mebendazole = Mebendazole_length,
                       Milbemycin = Milbemycin_length,
                       Moxidectin = Moxidectin_length,
                       Ricobendazole = Ricobendazole_length,
                       Selamectin = Selamectin_length,
                       Thiabendazole = Thiabendazole_length)
df.full.alt <- na.omit(df.full.alt.int)
rm(df.full.alt.int)

#perform spearman correlation, save matrices of correlations and p-values 
cor_matrix_csv <- function(df) {
  results <- rcorr(as.matrix(df), type = "spearman") #calculate spearman correlations and p-values
    correlations <- results$r #extract correlations
    p_values <- results$P #extract p-values
  
  cor_csv_name <- paste0("/projects/b1059/projects/Etta/BZML_correlations/data/cor_matrices/spearman.matrix.", 
                         gsub("^df\\.", "", deparse(substitute(df))), ".csv") #dynamic naming
  pval_csv_name <- paste0("/projects/b1059/projects/Etta/BZML_correlations/data/cor_matrices/spearman.p.values.",
                         gsub("^df\\.", "", deparse(substitute(df))), ".csv") #dynamic naming
  
  write.csv(correlations, file = cor_csv_name, row.names = TRUE) #save correlation matrix as csv
  write.csv(p_values, file = pval_csv_name, row.names = TRUE) #save p-values as csv
  
  assign(paste0("spearman.matrix.", gsub("^df\\.", "", deparse(substitute(df)))),
         correlations, envir = .GlobalEnv) #save correlation matrix as object
  assign(paste0("spearman.p.values.", gsub("^df\\.", "", deparse(substitute(df)))),
         p_values, envir = .GlobalEnv) #save p-values as object
  
  return(list(spearman_matrix = correlations, spearman_p_values = p_values)) #return correlations and p-values as a list

}

cor_matrix_csv(df.BZs) #make BZ matrices
cor_matrix_csv(df.MLs) #make ML matrices
cor_matrix_csv(df.MLs.alt) #make ML matrices, omit ivm 0.008
cor_matrix_csv(df.full) #make full matrices
cor_matrix_csv(df.full.alt) #make full matrices, omit ivm 0.008

#function that plots any correlation matrix, (matrix, method = '')
anyplot <- function(x, plot = '', method = '', gaps = T){
  if(plot == 'pheatmap' & gaps == T) {pheatmap(x, color = viridis(256, option = "D"), fontsize = 15,
                                   scale = "none", limits = c(-0.2, 1), breaks = seq(-0.2, 1, length.out = 256),
                                   angle_col = "45", cellwidth = 40, cellheight = 40, display_numbers = T, number_color = 'white',
                                   fontsize_number = 15,
                                   border_color = "white",
                                   cutree_rows = 2, cutree_cols = 2
                                 ) 
  }
  if(plot == 'pheatmap' & gaps == F){pheatmap(x, color = viridis(256, option = "D"), fontsize = 15,
                                                 scale = "none", limits = c(-0.2, 1), breaks = seq(-0.2, 1, length.out = 256),
                                                 angle_col = "45", cellwidth = 40, cellheight = 40, display_numbers = T, number_color = 'white',
                                                 fontsize_number = 15,
                                                 border_color = "white"
  ) 
  }
  
  if(plot == 'corrplot') {corrplot(x, method = method, order = "original", type = "lower", 
                                   diag = TRUE, addCoef.col = 1, number.cex = 1.8,
                                   col.lim = c(-1, 1),
                                   tl.col = 'black', tl.cex = 2.5, cl.cex = 2, col = viridis(256, option = "D"), 
                                   addgrid.col = "black")}
}

#function that generates a heatmap and saves it as a pdf or tiff
saveplot <- function(df, plot = '', method = '', type = '', gaps, width = '', height = ''){
  file.name <- paste(deparse(substitute(df)), ".", gsub("'", "", plot), ".", gsub("'", "", method), sep = "")
  if(type == 'pdf'){pdf(file = paste0("/projects/b1059/projects/Etta/BZML_correlations/figures/", file.name, "pdf"), width = width, height = height)}
  if(type == 'tiff'){tiff(file = paste0("/projects/b1059/projects/Etta/BZML_correlations/figures/", file.name, "tiff"), width = width, height = height, res = 500)}
  anyplot(df, plot, method, gaps)
  dev.off()
}

#generate plot, choose either 'pheatmap' for clustered heatmap with dendrogram
#or 'corrplot' for unclustered plot
saveplot(spearman.matrix.BZs, plot = 'pheatmap', type = 'pdf', width = 7, height = 6, gaps = F) 
saveplot(spearman.matrix.MLs, plot = 'pheatmap', type = 'pdf', width = 8.5, height = 7.25, gaps = F)
saveplot(spearman.matrix.MLs.alt, plot = 'pheatmap', type = 'pdf', width = 7.5, height = 6.5, gaps = F)
saveplot(spearman.matrix.full, plot = 'pheatmap', type = 'pdf', width = 12, height = 10.5, gaps = T)
saveplot(spearman.matrix.full.alt, plot = 'pheatmap', type = 'pdf', width = 12, height = 10.5, gaps = T)
saveplot(spearman.matrix.BZs, plot = 'pheatmap', type = 'tiff', width = 3500, height = 3250, gaps = F)
saveplot(spearman.matrix.MLs, plot = 'pheatmap', type = 'tiff', width = 4000, height = 3500, gaps = F)
saveplot(spearman.matrix.MLs.alt, plot = 'pheatmap', type = 'tiff', width = 3750, height = 3500, gaps = F)
saveplot(spearman.matrix.full, plot = 'pheatmap', type = 'tiff', width = 5750, height = 5500, gaps = T)
saveplot(spearman.matrix.full.alt, plot = 'pheatmap', type = 'tiff', width = 5500, height = 5000, gaps = T)
