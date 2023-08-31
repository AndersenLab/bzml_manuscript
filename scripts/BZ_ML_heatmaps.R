#BZ and ML correlations 
#By: ESS
#last updated 5/9/23, R version 4.1.1

#set the repository and wd
options(repos = "https://cloud.r-project.org/")
setwd("/projects/b1059/projects/Etta/BZML_correlations/")

#install required packages
#install.packages('data.table')
#install.packages('tidyverse')
#install.packages('corrplot')
#install.packages('viridis')
#install.packages('pheatmap')
#install.packages('grDevices')
#install.packages('Hmisc')
#install.packages('cowplot')

#load required packages
library(data.table) #fread
library(tidyverse) #reduce function
library(corrplot) #triangle plots
library(viridis) #colors
library(pheatmap) #clustered heatmaps
library(grDevices) #saving pdfs and pngs
library(Hmisc) #correlations and p values
library(cowplot) #comparing correlations from old and new trait files
library(ggplotify)

#####
#new trait files
#load new BZ dataframes, excluding cv traits
df.albendazole.new <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files_20230330/20230322_Albendazole_traitfile.tsv", drop = c("Albendazole_CV_length"))
df.benomyl.new <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files_20230330/20230322_Benomyl_traitfile.tsv", drop = c("Benomyl_CV_length")) 
df.fenbendazole.new <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files_20230330/20230322_Fenbendazole_traitfile.tsv", drop = c("Fenbendazole_CV_length"))
df.mebendazole.new <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files_20230330/20230322_Mebendazole_traitfile.tsv", drop = c("Mebendazole_CV_length"))
df.thiabendazole.new <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files_20230330/20230322_Thiabendazole_traitfile.tsv", drop = c("Thiabendazole_CV_length"))
df.ricobendazole.new <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files_20230330/20230322_Ricobendazole_traitfile.tsv", drop = c("Ricobendazole_CV_length"))

#load ML dataframes, excluding cv traits
df.abamectin.new <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files_20230330/20230322_Abamectin_traitfile.tsv", drop = c("Abamectin_CV_length"))
df.doramectin.new <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files_20230330/20230322_Doramectin_traitfile.tsv", drop = c("Doramectin_CV_length"))
df.eprinomectin.new <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files_20230330/20230322_Eprinomectin_traitfile.tsv", drop = c("Eprinomectin_CV_length"))
df.ivermectin_0.008.new <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files_20230330/20230322_Ivermectin_0.008_traitfile.tsv", drop = c("Ivermectin_0.008_CV_length"))
df.ivermectin.new <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files_20230330/20230322_Ivermectin_traitfile.tsv", drop = c("Ivermectin_CV_length"))
df.milbemycin.new <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files_20230330/20230322_Milbemycin_traitfile.tsv", drop = c("Milbemycin_CV_length"))
df.moxidectin.new <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files_20230330/20230322_Moxidectin_traitfile.tsv", drop = c("Moxidectin_CV_length"))
df.selamectin.new <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files_20230330/20230322_Selamectin_traitfile.tsv", drop = c("Selamectin_CV_length"))

#load control dataframes, excluding cv traits
#df.dmso.new <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files_20230330/20230322_DMSO_traitfile.tsv", drop = c("DMSO_CV_length"))
#df.water.new <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/trait_files_20230330/20230322_Water_traitfile.tsv", drop = c("Water_CV_length"))

#makes a data frame from whatever drugs you enter, with strains as row names 
merger <- function (df1, df2, ...) { 
  list.int  <-   list(df1, df2, ...)
  df.int <- list.int %>% reduce(full_join, by='strain')
  column_to_rownames(df.int, 'strain')
}

#merges data frames of interest
df.BZs.new <- merger(df.albendazole.new, df.benomyl.new, df.fenbendazole.new, df.mebendazole.new, 
                     df.thiabendazole.new,
                     df.ricobendazole.new) #BZs
df.MLs.new <- merger(df.abamectin.new, df.doramectin.new, df.eprinomectin.new, df.ivermectin.new, #all MLs
                     df.ivermectin_0.008.new, df.milbemycin.new, df.moxidectin.new, df.selamectin.new) 
df.MLs.alt.new <- merger(df.abamectin.new, df.doramectin.new, df.eprinomectin.new, #all MLs, omit second ivm dose
                         df.ivermectin.new, df.milbemycin.new, df.moxidectin.new, df.selamectin.new)

#df.conts.new <- merger(df.dmso.new, df.water.new) #controls 
df.full.new <- merger(df.albendazole.new, df.benomyl.new, df.fenbendazole.new, df.mebendazole.new, df.ricobendazole.new, 
                      df.thiabendazole.new, df.abamectin.new, df.doramectin.new, df.eprinomectin.new, df.ivermectin.new, #all MLs
                      df.ivermectin_0.008.new, df.milbemycin.new, df.moxidectin.new, df.selamectin.new)
df.full.alt.new <- merger(df.albendazole.new, df.benomyl.new, df.fenbendazole.new, df.mebendazole.new, df.ricobendazole.new,
                          df.thiabendazole.new, df.abamectin.new, df.doramectin.new, df.eprinomectin.new, df.ivermectin.new, #all MLs
                          df.milbemycin.new, df.moxidectin.new, df.selamectin.new) #all drugs, no controls, omit second ivm dose

#renames _length from BZ names, omits na
df.BZs.int.new  <- rename(df.BZs.new, 
                          Albendazole = Albendazole_length,
                          Benomyl = Benomyl_length,
                          Fendbendazole = Fenbendazole_length,
                          Mebendazole = Mebendazole_length,
                          Ricobendazole = Ricobendazole_length,
                          Thiabendazole = Thiabendazole_length)
df.BZs.new <- na.omit(df.BZs.int.new)
rm(df.BZs.int.new)

#renames _length from ML names, omits na
df.MLs.int.new  <- rename(df.MLs.new, 
                          Abamectin = Abamectin_length,
                          Doramectin = Doramectin_length,
                          Eprinomectin = Eprinomectin_length,
                          Ivermectin = Ivermectin_length,
                          Ivermectin_0.008 = Ivermectin_0.008_length,
                          Milbemycin = Milbemycin_length,
                          Moxidectin = Moxidectin_length,
                          Selamectin = Selamectin_length)
df.MLs.new <- na.omit(df.MLs.int.new)
rm(df.MLs.int.new)

#(no ivm 0.008), renames _length from ML names, omits na
df.MLs.alt.int.new  <- rename(df.MLs.alt.new, 
                              Abamectin = Abamectin_length,
                              Doramectin = Doramectin_length,
                              Eprinomectin = Eprinomectin_length,
                              Ivermectin = Ivermectin_length,
                              Milbemycin = Milbemycin_length,
                              Moxidectin = Moxidectin_length,
                              Selamectin = Selamectin_length)
df.MLs.alt.new <- na.omit(df.MLs.alt.int.new)
rm(df.MLs.alt.int.new)

#renames _length from drug names, omits na
df.full.int.new  <- rename(df.full.new, 
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
df.full.new <- na.omit(df.full.int.new)
rm(df.full.int.new)

#renames _length from drug names, omits na, omits ivm 0.008
df.full.alt.int.new  <- rename(df.full.alt.new, 
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
df.full.alt.new <- na.omit(df.full.alt.int.new)
rm(df.full.alt.int.new)

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

cor_matrix_csv(df.BZs.new) #make BZ matrices
cor_matrix_csv(df.MLs.new) #make ML matrices
cor_matrix_csv(df.MLs.alt.new) #make ML matrices, omit ivm 0.008
cor_matrix_csv(df.full.new) #make full matrices
cor_matrix_csv(df.full.alt.new) #make full matrices, omit ivm 0.008


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
  if(type == 'pdf'){pdf(file = paste0("/projects/b1059/projects/Etta/BZML_correlations/figures/heatmaps/", file.name, "pdf"), width = width, height = height)}
  if(type == 'tiff'){tiff(file = paste0("/projects/b1059/projects/Etta/BZML_correlations/figures/heatmaps/", file.name, "tiff"), width = width, height = height, res = 500)}
  anyplot(df, plot, method, gaps)
  dev.off()
}

#generate plot, choose either 'pheatmap' for clustered heatmap with dendrogram
#or 'corrplot' for unclustered plot
saveplot(spearman.matrix.BZs.new, plot = 'pheatmap', type = 'pdf', width = 7, height = 6, gaps = F) 
saveplot(spearman.matrix.MLs.new, plot = 'pheatmap', type = 'pdf', width = 8.5, height = 7.25, gaps = F)
saveplot(spearman.matrix.MLs.alt.new, plot = 'pheatmap', type = 'pdf', width = 7.5, height = 6.5, gaps = F)
saveplot(spearman.matrix.full.new, plot = 'pheatmap', type = 'pdf', width = 12, height = 10.5, gaps = T)
saveplot(spearman.matrix.full.alt.new, plot = 'pheatmap', type = 'pdf', width = 12, height = 10.5, gaps = T)
saveplot(spearman.matrix.BZs.new, plot = 'pheatmap', type = 'tiff', width = 3500, height = 3250, gaps = F)
saveplot(spearman.matrix.MLs.new, plot = 'pheatmap', type = 'tiff', width = 4000, height = 3500, gaps = F)
saveplot(spearman.matrix.MLs.alt.new, plot = 'pheatmap', type = 'tiff', width = 3750, height = 3500, gaps = F)
saveplot(spearman.matrix.full.new, plot = 'pheatmap', type = 'tiff', width = 5750, height = 5500, gaps = T)
saveplot(spearman.matrix.full.alt.new, plot = 'pheatmap', type = 'tiff', width = 5500, height = 5000, gaps = T)

#####
# old trait files #
##### 
#####
#load old BZ dataframes, excluding cv traits
df.albendazole.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/old_trait_files/20221103_Albendazole_traitfile.tsv", drop = c("Albendazole_CV_length"))
df.benomyl.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/old_trait_files/20221103_Benomyl_traitfile.tsv", drop = c("Benomyl_CV_length")) 
df.fenbendazole.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/old_trait_files/20221103_Fenbendazole_traitfile.tsv", drop = c("Fenbendazole_CV_length"))
df.mebendazole.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/old_trait_files/20221103_Mebendazole_traitfile.tsv", drop = c("Mebendazole_CV_length"))
df.thiabendazole.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/old_trait_files/20221103_Thiabendazole_traitfile.tsv", drop = c("Thiabendazole_CV_length"))
df.ricobendazole.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/old_trait_files/20221103_Ricobendazole_traitfile.tsv", drop = c("Ricobendazole_CV_length"))

#load ML dataframes, excluding cv traits
df.abamectin.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/old_trait_files/20221103_Abamectin_traitfile.tsv", drop = c("Abamectin_CV_length"))
df.doramectin.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/old_trait_files/20221103_Doramectin_traitfile.tsv", drop = c("Doramectin_CV_length"))
df.eprinomectin.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/old_trait_files/20221103_Eprinomectin_traitfile.tsv", drop = c("Eprinomectin_CV_length"))
df.ivermectin_0.008.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/old_trait_files/20221103_Ivermectin_0.008_traitfile.tsv", drop = c("Ivermectin_0.008_CV_length"))
df.ivermectin.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/old_trait_files/20221103_Ivermectin_traitfile.tsv", drop = c("Ivermectin_CV_length"))
df.milbemycin.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/old_trait_files/20221103_Milbemycin_traitfile.tsv", drop = c("Milbemycin_CV_length"))
df.moxidectin.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/old_trait_files/20221103_Moxidectin_traitfile.tsv", drop = c("Moxidectin_CV_length"))
df.selamectin.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/old_trait_files/20221103_Selamectin_traitfile.tsv", drop = c("Selamectin_CV_length"))

#load control dataframes, excluding cv traits
#df.dmso.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/old_trait_files/20221103_DMSO_traitfile.tsv", drop = c("DMSO_CV_length"))
#df.water.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/old_trait_files/20221103_Water_traitfile.tsv", drop = c("Water_CV_length"))

#makes a data frame from whatever drugs you enter, with strains as row names 
merger <- function (df1, df2, ...) { 
  list.int  <-   list(df1, df2, ...)
  df.int <- list.int %>% reduce(full_join, by='strain')
  column_to_rownames(df.int, 'strain')
}

#merges data frames of interest
df.BZs.old <- merger(df.albendazole.old, df.benomyl.old, df.fenbendazole.old, df.mebendazole.old, 
                     df.thiabendazole.old,
                     df.ricobendazole.old) #BZs
df.MLs.old <- merger(df.abamectin.old, df.doramectin.old, df.eprinomectin.old, df.ivermectin.old, #all MLs
                     df.ivermectin_0.008.old, df.milbemycin.old, df.moxidectin.old, df.selamectin.old) 
df.MLs.alt.old <- merger(df.abamectin.old, df.doramectin.old, df.eprinomectin.old, #all MLs, omit second ivm dose
                         df.ivermectin.old, df.milbemycin.old, df.moxidectin.old, df.selamectin.old)

#df.conts.old <- merger(df.dmso.old, df.water.old) #controls 
df.full.old <- merger(df.albendazole.old, df.benomyl.old, df.fenbendazole.old, df.mebendazole.old, df.ricobendazole.old, 
                      df.thiabendazole.old, df.abamectin.old, df.doramectin.old, df.eprinomectin.old, df.ivermectin.old, #all MLs
                      df.ivermectin_0.008.old, df.milbemycin.old, df.moxidectin.old, df.selamectin.old)
df.full.alt.old <- merger(df.albendazole.old, df.benomyl.old, df.fenbendazole.old, df.mebendazole.old, df.ricobendazole.old,
                          df.thiabendazole.old, df.abamectin.old, df.doramectin.old, df.eprinomectin.old, df.ivermectin.old, #all MLs
                          df.milbemycin.old, df.moxidectin.old, df.selamectin.old) #all drugs, no controls, omit second ivm dose

#renames _length from BZ names, omits na
df.BZs.int.old  <- rename(df.BZs.old, 
                          Albendazole = Albendazole_length,
                          Benomyl = Benomyl_length,
                          Fendbendazole = Fenbendazole_length,
                          Mebendazole = Mebendazole_length,
                          Ricobendazole = Ricobendazole_length,
                          Thiabendazole = Thiabendazole_length)
df.BZs.old <- na.omit(df.BZs.int.old)
rm(df.BZs.int.old)

#renames _length from ML names, omits na
df.MLs.int.old  <- rename(df.MLs.old, 
                          Abamectin = Abamectin_length,
                          Doramectin = Doramectin_length,
                          Eprinomectin = Eprinomectin_length,
                          Ivermectin = Ivermectin_length,
                          Ivermectin_0.008 = Ivermectin_0.008_length,
                          Milbemycin = Milbemycin_length,
                          Moxidectin = Moxidectin_length,
                          Selamectin = Selamectin_length)
df.MLs.old <- na.omit(df.MLs.int.old)
rm(df.MLs.int.old)

#(no ivm 0.008), renames _length from ML names, omits na
df.MLs.alt.int.old  <- rename(df.MLs.alt.old, 
                              Abamectin = Abamectin_length,
                              Doramectin = Doramectin_length,
                              Eprinomectin = Eprinomectin_length,
                              Ivermectin = Ivermectin_length,
                              Milbemycin = Milbemycin_length,
                              Moxidectin = Moxidectin_length,
                              Selamectin = Selamectin_length)
df.MLs.alt.old <- na.omit(df.MLs.alt.int.old)
rm(df.MLs.alt.int.old)

#renames _length from drug names, omits na
df.full.int.old  <- rename(df.full.old, 
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
df.full.old <- na.omit(df.full.int.old)
rm(df.full.int.old)

#renames _length from drug names, omits na, omits ivm 0.008
df.full.alt.int.old  <- rename(df.full.alt.old, 
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
df.full.alt.old <- na.omit(df.full.alt.int.old)
rm(df.full.alt.int.old)

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

cor_matrix_csv(df.BZs.old) #make BZ matrices
cor_matrix_csv(df.MLs.old) #make ML matrices
cor_matrix_csv(df.MLs.alt.old) #make ML matrices, omit ivm 0.008
cor_matrix_csv(df.full.old) #make full matrices
cor_matrix_csv(df.full.alt.old) #make full matrices, omit ivm 0.008


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
  if(type == 'pdf'){pdf(file = paste0("/projects/b1059/projects/Etta/BZML_correlations/figures/heatmaps/", file.name, "pdf"), width = width, height = height)}
  if(type == 'tiff'){tiff(file = paste0("/projects/b1059/projects/Etta/BZML_correlations/figures/heatmaps/", file.name, "tiff"), width = width, height = height, res = 500)}
  anyplot(df, plot, method, gaps)
  dev.off()
}

#generate plot, choose either 'pheatmap' for clustered heatmap with dendrogram
#or 'corrplot' for unclustered plot
saveplot(spearman.matrix.BZs.old, plot = 'pheatmap', type = 'pdf', width = 7, height = 6, gaps = F) 
saveplot(spearman.matrix.MLs.old, plot = 'pheatmap', type = 'pdf', width = 8.5, height = 7.25, gaps = F)
saveplot(spearman.matrix.MLs.alt.old, plot = 'pheatmap', type = 'pdf', width = 7.5, height = 6.5, gaps = F)
saveplot(spearman.matrix.full.old, plot = 'pheatmap', type = 'pdf', width = 12, height = 10.5, gaps = T)
saveplot(spearman.matrix.full.alt.old, plot = 'pheatmap', type = 'pdf', width = 12, height = 10.5, gaps = T)
saveplot(spearman.matrix.BZs.old, plot = 'pheatmap', type = 'tiff', width = 3500, height = 3250, gaps = F)
saveplot(spearman.matrix.MLs.old, plot = 'pheatmap', type = 'tiff', width = 4000, height = 3500, gaps = F)
saveplot(spearman.matrix.MLs.alt.old, plot = 'pheatmap', type = 'tiff', width = 3750, height = 3500, gaps = F)
saveplot(spearman.matrix.full.old, plot = 'pheatmap', type = 'tiff', width = 5750, height = 5500, gaps = T)
saveplot(spearman.matrix.full.alt.old, plot = 'pheatmap', type = 'tiff', width = 5500, height = 5000, gaps = T)




#####
#####
#going to experiment with a few ways of distinguishing drugs by subclass 
pheatmap(annotation_row =) 

drug_list <- colnames(df.full.old)
##### 
#old trait files 
#load old BZ dataframes, excluding cv traits
df.albendazole.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/trait_files/20221103_Albendazole_traitfile.tsv", drop = c("Albendazole_CV_length"))
df.benomyl.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/trait_files/20221103_Benomyl_traitfile.tsv", drop = c("Benomyl_CV_length")) 
df.fenbendazole.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/trait_files/20221103_Fenbendazole_traitfile.tsv", drop = c("Fenbendazole_CV_length"))
df.mebendazole.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/trait_files/20221103_Mebendazole_traitfile.tsv", drop = c("Mebendazole_CV_length"))
df.thiabendazole.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/trait_files/20221103_Thiabendazole_traitfile.tsv", drop = c("Thiabendazole_CV_length"))
df.ricobendazole.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/trait_files/20221103_Ricobendazole_traitfile.tsv", drop = c("Ricobendazole_CV_length"))

#load ML dataframes, excluding cv traits
df.abamectin.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/trait_files/20221103_Abamectin_traitfile.tsv", drop = c("Abamectin_CV_length"))
df.doramectin.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/trait_files/20221103_Doramectin_traitfile.tsv", drop = c("Doramectin_CV_length"))
df.eprinomectin.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/trait_files/20221103_Eprinomectin_traitfile.tsv", drop = c("Eprinomectin_CV_length"))
df.ivermectin_0.008.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/trait_files/20221103_Ivermectin_0.008_traitfile.tsv", drop = c("Ivermectin_0.008_CV_length"))
df.ivermectin.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/trait_files/20221103_Ivermectin_traitfile.tsv", drop = c("Ivermectin_CV_length"))
df.milbemycin.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/trait_files/20221103_Milbemycin_traitfile.tsv", drop = c("Milbemycin_CV_length"))
df.moxidectin.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/trait_files/20221103_Moxidectin_traitfile.tsv", drop = c("Moxidectin_CV_length"))
df.selamectin.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/trait_files/20221103_Selamectin_traitfile.tsv", drop = c("Selamectin_CV_length"))

#load control dataframes, excluding cv traits
df.dmso.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/trait_files/20221103_DMSO_traitfile.tsv", drop = c("DMSO_CV_length"))
df.water.old <- fread("/projects/b1059/projects/Etta/BZML_correlations/data/archive/trait_files/20221103_Water_traitfile.tsv", drop = c("Water_CV_length"))

#makes a data frame from whatever drugs you enter, with strains as row names 
merger <- function (df1, df2, ...) { 
list.int  <-   list(df1, df2, ...)
df.int <- list.int %>% reduce(full_join, by='strain')
column_to_rownames(df.int, 'strain')
}

#merges data frames of interest
df.BZs.old <- merger(df.albendazole.old, df.benomyl.old, df.fenbendazole.old, df.mebendazole.old, 
                 df.ricobendazole.old, df.thiabendazole.old) #BZs
df.MLs.old <- merger(df.abamectin.old, df.doramectin.old, df.eprinomectin.old, df.ivermectin.old, #all MLs
              df.ivermectin_0.008.old, df.milbemycin.old, df.moxidectin.old, df.selamectin.old) 
df.MLs.alt.old <- merger(df.abamectin.old, df.doramectin.old, df.eprinomectin.old, #all MLs, omit second ivm dose
                     df.ivermectin.old, df.milbemycin.old, df.moxidectin.old, df.selamectin.old)
df.conts.old <- merger(df.dmso.old, df.water.old) #controls
df.full.old <- merger(df.albendazole.old, df.benomyl.old, df.fenbendazole.old, df.mebendazole.old, 
                  df.ricobendazole.old, df.thiabendazole.old, df.abamectin.old, df.doramectin.old, df.eprinomectin.old, df.ivermectin.old, #all MLs
                  df.ivermectin_0.008.old, df.milbemycin.old, df.moxidectin.old, df.selamectin.old)
df.full.alt.old <- merger(df.albendazole.old, df.benomyl.old, df.fenbendazole.old, df.mebendazole.old, 
                      df.ricobendazole.old, df.thiabendazole.old, df.abamectin.old, df.doramectin.old, df.eprinomectin.old, df.ivermectin.old, #all MLs
                      df.milbemycin.old, df.moxidectin.old, df.selamectin.old) #all drugs, no controls, omit second ivm dose

#renames _length from BZ names, omits na
df.BZs.int.old  <- rename(df.BZs.old, 
                      Albendazole = Albendazole_length,
                      Benomyl = Benomyl_length,
                      Fendbendazole = Fenbendazole_length,
                      Mebendazole = Mebendazole_length,
                      Ricobendazole = Ricobendazole_length,
                      Thiabendazole = Thiabendazole_length)
df.BZs.old <- na.omit(df.BZs.int.old)
rm(df.BZs.int.old)

#renames _length from ML names, omits na
df.MLs.int.old  <- rename(df.MLs.old, 
                      Abamectin = Abamectin_length,
                      Doramectin = Doramectin_length,
                      Eprinomectin = Eprinomectin_length,
                      Ivermectin = Ivermectin_length,
                      Ivermectin_0.008 = Ivermectin_0.008_length,
                      Milbemycin = Milbemycin_length,
                      Moxidectin = Moxidectin_length,
                      Selamectin = Selamectin_length)
df.MLs.old <- na.omit(df.MLs.int.old)
rm(df.MLs.int.old)

#(no ivm 0.008), renames _length from ML names, omits na
df.MLs.alt.int.old  <- rename(df.MLs.alt.old, 
                          Abamectin = Abamectin_length,
                          Doramectin = Doramectin_length,
                          Eprinomectin = Eprinomectin_length,
                          Ivermectin = Ivermectin_length,
                          Milbemycin = Milbemycin_length,
                          Moxidectin = Moxidectin_length,
                          Selamectin = Selamectin_length)
df.MLs.alt.old <- na.omit(df.MLs.alt.int.old)
rm(df.MLs.alt.int.old)

#renames _length from drug names, omits na
df.full.int.old  <- rename(df.full.old, 
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
df.full <- na.omit(df.full.int.old)
rm(df.full.int.old)

#renames _length from drug names, omits na, omits ivm 0.008
df.full.alt.int.old  <- rename(df.full.alt.old, 
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
df.full.alt.old <- na.omit(df.full.alt.int.old)
rm(df.full.alt.int.old)

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

cor_matrix_csv(df.BZs.old) #make BZ matrices
cor_matrix_csv(df.MLs.old) #make ML matrices
cor_matrix_csv(df.MLs.alt.old) #make ML matrices, omit ivm 0.008
cor_matrix_csv(df.full.old) #make full matrices
cor_matrix_csv(df.full.alt.old) #make full matrices, omit ivm 0.008

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
  if(type == 'pdf'){pdf(file = paste0("/projects/b1059/projects/Etta/BZML_correlations/figures/archive/heatmaps/", file.name, "pdf"), width = width, height = height)}
  if(type == 'tiff'){tiff(file = paste0("/projects/b1059/projects/Etta/BZML_correlations/figures/archive/heatmaps/", file.name, "tiff"), width = width, height = height, res = 500)}
  anyplot(df, plot, method, gaps)
  dev.off()
}

#generate plot, choose either 'pheatmap' for clustered heatmap with dendrogram
#or 'corrplot' for unclustered plot
saveplot(spearman.matrix.BZs.old, plot = 'pheatmap', type = 'pdf', width = 7, height = 6, gaps = F) 
saveplot(spearman.matrix.MLs.old, plot = 'pheatmap', type = 'pdf', width = 8.5, height = 7.25, gaps = F)
saveplot(spearman.matrix.MLs.alt.old, plot = 'pheatmap', type = 'pdf', width = 7.5, height = 6.5, gaps = F)
saveplot(spearman.matrix.full.old, plot = 'pheatmap', type = 'pdf', width = 12, height = 10.5, gaps = T)
saveplot(spearman.matrix.full.alt.old, plot = 'pheatmap', type = 'pdf', width = 12, height = 10.5, gaps = T)
saveplot(spearman.matrix.BZs.old, plot = 'pheatmap', type = 'tiff', width = 3500, height = 3250, gaps = F)
saveplot(spearman.matrix.MLs.old, plot = 'pheatmap', type = 'tiff', width = 4000, height = 3500, gaps = F)
saveplot(spearman.matrix.MLs.alt.old, plot = 'pheatmap', type = 'tiff', width = 3750, height = 3500, gaps = F)
saveplot(spearman.matrix.full.old, plot = 'pheatmap', type = 'tiff', width = 5750, height = 5500, gaps = T)
saveplot(spearman.matrix.full.alt.old, plot = 'pheatmap', type = 'tiff', width = 5500, height = 5000, gaps = T)
R.Version()

#iterative heat maps
pdf(file = "/projects/b1059/projects/Etta/BZML_correlations/figures/heatmaps/iterative_heatmaps/full.old.pdf", width = 12, height = 10.5)
pheatmap(spearman.matrix.full.old, color = viridis(256, option = "D"), fontsize = 15,
         scale = "none", limits = c(-0.2, 1), breaks = seq(-0.2, 1, length.out = 256),
         angle_col = "45", cellwidth = 40, cellheight = 40, display_numbers = T, number_color = 'white',
         fontsize_number = 15,
         border_color = "white",
         cutree_cols = 2, cutree_rows = 2)
dev.off()

##### 