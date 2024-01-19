# BZML Competition Assays (v3) (DMSO, Albendazole) | Calculating competitive fitness 
# Updated 2024-01-12

library(readr)

# Functions 
comp_calc <- function(df, gen1_id = 1){
  
  calculate_gen_ci <- function(gen1_af, gen_af){
    ci <- log10(
      (gen1_af/gen_af - gen1_af)/( 1 - gen1_af)
    )
    return(ci)
  }
  
  #Calculates eq1 from Zhao for a strain, condition, and replicate
  gen1_af <- df %>%
    filter(generation == gen1_id) %>%
    .$allele_freq 
  
  generations <- c(1,3,5,7)
  #create a blank vector to store the ci values
  cis <- c()
  #print(cis)
  #print("starting loop")
  for(generation_id in generations){
    #print(generation_id)
    gen_af <- df %>%
      filter(generation == generation_id) %>%
      .$allele_freq
    #print(gen_af)
    #print(gen_af)    
    ci <- calculate_gen_ci(gen1_af, gen_af)
    #print(ci)
    #print(ci)
    #add the ci to list of cis
    cis <- c(cis, ci)
    
  }
  return(cis)
}

ci_lss <- function(ci_values){
  #takes a list of ci eq2 values and calculates ci 3 values
  # NOT SURE IF THE XS VALUES WILL ALWAYS BE TRUE
  xs <- c(0,2,4,6)
  xs_m <- as.matrix(xs) 
  reg <- lm.fit(xs_m, ci_values)
  val <- as.numeric(coef(reg))
  return(val) 
}

cal_fitness <- function(lm_coef) {
  fitness <- log2(
    1 / 10^lm_coef
  )
  return(as.numeric(fitness))
}


#Import df
v3 <- read.csv("data/competition_assays/v3_20240112.csv", header = T) 


v3 <- v3[v3$Accepted_Droplets >= 10000, ] #remove samples that have less than 10k accepted droplets


#pull out needed columns 
Shaver_v3_pruned_combined <- v3[, c("plate", "well", "strain", "replicate", "generation", "condition", "allele_freq")]


Shaver_v3_pruned_combined <- Shaver_v3_pruned_combined[complete.cases(Shaver_v3_pruned_combined[, c("replicate", "generation")]), ]

#identify redundant samples (samples ran twice due to low droplet counts - remove low droplet samples)
redundant_samples <- Shaver_v3_pruned_combined[duplicated(Shaver_v3_pruned_combined[, c("strain", "condition", "replicate", "generation")]) |
                                 duplicated(Shaver_v3_pruned_combined[, c("strain", "condition", "replicate", "generation")], fromLast = TRUE), ]


# remove redundant samples
Shaver_v3_pruned_combined <- Shaver_v3_pruned_combined[!(Shaver_v3_pruned_combined$strain == "ECA2690" &
                                                           Shaver_v3_pruned_combined$replicate == 10 &
                                                           Shaver_v3_pruned_combined$generation == 1 &
                                                           Shaver_v3_pruned_combined$condition == "albendazole" &
                                                           Shaver_v3_pruned_combined$plate == 3), ]

v3_total <- Shaver_v3_pruned_combined

write.csv(v3_total, "data/competition_assays/v3_total.csv", row.names = FALSE)

##########################################
#    N2 - DMSO and ABZ comp fitness df   #
##########################################
v3_n2_abz <- v3_total  %>% 
  dplyr::filter(strain == "N2") %>% 
  dplyr::mutate(generation = as.numeric(generation),
         replicate = as.numeric(replicate),
         allele_freq = as.numeric(allele_freq)) 

v3_n2_abz<- v3_n2_abz %>%
  dplyr::mutate(allele_freq = ifelse(allele_freq == 0, 0.000000001, allele_freq))

v3_n2_abz <- v3_n2_abz %>%
  dplyr::mutate(allele_freq = ifelse(allele_freq == 1, 0.999999999, allele_freq))

v3_n2_abz <- v3_n2_abz %>%
  dplyr::select(strain, condition, generation, replicate, allele_freq)

##########################################
# ECA2687 - DMSO and ABZ comp fitness df #
##########################################
v3_eca2687_abz <- Shaver_v3_pruned_combined %>% 
  dplyr::filter(strain == "ECA2687") %>% 
  dplyr::mutate(generation = as.numeric(generation),
         replicate = as.numeric(replicate),
         allele_freq = as.numeric(allele_freq)) 

v3_eca2687_abz <- v3_eca2687_abz %>%
  dplyr::mutate(allele_freq = ifelse(allele_freq == 0, 0.000000001, allele_freq))

v3_eca2687_abz <- v3_eca2687_abz %>%
  dplyr::mutate(allele_freq = ifelse(allele_freq == 1, 0.999999999, allele_freq))

##########################################
# ECA2688 - DMSO and ABZ comp fitness df #
##########################################
v3_eca2688_abz <- Shaver_v3_pruned_combined %>% 
  dplyr::filter(strain == "ECA2688") 
  
v3_eca2688_abz <- v3_eca2688_abz %>% 
  dplyr::mutate(generation = as.numeric(generation),
         replicate = as.numeric(replicate),
         allele_freq = as.numeric(allele_freq)) 

v3_eca2688_abz  <- v3_eca2688_abz  %>%
  dplyr::mutate(allele_freq = ifelse(allele_freq == 0, 0.000000001, allele_freq))

v3_eca2688_abz  <- v3_eca2688_abz  %>%
  dplyr::mutate(allele_freq = ifelse(allele_freq == 1, 0.999999999, allele_freq))

##########################################
# ECA2690 - DMSO and ABZ comp fitness df #
##########################################
v3_eca2690_abz <- Shaver_v3_pruned_combined %>% 
  dplyr::filter(strain == "ECA2690") %>% 
  dplyr::mutate(generation = as.numeric(generation),
         replicate = as.numeric(replicate),
         allele_freq = as.numeric(allele_freq)) 

v3_eca2690_abz  <- v3_eca2690_abz  %>%
  dplyr::mutate(allele_freq = ifelse(allele_freq == 0, 0.000000001, allele_freq))

v3_eca2690_abz  <- v3_eca2690_abz  %>%
  dplyr::mutate(allele_freq = ifelse(allele_freq == 1, 0.999999999, allele_freq))

##########################################
# ECA882 - DMSO and ABZ comp fitness df  #
##########################################
v3_eca882_abz <- Shaver_v3_pruned_combined %>% 
  dplyr::filter(strain == "ECA882") %>% 
  dplyr::mutate(generation = as.numeric(generation),
         replicate = as.numeric(replicate),
         allele_freq = as.numeric(allele_freq)) 

v3_eca882_abz <- v3_eca882_abz %>%
  dplyr::mutate(allele_freq = ifelse(allele_freq == 0, 0.000000001, allele_freq))

v3_eca882_abz <- v3_eca882_abz %>%
  dplyr::mutate(allele_freq = ifelse(allele_freq == 1, 0.999999999, allele_freq))

##########################################
# ECA2689 - DMSO and ABZ comp fitness df #
##########################################
v3_eca2689_abz <- Shaver_v3_pruned_combined %>% 
  dplyr::filter(strain == "ECA2689") %>% 
  dplyr::mutate(generation = as.numeric(generation),
         replicate = as.numeric(replicate),
         allele_freq = as.numeric(allele_freq)) 

v3_eca2689_abz <- v3_eca2689_abz %>%
  dplyr::mutate(allele_freq = ifelse(allele_freq == 0, 0.000000001, allele_freq))

v3_eca2689_abz <- v3_eca2689_abz %>%
  dplyr::mutate(allele_freq = ifelse(allele_freq == 1, 0.999999999, allele_freq))

##########################################
# PHX2481 - DMSO and ABZ comp fitness df #
##########################################
v3_phx2481_abz <- Shaver_v3_pruned_combined %>% 
  dplyr::filter(strain == "PHX2481") %>% 
  dplyr::mutate(generation = as.numeric(generation),
         replicate = as.numeric(replicate),
         allele_freq = as.numeric(allele_freq)) 

v3_phx2481_abz <- v3_phx2481_abz %>%
  dplyr::mutate(allele_freq = ifelse(allele_freq == 0, 0.000000001, allele_freq))

v3_phx2481_abz <- v3_phx2481_abz %>%
  dplyr::mutate(allele_freq = ifelse(allele_freq == 1, 0.999999999, allele_freq))

##########################################
# PHX2502 - DMSO and ABZ comp fitness df #
##########################################
v3_phx2502_abz <- Shaver_v3_pruned_combined %>% 
  dplyr::filter(strain == "PHX2502") %>% 
  dplyr::mutate(generation = as.numeric(generation),
         replicate = as.numeric(replicate),
         allele_freq = as.numeric(allele_freq)) 

v3_phx2502_abz  <- v3_phx2502_abz  %>%
  dplyr::mutate(allele_freq = ifelse(allele_freq == 0, 0.000000001, allele_freq))

v3_phx2502_abz  <- v3_phx2502_abz  %>%
  dplyr::mutate(allele_freq = ifelse(allele_freq == 1, 0.999999999, allele_freq))

##########################################
# PHX2647 - DMSO and ABZ comp fitness df #
##########################################
v3_phx2647_abz <- Shaver_v3_pruned_combined %>% 
  dplyr::filter(strain == "PHX2647") %>% 
  dplyr::mutate(generation = as.numeric(generation),
         replicate = as.numeric(replicate),
         allele_freq = as.numeric(allele_freq)) 

v3_phx2647_abz  <- v3_phx2647_abz  %>%
  dplyr::mutate(allele_freq = ifelse(allele_freq == 0, 0.000000001, allele_freq))

v3_phx2647_abz  <- v3_phx2647_abz  %>%
  dplyr::mutate(allele_freq = ifelse(allele_freq == 1, 0.999999999, allele_freq))

#####################################################
#                                                   #
#        Competitive Fitness  - N2 DMSO (v3)        #
#                                                   #
######################################################
# pull out one condition for one strain - N2 DMSO
v3_df_DMSO_N2_rep1 <- v3_n2_abz  %>% # Missing a gen
  dplyr::filter(strain == "N2")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 1)

# rounds df to 7 decimal places
# v3_df_DMSO_N2_rep1$allele_freq <- round(v3_df_DMSO_N2_rep1$allele_freq, 7)

v3_df_DMSO_N2_rep2 <- v3_n2_abz %>%
  dplyr::filter(strain == "N2")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 2)

v3_df_DMSO_N2_rep3 <- v3_n2_abz %>%
  dplyr::filter(strain == "N2")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 3)

v3_df_DMSO_N2_rep4 <- v3_n2_abz %>%
  dplyr::filter(strain == "N2")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 4)

v3_df_DMSO_N2_rep5 <- v3_n2_abz %>%
  dplyr::filter(strain == "N2")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 5)

v3_df_DMSO_N2_rep6 <- v3_n2_abz %>%
  dplyr::filter(strain == "N2")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 6)

v3_df_DMSO_N2_rep7 <- v3_n2_abz %>%
  dplyr::filter(strain == "N2")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 7)

v3_df_DMSO_N2_rep8 <- v3_n2_abz %>%
  dplyr::filter(strain == "N2")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 8)

v3_df_DMSO_N2_rep9 <- v3_n2_abz %>%
  dplyr::filter(strain == "N2")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 9)

v3_df_DMSO_N2_rep10 <- v3_n2_abz %>%
  dplyr::filter(strain == "N2")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 10)


#calculate comp fitness 
# rep 1 - DMSO N2
v3_grouped_DMSO_N2_rep1 <- v3_df_DMSO_N2_rep1 %>% 
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

v3_grouped_DMSO_N2_rep1_select <- v3_grouped_DMSO_N2_rep1 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 2 - DMSO N2
v3_grouped_DMSO_N2_rep2 <- v3_df_DMSO_N2_rep2 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_N2_rep2_select <- v3_grouped_DMSO_N2_rep2 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 3 - DMSO N2
grouped_DMSO_N2_rep3 <- v3_df_DMSO_N2_rep3 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_N2_rep3_select <- grouped_DMSO_N2_rep3 %>%
  dplyr::select(strain, condition, replicate, fitness)

#rep 4 - DMSO N2 
grouped_DMSO_N2_rep4 <- v3_df_DMSO_N2_rep4 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_N2_rep4_select <- grouped_DMSO_N2_rep4 %>%
  dplyr::select(strain, condition, replicate, fitness)

#rep 5 - DMSO N2 
grouped_DMSO_N2_rep5 <- v3_df_DMSO_N2_rep5 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_N2_rep5_select <- grouped_DMSO_N2_rep5 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 6 - DMSO N2 
grouped_DMSO_N2_rep6 <- v3_df_DMSO_N2_rep6 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_N2_rep6_select <- grouped_DMSO_N2_rep6 %>%
  select(strain, condition, replicate, fitness)

#rep 7 - DMSO N2
grouped_DMSO_N2_rep7 <- v3_df_DMSO_N2_rep7 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_N2_rep7_select <- grouped_DMSO_N2_rep7 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 8 - DMSO N2 
grouped_DMSO_N2_rep8 <- v3_df_DMSO_N2_rep8 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_N2_rep8_select <- grouped_DMSO_N2_rep8 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 9 - DMSO - N2 - Skip missing gen 7 
# grouped_DMSO_N2_rep9 <- v3_df_DMSO_N2_rep9 %>%
#   dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
#   dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
#   dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
#   dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))
# 
# grouped_DMSO_N2_rep9_select <- grouped_DMSO_N2_rep9 %>%
#   dplyr::select(strain, condition, replicate, fitness)

# rep 10 - DMSO - N2
grouped_DMSO_N2_rep10 <- v3_df_DMSO_N2_rep10 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_N2_rep10_select <- grouped_DMSO_N2_rep10 %>%
  dplyr::select(strain, condition, replicate, fitness)


fitness_df_n2_DMSO_v3 <- bind_rows(v3_grouped_DMSO_N2_rep1_select,
                                grouped_DMSO_N2_rep2_select,
                                grouped_DMSO_N2_rep3_select,
                                grouped_DMSO_N2_rep4_select,
                                grouped_DMSO_N2_rep5_select,
                                grouped_DMSO_N2_rep6_select,
                                grouped_DMSO_N2_rep7_select,
                                grouped_DMSO_N2_rep8_select,
                                # grouped_DMSO_N2_rep9_select, #missing gen
                                grouped_DMSO_N2_rep10_select)

######################################################
#                                                    #
#      Competitive Fitness  - N2 ABZ  (v3)           #
#                                                    #
######################################################
# pull out one condition for one strain - N2 ABZ
df_ABZ_N2_rep1 <- v3_n2_abz %>%
  dplyr::filter(strain == "N2")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 1)

df_ABZ_N2_rep2 <- v3_n2_abz %>%
  dplyr::filter(strain == "N2")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 2)

df_ABZ_N2_rep3 <- v3_n2_abz %>%
  dplyr::filter(strain == "N2")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 3)

df_ABZ_N2_rep4 <- v3_n2_abz %>%
  dplyr::filter(strain == "N2")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 4)

df_ABZ_N2_rep5 <- v3_n2_abz %>%
  dplyr::filter(strain == "N2")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 5)

df_ABZ_N2_rep6 <- v3_n2_abz %>%
  dplyr::filter(strain == "N2")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 6)

df_ABZ_N2_rep7 <- v3_n2_abz %>%
  dplyr::filter(strain == "N2")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 7)

df_ABZ_N2_rep8 <- v3_n2_abz %>%
  dplyr::filter(strain == "N2")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 8)

df_ABZ_N2_rep9 <- v3_n2_abz %>%
  dplyr::filter(strain == "N2")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 9)

df_ABZ_N2_rep10 <- v3_n2_abz %>%
  dplyr::filter(strain == "N2")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 10)

# rep 1 - N2 ABZ
grouped_ABZ_N2_rep1 <- df_ABZ_N2_rep1 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr:: mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_N2_rep1_select <- grouped_ABZ_N2_rep1 %>% 
  dplyr::select(strain, condition, replicate, fitness)

# rep 2 - N2 ABZ
grouped_ABZ_N2_rep2 <- df_ABZ_N2_rep2 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_N2_rep2_select <- grouped_ABZ_N2_rep2 %>% 
  dplyr::select(strain, condition, replicate, fitness)

# rep 3 - N2 ABZ
grouped_ABZ_N2_rep3 <- df_ABZ_N2_rep3 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_N2_rep3_select <- grouped_ABZ_N2_rep3 %>% 
  dplyr::select(strain, condition, replicate, fitness)

# rep 4 - N2 ABZ - missing gen 3
# grouped_ABZ_N2_rep4 <- df_ABZ_N2_rep4 %>%
#   dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
#   dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
#   dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
#   dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))
# 
# grouped_ABZ_N2_rep4_select <- grouped_ABZ_N2_rep4 %>% 
#   dplyr::select(strain, condition, replicate, fitness)

# rep 5 - N2 ABZ 
grouped_ABZ_N2_rep5 <- df_ABZ_N2_rep5 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_N2_rep5_select <- grouped_ABZ_N2_rep5 %>% 
  dplyr::select(strain, condition, replicate, fitness)

# rep 6 - N2 ABZ 
grouped_ABZ_N2_rep6 <- df_ABZ_N2_rep6 %>%
  group_nest(strain, condition, replicate, keep = TRUE) %>%
  mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_N2_rep6_select <- grouped_ABZ_N2_rep6 %>%
  select(strain, condition, replicate, fitness)

# rep 7 - N2 ABZ 
grouped_ABZ_N2_rep7 <- df_ABZ_N2_rep7 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_N2_rep7_select <- grouped_ABZ_N2_rep7 %>% 
  dplyr::select(strain, condition, replicate, fitness)

# rep 8 - N2 ABZ 
grouped_ABZ_N2_rep8 <- df_ABZ_N2_rep8 %>%
  group_nest(strain, condition, replicate, keep = TRUE) %>%
  mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_N2_rep8_select <- grouped_ABZ_N2_rep8 %>%
  select(strain, condition, replicate, fitness)

# rep 9 - N2 ABZ 
grouped_ABZ_N2_rep9 <- df_ABZ_N2_rep9 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_N2_rep9_select <- grouped_ABZ_N2_rep9 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 10 - N2 ABZ 
grouped_ABZ_N2_rep10 <- df_ABZ_N2_rep10 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_N2_rep10_select <- grouped_ABZ_N2_rep10 %>% 
  dplyr::select(strain, condition, replicate, fitness)


fitness_df_n2_ABZ <- bind_rows(grouped_ABZ_N2_rep1_select, 
                               grouped_ABZ_N2_rep2_select,
                               grouped_ABZ_N2_rep3_select,
                               # grouped_ABZ_N2_rep4_select, # missing gen
                               grouped_ABZ_N2_rep5_select,
                               grouped_ABZ_N2_rep6_select, 
                               grouped_ABZ_N2_rep7_select,
                               grouped_ABZ_N2_rep8_select, 
                               grouped_ABZ_N2_rep9_select,
                               grouped_ABZ_N2_rep10_select)

######################################################
#                                                    #
#          Competitive Fitness  - ECA882 ABZ         #
#                                                    #
######################################################
# pull out one condition for one strain - ECA882 ABZ
df_ABZ_ECA882_rep1 <- v3_eca882_abz %>%
  dplyr::filter(strain == "ECA882")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 1)

df_ABZ_ECA882_rep2 <- v3_eca882_abz  %>%
  dplyr::filter(strain == "ECA882")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 2)

df_ABZ_ECA882_rep3 <- v3_eca882_abz  %>%
  dplyr::filter(strain == "ECA882")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 3)

df_ABZ_ECA882_rep4 <- v3_eca882_abz  %>%
  dplyr::filter(strain == "ECA882")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 4)

df_ABZ_ECA882_rep5 <- v3_eca882_abz  %>%
  dplyr::filter(strain == "ECA882")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 5)

df_ABZ_ECA882_rep6 <- v3_eca882_abz  %>%
  dplyr::filter(strain == "ECA882")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 6)

df_ABZ_ECA882_rep7 <- v3_eca882_abz %>%
  dplyr::filter(strain == "ECA882")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 7)

df_ABZ_ECA882_rep8 <- v3_eca882_abz %>%
  dplyr::filter(strain == "ECA882")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 8)

df_ABZ_ECA882_rep9 <- v3_eca882_abz %>%
  dplyr::filter(strain == "ECA882")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 9)

df_ABZ_ECA882_rep10 <- v3_eca882_abz  %>%
  dplyr::filter(strain == "ECA882")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 10)

# df_ABZ_ECA882_rep1   <- df_ABZ_ECA882_rep1   %>%
#   dplyr::filter(!(allele_freq == 0.068586218))

# rep 1 - ECA882 ABZ
grouped_ABZ_ECA882_rep1 <- df_ABZ_ECA882_rep1 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA882_rep1_select <- grouped_ABZ_ECA882_rep1 %>% 
  dplyr::select(strain, condition, replicate, fitness)

# rep 2 - ECA882 ABZ
grouped_ABZ_ECA882_rep2 <- df_ABZ_ECA882_rep2 %>%
  group_nest(strain, condition, replicate, keep = TRUE) %>%
  mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA882_rep2_select <- grouped_ABZ_ECA882_rep2 %>%
  select(strain, condition, replicate, fitness)

# rep 3 - ECA882 ABZ
grouped_ABZ_ECA882_rep3 <- df_ABZ_ECA882_rep3 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA882_rep3_select <- grouped_ABZ_ECA882_rep3 %>% 
  dplyr::select(strain, condition, replicate, fitness)

# rep 4 - ECA882 ABZ
grouped_ABZ_ECA882_rep4 <- df_ABZ_ECA882_rep4 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA882_rep4_select <- grouped_ABZ_ECA882_rep4 %>% 
  dplyr::select(strain, condition, replicate, fitness)

# rep 5 - ECA882 ABZ
grouped_ABZ_ECA882_rep5 <- df_ABZ_ECA882_rep5 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA882_rep5_select <- grouped_ABZ_ECA882_rep5 %>% 
  dplyr::select(strain, condition, replicate, fitness)

# rep 6 - ECA882 ABZ
grouped_ABZ_ECA882_rep6 <- df_ABZ_ECA882_rep6 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA882_rep6_select <- grouped_ABZ_ECA882_rep6 %>% 
  dplyr::select(strain, condition, replicate, fitness)

# rep 7 - ECA882 ABZ
# df_ABZ_ECA882_rep7 <- df_ABZ_ECA882_rep7 %>%
#   dplyr::filter(!(allele_freq == 0.000000001))
grouped_ABZ_ECA882_rep7 <- df_ABZ_ECA882_rep7 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA882_rep7_select <- grouped_ABZ_ECA882_rep7 %>% 
  dplyr::select(strain, condition, replicate, fitness)

# rep 8 - ECA882 ABZ
grouped_ABZ_ECA882_rep8 <- df_ABZ_ECA882_rep8 %>%
  group_nest(strain, condition, replicate, keep = TRUE) %>%
  mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA882_rep8_select <- grouped_ABZ_ECA882_rep8 %>%
  select(strain, condition, replicate, fitness)

# rep 9 - ECA882 ABZ
#change gen 6 to 5 - typo in dataframe
df_ABZ_ECA882_rep9[df_ABZ_ECA882_rep9$generation == 6, "generation"] <- 5

grouped_ABZ_ECA882_rep9 <- df_ABZ_ECA882_rep9 %>%
  group_nest(strain, condition, replicate, keep = TRUE) %>%
  mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA882_rep9_select <- grouped_ABZ_ECA882_rep9 %>%
  select(strain, condition, replicate, fitness)

# rep 10 - ECA882 ABZ
grouped_ABZ_ECA882_rep10 <- df_ABZ_ECA882_rep10 %>%
  group_nest(strain, condition, replicate, keep = TRUE) %>%
  mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA882_rep10_select <- grouped_ABZ_ECA882_rep10 %>%
  select(strain, condition, replicate, fitness)


fitness_df_ECA882_ABZ <- bind_rows(grouped_ABZ_ECA882_rep1, 
                                   grouped_ABZ_ECA882_rep2,
                                   grouped_ABZ_ECA882_rep3,
                                   grouped_ABZ_ECA882_rep4,
                                   grouped_ABZ_ECA882_rep5,
                                   grouped_ABZ_ECA882_rep6,
                                   grouped_ABZ_ECA882_rep7,
                                   grouped_ABZ_ECA882_rep8,
                                   grouped_ABZ_ECA882_rep9,
                                   grouped_ABZ_ECA882_rep10)


# ######################################################
#                                                      #
#      Competitive Fitness  - ECA882 DMSO (v3)         #
#                                                      #
########################################################
# pull out one condition for one strain - ECA882 DMSO
df_DMSO_ECA882_rep1 <- v3_eca882_abz %>%
  dplyr::filter(strain == "ECA882") %>%
  dplyr::filter(condition == "DMSO") %>%
  dplyr::filter(replicate == 1)

df_DMSO_ECA882_rep2 <- v3_eca882_abz %>%
  filter(strain == "ECA882")%>%
  filter(condition == "DMSO")%>%
  filter(replicate == 2)

df_DMSO_ECA882_rep3 <- v3_eca882_abz %>%
  filter(strain == "ECA882")%>%
  filter(condition == "DMSO")%>%
  filter(replicate == 3)

df_DMSO_ECA882_rep4 <- v3_eca882_abz %>%
  filter(strain == "ECA882")%>%
  filter(condition == "DMSO")%>%
  filter(replicate == 4)

df_DMSO_ECA882_rep5 <- v3_eca882_abz %>%
  filter(strain == "ECA882")%>%
  filter(condition == "DMSO")%>%
  filter(replicate == 5)

df_DMSO_ECA882_rep6 <- v3_eca882_abz %>%
  filter(strain == "ECA882")%>%
  filter(condition == "DMSO")%>%
  filter(replicate == 6)

df_DMSO_ECA882_rep7 <- v3_eca882_abz %>%
  filter(strain == "ECA882")%>%
  filter(condition == "DMSO")%>%
  filter(replicate == 7)

df_DMSO_ECA882_rep8 <- v3_eca882_abz %>%
  filter(strain == "ECA882")%>%
  filter(condition == "DMSO")%>%
  filter(replicate == 8)

df_DMSO_ECA882_rep9 <- v3_eca882_abz %>%
  filter(strain == "ECA882")%>%
  filter(condition == "DMSO")%>%
  filter(replicate == 9)

df_DMSO_ECA882_rep10 <- v3_eca882_abz %>%
  filter(strain == "ECA882")%>%
  filter(condition == "DMSO")%>%
  filter(replicate == 10)

# rep 1 - DMSO - ECA882
grouped_DMSO_ECA882_rep1 <- df_DMSO_ECA882_rep1 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA882_rep1_select <- grouped_DMSO_ECA882_rep1 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 2 - DMSO - ECA882
grouped_DMSO_ECA882_rep2 <- df_DMSO_ECA882_rep2 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA882_rep2_select <- grouped_DMSO_ECA882_rep2 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 3 - DMSO - ECA882 0 missing gen 5 
# grouped_DMSO_ECA882_rep3 <- df_DMSO_ECA882_rep3 %>%
#   dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
#   dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
#   dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
#   dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))
# 
# grouped_DMSO_ECA882_rep3_select <- grouped_DMSO_ECA882_rep3 %>%
#   dplyr::select(strain, condition, replicate, fitness)

# rep 4 - DMSO - ECA882 
grouped_DMSO_ECA882_rep4 <- df_DMSO_ECA882_rep4 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA882_rep4_select <- grouped_DMSO_ECA882_rep4 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 5 - DMSO - ECA882
# df_DMSO_ECA882_rep5 <- df_DMSO_ECA882_rep5 %>%
#   filter(!(near(allele_freq, 0.08113632) & generation == 7))
grouped_DMSO_ECA882_rep5 <- df_DMSO_ECA882_rep5 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA882_rep5_select <- grouped_DMSO_ECA882_rep5 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 6 - DMSO - ECA882
# df_DMSO_ECA882_rep6 <- df_DMSO_ECA882_rep6 %>%
#   dplyr::filter(!(near(allele_freq, 0.374935445) & generation == 7))
grouped_DMSO_ECA882_rep6 <- df_DMSO_ECA882_rep6 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA882_rep6_select <- grouped_DMSO_ECA882_rep6 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 7 - DMSO - ECA882
grouped_DMSO_ECA882_rep7 <- df_DMSO_ECA882_rep7 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA882_rep7_select <- grouped_DMSO_ECA882_rep7 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 8 - DMSO - ECA882
grouped_DMSO_ECA882_rep8 <- df_DMSO_ECA882_rep8 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA882_rep8_select <- grouped_DMSO_ECA882_rep8 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 9 - DMSO - ECA882
# df_DMSO_ECA882_rep9 <- df_DMSO_ECA882_rep9 %>%
#   dplyr::filter(!(near(allele_freq, 0.923380282) & generation == 1))
grouped_DMSO_ECA882_rep9 <- df_DMSO_ECA882_rep9 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA882_rep9_select <- grouped_DMSO_ECA882_rep9 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 10 - DMSO - ECA882
grouped_DMSO_ECA882_rep10 <- df_DMSO_ECA882_rep10 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA882_rep10_select <- grouped_DMSO_ECA882_rep10 %>%
  dplyr::select(strain, condition, replicate, fitness)


fitness_df_ECA882_DMSO <- bind_rows(grouped_DMSO_ECA882_rep1_select,
                                    grouped_DMSO_ECA882_rep2_select,
                                    # grouped_DMSO_ECA882_rep3_select, # gen 5 missing
                                    grouped_DMSO_ECA882_rep4_select,
                                    grouped_DMSO_ECA882_rep5_select,
                                    grouped_DMSO_ECA882_rep6_select, 
                                    grouped_DMSO_ECA882_rep7_select,
                                    grouped_DMSO_ECA882_rep8_select,
                                    grouped_DMSO_ECA882_rep9_select,
                                    grouped_DMSO_ECA882_rep10_select)


######################################################
#                                                    #
#      Competitive Fitness  - PHX2647 DMSO  (V3)     #
#                                                    #
######################################################
# pull out one condition for one strain - PHX2647 DMSO
df_DMSO_PHX2647_rep1 <- v3_phx2647_abz %>%
  dplyr::filter(strain == "PHX2647")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 1)

df_DMSO_PHX2647_rep2 <- v3_phx2647_abz %>%
  dplyr::filter(strain == "PHX2647")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 2)

df_DMSO_PHX2647_rep3 <- v3_phx2647_abz %>%
  dplyr::filter(strain == "PHX2647")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 3)

df_DMSO_PHX2647_rep4 <- v3_phx2647_abz %>%
  dplyr::filter(strain == "PHX2647")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 4)

df_DMSO_PHX2647_rep5 <- v3_phx2647_abz %>%
  dplyr::filter(strain == "PHX2647")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 5)

df_DMSO_PHX2647_rep6 <- v3_phx2647_abz %>%
  dplyr::filter(strain == "PHX2647")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 6)

df_DMSO_PHX2647_rep7 <- v3_phx2647_abz %>%
  dplyr::filter(strain == "PHX2647")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 7)

df_DMSO_PHX2647_rep8 <- v3_phx2647_abz %>%
  dplyr::filter(strain == "PHX2647")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 8)

df_DMSO_PHX2647_rep9 <- v3_phx2647_abz %>%
  dplyr::filter(strain == "PHX2647")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 9)

df_DMSO_PHX2647_rep10 <- v3_phx2647_abz %>%
  dplyr::filter(strain == "PHX2647")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 10)

# rep 1 - DMSO - PHX2647 
# df_DMSO_PHX2647_rep1 <- df_DMSO_PHX2647_rep1 %>%
#   filter(!(near(allele_freq, 0.036006373) & generation == 7))
grouped_DMSO_PHX2647_rep1 <- df_DMSO_PHX2647_rep1 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_PHX2647_rep1_select <- grouped_DMSO_PHX2647_rep1 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 2 - DMSO - PHX2647
grouped_DMSO_PHX2647_rep2 <- df_DMSO_PHX2647_rep2 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_PHX2647_rep2_select <- grouped_DMSO_PHX2647_rep2 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 3 - DMSO - PHX2647
grouped_DMSO_PHX2647_rep3 <- df_DMSO_PHX2647_rep3 %>% 
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_PHX2647_rep3_select <- grouped_DMSO_PHX2647_rep3 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 4 - DMSO - PHX2647
grouped_DMSO_PHX2647_rep4 <- df_DMSO_PHX2647_rep4 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_PHX2647_rep4_select <- grouped_DMSO_PHX2647_rep4 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 5 - DMSO - PHX2647 
grouped_DMSO_PHX2647_rep5 <- df_DMSO_PHX2647_rep5 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_PHX2647_rep5_select <- grouped_DMSO_PHX2647_rep5 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 6 - DMSO - PHX2647
grouped_DMSO_PHX2647_rep6 <- df_DMSO_PHX2647_rep6 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_PHX2647_rep6_select <- grouped_DMSO_PHX2647_rep6 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 7 - DMSO - PHX2647 
grouped_DMSO_PHX2647_rep7 <- df_DMSO_PHX2647_rep7 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_PHX2647_rep7_select <- grouped_DMSO_PHX2647_rep7 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 8 - DMSO - PHX2647
grouped_DMSO_PHX2647_rep8 <- df_DMSO_PHX2647_rep8 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_PHX2647_rep8_select <- grouped_DMSO_PHX2647_rep8 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 9 - DMSO - PHX2647
grouped_DMSO_PHX2647_rep9 <- df_DMSO_PHX2647_rep9 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_PHX2647_rep9_select <- grouped_DMSO_PHX2647_rep9 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 10 - DMSO - PHX2647
grouped_DMSO_PHX2647_rep10 <- df_DMSO_PHX2647_rep10 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_PHX2647_rep10_select <- grouped_DMSO_PHX2647_rep10 %>%
  dplyr::select(strain, condition, replicate, fitness)

fitness_df_PHX2647_DMSO <- bind_rows(grouped_DMSO_PHX2647_rep1_select,
                                     grouped_DMSO_PHX2647_rep2_select,
                                     grouped_DMSO_PHX2647_rep3_select,
                                     grouped_DMSO_PHX2647_rep4_select,
                                     grouped_DMSO_PHX2647_rep5_select,
                                     grouped_DMSO_PHX2647_rep6_select,
                                     grouped_DMSO_PHX2647_rep7_select,
                                     grouped_DMSO_PHX2647_rep8_select,
                                     grouped_DMSO_PHX2647_rep9_select,
                                     grouped_DMSO_PHX2647_rep10_select)

######################################################
#                                                    #
#          Competitive Fitness  - PHX2647 ABZ        #
#                                                    #
######################################################
# pull out one condition for one strain - PHX2647 ABZ
df_ABZ_PHX2647_rep1 <- v3_phx2647_abz %>%
  dplyr::filter(strain == "PHX2647")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 1)

df_ABZ_PHX2647_rep2 <- v3_phx2647_abz %>%
  dplyr::filter(strain == "PHX2647")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 2)

df_ABZ_PHX2647_rep3 <- v3_phx2647_abz %>%
  dplyr::filter(strain == "PHX2647")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 3)

df_ABZ_PHX2647_rep4 <- v3_phx2647_abz %>%
  dplyr::filter(strain == "PHX2647")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 4)

df_ABZ_PHX2647_rep5 <- v3_phx2647_abz %>%
  dplyr::filter(strain == "PHX2647")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 5)

df_ABZ_PHX2647_rep6 <- v3_phx2647_abz %>%
  dplyr::filter(strain == "PHX2647")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 6)

df_ABZ_PHX2647_rep7 <- v3_phx2647_abz %>%
  dplyr::filter(strain == "PHX2647")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 7)

df_ABZ_PHX2647_rep8 <- v3_phx2647_abz %>%
  dplyr::filter(strain == "PHX2647")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 8)

df_ABZ_PHX2647_rep9 <- v3_phx2647_abz %>%
  dplyr::filter(strain == "PHX2647")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 9)

df_ABZ_PHX2647_rep10 <- v3_phx2647_abz%>%
  dplyr::filter(strain == "PHX2647")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 10)

# rep 1 - PHX2647 - ABZ 
grouped_ABZ_PHX2647_rep1 <- df_ABZ_PHX2647_rep1 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_PHX2647_rep1_select <- grouped_ABZ_PHX2647_rep1 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 2 - PHX2647 - ABZ 
grouped_ABZ_PHX2647_rep2 <- df_ABZ_PHX2647_rep2 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_PHX2647_rep2_select <- grouped_ABZ_PHX2647_rep2 %>% 
  dplyr::select(strain, condition, replicate, fitness)

# rep 3 - PHX2647 - ABZ
grouped_ABZ_PHX2647_rep3 <- df_ABZ_PHX2647_rep3 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_PHX2647_rep3_select <- grouped_ABZ_PHX2647_rep3 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 4 - PHX2647 - ABZ 
grouped_ABZ_PHX2647_rep4 <- df_ABZ_PHX2647_rep4 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_PHX2647_rep4_select <- grouped_ABZ_PHX2647_rep4 %>% 
  dplyr::select(strain, condition, replicate, fitness)

# rep 5- PHX2647 - ABZ 
grouped_ABZ_PHX2647_rep5 <- df_ABZ_PHX2647_rep5 %>%
  group_nest(strain, condition, replicate, keep = TRUE) %>%
  mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_PHX2647_rep5_select <- grouped_ABZ_PHX2647_rep5 %>%
  select(strain, condition, replicate, fitness)

# rep 6 - PHX2647 - ABZ
grouped_ABZ_PHX2647_rep6 <- df_ABZ_PHX2647_rep6 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_PHX2647_rep6_select <- grouped_ABZ_PHX2647_rep6 %>% 
  dplyr::select(strain, condition, replicate, fitness)

# rep 7 - PHX2647 - ABZ
grouped_ABZ_PHX2647_rep7 <- df_ABZ_PHX2647_rep7 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_PHX2647_rep7_select <- grouped_ABZ_PHX2647_rep7 %>% 
  dplyr::select(strain, condition, replicate, fitness)

# rep 8 - PHX2647 - ABZ
grouped_ABZ_PHX2647_rep8 <- df_ABZ_PHX2647_rep8 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_PHX2647_rep8_select <- grouped_ABZ_PHX2647_rep8 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 9 - PHX2647 - ABZ
grouped_ABZ_PHX2647_rep9 <- df_ABZ_PHX2647_rep9 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_PHX2647_rep9_select <- grouped_ABZ_PHX2647_rep9 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 10 - PHX2647 - ABZ
grouped_ABZ_PHX2647_rep10 <- df_ABZ_PHX2647_rep10 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_PHX2647_rep10_select <- grouped_ABZ_PHX2647_rep10 %>% 
  dplyr::select(strain, condition, replicate, fitness)

fitness_df_PHX2647_ABZ <- bind_rows(grouped_ABZ_PHX2647_rep1_select,
                                    grouped_ABZ_PHX2647_rep2_select,
                                    grouped_ABZ_PHX2647_rep3_select,
                                    grouped_ABZ_PHX2647_rep4_select,
                                    grouped_ABZ_PHX2647_rep5_select,
                                    grouped_ABZ_PHX2647_rep6_select,
                                    grouped_ABZ_PHX2647_rep7_select,
                                    grouped_ABZ_PHX2647_rep8_select,
                                    grouped_ABZ_PHX2647_rep9_select,
                                    grouped_ABZ_PHX2647_rep10_select)


######################################################
#                                                    #
#     Competitive Fitness  - PHX2502 DMSO  (v3)      #
#                                                    #
######################################################
# pull out one condition for one strain - PHX2502 DMSO
df_DMSO_PHX2502_rep1 <- v3_phx2502_abz %>%
  dplyr::filter(strain == "PHX2502")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 1)

df_DMSO_PHX2502_rep2 <- v3_phx2502_abz  %>%
  dplyr::filter(strain == "PHX2502")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 2)

df_DMSO_PHX2502_rep3 <- v3_phx2502_abz %>%
  filter(strain == "PHX2502")%>%
  filter(condition == "DMSO")%>%
  filter(replicate == 3)

df_DMSO_PHX2502_rep4 <- v3_phx2502_abz %>%
  filter(strain == "PHX2502")%>%
  filter(condition == "DMSO")%>%
  filter(replicate == 4)

df_DMSO_PHX2502_rep5 <- v3_phx2502_abz %>%
  filter(strain == "PHX2502")%>%
  filter(condition == "DMSO")%>%
  filter(replicate == 5)

df_DMSO_PHX2502_rep6 <- v3_phx2502_abz %>%
  filter(strain == "PHX2502")%>%
  filter(condition == "DMSO")%>%
  filter(replicate == 6)

df_DMSO_PHX2502_rep7 <- v3_phx2502_abz %>%
  filter(strain == "PHX2502")%>%
  filter(condition == "DMSO")%>%
  filter(replicate == 7)

df_DMSO_PHX2502_rep8 <- v3_phx2502_abz %>%
  filter(strain == "PHX2502")%>%
  filter(condition == "DMSO")%>%
  filter(replicate == 8)

df_DMSO_PHX2502_rep9 <- v3_phx2502_abz %>%
  filter(strain == "PHX2502")%>%
  filter(condition == "DMSO")%>%
  filter(replicate == 9)

df_DMSO_PHX2502_rep10 <- v3_phx2502_abz %>%
  filter(strain == "PHX2502")%>%
  filter(condition == "DMSO")%>%
  filter(replicate == 10)

# rep 1 - DMSO PHX2502 - missing gen 1 
# grouped_DMSO_PHX2502_rep1 <- df_DMSO_PHX2502_rep1 %>%
#   dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
#   dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
#   dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
#   dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))
# 
# grouped_DMSO_PHX2502_rep1_select <- grouped_DMSO_PHX2502_rep1 %>%
#   dplyr::select(strain, condition, replicate, fitness)

# rep 2 - DMSO PHX2502
grouped_DMSO_PHX2502_rep2 <- df_DMSO_PHX2502_rep2 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_PHX2502_rep2_select <- grouped_DMSO_PHX2502_rep2 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 3 - PHX2502 - DMSO 
grouped_DMSO_PHX2502_rep3 <- df_DMSO_PHX2502_rep3 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_PHX2502_rep3_select <- grouped_DMSO_PHX2502_rep3 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 4 - PHX2502 - DMSO
grouped_DMSO_PHX2502_rep4 <- df_DMSO_PHX2502_rep4 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_PHX2502_rep4_select <- grouped_DMSO_PHX2502_rep4 %>%
  dplyr::select(strain, condition, replicate, fitness)

#rep 5 - PHX2502 - DMSO
grouped_DMSO_PHX2502_rep5 <- df_DMSO_PHX2502_rep5 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_PHX2502_rep5_select <- grouped_DMSO_PHX2502_rep5 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 6 - DMSO - PHX2502
grouped_DMSO_PHX2502_rep6 <- df_DMSO_PHX2502_rep6 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_PHX2502_rep6_select <- grouped_DMSO_PHX2502_rep6 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 7 - DMSO - PHX2502
grouped_DMSO_PHX2502_rep7 <- df_DMSO_PHX2502_rep7 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_PHX2502_rep7_select <- grouped_DMSO_PHX2502_rep7 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 8 - PHX2502 - DMSO 
grouped_DMSO_PHX2502_rep8 <- df_DMSO_PHX2502_rep8 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_PHX2502_rep8_select <- grouped_DMSO_PHX2502_rep8 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 9 - PHX2502 - DMSO 
grouped_DMSO_PHX2502_rep9 <- df_DMSO_PHX2502_rep9 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_PHX2502_rep9_select <- grouped_DMSO_PHX2502_rep9 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 10 - PHX2502 - DMSO 
grouped_DMSO_PHX2502_rep10 <- df_DMSO_PHX2502_rep10 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_PHX2502_rep10_select <- grouped_DMSO_PHX2502_rep10 %>%
  dplyr::select(strain, condition, replicate, fitness)


fitness_df_PHX2502_DMSO <- bind_rows(#grouped_DMSO_PHX2502_rep1_select,
                                     grouped_DMSO_PHX2502_rep2_select,
                                     grouped_DMSO_PHX2502_rep3_select,
                                     grouped_DMSO_PHX2502_rep4_select,
                                     grouped_DMSO_PHX2502_rep5_select,
                                     grouped_DMSO_PHX2502_rep6_select,
                                     grouped_DMSO_PHX2502_rep7_select,
                                     grouped_DMSO_PHX2502_rep8_select,
                                     grouped_DMSO_PHX2502_rep9_select,
                                     grouped_DMSO_PHX2502_rep10_select)

######################################################
#                                                    #
#          Competitive Fitness  - PHX2502 ABZ        #
#                                                    #
######################################################
# pull out one condition for one strain - PHX2502 ABZ
df_ABZ_PHX2502_rep1 <- v3_phx2502_abz %>%
  dplyr::filter(strain == "PHX2502")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 1)

df_ABZ_PHX2502_rep2 <- v3_phx2502_abz  %>%
  dplyr::filter(strain == "PHX2502")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 2)

df_ABZ_PHX2502_rep3 <- v3_phx2502_abz  %>%
  dplyr::filter(strain == "PHX2502")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 3)

df_ABZ_PHX2502_rep4 <- v3_phx2502_abz  %>%
  dplyr::filter(strain == "PHX2502")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 4)

df_ABZ_PHX2502_rep5 <- v3_phx2502_abz  %>%
  dplyr::filter(strain == "PHX2502")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 5)

df_ABZ_PHX2502_rep6 <- v3_phx2502_abz  %>%
  dplyr::filter(strain == "PHX2502")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 6)

df_ABZ_PHX2502_rep7 <- v3_phx2502_abz  %>%
  dplyr::filter(strain == "PHX2502")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 7)

df_ABZ_PHX2502_rep8 <- v3_phx2502_abz  %>%
  dplyr::filter(strain == "PHX2502")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 8)

df_ABZ_PHX2502_rep9 <- v3_phx2502_abz  %>%
  dplyr::filter(strain == "PHX2502")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 9)

df_ABZ_PHX2502_rep10 <- v3_phx2502_abz  %>%
  dplyr::filter(strain == "PHX2502")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 10)


# df_ABZ_PHX2502_rep1  <- df_ABZ_PHX2502_rep1  %>%
#   filter(!(allele_freq == 0.000142898))

# rep 1 - PHX2502 - ABZ
grouped_ABZ_PHX2502_rep1 <- df_ABZ_PHX2502_rep1 %>% 
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_PHX2502_rep1_select <- grouped_ABZ_PHX2502_rep1 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 2 - PHX2502 - ABZ 
grouped_ABZ_PHX2502_rep2 <- df_ABZ_PHX2502_rep2 %>% #Issue 0 for gen 7
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_PHX2502_rep2_select <- grouped_ABZ_PHX2502_rep2 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 3 - PHX2502 - ABZ 
grouped_ABZ_PHX2502_rep3 <- df_ABZ_PHX2502_rep3 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_PHX2502_rep3_select <- grouped_ABZ_PHX2502_rep3 %>% 
  dplyr::select(strain, condition, replicate, fitness)

# rep 4 - PHX2502 - ABZ 
grouped_ABZ_PHX2502_rep4 <- df_ABZ_PHX2502_rep4 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_PHX2502_rep4_select <- grouped_ABZ_PHX2502_rep4 %>% 
  dplyr::select(strain, condition, replicate, fitness)

# rep 5 - PHX2502 - ABZ
grouped_ABZ_PHX2502_rep5 <- df_ABZ_PHX2502_rep5 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_PHX2502_rep5_select <- grouped_ABZ_PHX2502_rep5 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 6 - PHX2502 - ABZ 
# df_ABZ_PHX2502_rep6  <- df_ABZ_PHX2502_rep6  %>%
#   dplyr::filter(!(allele_freq == 0.255539568))

grouped_ABZ_PHX2502_rep6 <- df_ABZ_PHX2502_rep6 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_PHX2502_rep6_select <- grouped_ABZ_PHX2502_rep6 %>% 
  dplyr::select(strain, condition, replicate, fitness)

# rep 7 - PHX2502 - ABZ 
grouped_ABZ_PHX2502_rep7 <- df_ABZ_PHX2502_rep7 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_PHX2502_rep7_select <- grouped_ABZ_PHX2502_rep7 %>% 
  dplyr::select(strain, condition, replicate, fitness)

# rep 8 - PHX2502 - ABZ
grouped_ABZ_PHX2502_rep8 <- df_ABZ_PHX2502_rep8 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_PHX2502_rep8_select <- grouped_ABZ_PHX2502_rep8 %>% 
  dplyr::select(strain, condition, replicate, fitness)

# rep 9 - PHX2502 - ABZ 
grouped_ABZ_PHX2502_rep9 <- df_ABZ_PHX2502_rep9 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_PHX2502_rep9_select <- grouped_ABZ_PHX2502_rep9 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 10 - PHX2502 - ABZ
grouped_ABZ_PHX2502_rep10 <- df_ABZ_PHX2502_rep10 %>% 
  group_nest(strain, condition, replicate, keep = TRUE) %>%
  mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_PHX2502_rep10_select <- grouped_ABZ_PHX2502_rep10 %>%
  select(strain, condition, replicate, fitness)

fitness_df_PHX2502_ABZ <- bind_rows(grouped_ABZ_PHX2502_rep1_select,
                                    grouped_ABZ_PHX2502_rep2_select,
                                    grouped_ABZ_PHX2502_rep3_select,
                                    grouped_ABZ_PHX2502_rep4_select,
                                    grouped_ABZ_PHX2502_rep5_select,
                                    grouped_ABZ_PHX2502_rep6_select,
                                    grouped_ABZ_PHX2502_rep7_select,
                                    grouped_ABZ_PHX2502_rep8_select,
                                    grouped_ABZ_PHX2502_rep9_select,
                                    grouped_ABZ_PHX2502_rep10_select)



######################################################
#                                                    #
#          Competitive Fitness  - PHX2481 DMSO       #
#                                                    #
######################################################
# pull out one condition for one strain - PHX2481 DMSO
df_DMSO_PHX2481_rep1 <- v3_phx2481_abz %>%
  dplyr::filter(strain == "PHX2481")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 1)

df_DMSO_PHX2481_rep2 <- v3_phx2481_abz  %>%
  dplyr::filter(strain == "PHX2481")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 2)

df_DMSO_PHX2481_rep3 <- v3_phx2481_abz  %>%
  dplyr::filter(strain == "PHX2481")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 3)

df_DMSO_PHX2481_rep4 <- v3_phx2481_abz  %>%
  dplyr::filter(strain == "PHX2481")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 4)

df_DMSO_PHX2481_rep5 <- v3_phx2481_abz  %>%
  dplyr::filter(strain == "PHX2481")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 5)

df_DMSO_PHX2481_rep6 <- v3_phx2481_abz  %>%
  dplyr::filter(strain == "PHX2481")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 6)

df_DMSO_PHX2481_rep7 <- v3_phx2481_abz  %>%
  dplyr::filter(strain == "PHX2481")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 7)

df_DMSO_PHX2481_rep8 <- v3_phx2481_abz  %>%
  dplyr::filter(strain == "PHX2481")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 8)

df_DMSO_PHX2481_rep9 <- v3_phx2481_abz  %>%
  dplyr::filter(strain == "PHX2481")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 9)

df_DMSO_PHX2481_rep10 <- v3_phx2481_abz %>%
  dplyr::filter(strain == "PHX2481")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 10)

# rep 1 - PHX2481 - DMSO 
grouped_DMSO_PHX2481_rep1 <- df_DMSO_PHX2481_rep1 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_PHX2481_rep1_select <- grouped_DMSO_PHX2481_rep1 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 2 - PHX2481 - DMSO - Missing a gen 
grouped_DMSO_PHX2481_rep2 <- df_DMSO_PHX2481_rep2 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_PHX2481_rep2_select <- grouped_DMSO_PHX2481_rep2 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 3 - PHX2481 - DMSO 
grouped_DMSO_PHX2481_rep3 <- df_DMSO_PHX2481_rep3 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_PHX2481_rep3_select <- grouped_DMSO_PHX2481_rep3 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 4 - PHX2481 - DMSO 
# df_DMSO_PHX2481_rep4 <- df_DMSO_PHX2481_rep4 %>%
#   dplyr::filter(!(allele_freq == 0.99822127 & generation == 7))
grouped_DMSO_PHX2481_rep4 <- df_DMSO_PHX2481_rep4 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_PHX2481_rep4_select <- grouped_DMSO_PHX2481_rep4 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 5 - PHX2481 - DMSO 
grouped_DMSO_PHX2481_rep5 <- df_DMSO_PHX2481_rep5 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_PHX2481_rep5_select <- grouped_DMSO_PHX2481_rep5 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 6 - PHX2481 - DMSO - Missing gen 
grouped_DMSO_PHX2481_rep6 <- df_DMSO_PHX2481_rep6 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_PHX2481_rep6_select <- grouped_DMSO_PHX2481_rep6 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 7 - PHX2481 - DMSO 
grouped_DMSO_PHX2481_rep7 <- df_DMSO_PHX2481_rep7 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_PHX2481_rep7_select <- grouped_DMSO_PHX2481_rep7 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 8 - PHX2481 - DMSO 
grouped_DMSO_PHX2481_rep8 <- df_DMSO_PHX2481_rep8 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_PHX2481_rep8_select <- grouped_DMSO_PHX2481_rep8 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 9 - PHX2481 - DMSO 
# df_DMSO_PHX2481_rep9   <- df_DMSO_PHX2481_rep9   %>%
#   filter(!(allele_freq == 0.000132459))
grouped_DMSO_PHX2481_rep9 <- df_DMSO_PHX2481_rep9 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_PHX2481_rep9_select <- grouped_DMSO_PHX2481_rep9 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 10 - PHX2481 - DMSO 
grouped_DMSO_PHX2481_rep10 <- df_DMSO_PHX2481_rep10 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_PHX2481_rep10_select <- grouped_DMSO_PHX2481_rep10 %>%
  dplyr::select(strain, condition, replicate, fitness)


fitness_df_PHX2481_DMSO <- bind_rows(grouped_DMSO_PHX2481_rep1_select,
                                     grouped_DMSO_PHX2481_rep2_select,
                                     grouped_DMSO_PHX2481_rep3_select,
                                     grouped_DMSO_PHX2481_rep4_select,
                                     grouped_DMSO_PHX2481_rep5_select,
                                     grouped_DMSO_PHX2481_rep6_select,
                                     grouped_DMSO_PHX2481_rep7_select,
                                     grouped_DMSO_PHX2481_rep8_select,
                                     grouped_DMSO_PHX2481_rep9_select,
                                     grouped_DMSO_PHX2481_rep10_select)
View(fitness_df_PHX2481_DMSO)

######################################################
#                                                    #
#          Competitive Fitness  - PHX2481 ABZ        #
#                                                    #
######################################################
# pull out one condition for one strain - PHX2481 ABZ
df_ABZ_PHX2481_rep1 <- v3_phx2481_abz %>%
  dplyr::filter(strain == "PHX2481")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 1)

df_ABZ_PHX2481_rep2 <- v3_phx2481_abz  %>%
  dplyr::filter(strain == "PHX2481")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 2)

df_ABZ_PHX2481_rep3 <- v3_phx2481_abz  %>%
  dplyr::filter(strain == "PHX2481")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 3)

df_ABZ_PHX2481_rep4 <- v3_phx2481_abz  %>%
  dplyr::filter(strain == "PHX2481")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 4)

df_ABZ_PHX2481_rep5 <- v3_phx2481_abz %>%
  dplyr::filter(strain == "PHX2481")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 5)

df_ABZ_PHX2481_rep6 <- v3_phx2481_abz %>%
  dplyr::filter(strain == "PHX2481")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 6)

df_ABZ_PHX2481_rep7 <- v3_phx2481_abz  %>%
  dplyr::filter(strain == "PHX2481")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 7)

df_ABZ_PHX2481_rep8 <- v3_phx2481_abz  %>%
  dplyr::filter(strain == "PHX2481")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 8)

df_ABZ_PHX2481_rep9 <- v3_phx2481_abz  %>%
  dplyr::filter(strain == "PHX2481")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 9)

df_ABZ_PHX2481_rep10 <- v3_phx2481_abz  %>%
  dplyr::filter(strain == "PHX2481")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 10)

# rep 1 - PHX2481 - ABZ 
grouped_ABZ_PHX2481_rep1 <- df_ABZ_PHX2481_rep1 %>% 
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_PHX2481_rep1_select <- grouped_ABZ_PHX2481_rep1 %>% 
  dplyr::select(strain, condition, replicate, fitness)

# rep 2 - PHX2481 - ABZ 
grouped_ABZ_PHX2481_rep2 <- df_ABZ_PHX2481_rep2 %>% 
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_PHX2481_rep2_select <- grouped_ABZ_PHX2481_rep2 %>% 
  dplyr::select(strain, condition, replicate, fitness)

# rep 3 - ABZ - PHX2481
grouped_ABZ_PHX2481_rep3 <- df_ABZ_PHX2481_rep3 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_PHX2481_rep3_select <- grouped_ABZ_PHX2481_rep3 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 4 - ABZ - PHX2481 
# df_ABZ_PHX2481_rep4   <- df_ABZ_PHX2481_rep4   %>%
#   dplyr::filter(!(allele_freq == 0.007911553))
grouped_ABZ_PHX2481_rep4 <- df_ABZ_PHX2481_rep4 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_PHX2481_rep4_select <- grouped_ABZ_PHX2481_rep4 %>% 
  dplyr::select(strain, condition, replicate, fitness)

# rep 5 - ABZ - PHX2481 
grouped_ABZ_PHX2481_rep5 <- df_ABZ_PHX2481_rep5 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_PHX2481_rep5_select <- grouped_ABZ_PHX2481_rep5 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 6 - ABZ - PHX2481 
grouped_ABZ_PHX2481_rep6 <- df_ABZ_PHX2481_rep6 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_PHX2481_rep6_select <- grouped_ABZ_PHX2481_rep6 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 7 - ABZ - PHX2481
# df_ABZ_PHX2481_rep7   <- df_ABZ_PHX2481_rep7   %>%
#   dplyr::filter(!(allele_freq == 0.079226450))
grouped_ABZ_PHX2481_rep7 <- df_ABZ_PHX2481_rep7 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_PHX2481_rep7_select <- grouped_ABZ_PHX2481_rep7 %>% 
  dplyr::select(strain, condition, replicate, fitness)

# rep 8 - ABZ - PHX2481 
grouped_ABZ_PHX2481_rep8 <- df_ABZ_PHX2481_rep8 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_PHX2481_rep8_select <- grouped_ABZ_PHX2481_rep8 %>% 
  dplyr::select(strain, condition, replicate, fitness)

# rep 9 - ABZ - PHX2481 
# df_ABZ_PHX2481_rep9   <- df_ABZ_PHX2481_rep9   %>%
#   dplyr::filter(!(allele_freq == 0.0077868140))
# 
# df_ABZ_PHX2481_rep9   <- df_ABZ_PHX2481_rep9   %>%
#   dplyr::filter(!(allele_freq == 0.1459047620))
grouped_ABZ_PHX2481_rep9 <- df_ABZ_PHX2481_rep9 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_PHX2481_rep9_select <- grouped_ABZ_PHX2481_rep9 %>%
  dplyr::select(strain, condition, replicate, fitness)


# df_ABZ_PHX2481_rep10   <- df_ABZ_PHX2481_rep10  %>% skip gen 7 missing 
#   dplyr::filter(!(allele_freq == 0.1459047620))
# grouped_ABZ_PHX2481_rep10 <- df_ABZ_PHX2481_rep10 %>% 
#   dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
#   dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
#   dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
#   dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))
# 
# grouped_ABZ_PHX2481_rep10_select <- grouped_ABZ_PHX2481_rep10 %>%
#   dplyr::select(strain, condition, replicate, fitness)

fitness_df_PHX2481_ABZ <- bind_rows(grouped_ABZ_PHX2481_rep1_select,
                                    grouped_ABZ_PHX2481_rep2_select,
                                    grouped_ABZ_PHX2481_rep3_select,
                                    grouped_ABZ_PHX2481_rep4_select,
                                    grouped_ABZ_PHX2481_rep5_select,
                                    grouped_ABZ_PHX2481_rep6_select,
                                    grouped_ABZ_PHX2481_rep7_select,
                                    grouped_ABZ_PHX2481_rep8_select,
                                    grouped_ABZ_PHX2481_rep9_select)
                                    # grouped_ABZ_PHX2481_rep10_select)


######################################################
#                                                    #
#          Competitive Fitness  - ECA2687 DMSO       #
#                                                    #
######################################################
# pull out one condition for one strain - ECA2687 DMSO
df_DMSO_ECA2687_rep1 <- v3_eca2687_abz %>%
  dplyr::filter(strain == "ECA2687")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 1)

df_DMSO_ECA2687_rep2 <- v3_eca2687_abz %>%
  dplyr::filter(strain == "ECA2687")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 2)

df_DMSO_ECA2687_rep3 <- v3_eca2687_abz %>%
  dplyr::filter(strain == "ECA2687")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 3)

df_DMSO_ECA2687_rep4 <- v3_eca2687_abz %>%
  dplyr::filter(strain == "ECA2687")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 4)

df_DMSO_ECA2687_rep5 <- v3_eca2687_abz  %>%
  dplyr::filter(strain == "ECA2687")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 5)

df_DMSO_ECA2687_rep6 <- v3_eca2687_abz %>%
  dplyr::filter(strain == "ECA2687")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 6)

df_DMSO_ECA2687_rep7 <- v3_eca2687_abz %>%
  dplyr::filter(strain == "ECA2687")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 7)

df_DMSO_ECA2687_rep8 <- v3_eca2687_abz  %>%
  dplyr::filter(strain == "ECA2687")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 8)

df_DMSO_ECA2687_rep9 <- v3_eca2687_abz %>%
  dplyr::filter(strain == "ECA2687")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 9)

df_DMSO_ECA2687_rep10 <- v3_eca2687_abz %>%
  dplyr::filter(strain == "ECA2687")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 10)

# rep 1 - DMSO - ECA2687
grouped_DMSO_ECA2687_rep1 <- df_DMSO_ECA2687_rep1 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2687_rep1_select <- grouped_DMSO_ECA2687_rep1 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 2 - DMSO - ECA2687
grouped_DMSO_ECA2687_rep2 <- df_DMSO_ECA2687_rep2 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2687_rep2_select <- grouped_DMSO_ECA2687_rep2 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 3 - DMSO - ECA2687
grouped_DMSO_ECA2687_rep3 <- df_DMSO_ECA2687_rep3 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2687_rep3_select <- grouped_DMSO_ECA2687_rep3 %>%
  dplyr::select(strain, condition, replicate, fitness)


# rep 4 - DMSO - ECA2687
# View(df_DMSO_ECA2687_rep4)
# df_DMSO_ECA2687_rep4   <- df_DMSO_ECA2687_rep4   %>%
#   dplyr::filter(!(allele_freq == 0.226828546))
grouped_DMSO_ECA2687_rep4 <- df_DMSO_ECA2687_rep4 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2687_rep4_select <- grouped_DMSO_ECA2687_rep4 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 5 - ECA2687 - DMSO
grouped_DMSO_ECA2687_rep5 <- df_DMSO_ECA2687_rep5 %>% 
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2687_rep5_select <- grouped_DMSO_ECA2687_rep5 %>%
  select(strain, condition, replicate, fitness)

# rep 6 - ECA2687 - DMSO
grouped_DMSO_ECA2687_rep6 <- df_DMSO_ECA2687_rep6 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2687_rep6_select <- grouped_DMSO_ECA2687_rep6 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 7 - ECA2687 - DMSO 
grouped_DMSO_ECA2687_rep7 <- df_DMSO_ECA2687_rep7 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2687_rep7_select <- grouped_DMSO_ECA2687_rep7 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 8 - ECA2687 - DMSO
grouped_DMSO_ECA2687_rep8 <- df_DMSO_ECA2687_rep8 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2687_rep8_select <- grouped_DMSO_ECA2687_rep8 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep - 9 - ECA2687 - DMSO 
grouped_DMSO_ECA2687_rep9 <- df_DMSO_ECA2687_rep9 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2687_rep9_select <- grouped_DMSO_ECA2687_rep9 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 10 - ECA2687 - DMSO 
# df_DMSO_ECA2687_rep10   <- df_DMSO_ECA2687_rep10   %>%
#   dplyr::filter(!(allele_freq == 0.000000001 & generation == 5))
grouped_DMSO_ECA2687_rep10 <- df_DMSO_ECA2687_rep10 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2687_rep10_select <- grouped_DMSO_ECA2687_rep10 %>%
  dplyr::select(strain, condition, replicate, fitness)

fitness_df_ECA2687_DMSO <- bind_rows(grouped_DMSO_ECA2687_rep1_select,
                                     grouped_DMSO_ECA2687_rep2_select,
                                     grouped_DMSO_ECA2687_rep3_select,
                                     grouped_DMSO_ECA2687_rep4_select,
                                     grouped_DMSO_ECA2687_rep5_select,
                                     grouped_DMSO_ECA2687_rep6_select,
                                     grouped_DMSO_ECA2687_rep7_select,
                                     grouped_DMSO_ECA2687_rep8_select,
                                     grouped_DMSO_ECA2687_rep9_select,
                                     grouped_DMSO_ECA2687_rep10_select)


######################################################
#                                                    #
#          Competitive Fitness  - ECA2687 ABZ        #
#                                                    #
######################################################
# pull out one condition for one strain - ECA2687 ABZ
df_ABZ_ECA2687_rep1 <- v3_eca2687_abz %>%
  dplyr::filter(strain == "ECA2687")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 1)

df_ABZ_ECA2687_rep2 <- v3_eca2687_abz %>%
  dplyr::filter(strain == "ECA2687")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 2)

df_ABZ_ECA2687_rep3 <- v3_eca2687_abz  %>%
  dplyr::filter(strain == "ECA2687")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 3)

df_ABZ_ECA2687_rep4 <- v3_eca2687_abz %>%
  dplyr::filter(strain == "ECA2687")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 4)

df_ABZ_ECA2687_rep5 <- v3_eca2687_abz %>%
  dplyr::filter(strain == "ECA2687")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 5)

df_ABZ_ECA2687_rep6 <- v3_eca2687_abz %>%
  dplyr::filter(strain == "ECA2687")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 6)

df_ABZ_ECA2687_rep7 <- v3_eca2687_abz  %>%
  dplyr::filter(strain == "ECA2687")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 7)

df_ABZ_ECA2687_rep8 <- v3_eca2687_abz  %>%
  dplyr::filter(strain == "ECA2687")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 8)

df_ABZ_ECA2687_rep9 <- v3_eca2687_abz %>%
  dplyr::filter(strain == "ECA2687")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 9)

df_ABZ_ECA2687_rep10 <- v3_eca2687_abz %>%
  dplyr::filter(strain == "ECA2687")%>%
  dplyr::filter(condition == "albendazole")%>%
  dplyr::filter(replicate == 10)

# rep 1 - ECA2687 ABZ
# df_ABZ_ECA2687_rep1   <- df_ABZ_ECA2687_rep1   %>%
#   filter(!(allele_freq == 0.701568627))
grouped_ABZ_ECA2687_rep1 <- df_ABZ_ECA2687_rep1 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2687_rep1_select <- grouped_ABZ_ECA2687_rep1 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 2 - ECA2687 ABZ
grouped_ABZ_ECA2687_rep2 <- df_ABZ_ECA2687_rep2 %>% 
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2687_rep2_select <- grouped_ABZ_ECA2687_rep2 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 3 - ECA2687 ABZ
grouped_ABZ_ECA2687_rep3 <- df_ABZ_ECA2687_rep3 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2687_rep3_select <- grouped_ABZ_ECA2687_rep3 %>% 
  dplyr::select(strain, condition, replicate, fitness)

# rep 4 - ECA2687 ABZ
grouped_ABZ_ECA2687_rep4 <- df_ABZ_ECA2687_rep4 %>% 
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2687_rep4_select <- grouped_ABZ_ECA2687_rep4 %>% 
  dplyr::select(strain, condition, replicate, fitness)

# rep 5 - ECA2687 ABZ
grouped_ABZ_ECA2687_rep5 <- df_ABZ_ECA2687_rep5 %>% 
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2687_rep5_select <- grouped_ABZ_ECA2687_rep5 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 6 - ECA2687 ABZ
grouped_ABZ_ECA2687_rep6 <- df_ABZ_ECA2687_rep6 %>% #issue with 1 
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2687_rep6_select <- grouped_ABZ_ECA2687_rep6 %>% 
  dplyr::select(strain, condition, replicate, fitness)

# rep 7 - ECA2687 ABZ
grouped_ABZ_ECA2687_rep7 <- df_ABZ_ECA2687_rep7 %>%
  group_nest(strain, condition, replicate, keep = TRUE) %>%
  mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2687_rep7_select <- grouped_ABZ_ECA2687_rep7 %>% 
  dplyr::select(strain, condition, replicate, fitness)


# rep 8 - ECA2687 ABZ
# df_ABZ_ECA2687_rep8   <- df_ABZ_ECA2687_rep8  %>%
#   filter(!(allele_freq == 0.701568627))
grouped_ABZ_ECA2687_rep8 <- df_ABZ_ECA2687_rep8 %>%
  group_nest(strain, condition, replicate, keep = TRUE) %>%
  mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2687_rep8_select <- grouped_ABZ_ECA2687_rep8 %>%
  select(strain, condition, replicate, fitness)


# rep 9 - ECA2687 ABZ
# df_ABZ_ECA2687_rep9   <- df_ABZ_ECA2687_rep9  %>%
#   filter(!(allele_freq == 0.00005920 & generation == 7))
grouped_ABZ_ECA2687_rep9 <- df_ABZ_ECA2687_rep9 %>%
  group_nest(strain, condition, replicate, keep = TRUE) %>%
  mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2687_rep9_select <- grouped_ABZ_ECA2687_rep9 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 10 - ECA2687 ABZ
# df_ABZ_ECA2687_rep10   <- df_ABZ_ECA2687_rep10  %>%
#   filter(!(allele_freq == 0.000934361 & generation == 1))
grouped_ABZ_ECA2687_rep10 <- df_ABZ_ECA2687_rep10 %>% #Issue bc of 0 in rep
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2687_rep10_select <- grouped_ABZ_ECA2687_rep10 %>%
  dplyr::select(strain, condition, replicate, fitness)

fitness_df_ECA2687_ABZ <- bind_rows(grouped_ABZ_ECA2687_rep1,
                                    grouped_ABZ_ECA2687_rep2,
                                    grouped_ABZ_ECA2687_rep3,
                                    grouped_ABZ_ECA2687_rep4,
                                    grouped_ABZ_ECA2687_rep5,
                                    grouped_ABZ_ECA2687_rep6,
                                    grouped_ABZ_ECA2687_rep7,
                                    grouped_ABZ_ECA2687_rep8,
                                    grouped_ABZ_ECA2687_rep9,
                                    grouped_ABZ_ECA2687_rep10)

# ####################################################
#                                                    #
#          Competitive Fitness  - ECA2689 DMSO       #
#                                                    #
######################################################
# pull out one condition for one strain - ECA2689 DMSO
df_DMSO_ECA2689_rep1 <- v3_eca2689_abz %>%
  dplyr::filter(strain == "ECA2689")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 1)

df_DMSO_ECA2689_rep2 <- v3_eca2689_abz %>%
  dplyr::filter(strain == "ECA2689")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 2)

df_DMSO_ECA2689_rep3 <- v3_eca2689_abz  %>%
  dplyr::filter(strain == "ECA2689")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 3)

df_DMSO_ECA2689_rep4 <- v3_eca2689_abz %>%
  dplyr::filter(strain == "ECA2689")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 4)

df_DMSO_ECA2689_rep5 <- v3_eca2689_abz  %>%
  dplyr::filter(strain == "ECA2689")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 5)

df_DMSO_ECA2689_rep6 <- v3_eca2689_abz %>%
  dplyr::filter(strain == "ECA2689")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 6)

df_DMSO_ECA2689_rep7 <- v3_eca2689_abz %>%
  dplyr::filter(strain == "ECA2689")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 7)

df_DMSO_ECA2689_rep8 <- v3_eca2689_abz  %>%
  dplyr::filter(strain == "ECA2689")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 8)

df_DMSO_ECA2689_rep9 <- v3_eca2689_abz %>%
  dplyr::filter(strain == "ECA2689")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 9)

df_DMSO_ECA2689_rep10 <- v3_eca2689_abz %>%
  dplyr::filter(strain == "ECA2689")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 10)

# rep 1 - ECA2689 DMSO 
grouped_DMSO_ECA2689_rep1 <- df_DMSO_ECA2689_rep1 %>%
  group_nest(strain, condition, replicate, keep = TRUE) %>%
  mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2689_rep1_select <- grouped_DMSO_ECA2689_rep1 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 2 - ECA2689 DMSO 
# View(df_DMSO_ECA2689_rep2) - missing a gen
grouped_DMSO_ECA2689_rep2 <- df_DMSO_ECA2689_rep2 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2689_rep2_select <- grouped_DMSO_ECA2689_rep2 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 3 - ECA2689 DMSO 
# View(df_DMSO_ECA2689_rep3)
grouped_DMSO_ECA2689_rep3 <- df_DMSO_ECA2689_rep3 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2689_rep3_select <- grouped_DMSO_ECA2689_rep3 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 4 - ECA2689 DMSO 
grouped_DMSO_ECA2689_rep4 <- df_DMSO_ECA2689_rep4 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2689_rep4_select <- grouped_DMSO_ECA2689_rep4 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 5 - ECA2689 DMSO 
# gen 5 catg. as gen 6. Fix.
df_DMSO_ECA2689_rep5[df_DMSO_ECA2689_rep5$generation == 6, "generation"] <- 5
grouped_DMSO_ECA2689_rep5 <- df_DMSO_ECA2689_rep5 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2689_rep5_select <- grouped_DMSO_ECA2689_rep5 %>%
  dplyr::select(strain, condition, replicate, fitness)


# rep 6 - ECA2689 DMSO 
# df_DMSO_ECA2689_rep6   <- df_DMSO_ECA2689_rep6  %>%
#   filter(!(allele_freq == 0.316404128 & generation == 1))
grouped_DMSO_ECA2689_rep6 <- df_DMSO_ECA2689_rep6 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2689_rep6_select <- grouped_DMSO_ECA2689_rep6 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 7 - ECA2689 DMSO 
# gen 1 cat as gen 2 fix. 
df_DMSO_ECA2689_rep7[df_DMSO_ECA2689_rep7$generation == 2, "generation"] <- 1
grouped_DMSO_ECA2689_rep7 <- df_DMSO_ECA2689_rep7 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2689_rep7_select <- grouped_DMSO_ECA2689_rep7 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 8 - ECA2689 DMSO 
grouped_DMSO_ECA2689_rep8 <- df_DMSO_ECA2689_rep8 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2689_rep8_select <- grouped_DMSO_ECA2689_rep8 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 9 - ECA2689 DMSO 
df_DMSO_ECA2689_rep9[df_DMSO_ECA2689_rep9$generation == 2, "generation"] <- 1
grouped_DMSO_ECA2689_rep9 <- df_DMSO_ECA2689_rep9 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2689_rep9_select <- grouped_DMSO_ECA2689_rep9 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 10 - ECA2689 DMSO 
grouped_DMSO_ECA2689_rep10 <- df_DMSO_ECA2689_rep10 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2689_rep10_select <- grouped_DMSO_ECA2689_rep10 %>%
  dplyr::select(strain, condition, replicate, fitness)

fitness_df_ECA2689_DMSO <- bind_rows(grouped_DMSO_ECA2689_rep1_select,
                                     grouped_DMSO_ECA2689_rep2_select,
                                     grouped_DMSO_ECA2689_rep3_select,
                                     grouped_DMSO_ECA2689_rep4_select,
                                     grouped_DMSO_ECA2689_rep5_select,
                                     grouped_DMSO_ECA2689_rep6_select,
                                     grouped_DMSO_ECA2689_rep7_select,
                                     grouped_DMSO_ECA2689_rep8_select,
                                     grouped_DMSO_ECA2689_rep9_select,
                                     grouped_DMSO_ECA2689_rep10_select)


######################################################
#                                                    #
#          Competitive Fitness  - ECA2689 ABZ        #
#                                                    #
######################################################
# pull out one condition for one strain - ECA2689 ABZ
df_ABZ_ECA2689_rep1 <- v3_eca2689_abz %>%
  filter(strain == "ECA2689")%>%
  filter(condition == "albendazole")%>%
  filter(replicate == 1)

df_ABZ_ECA2689_rep2 <- v3_eca2689_abz %>%
  filter(strain == "ECA2689")%>%
  filter(condition == "albendazole")%>%
  filter(replicate == 2)

df_ABZ_ECA2689_rep3 <- v3_eca2689_abz  %>%
  filter(strain == "ECA2689")%>%
  filter(condition == "albendazole")%>%
  filter(replicate == 3)

df_ABZ_ECA2689_rep4 <- v3_eca2689_abz %>%
  filter(strain == "ECA2689")%>%
  filter(condition == "albendazole")%>%
  filter(replicate == 4)

df_ABZ_ECA2689_rep5 <- v3_eca2689_abz %>%
  filter(strain == "ECA2689")%>%
  filter(condition == "albendazole")%>%
  filter(replicate == 5)

df_ABZ_ECA2689_rep6 <- v3_eca2689_abz %>%
  filter(strain == "ECA2689")%>%
  filter(condition == "albendazole")%>%
  filter(replicate == 6)

df_ABZ_ECA2689_rep7 <- v3_eca2689_abz  %>%
  filter(strain == "ECA2689")%>%
  filter(condition == "albendazole")%>%
  filter(replicate == 7)

df_ABZ_ECA2689_rep8 <- v3_eca2689_abz  %>%
  filter(strain == "ECA2689")%>%
  filter(condition == "albendazole")%>%
  filter(replicate == 8)

df_ABZ_ECA2689_rep9 <- v3_eca2689_abz %>%
  filter(strain == "ECA2689")%>%
  filter(condition == "albendazole")%>%
  filter(replicate == 9)

df_ABZ_ECA2689_rep10 <- v3_eca2689_abz %>%
  filter(strain == "ECA2689")%>%
  filter(condition == "albendazole")%>%
  filter(replicate == 10)

# rep 1 - ECA2689 ABZ
grouped_ABZ_ECA2689_rep1 <- df_ABZ_ECA2689_rep1 %>%
  group_nest(strain, condition, replicate, keep = TRUE) %>%
  mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2689_rep1_select <- grouped_ABZ_ECA2689_rep1 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 2 - ECA2689 ABZ
grouped_ABZ_ECA2689_rep2 <- df_ABZ_ECA2689_rep2 %>% #issue with 1
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2689_rep2_select <- grouped_ABZ_ECA2689_rep2 %>%
  select(strain, condition, replicate, fitness)

# rep 3 - ECA2689 ABZ
grouped_ABZ_ECA2689_rep3 <- df_ABZ_ECA2689_rep3 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2689_rep3_select <- grouped_ABZ_ECA2689_rep3 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 4 - ECA2689 ABZ
grouped_ABZ_ECA2689_rep4 <- df_ABZ_ECA2689_rep4 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2689_rep4_select <- grouped_ABZ_ECA2689_rep4 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 5 - ECA2689 ABZ
grouped_ABZ_ECA2689_rep5 <- df_ABZ_ECA2689_rep5 %>% #issue with 1 
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2689_rep5_select <- grouped_ABZ_ECA2689_rep5 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 6 - ECA2689 ABZ
grouped_ABZ_ECA2689_rep6 <- df_ABZ_ECA2689_rep6 %>% #issue with 1 
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2689_rep6_select <- grouped_ABZ_ECA2689_rep6 %>% 
  dplyr::select(strain, condition, replicate, fitness)

# rep 7 - ECA2689 ABZ
grouped_ABZ_ECA2689_rep7 <- df_ABZ_ECA2689_rep7 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2689_rep7_select <- grouped_ABZ_ECA2689_rep7 %>% 
  dplyr::select(strain, condition, replicate, fitness)

# rep 8 - ECA2689 ABZ
grouped_ABZ_ECA2689_rep8 <- df_ABZ_ECA2689_rep8 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2689_rep8_select <- grouped_ABZ_ECA2689_rep8 %>% 
  dplyr::select(strain, condition, replicate, fitness)

# rep 9 - ECA2689 ABZ
df_ABZ_ECA2689_rep9[df_ABZ_ECA2689_rep9$generation == 4, "generation"] <- 3
grouped_ABZ_ECA2689_rep9 <- df_ABZ_ECA2689_rep9 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2689_rep9_select <- grouped_ABZ_ECA2689_rep9 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 10 - ECA2689 - ABZ
grouped_ABZ_ECA2689_rep10 <- df_ABZ_ECA2689_rep10 %>% #Issue bc of 0 in rep
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2689_rep10_select <- grouped_ABZ_ECA2689_rep10 %>%
  dplyr::select(strain, condition, replicate, fitness)

fitness_df_ECA2689_ABZ <- bind_rows(grouped_ABZ_ECA2689_rep1_select,
                                    grouped_ABZ_ECA2689_rep2_select,
                                    grouped_ABZ_ECA2689_rep3_select,
                                    grouped_ABZ_ECA2689_rep4_select,
                                    grouped_ABZ_ECA2689_rep5_select,
                                    grouped_ABZ_ECA2689_rep6_select,
                                    grouped_ABZ_ECA2689_rep7_select,
                                    grouped_ABZ_ECA2689_rep8_select,
                                    grouped_ABZ_ECA2689_rep9_select,
                                    grouped_ABZ_ECA2689_rep10_select)

# ####################################################
#                                                    #
#          Competitive Fitness  - ECA2688 DMSO       #
#                                                    #
######################################################
# pull out one condition for one strain - ECA2688 DMSO
df_DMSO_ECA2688_rep1 <- v3_eca2688_abz %>%
  dplyr::filter(strain == "ECA2688")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 1)

df_DMSO_ECA2688_rep2 <- v3_eca2688_abz %>%
  dplyr::filter(strain == "ECA2688")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 2)

df_DMSO_ECA2688_rep3 <- v3_eca2688_abz  %>%
  dplyr::filter(strain == "ECA2688")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 3)

df_DMSO_ECA2688_rep4 <- v3_eca2688_abz %>%
  dplyr::filter(strain == "ECA2688")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 4)

df_DMSO_ECA2688_rep5 <- v3_eca2688_abz %>%
  dplyr::filter(strain == "ECA2688")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 5)

df_DMSO_ECA2688_rep6 <- v3_eca2688_abz %>%
  dplyr::filter(strain == "ECA2688")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 6)

df_DMSO_ECA2688_rep7 <- v3_eca2688_abz %>%
  dplyr::filter(strain == "ECA2688")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 7)

df_DMSO_ECA2688_rep8 <- v3_eca2688_abz %>%
  dplyr::filter(strain == "ECA2688")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 8)

df_DMSO_ECA2688_rep9 <- v3_eca2688_abz %>%
  dplyr::filter(strain == "ECA2688")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 9)

df_DMSO_ECA2688_rep10 <- v3_eca2688_abz %>%
  dplyr::filter(strain == "ECA2688")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(replicate == 10)

# rep 1 - ECA2688 DMSO 
grouped_DMSO_ECA2688_rep1 <- df_DMSO_ECA2688_rep1 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2688_rep1_select <- grouped_DMSO_ECA2688_rep1 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 2 - ECA2688 DMSO 
grouped_DMSO_ECA2688_rep2 <- df_DMSO_ECA2688_rep2 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2688_rep2_select <- grouped_DMSO_ECA2688_rep2 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 3 - ECA2688 DMSO 
# df_DMSO_ECA2688_rep3   <- df_DMSO_ECA2688_rep3  %>%
#   filter(!(allele_freq == 0.002424242 & generation == 7))
grouped_DMSO_ECA2688_rep3 <- df_DMSO_ECA2688_rep3 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2688_rep3_select <- grouped_DMSO_ECA2688_rep3 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 4 - ECA2688 DMSO 
grouped_DMSO_ECA2688_rep4 <- df_DMSO_ECA2688_rep4 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2688_rep4_select <- grouped_DMSO_ECA2688_rep4 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 5 - ECA2688 DMSO 
grouped_DMSO_ECA2688_rep5 <- df_DMSO_ECA2688_rep5 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2688_rep5_select <- grouped_DMSO_ECA2688_rep5 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 6 - ECA2688 DMSO 
grouped_DMSO_ECA2688_rep6 <- df_DMSO_ECA2688_rep6 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2688_rep6_select <- grouped_DMSO_ECA2688_rep6 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 7 - ECA2688 DMSO 
grouped_DMSO_ECA2688_rep7 <- df_DMSO_ECA2688_rep7 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2688_rep7_select <- grouped_DMSO_ECA2688_rep7 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 8 - ECA2688 DMSO 
grouped_DMSO_ECA2688_rep8 <- df_DMSO_ECA2688_rep8 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2688_rep8_select <- grouped_DMSO_ECA2688_rep8 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 9 - ECA2688 DMSO 
grouped_DMSO_ECA2688_rep9 <- df_DMSO_ECA2688_rep9 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2688_rep9_select <- grouped_DMSO_ECA2688_rep9 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 10 - ECA2688 DMSO 
grouped_DMSO_ECA2688_rep10 <- df_DMSO_ECA2688_rep10 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2688_rep10_select <- grouped_DMSO_ECA2688_rep10 %>%
  dplyr::select(strain, condition, replicate, fitness)

fitness_df_ECA2688_DMSO <- bind_rows(grouped_DMSO_ECA2688_rep1_select,
                                     grouped_DMSO_ECA2688_rep2_select,
                                     grouped_DMSO_ECA2688_rep3_select,
                                     grouped_DMSO_ECA2688_rep4_select,
                                     grouped_DMSO_ECA2688_rep5_select,
                                     grouped_DMSO_ECA2688_rep6_select,
                                     grouped_DMSO_ECA2688_rep7_select,
                                     grouped_DMSO_ECA2688_rep8_select,
                                     grouped_DMSO_ECA2688_rep9_select,
                                     grouped_DMSO_ECA2688_rep10_select)


######################################################
#                                                    #
#          Competitive Fitness  - ECA2688 ABZ        #
#                                                    #
######################################################
# pull out one condition for one strain - ECA2688 ABZ
df_ABZ_ECA2688_rep1 <- v3_eca2688_abz %>%
  filter(strain == "ECA2688")%>%
  filter(condition == "albendazole")%>%
  filter(replicate == 1)

df_ABZ_ECA2688_rep2 <- v3_eca2688_abz  %>%
  filter(strain == "ECA2688")%>%
  filter(condition == "albendazole")%>%
  filter(replicate == 2)

df_ABZ_ECA2688_rep3 <- v3_eca2688_abz   %>%
  filter(strain == "ECA2688")%>%
  filter(condition == "albendazole")%>%
  filter(replicate == 3)

df_ABZ_ECA2688_rep4 <- v3_eca2688_abz %>%
  filter(strain == "ECA2688")%>%
  filter(condition == "albendazole")%>%
  filter(replicate == 4)

df_ABZ_ECA2688_rep5 <- v3_eca2688_abz  %>%
  filter(strain == "ECA2688")%>%
  filter(condition == "albendazole")%>%
  filter(replicate == 5)

df_ABZ_ECA2688_rep6 <- v3_eca2688_abz  %>%
  filter(strain == "ECA2688")%>%
  filter(condition == "albendazole")%>%
  filter(replicate == 6)

df_ABZ_ECA2688_rep7 <- v3_eca2688_abz   %>%
  filter(strain == "ECA2688")%>%
  filter(condition == "albendazole")%>%
  filter(replicate == 7)

df_ABZ_ECA2688_rep8 <- v3_eca2688_abz  %>%
  filter(strain == "ECA2688")%>%
  filter(condition == "albendazole")%>%
  filter(replicate == 8)

df_ABZ_ECA2688_rep9 <- v3_eca2688_abz  %>%
  filter(strain == "ECA2688")%>%
  filter(condition == "albendazole")%>%
  filter(replicate == 9)

df_ABZ_ECA2688_rep10 <- v3_eca2688_abz  %>%
  filter(strain == "ECA2688")%>%
  filter(condition == "albendazole")%>%
  filter(replicate == 10)

# rep 1 - ECA2688 ABZ
grouped_ABZ_ECA2688_rep1 <- df_ABZ_ECA2688_rep1 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2688_rep1_select <- grouped_ABZ_ECA2688_rep1 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 2 - ECA2688 ABZ
grouped_ABZ_ECA2688_rep2 <- df_ABZ_ECA2688_rep2 %>% 
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2688_rep2_select <- grouped_ABZ_ECA2688_rep2 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 3 - ECA2688 ABZ
grouped_ABZ_ECA2688_rep3 <- df_ABZ_ECA2688_rep3 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2688_rep3_select <- grouped_ABZ_ECA2688_rep3 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 4 - ECA2688 ABZ
grouped_ABZ_ECA2688_rep4 <- df_ABZ_ECA2688_rep4 %>% 
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2688_rep4_select <- grouped_ABZ_ECA2688_rep4 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 5 - ECA2688 ABZ
grouped_ABZ_ECA2688_rep5 <- df_ABZ_ECA2688_rep5 %>% 
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2688_rep5_select <- grouped_ABZ_ECA2688_rep5 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 6 - ECA2688 ABZ
# grouped_ABZ_ECA2688_rep6 <- df_ABZ_ECA2688_rep6 %>% #missing gen 7 
#   dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
#   dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
#   dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
#   dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))
# 
# grouped_ABZ_ECA2688_rep6_select <- grouped_ABZ_ECA2688_rep6 %>%
#   dplyr::select(strain, condition, replicate, fitness)

# rep 7 - ECA2688 ABZ
grouped_ABZ_ECA2688_rep7 <- df_ABZ_ECA2688_rep7 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2688_rep7_select <- grouped_ABZ_ECA2688_rep7 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 8 - ECA2688 ABZ
grouped_ABZ_ECA2688_rep8 <- df_ABZ_ECA2688_rep8 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2688_rep8_select <- grouped_ABZ_ECA2688_rep8 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 9 - ECA2688 ABZ
grouped_ABZ_ECA2688_rep9 <- df_ABZ_ECA2688_rep9 %>%
  group_nest(strain, condition, replicate, keep = TRUE) %>%
  mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2688_rep9_select <- grouped_ABZ_ECA2688_rep9 %>%
  select(strain, condition, replicate, fitness)

# rep 10 - ECA2688 ABZ
grouped_ABZ_ECA2688_rep10 <- df_ABZ_ECA2688_rep10 %>% #Issue bc of 0 in rep
  group_nest(strain, condition, replicate, keep = TRUE) %>%
  mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2688_rep10_select <- grouped_ABZ_ECA2688_rep10 %>%
  select(strain, condition, replicate, fitness)

fitness_df_ECA2688_ABZ <- bind_rows(grouped_ABZ_ECA2688_rep1_select,
                                    grouped_ABZ_ECA2688_rep2_select,
                                    grouped_ABZ_ECA2688_rep3_select,
                                    grouped_ABZ_ECA2688_rep4_select,
                                    grouped_ABZ_ECA2688_rep5_select,
                                    # grouped_ABZ_ECA2688_rep6_select,
                                    grouped_ABZ_ECA2688_rep7_select,
                                    grouped_ABZ_ECA2688_rep8_select,
                                    grouped_ABZ_ECA2688_rep9_select,
                                    grouped_ABZ_ECA2688_rep10_select)


######################################################
#                                                    #
#          Competitive Fitness  - ECA2690 DMSO       #
#                                                    #
######################################################
# pull out one condition for one strain - ECA2690 DMSO
df_DMSO_ECA2690_rep1 <- v3_eca2690_abz %>%
  filter(strain == "ECA2690")%>%
  filter(condition == "DMSO")%>%
  filter(replicate == 1)

df_DMSO_ECA2690_rep2 <- v3_eca2690_abz %>%
  filter(strain == "ECA2690")%>%
  filter(condition == "DMSO")%>%
  filter(replicate == 2)

df_DMSO_ECA2690_rep3 <- v3_eca2690_abz  %>%
  filter(strain == "ECA2690")%>%
  filter(condition == "DMSO")%>%
  filter(replicate == 3)

df_DMSO_ECA2690_rep4 <- v3_eca2690_abz %>%
  filter(strain == "ECA2690")%>%
  filter(condition == "DMSO")%>%
  filter(replicate == 4)

df_DMSO_ECA2690_rep5 <- v3_eca2690_abz %>%
  filter(strain == "ECA2690")%>%
  filter(condition == "DMSO")%>%
  filter(replicate == 5)

df_DMSO_ECA2690_rep6 <- v3_eca2690_abz %>%
  filter(strain == "ECA2690")%>%
  filter(condition == "DMSO")%>%
  filter(replicate == 6)

df_DMSO_ECA2690_rep7 <- v3_eca2690_abz %>%
  filter(strain == "ECA2690")%>%
  filter(condition == "DMSO")%>%
  filter(replicate == 7)

df_DMSO_ECA2690_rep8 <- v3_eca2690_abz %>%
  filter(strain == "ECA2690")%>%
  filter(condition == "DMSO")%>%
  filter(replicate == 8)

df_DMSO_ECA2690_rep9 <- v3_eca2690_abz %>%
  filter(strain == "ECA2690")%>%
  filter(condition == "DMSO")%>%
  filter(replicate == 9)

df_DMSO_ECA2690_rep10 <- v3_eca2690_abz %>%
  filter(strain == "ECA2690")%>%
  filter(condition == "DMSO")%>%
  filter(replicate == 10)

# rep 1 - ECA2690 DMSO
grouped_DMSO_ECA2690_rep1 <- df_DMSO_ECA2690_rep1 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2690_rep1_select <- grouped_DMSO_ECA2690_rep1 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 2 - ECA2690 DMSO
grouped_DMSO_ECA2690_rep2 <- df_DMSO_ECA2690_rep2 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2690_rep2_select <- grouped_DMSO_ECA2690_rep2 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 3 - ECA2690 DMSO
grouped_DMSO_ECA2690_rep3 <- df_DMSO_ECA2690_rep3 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2690_rep3_select <- grouped_DMSO_ECA2690_rep3 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 4 - ECA2690 DMSO
grouped_DMSO_ECA2690_rep4 <- df_DMSO_ECA2690_rep4 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2690_rep4_select <- grouped_DMSO_ECA2690_rep4 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 5 - ECA2690 DMSO
grouped_DMSO_ECA2690_rep5 <- df_DMSO_ECA2690_rep5 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2690_rep5_select <- grouped_DMSO_ECA2690_rep5 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 6 - ECA2690 DMSO
grouped_DMSO_ECA2690_rep6 <- df_DMSO_ECA2690_rep6 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2690_rep6_select <- grouped_DMSO_ECA2690_rep6 %>%
  dplyr::select(strain, condition, replicate, fitness)


# rep 7 - ECA2690 DMSO
# df_DMSO_ECA2690_rep7   <- df_DMSO_ECA2690_rep7  %>%
#   filter(!(allele_freq == 0.000121205 & generation == 3))

grouped_DMSO_ECA2690_rep7 <- df_DMSO_ECA2690_rep7 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2690_rep7_select <- grouped_DMSO_ECA2690_rep7 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 8 - ECA2690 DMSO
grouped_DMSO_ECA2690_rep8 <- df_DMSO_ECA2690_rep8 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2690_rep8_select <- grouped_DMSO_ECA2690_rep8 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 9 - ECA2690 DMSO
grouped_DMSO_ECA2690_rep9 <- df_DMSO_ECA2690_rep9 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2690_rep9_select <- grouped_DMSO_ECA2690_rep9 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 10 - ECA2690 DMSO
# View(df_DMSO_ECA2690_rep10)
# df_DMSO_ECA2690_rep10   <- df_DMSO_ECA2690_rep10  %>%
#   filter(!(allele_freq == 0.009371043 & generation == 3))
# df_DMSO_ECA2690_rep10   <- df_DMSO_ECA2690_rep10  %>%
#   filter(!(allele_freq == 0.003306635 & generation == 5))
grouped_DMSO_ECA2690_rep10 <- df_DMSO_ECA2690_rep10 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_DMSO_ECA2690_rep10_select <- grouped_DMSO_ECA2690_rep10 %>%
  dplyr::select(strain, condition, replicate, fitness)

fitness_df_ECA2690_DMSO <- bind_rows(grouped_DMSO_ECA2690_rep1_select,
                                     grouped_DMSO_ECA2690_rep2_select,
                                     grouped_DMSO_ECA2690_rep3_select,
                                     grouped_DMSO_ECA2690_rep4_select,
                                     grouped_DMSO_ECA2690_rep5_select,
                                     grouped_DMSO_ECA2690_rep6_select,
                                     grouped_DMSO_ECA2690_rep7_select,
                                     grouped_DMSO_ECA2690_rep8_select,
                                     grouped_DMSO_ECA2690_rep9_select,
                                     grouped_DMSO_ECA2690_rep10_select)


######################################################
#                                                    #
#          Competitive Fitness  - ECA2690 ABZ        #
#                                                    #
######################################################
# pull out one condition for one strain - ECA2690 ABZ
df_ABZ_ECA2690_rep1 <- v3_eca2690_abz %>%
  filter(strain == "ECA2690")%>%
  filter(condition == "albendazole")%>%
  filter(replicate == 1)

df_ABZ_ECA2690_rep2 <- v3_eca2690_abz  %>%
  filter(strain == "ECA2690")%>%
  filter(condition == "albendazole")%>%
  filter(replicate == 2)

df_ABZ_ECA2690_rep3 <- v3_eca2690_abz  %>%
  filter(strain == "ECA2690")%>%
  filter(condition == "albendazole")%>%
  filter(replicate == 3)

df_ABZ_ECA2690_rep4 <- v3_eca2690_abz %>%
  filter(strain == "ECA2690")%>%
  filter(condition == "albendazole")%>%
  filter(replicate == 4)

df_ABZ_ECA2690_rep5 <- v3_eca2690_abz  %>%
  filter(strain == "ECA2690")%>%
  filter(condition == "albendazole")%>%
  filter(replicate == 5)

df_ABZ_ECA2690_rep6 <- v3_eca2690_abz  %>%
  filter(strain == "ECA2690")%>%
  filter(condition == "albendazole")%>%
  filter(replicate == 6)

df_ABZ_ECA2690_rep7 <- v3_eca2690_abz   %>%
  filter(strain == "ECA2690")%>%
  filter(condition == "albendazole")%>%
  filter(replicate == 7)

df_ABZ_ECA2690_rep8 <- v3_eca2690_abz  %>%
  filter(strain == "ECA2690")%>%
  filter(condition == "albendazole")%>%
  filter(replicate == 8)

df_ABZ_ECA2690_rep9 <- v3_eca2690_abz  %>%
  filter(strain == "ECA2690")%>%
  filter(condition == "albendazole")%>%
  filter(replicate == 9)

df_ABZ_ECA2690_rep10 <- v3_eca2690_abz  %>%
  filter(strain == "ECA2690")%>%
  filter(condition == "albendazole")%>%
  filter(replicate == 10)

# rep 1 - ECA2690 ABZ
grouped_ABZ_ECA2690_rep1 <- df_ABZ_ECA2690_rep1 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2690_rep1_select <- grouped_ABZ_ECA2690_rep1 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 2 - ECA2690 ABZ
grouped_ABZ_ECA2690_rep2 <- df_ABZ_ECA2690_rep2 %>% #issue with 1
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2690_rep2_select <- grouped_ABZ_ECA2690_rep2 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 3 - ECA2690 ABZ
grouped_ABZ_ECA2690_rep3 <- df_ABZ_ECA2690_rep3 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2690_rep3_select <- grouped_ABZ_ECA2690_rep3 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 4 - ECA2690 ABZ
# df_ABZ_ECA2690_rep4 <- df_ABZ_ECA2690_rep4 %>% 
#   filter(!(allele_freq == 0.002374894 & generation ==3))
grouped_ABZ_ECA2690_rep4 <- df_ABZ_ECA2690_rep4 %>% #issue with 1
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2690_rep4_select <- grouped_ABZ_ECA2690_rep4 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 5 - ECA2690 ABZ
# df_ABZ_ECA2690_rep5 <- df_ABZ_ECA2690_rep5 %>% 
#   filter(!(allele_freq == 0.062372517))
# df_ABZ_ECA2690_rep5 <- df_ABZ_ECA2690_rep5 %>% 
#   filter(!(allele_freq == 0.004906364))
grouped_ABZ_ECA2690_rep5 <- df_ABZ_ECA2690_rep5 %>% #issue with 1
  group_nest(strain, condition, replicate, keep = TRUE) %>%
  mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2690_rep5_select <- grouped_ABZ_ECA2690_rep5 %>%
  select(strain, condition, replicate, fitness)

# rep 6 - ECA2690 ABZ
#typo with gen fix 2 > 1 
df_ABZ_ECA2690_rep6[df_ABZ_ECA2690_rep6$generation == 2, "generation"] <- 1

grouped_ABZ_ECA2690_rep6 <- df_ABZ_ECA2690_rep6 %>% #issue with 1
  group_nest(strain, condition, replicate, keep = TRUE) %>%
  mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2690_rep6_select <- grouped_ABZ_ECA2690_rep6 %>%
  select(strain, condition, replicate, fitness)

# rep 7 - ECA2690 ABZ
grouped_ABZ_ECA2690_rep7 <- df_ABZ_ECA2690_rep7 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2690_rep7_select <- grouped_ABZ_ECA2690_rep7 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 8 - ECA2690 ABZ
grouped_ABZ_ECA2690_rep8 <- df_ABZ_ECA2690_rep8 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2690_rep8_select <- grouped_ABZ_ECA2690_rep8 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 9 - ECA2690 ABZ
# df_ABZ_ECA2690_rep9 <- df_ABZ_ECA2690_rep9 %>% 
#   filter(!(allele_freq == 0.001837485))
grouped_ABZ_ECA2690_rep9 <- df_ABZ_ECA2690_rep9 %>%
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2690_rep9_select <- grouped_ABZ_ECA2690_rep9 %>%
  dplyr::select(strain, condition, replicate, fitness)

# rep 10 - ECA2690 ABZ
grouped_ABZ_ECA2690_rep10 <- df_ABZ_ECA2690_rep10 %>% #Issue bc of 0 in rep
  dplyr::group_nest(strain, condition, replicate, keep = TRUE) %>%
  dplyr::mutate(ci_values = map(data, comp_calc, gen1_id = 1))%>%
  dplyr::mutate(lm_coef = map_dbl(ci_values, ci_lss))%>%
  dplyr::mutate(fitness = map_dbl(lm_coef, cal_fitness))

grouped_ABZ_ECA2690_rep10_select <- grouped_ABZ_ECA2690_rep10 %>%
  dplyr::select(strain, condition, replicate, fitness)

fitness_df_ECA2690_ABZ <- bind_rows(grouped_ABZ_ECA2690_rep1_select,
                                    grouped_ABZ_ECA2690_rep2_select,
                                    grouped_ABZ_ECA2690_rep3_select,
                                    grouped_ABZ_ECA2690_rep4_select,
                                    grouped_ABZ_ECA2690_rep5_select,
                                    grouped_ABZ_ECA2690_rep6_select,
                                    grouped_ABZ_ECA2690_rep7_select,
                                    grouped_ABZ_ECA2690_rep8_select,
                                    grouped_ABZ_ECA2690_rep9_select,
                                    grouped_ABZ_ECA2690_rep10_select)

# Total fitness 
BZML_total_abz_fitness <- bind_rows(fitness_df_ECA882_ABZ,
                                fitness_df_n2_ABZ,
                                fitness_df_PHX2647_ABZ,
                                fitness_df_PHX2502_ABZ,
                                fitness_df_PHX2481_ABZ,
                                fitness_df_ECA2687_ABZ,
                                fitness_df_ECA2689_ABZ,
                                fitness_df_ECA2688_ABZ,
                                fitness_df_ECA2690_ABZ)

save(BZML_total_abz_fitness, file = "data/competition_assays/BZML_total_abz_fitness_20240113.csv")


BZML_total_dsmo_fitness_v3 <- bind_rows(fitness_df_ECA882_DMSO,
                                        fitness_df_n2_DMSO_v3,
                                    fitness_df_PHX2647_DMSO,
                                    fitness_df_PHX2502_DMSO,
                                    fitness_df_PHX2481_DMSO,
                                    fitness_df_ECA2687_DMSO,
                                    fitness_df_ECA2689_DMSO,
                                    fitness_df_ECA2688_DMSO,
                                    fitness_df_ECA2690_DMSO)

save(BZML_total_dsmo_fitness_v3, file = "data/competition_assays/BZML_total_dsmo_fitness_v3_20240113.csv")
