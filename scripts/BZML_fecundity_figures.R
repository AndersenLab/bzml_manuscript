# BZML manuscript | fecundity analysis and figures 

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/"))
setwd("../")

### load R packages####
install.packages("multcomp")
install.packages("agricolae")
install.packages("ggpubr")
library(ggpubr)
library(ggplot2)
library(agricolae)
library(tidyverse)
library(ggpubr)
library(ggtree) #issue with R.4.1.1   
library(cowplot)
library(ggalt) #issue with R.4.1.1   
library(ggthemes) 
library(ggrepel)
library(ggplot2)
library(ggsignif)
library(dplyr)
library(stats)
library(multcomp)
library(agricolae)
library(plyr)

##########################################
#           Figure X                     #
#         Fecundity Assays               #
##########################################
## Lifetime fecundity bar plot

#Set strain colors 
col_comp <- c("N2" = "#FFA500", #orange
              "ECA882" = "#FFE400", #yellow
              "PHX2647" = "#9DDD00", #lime green
              "PHX2502" = "#00B558", #medium green
              "PHX2481" = "#008878", #teal
              "ECA2687" = "#14517F", #lighest dark blue
              "ECA2688" = "#00C3D7", #medium blue
              "ECA2689" = "#2D52E4", #dark blue
              "ECA2690" = "#3A0045") #dark blue

#Assign strain labels
strain_labels <- c("N2" = "N2", "ECA882" = "ben-1", "PHX2647" = "avr-14",
                   "PHX2502" = "avr-15", "PHX2481" = "glc-1", "ECA2687" = "avr-14, avr-15",
                   "ECA2688" = "avr-14, glc-1", "ECA2689" = "avr-15 glc-1", "ECA2690" = "avr-14 avr-15,glc-1")

strain_order <- c("N2", "ECA882", "PHX2647", "PHX2502", "PHX2481", "ECA2687", "ECA2688", "ECA2689", "ECA2690")

#############################################
#           Figure X A                      #
#  Lifetime Fecundity in DMSO (Sets 1-5)    #
#############################################
sets1_5_mean_sd_DMSO  <- mean_sd_DMSO %>% #Filter out sets 6 and 7
  filter(!(set %in% c(6, 7)))

DMSO_outliers <- sets1_5_mean_sd_DMSO %>% #DMSO outliers
  filter(set %in% c("1", "2") & avg_total_offspring < 5)

#Remove the following samples from the dataframe due to less than 5 offspring in set 1. 
#ECA2687 rep 3 and 7 (Divide by 8)
#ECA2688 rep 2, 6, 10 (Divide by 7)
#ECA2689 rep 3 (Divide by 9)
#ECA2690 rep 1 (Divide by 9)
#ECA882 rep 2 (Divide by 9)
#PHX2481 rep 9 (Divide by 9)
# PHX2502 rep 3 (Divide by 9)
sets1_5_mean_sd_DMSO_outliers_removed <- sets1_5_mean_sd_DMSO %>%
  filter(!(strain == "ECA2687" & (replicate == 3 | replicate == 7))) %>% 
  filter(!(strain == "ECA2688" & (replicate == 2 | replicate == 6 | replicate == 10))) %>% 
  filter(!(strain == "ECA2689" & (replicate == 3))) %>% 
  filter(!(strain == "ECA2690" & (replicate == 1))) %>%           
  filter(!(strain == "ECA882" & (replicate == 2))) %>%  
  filter(!(strain == "PHX2481" & (replicate == 9))) %>%   
  filter(!(strain == "PHX2502" & (replicate == 3)))   

# strain_order <- c("N2", "ECA882", "PHX2647", "PHX2502", "PHX2481", "ECA2687", "ECA2688", "ECA2689", "ECA2690")

write.csv(sets1_5_mean_sd_DMSO_outliers_removed, "sets1_5_mean_sd_DMSO_outliers_removed.csv", row.names = FALSE)

#Import dataframe 
sets1_5_mean_sd_DMSO_outliers_removed <- read.csv("sets1_5_mean_sd_DMSO_outliers_removed.csv")

#Pull out what we need for plotting - WHY AREN'T YOU WORKING 
# sets1_5_mean_sd_DMSO_outliers_removed <- sets1_5_mean_sd_DMSO_outliers_removed %>%
#   select(- scorer_initials.x,
#          - scorer_initials.y,
#          - scorer_initials.x.x,
#          - scorer_initials.y.y,
#          - adult_questionable,
#          - plate_notes,
#          - plate_created,
#          - mean.y,
#          - sd.y)
         

#ECA2687 df - divide by 8
DMSO_ECA2687 <- sets1_5_mean_sd_DMSO_outliers_removed  %>%
  filter(strain == "ECA2687")

T_mean_sd_DMSO_ECA2687 <- DMSO_ECA2687 %>%
  group_by(strain, condition) %>%
  summarize(
    lifetime_1 = sum(total_offspring_1) / 8,
    lifetime_2 = sum(total_offspring_2) / 8,
    lifetime_3 = sum(total_offspring_3) / 8,
    lifetime_4 = sum(total_offspring_4) / 8,
    lifetime_mean = mean(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)),
    lifetime_sd = sd(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)))

T_mean_sd_DMSO_ECA2687 <- T_mean_sd_DMSO_ECA2687 %>%
  mutate(condition = "DMSO",
         strain = "ECA2687")


#ECA2688 - divide by 7 
DMSO_ECA2688 <- sets1_5_mean_sd_DMSO_outliers_removed  %>%
  filter(strain == "ECA2688")

T_mean_sd_DMSO_ECA2688 <- DMSO_ECA2688 %>%
  group_by(strain, condition) %>%
  summarize(
    lifetime_1 = sum(total_offspring_1) / 7,
    lifetime_2 = sum(total_offspring_2) / 7,
    lifetime_3 = sum(total_offspring_3) / 7,
    lifetime_4 = sum(total_offspring_4) / 7,
    lifetime_mean = mean(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)),
    lifetime_sd = sd(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)))

T_mean_sd_DMSO_ECA2688 <- T_mean_sd_DMSO_ECA2688 %>%
  mutate(condition = "DMSO",
         strain = "ECA2688")

#ECA2689 - divide by 9
DMSO_ECA2689 <- sets1_5_mean_sd_DMSO_outliers_removed %>%
  filter(strain == "ECA2689")

T_mean_sd_DMSO_ECA2689 <- DMSO_ECA2689 %>%
  group_by(strain, condition) %>%
  summarize(
    lifetime_1 = sum(total_offspring_1) / 9,
    lifetime_2 = sum(total_offspring_2) / 9,
    lifetime_3 = sum(total_offspring_3) / 9,
    lifetime_4 = sum(total_offspring_4) / 9,
    lifetime_mean = mean(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)),
    lifetime_sd = sd(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)))

T_mean_sd_DMSO_ECA2689 <- T_mean_sd_DMSO_ECA2689 %>%
  mutate(condition = "DMSO",
         strain = "ECA2689")

#ECA2690 - divide by 9
DMSO_ECA2690 <- sets1_5_mean_sd_DMSO_outliers_removed  %>%
  filter(strain == "ECA2690")

T_mean_sd_DMSO_ECA2690 <- DMSO_ECA2690 %>%
  group_by(strain, condition) %>%
  summarize(
    lifetime_1 = sum(total_offspring_1) / 9,
    lifetime_2 = sum(total_offspring_2) / 9,
    lifetime_3 = sum(total_offspring_3) / 9,
    lifetime_4 = sum(total_offspring_4) / 9,
    lifetime_mean = mean(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)),
    lifetime_sd = sd(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)))

T_mean_sd_DMSO_ECA2690 <- T_mean_sd_DMSO_ECA2690 %>%
  mutate(condition = "DMSO",
         strain = "ECA2690")

#ECA882 - divide by 9
DMSO_ECA882 <- sets1_5_mean_sd_DMSO_outliers_removed  %>%
  filter(strain == "ECA882")

T_mean_sd_DMSO_ECA882 <- DMSO_ECA882 %>%
  group_by(strain, condition) %>%
  summarize(
    lifetime_1 = sum(total_offspring_1) / 9,
    lifetime_2 = sum(total_offspring_2) / 9,
    lifetime_3 = sum(total_offspring_3) / 9,
    lifetime_4 = sum(total_offspring_4) / 9,
    lifetime_mean = mean(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)),
    lifetime_sd = sd(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)))

T_mean_sd_DMSO_ECA882 <- T_mean_sd_DMSO_ECA882 %>%
  mutate(condition = "DMSO",
         strain = "ECA882")

#PHX2481 - divide by 9
DMSO_PHX2481 <- sets1_5_mean_sd_DMSO_outliers_removed  %>%
  filter(strain == "PHX2481")

T_mean_sd_DMSO_PHX2481 <- DMSO_PHX2481 %>%
  group_by(strain, condition) %>%
  summarize(
    lifetime_1 = sum(total_offspring_1) / 9,
    lifetime_2 = sum(total_offspring_2) / 9,
    lifetime_3 = sum(total_offspring_3) / 9,
    lifetime_4 = sum(total_offspring_4) / 9,
    lifetime_mean = mean(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)),
    lifetime_sd = sd(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)))

T_mean_sd_DMSO_PHX2481 <- T_mean_sd_DMSO_PHX2481 %>%
  mutate(condition = "DMSO",
         strain = "PHX2481")

# PHX2502 - divide by 9
DMSO_PHX2502 <- sets1_5_mean_sd_DMSO_outliers_removed  %>%
  filter(strain == "PHX2502")

T_mean_sd_DMSO_PHX2502 <- DMSO_PHX2502 %>%
  group_by(strain, condition) %>%
  summarize(
    lifetime_1 = sum(total_offspring_1) / 9,
    lifetime_2 = sum(total_offspring_2) / 9,
    lifetime_3 = sum(total_offspring_3) / 9,
    lifetime_4 = sum(total_offspring_4) / 9,
    lifetime_mean = mean(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)),
    lifetime_sd = sd(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)))

T_mean_sd_DMSO_PHX2502 <- T_mean_sd_DMSO_PHX2502 %>%
  mutate(condition = "DMSO",
         strain = "PHX2502")
 
# N2 - divide by 10
DMSO_N2 <- sets1_5_mean_sd_DMSO_outliers_removed  %>%
  filter(strain == "N2")

T_mean_sd_DMSO_N2 <- DMSO_N2 %>%
  group_by(strain, condition) %>%
  summarize(
    lifetime_1 = sum(total_offspring_1) / 10,
    lifetime_2 = sum(total_offspring_2) / 10,
    lifetime_3 = sum(total_offspring_3) / 10,
    lifetime_4 = sum(total_offspring_4) / 10,
    lifetime_mean = mean(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)),
    lifetime_sd = sd(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)))

T_mean_sd_DMSO_N2 <- T_mean_sd_DMSO_N2 %>%
  mutate(condition = "DMSO",
         strain = "N2")

# PHX2647 - divide by 10
DMSO_PHX2647 <- sets1_5_mean_sd_DMSO_outliers_removed %>%
  filter(strain == "PHX2647")

T_mean_sd_DMSO_PHX2647 <- DMSO_PHX2647 %>%
  group_by(strain, condition) %>%
  summarize(
    lifetime_1 = sum(total_offspring_1) / 10,
    lifetime_2 = sum(total_offspring_2) / 10,
    lifetime_3 = sum(total_offspring_3) / 10,
    lifetime_4 = sum(total_offspring_4) / 10,
    lifetime_mean = mean(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)),
    lifetime_sd = sd(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)))

T_mean_sd_DMSO_PHX2647 <- T_mean_sd_DMSO_PHX2647 %>%
  mutate(condition = "DMSO",
         strain = "PHX2647")

merged_lifetime_DMSO <- bind_rows(T_mean_sd_DMSO_N2, 
                                  T_mean_sd_DMSO_ECA882, 
                                  T_mean_sd_DMSO_PHX2647,
                                  T_mean_sd_DMSO_PHX2502,
                                  T_mean_sd_DMSO_PHX2481,
                                  T_mean_sd_DMSO_ECA2687,
                                  T_mean_sd_DMSO_ECA2688,
                                  T_mean_sd_DMSO_ECA2689,
                                  T_mean_sd_DMSO_ECA2690)
View(merged_lifetime_DMSO)


#Need to transform df - already done - import 
merged_lifetime <- read.csv("data/fecundity/lifetime_fecundity_allconditions_20231219.csv")

merged_lifetime_DMSO <- read.csv("data/fecundity/lifetime_fecundity_DMSO_20231219.csv")


# Perform ANOVA
merged_lifetime_DMSO_AOV <- merged_lifetime %>%
  dplyr::filter(condition == "DMSO")%>%
  aov(lifetime ~ strain, data = .)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(group1 == 'N2' | group2 == 'N2')

View(merged_lifetime_DMSO_AOV)

# # Check to see if the double (avr-14, avr-15) and triple are sig. different from each other 
# merged_lifetime_DMSOdoubletriple_AOV <- merged_lifetime %>%
#   dplyr::filter(condition == "DMSO")%>%
#   aov(lifetime ~ strain, data = .)%>%
#   rstatix::tukey_hsd()%>%
#   dplyr::filter(group1 == 'ECA2690' | group2 == 'ECA2690')
# 
# View(merged_lifetime_DMSOdoubletriple_AOV) 

# Flip to always have N2 represented in the group 2 column
merged_lifetime_DMSO_AOV <- merged_lifetime_DMSO_AOV %>%
  mutate(
    temp = group1,
    group1 = ifelse(group1 == 'N2', group2, group1),
    group2 = ifelse(temp == 'N2', temp, group2)
  ) %>%
  dplyr::select(-temp)


merged_lifetime_DMSO_AOV <- merged_lifetime_DMSO_AOV %>%
  mutate(
    group2 = case_when(
      group2 %in% c("PHX2481", "PHX2502", "PHX2647") ~ "N2",
      TRUE ~ group2))

merged_lifetime_DMSO_AOV  <- merged_lifetime_DMSO_AOV %>% 
  dplyr::mutate(strain= group1)

# ggsave("figures/fecundity/DMSO_lifetime_sets1_5 _20231218.png", DMSO_lifetime_sets1_5, width = 8, height = 6, dpi = 300)

write.csv(merged_lifetime_DMSO, "data/fecundity/lifetime_fecundity_DMSO_20231219.csv", row.names = FALSE)

#############################################
#           Figure X B                      #
#  Lifetime Fecundity in ABZ (Sets 1-5)     #
#############################################
sets1_5_mean_sd_ABZ  <- mean_sd_ABZ %>% #Filter out sets 6 and 7
  filter(!(set %in% c(6, 7)))

ABZ_outliers <- sets1_5_mean_sd_ABZ %>% #DMSO outliers
  filter(set %in% c("1", "2") & avg_total_offspring < 5)

View(ABZ_outliers)

#Remove the following samples from the dataframe due to less than 5 offspring in set 1. 
#ECA2687 rep 1, 3, and 6 (Divide by 7)
#ECA2690 rep 2 and 4 (Divide by 8)
#PHX2481 reps 1, 2, and 4 (Divide by 7)
#PHX2502 rep 1 (Divide by 9)
sets1_5_mean_sd_ABZ_outliers_removed <- sets1_5_mean_sd_ABZ %>%
  filter(!(strain == "ECA2687" & (replicate == 1 | replicate == 3 | replicate == 6))) %>% 
  filter(!(strain == "ECA2690" & (replicate == 2 | replicate == 4 ))) %>%           
  filter(!(strain == "PHX2481" & (replicate == 1 | replicate == 2 | replicate == 4))) %>%  
  filter(!(strain == "PHX2502" & (replicate == 1)))   

write.csv(sets1_5_mean_sd_ABZ_outliers_removed, "sets1_5_mean_sd_ABZ_outliers_removed.csv", row.names = FALSE)

sets1_5_mean_sd_ABZ_outliers_removed <- read.csv("sets1_5_mean_sd_ABZ_outliers_removed.csv")

strain_order <- c("N2", "ECA882", "PHX2647", "PHX2502", "PHX2481", "ECA2687", "ECA2688", "ECA2689", "ECA2690")

#ECA2687 df - divide by 7 - ABZ
ABZ_ECA2687 <- sets1_5_mean_sd_ABZ_outliers_removed  %>%
  filter(strain == "ECA2687")

T_mean_sd_ABZ_ECA2687 <- ABZ_ECA2687 %>%
  group_by(strain, condition) %>%
  summarize(
    lifetime_1 = sum(total_offspring_1) / 7,
    lifetime_2 = sum(total_offspring_2) / 7,
    lifetime_3 = sum(total_offspring_3) / 7,
    lifetime_4 = sum(total_offspring_4) / 7,
    lifetime_mean = mean(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)),
    lifetime_sd = sd(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)))

T_mean_sd_ABZ_ECA2687 <- T_mean_sd_ABZ_ECA2687 %>%
  mutate(condition = "albendazole",
         strain = "ECA2687")

#ECA2688 - divide by 10 - ABZ
ABZ_ECA2688 <- sets1_5_mean_sd_ABZ_outliers_removed  %>%
  filter(strain == "ECA2688")

T_mean_sd_ABZ_ECA2688 <- ABZ_ECA2688 %>%
  group_by(strain, condition) %>%
  summarize(
    lifetime_1 = sum(total_offspring_1) / 10,
    lifetime_2 = sum(total_offspring_2) / 10,
    lifetime_3 = sum(total_offspring_3) / 10,
    lifetime_4 = sum(total_offspring_4) / 10,
    lifetime_mean = mean(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)),
    lifetime_sd = sd(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)))

T_mean_sd_ABZ_ECA2688 <- T_mean_sd_ABZ_ECA2688 %>%
  mutate(condition = "albendazole",
         strain = "ECA2688")

#ECA2689 - divide by 10 - ABZ
ABZ_ECA2689 <- sets1_5_mean_sd_ABZ_outliers_removed  %>%
  filter(strain == "ECA2689")

T_mean_sd_ABZ_ECA2689 <- ABZ_ECA2689 %>%
  group_by(strain, condition) %>%
  summarize(
    lifetime_1 = sum(total_offspring_1) / 10,
    lifetime_2 = sum(total_offspring_2) / 10,
    lifetime_3 = sum(total_offspring_3) / 10,
    lifetime_4 = sum(total_offspring_4) / 10,
    lifetime_mean = mean(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)),
    lifetime_sd = sd(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)))

T_mean_sd_ABZ_ECA2689 <- T_mean_sd_ABZ_ECA2689 %>%
  mutate(condition = "albendazole",
         strain = "ECA2689")

#ECA2690 - divide by 8 - ABZ
ABZ_ECA2690 <- sets1_5_mean_sd_ABZ_outliers_removed  %>%
  filter(strain == "ECA2690")

T_mean_sd_ABZ_ECA2690 <- ABZ_ECA2690 %>%
  group_by(strain, condition) %>%
  summarize(
    lifetime_1 = sum(total_offspring_1) / 8,
    lifetime_2 = sum(total_offspring_2) / 8,
    lifetime_3 = sum(total_offspring_3) / 8,
    lifetime_4 = sum(total_offspring_4) / 8,
    lifetime_mean = mean(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)),
    lifetime_sd = sd(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)))

T_mean_sd_ABZ_ECA2690 <- T_mean_sd_ABZ_ECA2690 %>%
  mutate(condition = "albendazole",
         strain = "ECA2690")

#ECA882 - divide by 10 - ABZ
ABZ_ECA882 <- sets1_5_mean_sd_ABZ_outliers_removed %>%
  filter(strain == "ECA882")

T_mean_sd_ABZ_ECA882 <- ABZ_ECA882 %>%
  group_by(strain, condition) %>%
  summarize(
    lifetime_1 = sum(total_offspring_1) / 10,
    lifetime_2 = sum(total_offspring_2) / 10,
    lifetime_3 = sum(total_offspring_3) / 10,
    lifetime_4 = sum(total_offspring_4) / 10,
    lifetime_mean = mean(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)),
    lifetime_sd = sd(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)))

T_mean_sd_ABZ_ECA882 <- T_mean_sd_ABZ_ECA882 %>%
  mutate(condition = "albendazole",
         strain = "ECA882")

#PHX2481 - divide by 7 - ABZ
ABZ_PHX2481 <- sets1_5_mean_sd_ABZ_outliers_removed  %>%
  filter(strain == "PHX2481")

T_mean_sd_ABZ_PHX2481 <- ABZ_PHX2481 %>%
  group_by(strain, condition) %>%
  summarize(
    lifetime_1 = sum(total_offspring_1) / 7,
    lifetime_2 = sum(total_offspring_2) / 7,
    lifetime_3 = sum(total_offspring_3) / 7,
    lifetime_4 = sum(total_offspring_4) / 7,
    lifetime_mean = mean(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)),
    lifetime_sd = sd(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)))

T_mean_sd_ABZ_PHX2481 <- T_mean_sd_ABZ_PHX2481 %>%
  mutate(condition = "albendazole",
         strain = "PHX2481")

# PHX2502 - divide by 9 - ABZ
ABZ_PHX2502 <- sets1_5_mean_sd_ABZ_outliers_removed  %>%
  filter(strain == "PHX2502")

T_mean_sd_ABZ_PHX2502 <- ABZ_PHX2502 %>%
  group_by(strain, condition) %>%
  summarize(
    lifetime_1 = sum(total_offspring_1) / 9,
    lifetime_2 = sum(total_offspring_2) / 9,
    lifetime_3 = sum(total_offspring_3) / 9,
    lifetime_4 = sum(total_offspring_4) / 9,
    lifetime_mean = mean(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)),
    lifetime_sd = sd(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)))

T_mean_sd_ABZ_PHX2502 <- T_mean_sd_ABZ_PHX2502 %>%
  mutate(condition = "albendazole",
         strain = "PHX2502")

# N2 - divide by 10 - ABZ
ABZ_N2 <- sets1_5_mean_sd_ABZ_outliers_removed  %>%
  filter(strain == "N2")

T_mean_sd_ABZ_N2 <- ABZ_N2 %>%
  group_by(strain, condition) %>%
  summarize(
    lifetime_1 = sum(total_offspring_1) / 10,
    lifetime_2 = sum(total_offspring_2) / 10,
    lifetime_3 = sum(total_offspring_3) / 10,
    lifetime_4 = sum(total_offspring_4) / 10,
    lifetime_mean = mean(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)),
    lifetime_sd = sd(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)))

T_mean_sd_ABZ_N2 <- T_mean_sd_ABZ_N2 %>%
  mutate(condition = "albendazole",
         strain = "N2")

# PHX2647 - divide by 10 - ABZ
ABZ_PHX2647 <- sets1_5_mean_sd_ABZ_outliers_removed  %>%
  filter(strain == "PHX2647")

T_mean_sd_ABZ_PHX2647 <- ABZ_PHX2647 %>%
  group_by(strain, condition) %>%
  summarize(
    lifetime_1 = sum(total_offspring_1) / 10,
    lifetime_2 = sum(total_offspring_2) / 10,
    lifetime_3 = sum(total_offspring_3) / 10,
    lifetime_4 = sum(total_offspring_4) / 10,
    lifetime_mean = mean(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)),
    lifetime_sd = sd(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)))

T_mean_sd_ABZ_PHX2647 <- T_mean_sd_ABZ_PHX2647 %>%
  mutate(condition = "albendazole",
         strain = "PHX2647")

merged_lifetime_ABZ <- bind_rows(T_mean_sd_ABZ_N2, 
                                  T_mean_sd_ABZ_ECA882, 
                                  T_mean_sd_ABZ_PHX2647,
                                  T_mean_sd_ABZ_PHX2502,
                                  T_mean_sd_ABZ_PHX2481,
                                  T_mean_sd_ABZ_ECA2687,
                                  T_mean_sd_ABZ_ECA2688,
                                  T_mean_sd_ABZ_ECA2689,
                                  T_mean_sd_ABZ_ECA2690)
View(merged_lifetime_ABZ)

merged_lifetime_ABZ<- read.csv("data/fecundity/lifetime_fecundity_ABZ_20231219.csv")

# ggsave("figures/fecundity/ABZ_lifetime_sets1_5 _20231218.png", ABZ_lifetime_sets1_5, width = 8, height = 6, dpi = 300)

write.csv(merged_lifetime_ABZ, "data/fecundity/lifetime_fecundity_ABZ_20231219.csv", row.names = FALSE)

merged_lifetime <- read.csv("data/fecundity/lifetime_fecundity_allconditions_20231219.csv")

View(merged_lifetime)

# Perform ANOVA
merged_lifetime_ABZ_AOV <- merged_lifetime %>%
  dplyr::filter(condition == "albendazole")%>%
  aov(lifetime ~ strain, data = .)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(group1 == 'N2' | group2 == 'N2')


# Check to see if the double (avr-14, avr-15) and triple are sig. different from each other 
merged_lifetime_ABZdoubletriple_AOV <- merged_lifetime %>%
  dplyr::filter(condition == "albendazole")%>%
  aov(lifetime ~ strain, data = .)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(group1 == 'ECA2690' | group2 == 'ECA2690')

View(merged_lifetime_ABZdoubletriple_AOV) 

# Flip to always have N2 represented in the group 2 column
merged_lifetime_ABZ_AOV <- merged_lifetime_ABZ_AOV %>%
  mutate(
    temp = group1,
    group1 = ifelse(group1 == 'N2', group2, group1),
    group2 = ifelse(temp == 'N2', temp, group2)
  ) %>%
  dplyr::select(-temp)


merged_lifetime_ABZ_AOV <- merged_lifetime_ABZ_AOV %>%
  mutate(
    group2 = case_when(
      group2 %in% c("PHX2481", "PHX2502", "PHX2647") ~ "N2",
      TRUE ~ group2))

merged_lifetime_ABZ_AOV  <- merged_lifetime_ABZ_AOV %>% 
  dplyr::mutate(strain= group1)

View(merged_lifetime_ABZ_AOV)

write.csv(merged_lifetime_ABZ, "data/fecundity/lifetime_fecundity_ABZ_20231219.csv", row.names = FALSE)

#############################################
#           Figure X C                      #
#  Lifetime Fecundity in IVM (Sets 1-5)     #
#############################################
sets1_5_mean_sd_IVM  <- mean_sd_IVM %>% #Filter out sets 6 and 7
  filter(!(set %in% c(6, 7)))

IVM_outliers <- sets1_5_mean_sd_IVM %>% #DMSO outliers
  filter(set %in% c("1", "2") & avg_total_offspring < 5)

View(IVM_outliers)

#Remove the following samples from the dataframe due to less than 5 offspring in set 1. 
#ECA2689 rep 3 (Divide by 9)
#ECA2690 rep 2 (Divide by 9)
#N2 rep 7 (Divide by 9)
#PHX2481 reps 6 (Divide by 9)
#PHX2502 rep 5, 8, and 9 (Divide by 7)

sets1_5_mean_sd_IVM_outliers_removed <- sets1_5_mean_sd_IVM %>%
  filter(!(strain == "ECA2689" & (replicate == 3))) %>% 
  filter(!(strain == "ECA2690" & (replicate == 2))) %>%   
  filter(!(strain == "N2" & (replicate == 7))) %>%   
  filter(!(strain == "PHX2481" & (replicate == 6))) %>%  
  filter(!(strain == "PHX2502" & (replicate == 5 | replicate == 8 | replicate == 9)))   

write.csv(sets1_5_mean_sd_IVM_outliers_removed, "sets1_5_mean_sd_IVM_outliers_removed.csv", row.names = FALSE)

sets1_5_mean_sd_IVM_outliers_removed <- read.csv("sets1_5_mean_sd_IVM_outliers_removed.csv")

# strain_order <- c("N2", "ECA882", "PHX2647", "PHX2502", "PHX2481", "ECA2687", "ECA2688", "ECA2689", "ECA2690")

#ECA2687 df - divide by 10 - IVM DONE
IVM_ECA2687 <- mean_sd_IVM  %>%
  filter(strain == "ECA2687")

T_mean_sd_IVM_ECA2687 <- IVM_ECA2687 %>%
  group_by(strain, condition) %>%
  summarize(
    lifetime_1 = sum(total_offspring_1) / 10,
    lifetime_2 = sum(total_offspring_2) / 10,
    lifetime_3 = sum(total_offspring_3) / 10,
    lifetime_4 = sum(total_offspring_4) / 10,
    lifetime_mean = mean(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)),
    lifetime_sd = sd(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)))

T_mean_sd_IVM_ECA2687 <- T_mean_sd_IVM_ECA2687 %>%
  mutate(condition = "ivermectin",
         strain = "ECA2687")

#ECA2688 - divide by 10 - IVM - DONE
IVM_ECA2688 <- mean_sd_IVM  %>%
  filter(strain == "ECA2688")

T_mean_sd_IVM_ECA2688 <- IVM_ECA2688 %>%
  group_by(strain, condition) %>%
  summarize(
    lifetime_1 = sum(total_offspring_1) / 10,
    lifetime_2 = sum(total_offspring_2) / 10,
    lifetime_3 = sum(total_offspring_3) / 10,
    lifetime_4 = sum(total_offspring_4) / 10,
    lifetime_mean = mean(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)),
    lifetime_sd = sd(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)))

T_mean_sd_IVM_ECA2688 <- T_mean_sd_IVM_ECA2688 %>%
  mutate(condition = "ivermectin",
         strain = "ECA2688")

#ECA2689 - divide by 9 - IVM - DONE
IVM_ECA2689 <- mean_sd_IVM  %>%
  filter(strain == "ECA2689")

T_mean_sd_IVM_ECA2689 <- IVM_ECA2689 %>%
  group_by(strain, condition) %>%
  summarize(
    lifetime_1 = sum(total_offspring_1) / 9,
    lifetime_2 = sum(total_offspring_2) / 9,
    lifetime_3 = sum(total_offspring_3) / 9,
    lifetime_4 = sum(total_offspring_4) / 9,
    lifetime_mean = mean(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)),
    lifetime_sd = sd(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)))

T_mean_sd_IVM_ECA2689 <- T_mean_sd_IVM_ECA2689 %>%
  mutate(condition = "ivermectin",
         strain = "ECA2689")

#ECA2690 - divide by 9 - IVM - DONE
IVM_ECA2690 <- mean_sd_IVM  %>%
  filter(strain == "ECA2690")

T_mean_sd_IVM_ECA2690 <- IVM_ECA2690 %>%
  group_by(strain, condition) %>%
  summarize(
    lifetime_1 = sum(total_offspring_1) / 9,
    lifetime_2 = sum(total_offspring_2) / 9,
    lifetime_3 = sum(total_offspring_3) / 9,
    lifetime_4 = sum(total_offspring_4) / 9,
    lifetime_mean = mean(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)),
    lifetime_sd = sd(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)))

T_mean_sd_IVM_ECA2690 <- T_mean_sd_IVM_ECA2690 %>%
  mutate(condition = "ivermectin",
         strain = "ECA2690")

#ECA882 - divide by 10 - IVM
IVM_ECA882 <- mean_sd_IVM  %>%
  filter(strain == "ECA882")

T_mean_sd_IVM_ECA882 <- IVM_ECA882 %>%
  group_by(strain, condition) %>%
  summarize(
    lifetime_1 = sum(total_offspring_1) / 10,
    lifetime_2 = sum(total_offspring_2) / 10,
    lifetime_3 = sum(total_offspring_3) / 10,
    lifetime_4 = sum(total_offspring_4) / 10,
    lifetime_mean = mean(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)),
    lifetime_sd = sd(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)))

T_mean_sd_IVM_ECA882 <- T_mean_sd_IVM_ECA882 %>%
  mutate(condition = "ivermectin",
         strain = "ECA882")

#PHX2481 - divide by 9 - IVM - DONE
IVM_PHX2481 <- mean_sd_IVM  %>%
  filter(strain == "PHX2481")

T_mean_sd_IVM_PHX2481 <- IVM_PHX2481 %>%
  group_by(strain, condition) %>%
  summarize(
    lifetime_1 = sum(total_offspring_1) / 9,
    lifetime_2 = sum(total_offspring_2) / 9,
    lifetime_3 = sum(total_offspring_3) / 9,
    lifetime_4 = sum(total_offspring_4) / 9,
    lifetime_mean = mean(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)),
    lifetime_sd = sd(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)))

T_mean_sd_IVM_PHX2481 <- T_mean_sd_IVM_PHX2481 %>%
  mutate(condition = "ivermectin",
         strain = "PHX2481")

# PHX2502 - divide by 7 - IVM
IVM_PHX2502 <- mean_sd_IVM  %>%
  filter(strain == "PHX2502")

T_mean_sd_IVM_PHX2502 <- IVM_PHX2502 %>%
  group_by(strain, condition) %>%
  summarize(
    lifetime_1 = sum(total_offspring_1) / 7,
    lifetime_2 = sum(total_offspring_2) / 7,
    lifetime_3 = sum(total_offspring_3) / 7,
    lifetime_4 = sum(total_offspring_4) / 7,
    lifetime_mean = mean(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)),
    lifetime_sd = sd(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)))

T_mean_sd_IVM_PHX2502 <- T_mean_sd_IVM_PHX2502 %>%
  mutate(condition = "ivermectin",
         strain = "PHX2502")

# N2 - divide by 9 - IVM
IVM_N2 <- mean_sd_IVM  %>%
  filter(strain == "N2")

T_mean_sd_IVM_N2 <- IVM_N2 %>%
  group_by(strain, condition) %>%
  summarize(
    lifetime_1 = sum(total_offspring_1) / 9,
    lifetime_2 = sum(total_offspring_2) / 9,
    lifetime_3 = sum(total_offspring_3) / 9,
    lifetime_4 = sum(total_offspring_4) / 9,
    lifetime_mean = mean(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)),
    lifetime_sd = sd(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)))

T_mean_sd_IVM_N2 <- T_mean_sd_IVM_N2 %>%
  mutate(condition = "ivermectin",
         strain = "N2")

# PHX2647 - divide by 10 - IVM
IVM_PHX2647 <- mean_sd_IVM  %>%
  filter(strain == "PHX2647")

T_mean_sd_IVM_PHX2647 <- IVM_PHX2647 %>%
  group_by(strain, condition) %>%
  summarize(
    lifetime_1 = sum(total_offspring_1) / 10,
    lifetime_2 = sum(total_offspring_2) / 10,
    lifetime_3 = sum(total_offspring_3) / 10,
    lifetime_4 = sum(total_offspring_4) / 10,
    lifetime_mean = mean(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)),
    lifetime_sd = sd(c(lifetime_1, lifetime_2, lifetime_3, lifetime_4)))

T_mean_sd_IVM_PHX2647 <- T_mean_sd_IVM_PHX2647 %>%
  mutate(condition = "ivermectin",
         strain = "PHX2647")

merged_lifetime_IVM <- bind_rows(T_mean_sd_IVM_N2, 
                                 T_mean_sd_IVM_ECA882, 
                                 T_mean_sd_IVM_PHX2647,
                                 T_mean_sd_IVM_PHX2502,
                                 T_mean_sd_IVM_PHX2481,
                                 T_mean_sd_IVM_ECA2687,
                                 T_mean_sd_IVM_ECA2688,
                                 T_mean_sd_IVM_ECA2689,
                                 T_mean_sd_IVM_ECA2690)
View(merged_lifetime_IVM)

merged_lifetime_IVM <- read.csv("data/fecundity/lifetime_fecundity_IVM_20231219.csv")


merged_lifetime_DMSO <- merged_lifetime %>%
  dplyr::filter(condition == "DMSO")


write.csv(merged_lifetime_IVM, "data/fecundity/lifetime_fecundity_IVM_20231219.csv", row.names = FALSE)

merged_lifetime <- read.csv("data/fecundity/lifetime_fecundity_allconditions_20231219.csv")

# Perform ANOVA
merged_lifetime_IVM_AOV <- merged_lifetime %>%
  dplyr::filter(condition == "ivermectin")%>%
  aov(lifetime ~ strain, data = .)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(group1 == 'N2' | group2 == 'N2')

View(merged_lifetime_IVM_AOV)

# Flip to always have N2 represented in the group 2 column
merged_lifetime_IVM_AOV <- merged_lifetime_IVM_AOV %>%
  mutate(
    temp = group1,
    group1 = ifelse(group1 == 'N2', group2, group1),
    group2 = ifelse(temp == 'N2', temp, group2)
  ) %>%
  dplyr::select(-temp)


merged_lifetime_IVM_AOV <- merged_lifetime_IVM_AOV %>%
  mutate(
    group2 = case_when(
      group2 %in% c("PHX2481", "PHX2502", "PHX2647") ~ "N2",
      TRUE ~ group2))

merged_lifetime_IVM_AOV  <- merged_lifetime_IVM_AOV %>% 
  dplyr::mutate(strain= group1)

View(merged_lifetime_IVM_AOV)

write.csv(merged_lifetime_IVM, "data/fecundity/lifetime_fecundity_IVM_20231219.csv", row.names = FALSE)

#############################################
#                                           #
#       Lifetime Fecundity Figures          #
#                                           #
#############################################
# Lifetime fecundity figures 
#Plot DMSO lifetime fecundity
DMSO_lifetime <- merged_lifetime_DMSO  %>%
  mutate(strain = factor(strain, levels = strain_order)) %>%  
  ggplot(aes(x = strain, y = lifetime_mean, fill = strain)) +
  geom_bar(stat = 'identity', color = 'black', width = 0.5, position = position_dodge(0.2), alpha = 0.7) +
  geom_errorbar(aes(ymin = lifetime_mean - lifetime_sd, ymax = lifetime_mean + lifetime_sd),
                width = 0.25, position = position_dodge(0.9)) +  # Add error bars
  geom_point(data = merged_lifetime_DMSO  , aes(y = lifetime_1), position = position_dodge(width = 0.9), color = "black") + #shape = 1, 
  geom_point(data = merged_lifetime_DMSO  , aes(y = lifetime_2), position = position_dodge(width = 0.9), color = "black") + #shape = 2, 
  geom_point(data = merged_lifetime_DMSO  , aes(y = lifetime_3), position = position_dodge(width = 0.9), color = "black") + #shape = 3, 
  geom_point(data = merged_lifetime_DMSO  , aes(y = lifetime_4), position = position_dodge(width = 0.9), color = "black") + #shape = 4, 
  geom_text(data = merged_lifetime_DMSO  , aes(x = strain, y = 350, label = round(lifetime_mean, 2), size = 18),
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3, color = "black") +
  labs(x = NULL, y = NULL, fill = "Strain") +
  scale_fill_manual(values = col_comp, breaks = names(col_comp), labels = strain_labels) +
  stat_pvalue_manual(merged_lifetime_DMSO_AOV, label = "p.adj.signif",y.position = c(380), xmax="group1", remove.bracket = TRUE) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 400), sec.axis = dup_axis(name = "DMSO")) +  # Set y-axis limits
  theme_minimal() +
  theme(axis.text.x = element_blank(),        
        axis.text.y = element_text(size = 12, color = "black"),
        axis.line = element_line(color = "black"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        legend.position = "none",
        legend.box.background = element_blank(),
        legend.key.size = unit(1.5, "lines"),
        plot.margin = margin(5, 10, 10, 10, "pt"),
        plot.title = element_text(hjust = 0, vjust = 2 , margin = margin(10, 0, 20, 0)))

DMSO_lifetime

# ABZ Lifetime
ABZ_lifetime <- merged_lifetime_ABZ  %>%
  mutate(strain = factor(strain, levels = strain_order)) %>%  
  ggplot(aes(x = strain, y = lifetime_mean, fill = strain)) +
  geom_bar(stat = 'identity', color = 'black', width = 0.5, position = position_dodge(0.2), alpha = 0.7) +
  geom_errorbar(aes(ymin = lifetime_mean - lifetime_sd, ymax = lifetime_mean + lifetime_sd),
                width = 0.25, position = position_dodge(0.9)) + 
  geom_point(data = merged_lifetime_ABZ  , aes(y = lifetime_1), position = position_dodge(width = 0.9), color = "black") + #shape = 1, 
  geom_point(data = merged_lifetime_ABZ  , aes(y = lifetime_2), position = position_dodge(width = 0.9), color = "black") + #shape = 2, 
  geom_point(data = merged_lifetime_ABZ  , aes(y = lifetime_3), position = position_dodge(width = 0.9), color = "black") + #shape = 3, 
  geom_point(data = merged_lifetime_ABZ  , aes(y = lifetime_4), position = position_dodge(width = 0.9), color = "black") + #shape = 4, 
  geom_text(data = merged_lifetime_ABZ  , aes(x = strain, y = 350, label = round(lifetime_mean, 2), size = 18),
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3, color = "black") +
  labs(x = NULL, y = "Lifetime fecundity", fill = "Strain") +
  # ggtitle("Albendazole") +
  scale_fill_manual(values = col_comp, breaks = names(col_comp), labels = strain_labels) +
  # scale_x_discrete(breaks = c("N2","ECA882", "PHX2647","PHX2502","PHX2481","ECA2687","ECA2688","ECA2689","ECA2690"),
  #                  labels = c("N2\n \n ", "ben-1\n \n ","avr-14\n \n ","avr-15\n \n ","glc-1\n \n ","avr-14\navr-15\n ","avr-14\nglc-1\n ","avr-15\nglc-1\n ","avr-14\navr-15\nglc-1"))+
  stat_pvalue_manual(merged_lifetime_ABZ_AOV, label = "p.adj.signif",y.position = c(380), xmax="group1", remove.bracket = TRUE) +
  # scale_y_continuous(expand = c(0, 0), limits = c(0, 390)) +  # Set y-axis limits
  scale_y_continuous(expand = c(0, 0), limits = c(0, 400), sec.axis = dup_axis(name = "Albendazole")) +  # Set y-axis limits
  # scale_y_continuous(expand = c(0, 0), limits = c(0, max(merged_lifetime_ABZ $lifetime_mean) + max(merged_lifetime_ABZ$lifetime_sd) + 10)) +  # Adjust y-axis limits
  theme_minimal() +
  theme(axis.text.x = element_blank(),        
        axis.text.y = element_text(size = 12, color = "black"),
        axis.line = element_line(color = "black"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        legend.position = "none",
        legend.box.background = element_blank(),
        legend.key.size = unit(1.5, "lines"),
        plot.margin = margin(5, 10, 10, 10, "pt"),
        plot.title = element_text(hjust = 0, vjust = 2 , margin = margin(10, 0, 20, 0)))

# IVM Lifetime
IVM_lifetime <- merged_lifetime_IVM  %>%
  mutate(strain = factor(strain, levels = strain_order)) %>%  
  ggplot(aes(x = strain, y = lifetime_mean, fill = strain)) +
  geom_bar(stat = 'identity', color = 'black', width = 0.5, position = position_dodge(0.2), alpha = 0.7) +
  geom_errorbar(aes(ymin = lifetime_mean - lifetime_sd, ymax = lifetime_mean + lifetime_sd),
                width = 0.25, position = position_dodge(0.9)) +  
  geom_point(data = merged_lifetime_IVM  , aes(y = lifetime_1), position = position_dodge(width = 0.9), color = "black") + #shape = 1, 
  geom_point(data = merged_lifetime_IVM  , aes(y = lifetime_2), position = position_dodge(width = 0.9), color = "black") + #shape = 2, 
  geom_point(data = merged_lifetime_IVM  , aes(y = lifetime_3), position = position_dodge(width = 0.9), color = "black") + #shape = 3, 
  geom_point(data = merged_lifetime_IVM  , aes(y = lifetime_4), position = position_dodge(width = 0.9), color = "black") + #shape = 4, 
  geom_text(data = merged_lifetime_IVM  , aes(x = strain, y = 350, label = round(lifetime_mean, 2), size = 18),
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3, color = "black") +
  labs(x = NULL, y = NULL, fill = "Strain") +
  scale_fill_manual(values = col_comp, breaks = names(col_comp), labels = strain_labels) +
  scale_x_discrete(breaks = c("N2","ECA882", "PHX2647","PHX2502","PHX2481","ECA2687","ECA2688","ECA2689","ECA2690"),
                   labels = c("N2\n \n ", "ben-1\n \n ","avr-14\n \n ","avr-15\n \n ","glc-1\n \n ","avr-14\navr-15\n ","avr-14\nglc-1\n ","avr-15\nglc-1\n ","avr-14\navr-15\nglc-1"))+
  stat_pvalue_manual(merged_lifetime_IVM_AOV, label = "p.adj.signif",y.position = c(380), xmax="group1", remove.bracket = TRUE) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 410), sec.axis = dup_axis(name = "Ivermectin")) +  # Set y-axis limits
  # scale_y_continuous(expand = c(0, 0), limits = c(0, max(merged_lifetime_ABZ $lifetime_mean) + max(merged_lifetime_ABZ$lifetime_sd) + 10)) +  # Adjust y-axis limits
  theme_minimal() +
        theme(axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 1, margin = unit(c(0, 0, 0, 0), units = "in"), face = "italic"), axis.title = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.line = element_line(color = "black"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        legend.position = "none",
        legend.box.background = element_blank(),
        legend.key.size = unit(1.5, "lines"),
        plot.margin = margin(5, 10, 10, 10, "pt"),
        plot.title = element_text(hjust = 0, vjust = 2 , margin = margin(10, 0, 20, 0)))

IVM_lifetime 




# ggsave("figures/fecundity/IVM_lifetime_sets1_5 _20231218.png", IVM_lifetime_sets1_5, width = 8, height = 6, dpi = 300)
#Combine all plots for lifetimefecundity figure
comb_lifetimefecundity <- cowplot::plot_grid(DMSO_lifetime, 
                                             ABZ_lifetime, 
                                             IVM_lifetime,
                                             ncol = 1,
                                             nrow = 3,
                                             labels = c("A","B","C"),
                                             align = "vh",
                                             axis = "lrbt",
                                             label_size = 12,
                                             label_fontfamily = "Arial",
                                             rel_widths = c(1,1,1),
                                             rel_heights = c(0.8,0.8,0.8))

comb_lifetimefecundity

ggsave(filename = "figures/comb_lifetimefecundity_20240105.eps",plot = comb_lifetimefecundity ,width = 5,height = 6,units = "in")
ggsave(filename = "figures/comb_lifetimefecundity_20240105.png",plot = comb_lifetimefecundity ,width = 5,height = 6,units = "in")
ggsave(filename = "figures/comb_lifetimefecundity_20240105.svg",plot = comb_lifetimefecundity ,width = 5,height = 6,units = "in")



##########################################
#                                        #
#           DAILY FECUNDITY              #
#                                        #
##########################################
#############################################
#           Figure X                        #
#  Daily Fecundity in DMSO (Sets 1-5)       #
#############################################
#Combine all outlier removed datasets for daily fecundity 

#DMSO ANOVA
DMSO_set1 <- sets1_5_mean_sd_DMSO_outliers_removed %>% 
  dplyr::filter(set == "1")
View(DMSO_set1)

DMSO_set2 <- sets1_5_mean_sd_DMSO_outliers_removed %>% 
  dplyr::filter(set == "2")

DMSO_set3 <- sets1_5_mean_sd_DMSO_outliers_removed %>% 
  dplyr::filter(set == "3")

DMSO_set4 <- sets1_5_mean_sd_DMSO_outliers_removed %>% 
  dplyr::filter(set == "4")

DMSO_set5 <- sets1_5_mean_sd_DMSO_outliers_removed %>% 
  dplyr::filter(set == "5")

# Perform ANOVA - Set 1 
daily_fecundity_DMSO_AOV_set1 <- DMSO_set1 %>%
  dplyr::filter(condition == "DMSO")%>%
  aov(avg_total_offspring ~ strain, data = .)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(group1 == 'N2' | group2 == 'N2')

View(DMSO_set1)
View(daily_fecundity_DMSO_AOV_set1)

# Flip to always have N2 represented in the group 2 column
daily_fecundity_DMSO_AOV_set1 <- daily_fecundity_DMSO_AOV_set1 %>%
  mutate(
    temp = group1,
    group1 = ifelse(group1 == 'N2', group2, group1),
    group2 = ifelse(temp == 'N2', temp, group2)
  ) %>%
  dplyr::select(-temp)


daily_fecundity_DMSO_AOV_set1 <- daily_fecundity_DMSO_AOV_set1 %>%
  mutate(
    group2 = case_when(
      group2 %in% c("PHX2481", "PHX2502", "PHX2647") ~ "N2",
      TRUE ~ group2))

daily_fecundity_DMSO_AOV_set1  <- daily_fecundity_DMSO_AOV_set1 %>% 
  dplyr::mutate(strain= group1)

View(daily_fecundity_DMSO_AOV_set1)

# Perform ANOVA - Set 2 
daily_fecundity_DMSO_AOV_set2 <- DMSO_set2 %>%
  dplyr::filter(condition == "DMSO")%>%
  aov(avg_total_offspring ~ strain, data = .)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(group1 == 'N2' | group2 == 'N2')

# Flip to always have N2 represented in the group 2 column
daily_fecundity_DMSO_AOV_set2 <- daily_fecundity_DMSO_AOV_set2 %>%
  mutate(
    temp = group1,
    group1 = ifelse(group1 == 'N2', group2, group1),
    group2 = ifelse(temp == 'N2', temp, group2)
  ) %>%
  dplyr::select(-temp)


daily_fecundity_DMSO_AOV_set2 <- daily_fecundity_DMSO_AOV_set2 %>%
  mutate(
    group2 = case_when(
      group2 %in% c("PHX2481", "PHX2502", "PHX2647") ~ "N2",
      TRUE ~ group2))

daily_fecundity_DMSO_AOV_set2  <- daily_fecundity_DMSO_AOV_set2 %>% 
  dplyr::mutate(strain= group1)
View(daily_fecundity_DMSO_AOV_set2)

# Perform ANOVA - Set 3 
daily_fecundity_DMSO_AOV_set3 <- DMSO_set3 %>%
  dplyr::filter(condition == "DMSO")%>%
  aov(avg_total_offspring ~ strain, data = .)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(group1 == 'N2' | group2 == 'N2')

# Flip to always have N2 represented in the group 2 column
daily_fecundity_DMSO_AOV_set3 <- daily_fecundity_DMSO_AOV_set3 %>%
  mutate(
    temp = group1,
    group1 = ifelse(group1 == 'N2', group2, group1),
    group2 = ifelse(temp == 'N2', temp, group2)
  ) %>%
  dplyr::select(-temp)


daily_fecundity_DMSO_AOV_set3 <- daily_fecundity_DMSO_AOV_set3 %>%
  mutate(
    group2 = case_when(
      group2 %in% c("PHX2481", "PHX2502", "PHX2647") ~ "N2",
      TRUE ~ group2))

daily_fecundity_DMSO_AOV_set3  <- daily_fecundity_DMSO_AOV_set3 %>% 
  dplyr::mutate(strain= group1)
View(daily_fecundity_DMSO_AOV_set3)

# Perform ANOVA - Set 4 
daily_fecundity_DMSO_AOV_set4 <- DMSO_set4 %>%
  dplyr::filter(condition == "DMSO")%>%
  aov(avg_total_offspring ~ strain, data = .)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(group1 == 'N2' | group2 == 'N2')

# Flip to always have N2 represented in the group 2 column
daily_fecundity_DMSO_AOV_set4 <- daily_fecundity_DMSO_AOV_set4 %>%
  mutate(
    temp = group1,
    group1 = ifelse(group1 == 'N2', group2, group1),
    group2 = ifelse(temp == 'N2', temp, group2)
  ) %>%
  dplyr::select(-temp)


daily_fecundity_DMSO_AOV_set4 <- daily_fecundity_DMSO_AOV_set4 %>%
  mutate(
    group2 = case_when(
      group2 %in% c("PHX2481", "PHX2502", "PHX2647") ~ "N2",
      TRUE ~ group2))

daily_fecundity_DMSO_AOV_set4  <- daily_fecundity_DMSO_AOV_set4 %>% 
  dplyr::mutate(strain= group1)
View(daily_fecundity_DMSO_AOV_set4)

# Perform ANOVA - Set 5 
daily_fecundity_DMSO_AOV_set5 <- DMSO_set5 %>%
  dplyr::filter(condition == "DMSO")%>%
  aov(avg_total_offspring ~ strain, data = .)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(group1 == 'N2' | group2 == 'N2')

# Flip to always have N2 represented in the group 2 column
daily_fecundity_DMSO_AOV_set5 <- daily_fecundity_DMSO_AOV_set5 %>%
  mutate(
    temp = group1,
    group1 = ifelse(group1 == 'N2', group2, group1),
    group2 = ifelse(temp == 'N2', temp, group2)
  ) %>%
  dplyr::select(-temp)


daily_fecundity_DMSO_AOV_set5 <- daily_fecundity_DMSO_AOV_set5 %>%
  mutate(
    group2 = case_when(
      group2 %in% c("PHX2481", "PHX2502", "PHX2647") ~ "N2",
      TRUE ~ group2))

daily_fecundity_DMSO_AOV_set5  <- daily_fecundity_DMSO_AOV_set5 %>% 
  dplyr::mutate(strain= group1)
View(daily_fecundity_DMSO_AOV_set5)

#clean up dataframes 
write.csv(DMSO_set1, "data/fecundity/DMSO_set1.csv", row.names = FALSE, col.names = FALSE)
write.csv(DMSO_set2, "data/fecundity/DMSO_set2.csv", row.names = FALSE, col.names = FALSE)
write.csv(DMSO_set3, "data/fecundity/DMSO_set3.csv", row.names = FALSE, col.names = FALSE)
write.csv(DMSO_set4, "data/fecundity/DMSO_set4.csv", row.names = FALSE, col.names = FALSE)
write.csv(DMSO_set5, "data/fecundity/DMSO_set5.csv", row.names = FALSE, col.names = FALSE)
write.csv(daily_fecundity_DMSO_AOV_set1, "data/fecundity/daily_fecundity_DMSO_AOV_set1.csv", row.names = FALSE, col.names = FALSE)
write.csv(daily_fecundity_DMSO_AOV_set2, "data/fecundity/daily_fecundity_DMSO_AOV_set2.csv", row.names = FALSE, col.names = FALSE)
write.csv(daily_fecundity_DMSO_AOV_set3, "data/fecundity/daily_fecundity_DMSO_AOV_set3.csv", row.names = FALSE, col.names = FALSE)
write.csv(daily_fecundity_DMSO_AOV_set4, "data/fecundity/daily_fecundity_DMSO_AOV_set4.csv", row.names = FALSE, col.names = FALSE)
write.csv(daily_fecundity_DMSO_AOV_set5, "data/fecundity/daily_fecundity_DMSO_AOV_set5.csv", row.names = FALSE, col.names = FALSE)


merged_daily <- read.csv("data/fecundity/DMSO_merged_daily.csv")
View(merged_daily)


#Set custom order for how strains should be displayed in boxplot 
custom_order <- c("N2", "ECA882", "PHX2647", "PHX2502", "PHX2481", "ECA2687", "ECA2688", "ECA2689", "ECA2690")


# Create the box plot with jittered points and set boxes
DMSO_daily <- ggplot(merged_daily, aes(x = factor(strain, levels = custom_order), y = avg_total_offspring, fill = factor(strain, levels = custom_order))) +
  geom_boxplot(width = 0.5, color = "black", outlier.shape = NA) +
  geom_jitter(shape = 21, position = position_jitter(0.2), size = 0.5, alpha = 0.5, fill = "black") +
  labs(x = NULL, y = NULL) +
  scale_fill_manual(values = col_comp, breaks = custom_order, labels = strain_labels) +
  stat_pvalue_manual(merged_daily, label = "p.adj.signif", y.position = c(140), xmax = "group1", remove.bracket = TRUE) +
  theme_bw() +   
  facet_grid(. ~ set, switch = "y", labeller = labeller(set = c("Day 1", "Day 2", "Day 3", "Day 4", "Days 5 - 7"))) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.line = element_line(color = "black"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 12, vjust = 1, color = "black"),
        strip.background = element_blank(), 
        legend.position = "none",
        legend.box.background = element_blank(),
        legend.key.size = unit(1.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(10, 10, 20, 10, "pt")) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 25), oob = scales::squish, sec.axis = dup_axis(name = "DMSO"))

# Print the plot
print(DMSO_daily)


# TEST - DMSO Daily as its own supplemental figure
DMSO_daily_individual <- ggplot(merged_daily, aes(x = factor(strain, levels = custom_order), y = avg_total_offspring, fill = factor(strain, levels = custom_order))) +
  geom_boxplot(width = 0.5, color = "black", outlier.shape = NA) +
  geom_jitter(shape = 21, position = position_jitter(0.2), size = 0.5, alpha = 0.5, fill = "black") +
  labs(x = "Strain", y = "Fecundity") +
  scale_fill_manual(values = col_comp, breaks = custom_order, labels = strain_labels) +
  stat_pvalue_manual(merged_daily, label = "p.adj.signif", y.position = c(140), xmax = "group1", remove.bracket = TRUE) +
  theme_bw() +   
  facet_wrap(.~ set, ncol = 3) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.line = element_line(color = "black"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        legend.text = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 12, vjust = 1, color = "black"),
        strip.background = element_blank(), 
        legend.position = "none",
        legend.box.background = element_blank(),
        legend.key.size = unit(1.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(10, 10, 20, 10, "pt")) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 25), oob = scales::squish, sec.axis = dup_axis(name = "DMSO"))

legend_plot <- ggplot(merged_daily) +
  geom_point(aes(x = strain, y = 1, color = factor(col_comp, levels = custom_order)), size = 4) +
  theme_void() +
  theme(legend.position = "none")

legend_plot

# Print the plot
print(DMSO_daily_individual)

ggsave(filename = "figures/DMSO_daily_individual_20231221.svg",plot = DMSO_daily_individual ,width = 10,height = 6,units = "in")


#############################################
#           Figure X                        #
#  Daily Fecundity in ABZ (Sets 1-5)        #
#############################################

#ABZ ANOVA
ABZ_set1 <- sets1_5_mean_sd_ABZ_outliers_removed %>% 
  dplyr::filter(set == "1")
View(ABZ_set1)

ABZ_set2 <- sets1_5_mean_sd_ABZ_outliers_removed %>% 
  dplyr::filter(set == "2")

ABZ_set3 <- sets1_5_mean_sd_ABZ_outliers_removed %>% 
  dplyr::filter(set == "3")

ABZ_set4 <- sets1_5_mean_sd_ABZ_outliers_removed %>% 
  dplyr::filter(set == "4")

ABZ_set5 <- sets1_5_mean_sd_ABZ_outliers_removed %>% 
  dplyr::filter(set == "5")

# Perform ANOVA - Set 1 
daily_fecundity_ABZ_AOV_set1 <- ABZ_set1 %>%
  dplyr::filter(condition == "albendazole")%>%
  aov(avg_total_offspring ~ strain, data = .)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(group1 == 'N2' | group2 == 'N2')

View(ABZ_set1)
View(daily_fecundity_ABZ_AOV_set1)

# Flip to always have N2 represented in the group 2 column
daily_fecundity_ABZ_AOV_set1 <- daily_fecundity_ABZ_AOV_set1 %>%
  mutate(
    temp = group1,
    group1 = ifelse(group1 == 'N2', group2, group1),
    group2 = ifelse(temp == 'N2', temp, group2)
  ) %>%
  dplyr::select(-temp)


daily_fecundity_ABZ_AOV_set1 <- daily_fecundity_ABZ_AOV_set1 %>%
  mutate(
    group2 = case_when(
      group2 %in% c("PHX2481", "PHX2502", "PHX2647") ~ "N2",
      TRUE ~ group2))

daily_fecundity_ABZ_AOV_set1  <- daily_fecundity_ABZ_AOV_set1 %>% 
  dplyr::mutate(strain= group1)

View(daily_fecundity_ABZ_AOV_set1)

# Perform ANOVA - Set 2 
daily_fecundity_ABZ_AOV_set2 <- ABZ_set2 %>%
  dplyr::filter(condition == "albendazole")%>%
  aov(avg_total_offspring ~ strain, data = .)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(group1 == 'N2' | group2 == 'N2')

# Flip to always have N2 represented in the group 2 column
daily_fecundity_ABZ_AOV_set2 <- daily_fecundity_ABZ_AOV_set2 %>%
  mutate(
    temp = group1,
    group1 = ifelse(group1 == 'N2', group2, group1),
    group2 = ifelse(temp == 'N2', temp, group2)
  ) %>%
  dplyr::select(-temp)


daily_fecundity_ABZ_AOV_set2 <- daily_fecundity_ABZ_AOV_set2 %>%
  mutate(
    group2 = case_when(
      group2 %in% c("PHX2481", "PHX2502", "PHX2647") ~ "N2",
      TRUE ~ group2))

daily_fecundity_ABZ_AOV_set2  <- daily_fecundity_ABZ_AOV_set2 %>% 
  dplyr::mutate(strain= group1)
View(daily_fecundity_ABZ_AOV_set2)

# Perform ANOVA - Set 3 
daily_fecundity_ABZ_AOV_set3 <- ABZ_set3 %>%
  dplyr::filter(condition == "albendazole")%>%
  aov(avg_total_offspring ~ strain, data = .)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(group1 == 'N2' | group2 == 'N2')

# Flip to always have N2 represented in the group 2 column
daily_fecundity_ABZ_AOV_set3 <- daily_fecundity_ABZ_AOV_set3 %>%
  mutate(
    temp = group1,
    group1 = ifelse(group1 == 'N2', group2, group1),
    group2 = ifelse(temp == 'N2', temp, group2)
  ) %>%
  dplyr::select(-temp)


daily_fecundity_ABZ_AOV_set3 <- daily_fecundity_ABZ_AOV_set3 %>%
  mutate(
    group2 = case_when(
      group2 %in% c("PHX2481", "PHX2502", "PHX2647") ~ "N2",
      TRUE ~ group2))

daily_fecundity_ABZ_AOV_set3  <- daily_fecundity_ABZ_AOV_set3 %>% 
  dplyr::mutate(strain= group1)
View(daily_fecundity_ABZ_AOV_set3)

# Perform ANOVA - Set 4 
daily_fecundity_ABZ_AOV_set4 <- ABZ_set4 %>%
  dplyr::filter(condition == "albendazole")%>%
  aov(avg_total_offspring ~ strain, data = .)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(group1 == 'N2' | group2 == 'N2')

# Flip to always have N2 represented in the group 2 column
daily_fecundity_ABZ_AOV_set4 <- daily_fecundity_ABZ_AOV_set4 %>%
  mutate(
    temp = group1,
    group1 = ifelse(group1 == 'N2', group2, group1),
    group2 = ifelse(temp == 'N2', temp, group2)
  ) %>%
  dplyr::select(-temp)


daily_fecundity_ABZ_AOV_set4 <- daily_fecundity_ABZ_AOV_set4 %>%
  mutate(
    group2 = case_when(
      group2 %in% c("PHX2481", "PHX2502", "PHX2647") ~ "N2",
      TRUE ~ group2))

daily_fecundity_ABZ_AOV_set4  <- daily_fecundity_ABZ_AOV_set4 %>% 
  dplyr::mutate(strain= group1)
View(daily_fecundity_ABZ_AOV_set4)

# Perform ANOVA - Set 5 
daily_fecundity_ABZ_AOV_set5 <- ABZ_set5 %>%
  dplyr::filter(condition == "albendazole")%>%
  aov(avg_total_offspring ~ strain, data = .)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(group1 == 'N2' | group2 == 'N2')

# Flip to always have N2 represented in the group 2 column
daily_fecundity_ABZ_AOV_set5 <- daily_fecundity_ABZ_AOV_set5 %>%
  mutate(
    temp = group1,
    group1 = ifelse(group1 == 'N2', group2, group1),
    group2 = ifelse(temp == 'N2', temp, group2)
  ) %>%
  dplyr::select(-temp)


daily_fecundity_ABZ_AOV_set5 <- daily_fecundity_ABZ_AOV_set5 %>%
  mutate(
    group2 = case_when(
      group2 %in% c("PHX2481", "PHX2502", "PHX2647") ~ "N2",
      TRUE ~ group2))

daily_fecundity_ABZ_AOV_set5  <- daily_fecundity_ABZ_AOV_set5 %>% 
  dplyr::mutate(strain= group1)
View(daily_fecundity_ABZ_AOV_set5)



#clean up dataframes 
write.csv(ABZ_set1, "data/fecundity/ABZ_set1.csv", row.names = FALSE, col.names = FALSE)
write.csv(ABZ_set2, "data/fecundity/ABZ_set2.csv", row.names = FALSE, col.names = FALSE)
write.csv(ABZ_set3, "data/fecundity/ABZ_set3.csv", row.names = FALSE, col.names = FALSE)
write.csv(ABZ_set4, "data/fecundity/ABZ_set4.csv", row.names = FALSE, col.names = FALSE)
write.csv(ABZ_set5, "data/fecundity/ABZ_set5.csv", row.names = FALSE, col.names = FALSE)
write.csv(daily_fecundity_ABZ_AOV_set1, "data/fecundity/daily_fecundity_ABZ_AOV_set1.csv", row.names = FALSE, col.names = FALSE)
write.csv(daily_fecundity_ABZ_AOV_set2, "data/fecundity/daily_fecundity_ABZ_AOV_set2.csv", row.names = FALSE, col.names = FALSE)
write.csv(daily_fecundity_ABZ_AOV_set3, "data/fecundity/daily_fecundity_ABZ_AOV_set3.csv", row.names = FALSE, col.names = FALSE)
write.csv(daily_fecundity_ABZ_AOV_set4, "data/fecundity/daily_fecundity_ABZ_AOV_set4.csv", row.names = FALSE, col.names = FALSE)
write.csv(daily_fecundity_ABZ_AOV_set5, "data/fecundity/daily_fecundity_ABZ_AOV_set5.csv", row.names = FALSE, col.names = FALSE)


merged_daily_ABZ <- read.csv("data/fecundity/ABZ_merged_daily.csv") #UPLOAD THIS FILE TO GITHUB
View(merged_daily_ABZ)

# calculate intrinsic growth rate (ABZ)
ABZ_merged_daily <- mutate(ABZ_merged_daily, #create age column with appropriate ages
                             age = case_when(
                               set == 1 ~ 3,
                               set == 2 ~ 4,
                               set == 3 ~ 5,
                               set == 4 ~ 6,
                               set == 5 ~ 7,
                               TRUE ~ NA_integer_
                             ))

ABZ_merged_daily <- mutate(ABZ_merged_daily,
                           intrinsic_growth_rate = log(avg_total_offspring/ age))



#Set custom order for how strains should be displayed in boxplot 
custom_order <- c("N2", "ECA882", "PHX2647", "PHX2502", "PHX2481", "ECA2687", "ECA2688", "ECA2689", "ECA2690")


# Create the box plot with jittered points and set boxes
ABZ_daily <- ggplot(ABZ_merged_daily, aes(x = factor(strain, levels = custom_order), y = avg_total_offspring, fill = factor(strain, levels = custom_order))) +
  geom_boxplot(width = 0.5, color = "black", outlier.shape = NA) +
  geom_jitter(shape = 21, position = position_jitter(0.2), size = 0.5, alpha = 0.5, fill = "black") +
  labs(x = NULL, y = "Fecundity") +
  scale_fill_manual(values = col_comp, breaks = custom_order, labels = strain_labels) +
  stat_pvalue_manual(ABZ_merged_daily, label = "p.adj.signif", y.position = c(140), xmax = "group1", remove.bracket = TRUE) +
  theme_bw() +   
  facet_grid(. ~ set, switch = "y") +  # Add switch = "y" to remove labels at the top
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        # axis.title = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.line = element_line(color = "black"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        legend.text = element_text(size = 12, color = "black"),
        # legend.title = element_text(size = 12, color = "black"),
        strip.text.x = element_blank(),  # Remove facet labels
        strip.background = element_blank(), 
        legend.position = "none",
        legend.box.background = element_blank(),
        legend.key.size = unit(1.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(10, 10, 20, 10, "pt")) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 25), oob = scales::squish, sec.axis = dup_axis(name = "Albendazole"))

# Print the plot
print(ABZ_daily)

ABZ_daily_individual <- ggplot(ABZ_merged_daily, aes(x = factor(strain, levels = custom_order), y = avg_total_offspring, fill = factor(strain, levels = custom_order))) +
  geom_boxplot(width = 0.5, color = "black", outlier.shape = NA) +
  geom_jitter(shape = 21, position = position_jitter(0.2), size = 0.5, alpha = 0.5, fill = "black") +
  labs(x = "Strain", y = "Fecundity") +
  scale_fill_manual(values = col_comp, breaks = custom_order, labels = strain_labels) +
  stat_pvalue_manual(ABZ_merged_daily, label = "p.adj.signif", y.position = c(150), xmax = "group1", remove.bracket = TRUE) +
  theme_bw() +   
  facet_wrap(.~ set, ncol = 3) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.line = element_line(color = "black"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        legend.text = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 12, vjust = 1, color = "black"),
        strip.background = element_blank(), 
        legend.position = "right",
        legend.box.background = element_blank(),
        legend.key.size = unit(1.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(10, 10, 20, 10, "pt")) +
  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, 50), oob = scales::squish, sec.axis = dup_axis(name = "Albendazole"))

legend_plot <- ggplot(ABZ_merged_daily) +
  geom_point(aes(x = strain, y = 1, color = factor(col_comp, levels = custom_order)), size = 4) +
  theme_void() +
  theme(legend.position = "none")

legend_plot

# Print the plot
print(ABZ_daily_individual)

ggsave(filename = "figures/ABZ_daily_individual_20231221_3.svg",plot = ABZ_daily_individual ,width = 10,height = 6,units = "in")



#############################################
#           Figure X                        #
#  Daily Fecundity in IVM (Sets 1-5)        #
#############################################

#IVM ANOVA
IVM_set1 <- sets1_5_mean_sd_IVM_outliers_removed %>% 
  dplyr::filter(set == "1")
View(IVM_set1)

IVM_set2 <- sets1_5_mean_sd_IVM_outliers_removed %>% 
  dplyr::filter(set == "2")

IVM_set3 <- sets1_5_mean_sd_IVM_outliers_removed %>% 
  dplyr::filter(set == "3")

IVM_set4 <- sets1_5_mean_sd_IVM_outliers_removed %>% 
  dplyr::filter(set == "4")

IVM_set5 <- sets1_5_mean_sd_IVM_outliers_removed %>% 
  dplyr::filter(set == "5")

# Perform ANOVA - Set 1 
daily_fecundity_IVM_AOV_set1 <- IVM_set1 %>%
  dplyr::filter(condition == "ivermectin")%>%
  aov(avg_total_offspring ~ strain, data = .)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(group1 == 'N2' | group2 == 'N2')

View(IVM_set1)
View(daily_fecundity_IVM_AOV_set1)

# Flip to always have N2 represented in the group 2 column
daily_fecundity_IVM_AOV_set1 <- daily_fecundity_IVM_AOV_set1 %>%
  mutate(
    temp = group1,
    group1 = ifelse(group1 == 'N2', group2, group1),
    group2 = ifelse(temp == 'N2', temp, group2)
  ) %>%
  dplyr::select(-temp)


daily_fecundity_IVM_AOV_set1 <- daily_fecundity_IVM_AOV_set1 %>%
  mutate(
    group2 = case_when(
      group2 %in% c("PHX2481", "PHX2502", "PHX2647") ~ "N2",
      TRUE ~ group2))

daily_fecundity_IVM_AOV_set1  <- daily_fecundity_IVM_AOV_set1 %>% 
  dplyr::mutate(strain= group1)

View(daily_fecundity_IVM_AOV_set1)

# Perform ANOVA - Set 2 
daily_fecundity_IVM_AOV_set2 <- IVM_set2 %>%
  dplyr::filter(condition == "ivermectin")%>%
  aov(avg_total_offspring ~ strain, data = .)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(group1 == 'N2' | group2 == 'N2')

# Flip to always have N2 represented in the group 2 column
daily_fecundity_IVM_AOV_set2 <- daily_fecundity_IVM_AOV_set2 %>%
  mutate(
    temp = group1,
    group1 = ifelse(group1 == 'N2', group2, group1),
    group2 = ifelse(temp == 'N2', temp, group2)
  ) %>%
  dplyr::select(-temp)


daily_fecundity_IVM_AOV_set2 <- daily_fecundity_IVM_AOV_set2 %>%
  mutate(
    group2 = case_when(
      group2 %in% c("PHX2481", "PHX2502", "PHX2647") ~ "N2",
      TRUE ~ group2))

daily_fecundity_IVM_AOV_set2  <- daily_fecundity_IVM_AOV_set2 %>% 
  dplyr::mutate(strain= group1)
View(daily_fecundity_IVM_AOV_set2)

# Perform ANOVA - Set 3 
daily_fecundity_IVM_AOV_set3 <- IVM_set3 %>%
  dplyr::filter(condition == "ivermectin")%>%
  aov(avg_total_offspring ~ strain, data = .)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(group1 == 'N2' | group2 == 'N2')

# Flip to always have N2 represented in the group 2 column
daily_fecundity_IVM_AOV_set3 <- daily_fecundity_IVM_AOV_set3 %>%
  mutate(
    temp = group1,
    group1 = ifelse(group1 == 'N2', group2, group1),
    group2 = ifelse(temp == 'N2', temp, group2)
  ) %>%
  dplyr::select(-temp)


daily_fecundity_IVM_AOV_set3 <- daily_fecundity_IVM_AOV_set3 %>%
  mutate(
    group2 = case_when(
      group2 %in% c("PHX2481", "PHX2502", "PHX2647") ~ "N2",
      TRUE ~ group2))

daily_fecundity_IVM_AOV_set3  <- daily_fecundity_IVM_AOV_set3 %>% 
  dplyr::mutate(strain= group1)
View(daily_fecundity_IVM_AOV_set3)

# Perform ANOVA - Set 4 
daily_fecundity_IVM_AOV_set4 <- IVM_set4 %>%
  dplyr::filter(condition == "ivermectin")%>%
  aov(avg_total_offspring ~ strain, data = .)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(group1 == 'N2' | group2 == 'N2')

# Flip to always have N2 represented in the group 2 column
daily_fecundity_IVM_AOV_set4 <- daily_fecundity_IVM_AOV_set4 %>%
  mutate(
    temp = group1,
    group1 = ifelse(group1 == 'N2', group2, group1),
    group2 = ifelse(temp == 'N2', temp, group2)
  ) %>%
  dplyr::select(-temp)


daily_fecundity_IVM_AOV_set4 <- daily_fecundity_IVM_AOV_set4 %>%
  mutate(
    group2 = case_when(
      group2 %in% c("PHX2481", "PHX2502", "PHX2647") ~ "N2",
      TRUE ~ group2))

daily_fecundity_IVM_AOV_set4  <- daily_fecundity_IVM_AOV_set4 %>% 
  dplyr::mutate(strain= group1)
View(daily_fecundity_IVM_AOV_set4)

# Perform ANOVA - Set 5 
daily_fecundity_IVM_AOV_set5 <- IVM_set5 %>%
  dplyr::filter(condition == "ivermectin")%>%
  aov(avg_total_offspring ~ strain, data = .)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(group1 == 'N2' | group2 == 'N2')

# Flip to always have N2 represented in the group 2 column
daily_fecundity_IVM_AOV_set5 <- daily_fecundity_IVM_AOV_set5 %>%
  mutate(
    temp = group1,
    group1 = ifelse(group1 == 'N2', group2, group1),
    group2 = ifelse(temp == 'N2', temp, group2)
  ) %>%
  dplyr::select(-temp)


daily_fecundity_IVM_AOV_set5 <- daily_fecundity_IVM_AOV_set5 %>%
  mutate(
    group2 = case_when(
      group2 %in% c("PHX2481", "PHX2502", "PHX2647") ~ "N2",
      TRUE ~ group2))

daily_fecundity_IVM_AOV_set5  <- daily_fecundity_IVM_AOV_set5 %>% 
  dplyr::mutate(strain= group1)
View(daily_fecundity_IVM_AOV_set5)



#clean up dataframes 
write.csv(IVM_set1, "data/fecundity/IVM_set1.csv", row.names = FALSE, col.names = FALSE)
write.csv(IVM_set2, "data/fecundity/IVM_set2.csv", row.names = FALSE, col.names = FALSE)
write.csv(IVM_set3, "data/fecundity/IVM_set3.csv", row.names = FALSE, col.names = FALSE)
write.csv(IVM_set4, "data/fecundity/IVM_set4.csv", row.names = FALSE, col.names = FALSE)
write.csv(IVM_set5, "data/fecundity/IVM_set5.csv", row.names = FALSE, col.names = FALSE)
write.csv(daily_fecundity_IVM_AOV_set1, "data/fecundity/daily_fecundity_IVM_AOV_set1.csv", row.names = FALSE, col.names = FALSE)
write.csv(daily_fecundity_IVM_AOV_set2, "data/fecundity/daily_fecundity_IVM_AOV_set2.csv", row.names = FALSE, col.names = FALSE)
write.csv(daily_fecundity_IVM_AOV_set3, "data/fecundity/daily_fecundity_IVM_AOV_set3.csv", row.names = FALSE, col.names = FALSE)
write.csv(daily_fecundity_IVM_AOV_set4, "data/fecundity/daily_fecundity_IVM_AOV_set4.csv", row.names = FALSE, col.names = FALSE)
write.csv(daily_fecundity_IVM_AOV_set5, "data/fecundity/daily_fecundity_IVM_AOV_set5.csv", row.names = FALSE, col.names = FALSE)



merged_daily_IVM <- read.csv("data/fecundity/IVM_combo.csv") #UPLOAD THIS FILE TO GITHUB
View(merged_daily_IVM)

merged_daily_IVM  <- read.csv("data/fecundity/IVM_combo.csv")

#Set custom order for how strains should be displayed in boxplot 
custom_order <- c("N2", "ECA882", "PHX2647", "PHX2502", "PHX2481", "ECA2687", "ECA2688", "ECA2689", "ECA2690")


# Create the box plot with jittered points and set boxes
IVM_daily <- ggplot(merged_daily_IVM, aes(x = factor(strain, levels = custom_order), y = avg_total_offspring, fill = factor(strain, levels = custom_order))) +
  geom_boxplot(width = 0.5, color = "black", outlier.shape = NA) +
  geom_jitter(shape = 21, position = position_jitter(0.2), size = 0.5, alpha = 0.5, fill = "black") +
  labs(x = "Strain", y = NULL) +
  scale_fill_manual(values = col_comp, breaks = custom_order, labels = strain_labels) +
  stat_pvalue_manual(merged_daily_IVM, label = "p.adj.signif", y.position = c(140), xmax = "group1", remove.bracket = TRUE) +
  theme_bw() +   
  facet_grid(. ~ set) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        # axis.title = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.line = element_line(color = "black"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        legend.text = element_text(size = 12, color = "black"),
        # legend.title = element_text(size = 12, color = "black"),
        strip.text.x = element_blank(),  # Remove facet labels
        strip.background = element_blank(), 
        legend.position = "none",
        legend.box.background = element_blank(),
        legend.key.size = unit(1.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(10, 10, 20, 10, "pt")) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 25), oob = scales::squish, sec.axis = dup_axis(name = "Ivermectin"))

# Print the plot
print(IVM_daily)

IVM_daily_individual <- ggplot(merged_daily_IVM, aes(x = factor(strain, levels = custom_order), y = avg_total_offspring, fill = factor(strain, levels = custom_order))) +
  geom_boxplot(width = 0.5, color = "black", outlier.shape = NA) +
  geom_jitter(shape = 21, position = position_jitter(0.2), size = 0.5, alpha = 0.5, fill = "black") +
  labs(x = "Strain", y = "Fecundity") +
  scale_fill_manual(values = col_comp, breaks = custom_order, labels = strain_labels) +
  stat_pvalue_manual(merged_daily_IVM, label = "p.adj.signif", y.position = c(140), xmax = "group1", remove.bracket = TRUE) +
  theme_bw() +   
  facet_wrap(.~ set, ncol = 3) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.line = element_line(color = "black"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        legend.text = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 12, vjust = 1, color = "black"),
        strip.background = element_blank(), 
        legend.position = "none",
        legend.box.background = element_blank(),
        legend.key.size = unit(1.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(10, 10, 20, 10, "pt")) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 25), oob = scales::squish, sec.axis = dup_axis(name = "Ivermectin"))

legend_plot <- ggplot(merged_daily) +
  geom_point(aes(x = strain, y = 1, color = factor(col_comp, levels = custom_order)), size = 4) +
  theme_void() +
  theme(legend.position = "none")

legend_plot

# Print the plot
print(IVM_daily_individual)

ggsave(filename = "figures/IVM_daily_individual_20231221.svg",plot = IVM_daily_individual ,width = 10,height = 6,units = "in")










#Combine all plots for lifetimefecundity figure
comb_daily_fecundity <- cowplot::plot_grid(DMSO_daily, 
                                             ABZ_daily, 
                                             IVM_daily,
                                             ncol = 1,
                                             nrow = 3,
                                             labels = c("A","B","C"),
                                             align = "vh",
                                             axis = "lrbt",
                                             label_size = 12,
                                             label_fontfamily = "Arial",
                                             rel_widths = c(1,1,1),
                                             rel_heights = c(1,1,1))

comb_daily_fecundity 

ggsave(filename = "figures/comb_daily_fecundity_20231220.jpeg",plot = comb_daily_fecundity ,width = 12,height = 10,units = "in")



















DMSO_daily <- ggplot(combined_daily_allsets, aes(x = factor(strain, levels = custom_order), y = avg_total_offspring, fill = factor(strain, levels = custom_order))) +
  geom_boxplot(width = 0.5, color = "black", outlier.shape = NA) +
  geom_jitter(shape = 21, position = position_jitter(0.2), size = 0.5, alpha = 0.5, fill = "black") +
  labs(x = NULL, y = "Fecundity") +
  scale_fill_manual(values = col_comp, breaks = custom_order, labels = strain_labels) +
  # scale_x_discrete(labels = strain_labels) +  # Set custom strain labels 
  theme_bw() +   facet_grid(. ~ set) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.line = element_line(color = "black"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 12, vjust = 1, color = "black"),
        strip.background = element_blank(), 
        legend.position = "none",
        legend.box.background = element_blank(),
        legend.key.size = unit(1.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(10, 10, 20, 10, "pt")) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 25), oob = scales::squish, sec.axis = dup_axis(name = "DMSO")+
                       geom_signif(comparisons = combined_daily_allsets, annotations = "p.adj.signif", y_position = 140, step_increase = 20))  # Add adjusted p-values






# Add adjusted p-values using geom_signif
DMSO_daily <- DMSO_daily +
  geom_signif(comparisons = your_comparisons_data_frame, annotations = "p_values", y_position = 140, step_increase = 20)













                    

  

DMSO_daily
ggsave("figures/fecundity/DMSO_daily_fecundity_fulldataset_20231218.png", p, width = 10, height = 6, dpi = 300)


# #Testing cumulative density function plots 
# 
# DMSO_daily <- ggplot(sets1_5_mean_sd_DMSO_outliers_removed, aes(x = factor(set), y = avg_total_offspring, fill = factor(strain))) + 
#   stat_ecdf(geom = "point") +
#   labs(y = "Fecundity", x = "set")+
#   theme_classic()
# 
# DMSO_daily











#### old 

















##########################################
#              Figure 1B                 #
#      Daily fecundity in DMSO           #
##########################################
# mean_sd_DMSO_filtered <- mean_sd_DMSO_filtered %>% 
#   select(strain, set, avg_total_offspring, mean, sd, lifetime_mean)

# stat_box_data <- function(y) {
#   upper_limit <- max(na.omit(y)) * 1.15
#   return( 
#     data.frame(
#       y = 0.98 * upper_limit,
#       label = paste(length(y), '\n')))}



# Create the violin plot with jittered points
# p <- ggplot(mean_sd_DMSO_filtered, aes(x = reorder_within(strain, avg_total_offspring, custom_order), y = avg_total_offspring, fill = strain)) +
#   facet_grid(. ~ set) +
#   geom_violin(scale = "width", trim = FALSE) +
#   geom_jitter(position = position_jitter(0.2), alpha = 0.5, color = "black") +
#   geom_boxplot(width = 0.1, fill = "white", color = "black") +
#   labs(x = "Strain", y = "Fecundity") +
#   scale_fill_manual(values = col_comp, breaks = custom_order) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         axis.title = element_text(size = 12, color = "black"),
#         axis.text.y = element_text(size = 12, color = "black"),
#         axis.line = element_line(color = "black"),
#         legend.text = element_text(size = 12, color = "black"),
#         legend.title = element_text(size = 12, color = "black"),
#         strip.text = element_text(size = 12, vjust = 1, color = "black"),
#         strip.background = element_blank(),
#         legend.position = "right",
#         legend.box.background = element_blank(),
#         legend.key.size = unit(1.5, "lines"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         plot.margin = margin(10, 10, 20, 10, "pt")) +
#   scale_y_continuous(limits = c(0, 175), breaks = seq(0, 175, 25), oob = scales::squish)
# 

# Add significance annotations manually using custom symbols
# comparison_annotations <- map_p_values(anova_p_values)
# annotation_data <- data.frame(set = unique(mean_sd_DMSO$set), comparison_annotations)
# 
# p + geom_signif(comparisons = list(c("N2", custom_order[-1])),
#                 aes(y = 160, ymax = 160, annotations = comparison_annotations),
#                 data = annotation_data,
#                 map_signif_level = FALSE,
#                 textsize = 3, vjust = 0.5, step_increase = 0.05)
  

##########################################
#              Figure XX                 #
# Daily fecundity in DMSO - Full dataset #
##########################################              
#Works! 
# Create the box plot with jittered points and set boxes
p <- ggplot(mean_sd_DMSO, aes(x = factor(strain, levels = custom_order), y = avg_total_offspring, fill = factor(strain, levels = custom_order))) +
  facet_grid(. ~ set) + 
  geom_boxplot(width = 0.3, color = "black", outlier.shape = NA) +
  geom_jitter(position = position_jitter(0.2), alpha = 0.5, color = "black") +
  labs(x = "Strain", y = "Fecundity") +
  scale_fill_manual(values = col_comp, breaks = custom_order, labels = strain_labels) +
  scale_x_discrete(labels = strain_labels) +  # Set custom strain labels 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 12, vjust = 1, color = "black"),
        strip.background = element_blank(), 
        legend.position = "none",
        legend.box.background = element_blank(),
        legend.key.size = unit(1.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(10, 10, 20, 10, "pt")) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 25), oob = scales::squish) +
  ggtitle("DMSO")  # Add the title "DMSO"

p
ggsave("figures/fecundity/DMSO_daily_fecundity_fulldataset_20231218.png", p, width = 10, height = 6, dpi = 300)

##########################################
#              Figure XX                 #
# Daily fecundity in DMSO - Sets 1 - 5   #
##########################################              
#Works! 
# Create the box plot with jittered points and set boxes
D_sets1_5 <- ggplot(sets1_5_mean_sd_DMSO, aes(x = factor(strain, levels = custom_order), y = avg_total_offspring, fill = factor(strain, levels = custom_order))) +
  facet_grid(. ~ set) + 
  geom_boxplot(width = 0.3, color = "black", outlier.shape = NA) +
  geom_jitter(position = position_jitter(0.2), alpha = 0.5, color = "black") +
  labs(x = "Strain", y = "Fecundity") +
  scale_fill_manual(values = col_comp, breaks = custom_order, labels = strain_labels) +
  scale_x_discrete(labels = strain_labels) +  # Set custom strain labels 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 12, vjust = 1, color = "black"),
        strip.background = element_blank(), 
        legend.position = "none",
        legend.box.background = element_blank(),
        legend.key.size = unit(1.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(10, 10, 20, 10, "pt")) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 25), oob = scales::squish) +
  ggtitle("DMSO - Sets 1 - 5") 

D_sets1_5
ggsave("figures/fecundity/DMSO_daily_fecundity_sets1_5_20231218.png", D_sets1_5, width = 10, height = 6, dpi = 300)


##########################################
#              Figure XX                 #
# Daily fecundity in ABZ - Full dataset  #
##########################################

A <- ggplot(mean_sd_ABZ, aes(x = factor(strain, levels = custom_order), y = avg_total_offspring, fill = factor(strain, levels = custom_order))) +
  facet_grid(. ~ set) + 
  geom_boxplot(width = 0.3, color = "black", outlier.shape = NA) +
  geom_jitter(position = position_jitter(0.2), alpha = 0.5, color = "black") +
  labs(x = "Strain", y = "Fecundity") +
  scale_fill_manual(values = col_comp, breaks = custom_order, labels = strain_labels) +
  scale_x_discrete(labels = strain_labels) +  # Set custom strain labels 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 12, vjust = 1, color = "black"),
        strip.background = element_blank(), 
        legend.position = "none",
        legend.box.background = element_blank(),
        legend.key.size = unit(1.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(10, 10, 20, 10, "pt")) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 25), oob = scales::squish) +
  ggtitle("Albendazole - Full Dataset")  # Add the title "DMSO"

A
ggsave("figures/fecundity/ABZ_daily_fecundity_fulldataset_20231218.png", A, width = 10, height = 6, dpi = 300)

##########################################
#              Figure XX                 #
# Daily fecundity in ABZ - Sets 1 - 5    #
##########################################

A_sets1_5 <- ggplot(sets1_5_mean_sd_ABZ, aes(x = factor(strain, levels = custom_order), y = avg_total_offspring, fill = factor(strain, levels = custom_order))) +
  facet_grid(. ~ set) + 
  geom_boxplot(width = 0.3, color = "black", outlier.shape = NA) +
  geom_jitter(position = position_jitter(0.2), alpha = 0.5, color = "black") +
  labs(x = "Strain", y = "Fecundity") +
  scale_fill_manual(values = col_comp, breaks = custom_order, labels = strain_labels) +
  scale_x_discrete(labels = strain_labels) +  # Set custom strain labels 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 12, vjust = 1, color = "black"),
        strip.background = element_blank(), 
        legend.position = "none",
        legend.box.background = element_blank(),
        legend.key.size = unit(1.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(10, 10, 20, 10, "pt")) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 25), oob = scales::squish) +
  ggtitle("Albendazole - Sets 1 - 5")  # Add the title "DMSO"

A_sets1_5
ggsave("figures/fecundity/ABZ_daily_fecundity_sets_1_5_20231218.png", A_sets1_5, width = 10, height = 6, dpi = 300)


##########################################
#              Figure XX                 #
# Daily fecundity in IVM - Full dataset  #
##########################################
I <- ggplot(mean_sd_IVM, aes(x = factor(strain, levels = custom_order), y = avg_total_offspring, fill = factor(strain, levels = custom_order))) +
  facet_grid(. ~ set) + 
  geom_boxplot(width = 0.3, color = "black", outlier.shape = NA) +
  geom_jitter(position = position_jitter(0.2), alpha = 0.5, color = "black") +
  labs(x = "Strain", y = "Fecundity") +
  scale_fill_manual(values = col_comp, breaks = custom_order, labels = strain_labels) +
  scale_x_discrete(labels = strain_labels) +  # Set custom strain labels 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 12, vjust = 1, color = "black"),
        strip.background = element_blank(), 
        legend.position = "none",
        legend.box.background = element_blank(),
        legend.key.size = unit(1.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(10, 10, 20, 10, "pt")) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 25), oob = scales::squish) +
  ggtitle("Ivermectin - Full dataset")  

I  
ggsave("figures/fecundity/IVM_daily_fecundity_fulldataset_20231218.png", I, width = 10, height = 6, dpi = 300)

##########################################
#              Figure XX                 #
# Daily fecundity in IVM - Sets 1 - 5    #
##########################################
I_sets1_5 <- ggplot(sets1_5_mean_sd_IVM, aes(x = factor(strain, levels = custom_order), y = avg_total_offspring, fill = factor(strain, levels = custom_order))) +
  facet_grid(. ~ set) + 
  geom_boxplot(width = 0.3, color = "black", outlier.shape = NA) +
  geom_jitter(position = position_jitter(0.2), alpha = 0.5, color = "black") +
  labs(x = "Strain", y = "Fecundity") +
  scale_fill_manual(values = col_comp, breaks = custom_order, labels = strain_labels) +
  scale_x_discrete(labels = strain_labels) +  # Set custom strain labels 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 12, vjust = 1, color = "black"),
        strip.background = element_blank(), 
        legend.position = "none",
        legend.box.background = element_blank(),
        legend.key.size = unit(1.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(10, 10, 20, 10, "pt")) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 25), oob = scales::squish) +
  ggtitle("Ivermectin - Sets 1 - 5")  

I_sets1_5
ggsave("figures/fecundity/IVM_daily_fecundity_sets1_5_20231218.png", I_sets1_5, width = 10, height = 6, dpi = 300)



















        
        
        
    

    















IVM_lifetime <- mean_sd_IVM_TOTAL %>%
  ggplot(aes(x = strain, y = total_mean, fill = strain)) +
  geom_bar(stat = 'identity', color = 'gray51', width = 0.7, position = position_dodge(0.9)) +
  labs(x = "Strain", y = "Lifetime Fecundity", fill = "Strain") +
  ggtitle("Ivermectin")+
  scale_fill_manual(values = col_comp, breaks = names(col_comp), labels = strain_labels) +
  scale_x_discrete(labels = strain_labels) + #set custom strain labels 
  scale_y_continuous(expand = c(0, 0)) +
  theme_minimal() +

  # Assuming you have a data frame named "mean_sd_DMSO" with columns "strain", "total_mean", "set", and "Brood_size"
  
  # Remove rows with missing values but retain rows with zeros
  mean_sd_DMSO <- mean_sd_DMSO[complete.cases(mean_sd_DMSO), ]

# Create the plot using ggplot
fig2b_r <- ggplot(mean_sd_DMSO, aes(x = strain, y = total_mean, fill = strain)) +
  # facet_grid(. ~ set) +
  geom_jitter(shape = 21, position = position_jitter(0.2), size = 0.5, fill = "black") +
  geom_boxplot(aes(fill = strain), outlier.shape = NA, alpha = 0.5) +
  scale_fill_manual(values = col_comp) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 12, vjust = 1, color = "black"),
        strip.background = element_blank(),
        legend.position = "none",
        panel.grid = ggplot2::element_blank(),
        legend.spacing.x = unit(1, 'mm'),
        legend.background = element_rect(fill = alpha('white', 0)),
        text = element_text(family = "Helvetica")) +
  ylab("Brood size")

# Print or continue with further operations on fig2b_r
print(fig2b_r)

  # scale_y_continuous(limits =c(0,210),breaks=seq(0,50,100,150,200), oob = scales::squish) #+
  # stat_compare_means(comparisons = list(c("swept","divergent")),label.y = c(201),label = "p.signif") +
  # labs(color="Genotype",fill="Genotype")





mean_sd_DMSO <- mean_sd_DMSO %>% ungroup()

ggplot(mean_sd_DMSO, aes(x = strain, y = total_mean, fill = strain)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +  # Box plot without outliers and slightly opaque
  geom_jitter(shape = 21, position = position_jitter(0.2), size = 2, alpha = 0.7) +  # Jitter points with fill "black" and slightly opaque
  facet_grid(. ~ set) +  # Facet by the "set" variable
  scale_fill_manual(values = col_comp) +  # Set the colors for each strain
  theme_bw() +  # Use a white background theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 12, vjust = 1, color = "black"),
        strip.background = element_blank(),
        legend.position = "none",  # Remove the legend
        panel.grid = ggplot2::element_blank(),
        legend.spacing.x = unit(1, 'mm'),
        legend.background = element_rect(fill = alpha('white', 0)),
        text = element_text(family = "Helvetica")) +  
  ylab("Brood size") + 
  scale_y_continuous(limits = c(0, 210), breaks = seq(0, 210, 50), oob = scales::squish)

### cow fig2 



fig2b <- cowplot::plot_grid(fig2b_l,fig2b_r,
                            nrow = 1,
                            ncol = 2,
                            label_size = 12, 
                            label_fontfamily="Helvetica",
                            axis = "tb",
                            rel_widths = c(1.5,5))



fig_2 <- cowplot::plot_grid(fig2a, fig2b,
                            labels = c("A","B"), 
                            nrow = 2,
                            ncol = 1,
                            label_size = 12, 
                            label_fontfamily="Helvetica",
                            axis = "lr",
                            rel_heights =  c(1.1,2))

ggsave(fig_2, filename = paste( "figures/Fig_2.pdf",sep = ""), device = cairo_pdf, units = "mm",height = 140, width = 170)

##########################################
#              Figure ??                 #
#      Daily Starvation in DMSO          #
##########################################

starvation_alldrugs <- BZML_Fecundity_Counts_Amanda_BZML_fecundity_counts_amanda %>%
  filter(set >= 1 & set <= 5)

#To assess how long it took each strain to starve out in each condition to ensure 
#Drugs were selective 

#Filter to condition - DMSO 
starvation_DMSO <- starvation_alldrugs %>% 
  dplyr::filter(condition == "DMSO")

#Make new df w/ counts column and count the number of each strain in each set 
counts_starv_DMSO <- starvation_DMSO %>%
  count(strain, days_to_starvation_edited, name = "counts")

# Convert the "strain" variable to a factor with the desired order
counts_starv_DMSO$strain <- factor(counts_starv_DMSO$strain, levels = strain_order)

#Need to make zeros for strains that don't have data on a given dataset 
new_DMSO <- data.frame(
  days_to_starvation_edited = c(12, 13),
  strain = c("ECA2687", "ECA2687"),
  counts = c(0, 0))

ECA2688_DMSO <- data.frame(
  days_to_starvation_edited = c(11, 12, 13),
  strain = c("ECA2688", "ECA2688", "ECA2688"),
  counts = c(0, 0, 0))

ECA2689_DMSO <- data.frame(
  days_to_starvation_edited = c(6, 10, 11, 12, 13),
  strain = c("ECA2689", "ECA2689", "ECA2689", "ECA2689", "ECA2689"),
  counts = c(0, 0, 0, 0, 0))

ECA2690_DMSO <- data.frame(
  days_to_starvation_edited = c(6, 10, 11, 12, 13),
  strain = c("ECA2690", "ECA2690", "ECA2690", "ECA2690", "ECA2690"),
  counts = c(0, 0, 0, 0, 0))

ECA882_DMSO <- data.frame(
  days_to_starvation_edited = c(11, 12, 13),
  strain = c("ECA882", "ECA882", "ECA882"),
  counts = c(0, 0, 0))

N2_DMSO <- data.frame(
  days_to_starvation_edited = c(11, 12, 13),
  strain = c("N2", "N2", "N2"),
  counts = c(0, 0, 0))

PHX2481_DMSO <- data.frame(
  days_to_starvation_edited = c(10, 11, 12, 13),
  strain = c("PHX2481", "PHX2481", "PHX2481", "PHX2481"),
  counts = c(0, 0, 0, 0))

PHX2502_DMSO <- data.frame(
  days_to_starvation_edited = c(6, 10, 11, 12, 13),
  strain = c("PHX2502", "PHX2502", "PHX2502", "PHX2502", "PHX2502"),
  counts = c(0, 0, 0, 0, 0))

PHX2647_DMSO <- data.frame(
  days_to_starvation_edited = c(11, 12, 13),
  strain = c("PHX2647", "PHX2647", "PHX2647"),
  counts = c(0, 0, 0, 0, 0))

#combine all new dfs
combined_starved <- rbind(new_DMSO, ECA2688_DMSO, ECA2689_DMSO,
                          ECA2690_DMSO, ECA882_DMSO, N2_DMSO,
                          PHX2481_DMSO, PHX2502_DMSO, PHX2647_DMSO)

combined_starved$days_to_starvation_edited <- as.numeric(as.character(combined_starved$days_to_starvation_edited))
View(combined_starved)

counts_starv_DMSO$days_to_starvation_edited <- as.numeric(as.character(counts_starv_DMSO$days_to_starvation_edited))
counts_starv_DMSO$strain <- as.character(counts_starv_DMSO$strain)
combined_starved$days_to_starvation <- NULL

#combine new dfs with data
NEW_counts_starv_DMSO <- rbind(counts_starv_DMSO, combined_starved)

#remove NAs
NEW_counts_starv_DMSO <- na.omit(NEW_counts_starv_DMSO)


fig_starv_dmso <- NEW_counts_starv_DMSO %>%
  ggplot(aes(x = factor(days_to_starvation_edited), y = counts, fill = strain)) +
  geom_bar(stat = "identity", position = position_dodge2(width = 0.8), width = 0.7, color = "black", alpha = 0.7) +
  scale_fill_manual(values = col_comp, name = "Strain") +
  theme_minimal() +
  labs(x = "Days to Starvation", y = "Number of Starved Plates", title = "DMSO") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Set the font size for x-axis labels
        axis.title = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.line = element_line(color = "black"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),  # Remove the default panel border
        panel.background = element_blank(),
        strip.text = element_text(size = 12, vjust = 1.5),
        panel.spacing = unit(0.2, "lines")) +  # Add space between facets
  facet_wrap(~ strain, ncol = 3, scales = "free_x") +
  geom_text(aes(label = counts), position = position_dodge(width = 0.7), vjust = -0.5, size = 3) +
  scale_x_discrete(labels = function(x) ifelse(x == "14", "14 +", x), breaks = unique(factor(NEW_counts_starv_DMSO$days_to_starvation_edited)))

fig_starv_dmso
#Save figure 
ggsave(fig_starv_dmso, filename = "figures/Fig_starv_dmso.jpg", bg = "white", units = "in", height = 7, width = 8.5)

##########################################
#              Figure ??                 #
#      Daily Starvation in ABZ           #
##########################################

#Filter to condition - ABZ
starvation_ABZ <- starvation_alldrugs %>% 
  dplyr::filter(condition == "albendazole")

#Make new df w/ counts column and count the number of each strain in each set 
counts_starv_ABZ <- starvation_ABZ %>%
  count(strain, days_to_starvation_edited, name = "counts")

# Convert the "strain" variable to a factor with the desired order
counts_starv_ABZ$strain <- factor(counts_starv_ABZ$strain, levels = strain_order)

#Remove unncessary zeros 
data_filtered_ABZ <- counts_starv_ABZ %>%
  filter(!(counts == 0 & strain == "ECA882" & days_to_starvation_edited == "11"))
data_filtered_ABZ <- counts_starv_ABZ %>%
  filter(!(counts == 0 & strain == "ECA882" & days_to_starvation_edited == "12"))
data_filtered_ABZ <- counts_starv_ABZ %>%
  filter(!(counts == 0 & strain == "ECA882" & days_to_starvation_edited == "13"))
data_filtered_ABZ <- counts_starv_ABZ %>%
  filter(!(counts == 0 & strain == "PHX2481" & days_to_starvation_edited == "10"))
data_filtered_ABZ <- counts_starv_ABZ %>%
  filter(!(counts == 0 & strain == "PHX2481" & days_to_starvation_edited == "11"))
data_filtered_ABZ <- counts_starv_ABZ %>%
  filter(!(counts == 0 & strain == "ECA2688" & days_to_starvation_edited == "11"))
data_filtered_ABZ <- counts_starv_ABZ %>%
  filter(!(counts == 0 & strain == "ECA2689" & days_to_starvation_edited == "10"))
data_filtered_ABZ <- counts_starv_ABZ %>%
  filter(!(counts == 0 & strain == "ECA2690" & days_to_starvation_edited == "10"))
data_filtered_ABZ <- counts_starv_ABZ %>%
  filter(!(counts == 0 & strain == "ECA2690" & days_to_starvation_edited == "11"))
data_filtered_ABZ <- counts_starv_ABZ %>%
  filter(!(counts == 0 & strain == "ECA2690" & days_to_starvation_edited == "12"))

#combine new dfs with data
NEW_counts_starv_ABZ <- rbind(counts_starv_ABZ, combined_starved)

#remove NAs
NEW_counts_starv_ABZ <- na.omit(NEW_counts_starv_ABZ)


filtered_counts_starv_ABZ <- NEW_counts_starv_ABZ %>%
  group_by(strain, days_to_starvation_edited) %>%
  filter(!(counts == 0 & any(counts != 0)))

fig_starv_abz <- filtered_counts_starv_ABZ %>%
  ggplot(aes(x = factor(days_to_starvation_edited), y = counts, fill = strain)) +
  geom_bar(stat = "identity", position = position_dodge2(width = 0.8), width = 0.7, color = "black", alpha = 0.7) +
  scale_fill_manual(values = col_comp, name = "Strain") +
  theme_minimal() +
  labs(x = "Days to Starvation", y = "Number of Starved Plates", title = "Albendazole") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Set the font size for x-axis labels
        axis.title = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.line = element_line(color = "black"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),  # Remove the default panel border
        panel.background = element_blank(),
        strip.text = element_text(size = 12, vjust = 1.5),
        panel.spacing = unit(0.2, "lines")) +  # Add space between facets
  facet_wrap(~ strain, ncol = 3, scales = "free_x") +
  geom_text(aes(label = counts), position = position_dodge(width = 0.7), vjust = -0.5, size = 3) +
  scale_x_discrete(labels = function(x) ifelse(x == "14", "14 +", x), breaks = unique(factor(filtered_counts_starv_ABZ$days_to_starvation_edited)))

fig_starv_abz

#Save figure 
ggsave(fig_starv_abz, filename = "figures/fig_starv_abz.jpg", bg = "white", units = "in", height = 7, width = 8.5)

##########################################
#              Figure ??                 #
#      Daily Starvation in IVM           #
##########################################

#Filter to condition - IVM
starvation_IVM <- starvation_alldrugs %>% 
  dplyr::filter(condition == "ivermectin")

#Make new df w/ counts column and count the number of each strain in each set 
counts_starv_IVM <- starvation_IVM %>%
  count(strain, days_to_starvation_edited, name = "counts")



# Convert the "strain" variable to a factor with the desired order
counts_starv_IVM$strain <- factor(counts_starv_IVM$strain, levels = strain_order)

#Remove unncessary zeros 
data_filtered_IVM <- counts_starv_IVM %>%
  filter(!(counts == 0 & strain == "N2" & days_to_starvation_edited == "11"))
data_filtered_IVM <- counts_starv_IVM %>%
  filter(!(counts == 0 & strain == "ECA882" & days_to_starvation_edited == "11"))
data_filtered_IVM <- counts_starv_IVM %>%
  filter(!(counts == 0 & strain == "PX2647" & days_to_starvation_edited == "11"))
data_filtered_IVM <- counts_starv_IVM %>%
  filter(!(counts == 0 & strain == "PHX2502" & days_to_starvation_edited == "10"))
data_filtered_IVM <- counts_starv_IVM %>%
  filter(!(counts == 0 & strain == "PHX2502" & days_to_starvation_edited == "11"))
data_filtered_IVM <- counts_starv_IVM %>%
  filter(!(counts == 0 & strain == "PHX2502" & days_to_starvation_edited == "12"))
data_filtered_IVM <- counts_starv_IVM %>%
  filter(!(counts == 0 & strain == "PHX2502" & days_to_starvation_edited == "13"))
data_filtered_IVM <- counts_starv_IVM %>%
  filter(!(counts == 0 & strain == "PHX2481" & days_to_starvation_edited == "13"))
data_filtered_IVM <- counts_starv_IVM %>%
  filter(!(counts == 0 & strain == "ECA2689" & days_to_starvation_edited == "10"))
data_filtered_IVM <- counts_starv_IVM %>%
  filter(!(counts == 0 & strain == "ECA2689" & days_to_starvation_edited == "11"))
data_filtered_IVM <- counts_starv_IVM %>%
  filter(!(counts == 0 & strain == "ECA2689" & days_to_starvation_edited == "12"))
data_filtered_IVM <- counts_starv_IVM %>%
  filter(!(counts == 0 & strain == "ECA2689" & days_to_starvation_edited == "13"))
data_filtered_IVM <- counts_starv_IVM %>%
  filter(!(counts == 0 & strain == "ECA2690" & days_to_starvation_edited == "6"))
data_filtered_IVM <- counts_starv_IVM %>%
  filter(!(counts == 0 & strain == "ECA2690" & days_to_starvation_edited == "10"))


#combine new dfs with data
NEW_counts_starv_IVM <- rbind(data_filtered_IVM, combined_starved)

#remove NAs
NEW_counts_starv_IVM <- na.omit(NEW_counts_starv_IVM)


filtered_counts_starv_IVM <- NEW_counts_starv_IVM  %>%
  group_by(strain, days_to_starvation_edited) %>%
  filter(!(counts == 0 & any(counts != 0)))

fig_starv_ivm <- filtered_counts_starv_IVM  %>%
  ggplot(aes(x = factor(days_to_starvation_edited), y = counts, fill = strain)) +
  geom_bar(stat = "identity", position = position_dodge2(width = 0.8), width = 0.7, color = "black", alpha = 0.7) +
  scale_fill_manual(values = col_comp, name = "Strain") +
  theme_minimal() +
  labs(x = "Days to Starvation", y = "Number of Starved Plates", title = "Ivermectin") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Set the font size for x-axis labels
        axis.title = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.line = element_line(color = "black"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),  # Remove the default panel border
        panel.background = element_blank(),
        strip.text = element_text(size = 12, vjust = 1.5),
        panel.spacing = unit(0.2, "lines")) +  # Add space between facets
  facet_wrap(~ strain, ncol = 3, scales = "free_x") +
  geom_text(aes(label = counts), position = position_dodge(width = 0.7), vjust = -0.5, size = 3) +
  scale_x_discrete(labels = function(x) ifelse(x == "14", "14 +", x), breaks = unique(factor(filtered_counts_starv_IVM $days_to_starvation_edited)))

fig_starv_ivm

#Save figure 
ggsave(fig_starv_ivm, filename = "figures/fig_starv_ivm.jpg", bg = "white", units = "in", height = 7, width = 8.5)

