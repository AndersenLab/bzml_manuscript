# Fig 1 plots - allele frequencies + competitive fitness - DMSO, ABZ, and IVM
# A = DMSO allele freq
# B = DMSO comp fit 
# C = ABZ allele freq 
# D = ABZ comp fit 
# E = IVM allele freq 
# F = IVM comp fit 

# Fig SX - Extra v3 DMSO plot 
# A = DMSO allele freq
# B = DMSO comp fit

#load libraries 
library(ggstatsplot)
library(dplyr)
library(ggpubr)
library(rstatix)

#strain colors
col_comp <- c("N2" = "#FFA500", #orange
              "ECA882" = "#FFE400", #yellow
              "PHX2647" = "#9DDD00", #lime green
              "PHX2502" = "#00B558", #medium green
              "PHX2481" = "#008878", #teal
              "ECA2687" = "#14517F", #lighest dark blue
              "ECA2688" = "#00C3D7", #medium blue
              "ECA2689" = "#2D52E4", #dark blue
              "ECA2690" = "#3A0045") #dark blue

###################################
#                                 #
#       Fig 1A + B  v2 DMSO       #
#                                 #
###################################

# Import v2 dataframe - DMSO and IVM
comp_assays_ivm_mean_sd_df <- read_csv("data/competition_assays/comp_assays_ivm_mean_sd_df.csv")

# Fig 1 A - DMSO allele freq
allelefreq_DMSO_v2 <- BZML_total_dsmo_fitness_v2  %>%
  dplyr::filter(condition == "DMSO") %>%
  ggplot() +
  aes(x = generation, y = mean, color = strain) +
  geom_point(size = 3) +
  ylab("Relative allele frequency") +
  geom_errorbar(
    aes(ymin = (mean - sd), ymax = (mean + sd)),
    width = 0.25,
    size = 0.5
  ) +
  geom_line(size = 1) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), limits = c(0, 1.0)) +
  scale_color_manual(
    name = "Strain",
    labels = c(
      "N2" = "N2",
      "PHX2647" = "avr-14",
      "PHX2502" = "avr-15",
      "PHX2481" = "glc-1",
      "ECA2687" = "avr-14;avr-15",
      "ECA2688" = "avr-14;glc-1",
      "ECA2689" = "avr-15,glc-1",
      "ECA2690" = "avr-14;avr-15,glc-1"
    ),
    values = col_comp
  ) +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 7)) +
  cowplot::theme_cowplot(12) +
theme(legend.position = "none",
      # axis.text.y.right = element_blank(),
      # axis.ticks.y.right = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      plot.background = element_rect(fill = "white"))

allelefreq_DMSO

# Fig 1 B - DMSO comp fitness

# Calculate v3 DMSO stats for comp fitness assay
BZML_stats_v2 <- BZML_total_dsmo_fitness_v2 %>%
  dplyr::filter(condition == "DMSO")%>%
  aov(fitness ~ strain, data = .)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(group1 == 'N2' | group2 == 'N2')

# Flip to always have N2 represented in the group 2 column
BZML_stats_v2  <- BZML_stats_v2  %>%
  mutate(
    temp = group1,
    group1 = ifelse(group1 == 'N2', group2, group1),
    group2 = ifelse(temp == 'N2', temp, group2)
  ) %>%
  dplyr::select(-temp)

View(BZML_stats_v2 )

BZML_stats_v2   <- BZML_stats_v2   %>%
  mutate(
    group2 = case_when(
      group2 %in% c("PHX2481", "PHX2502", "PHX2647") ~ "N2",
      TRUE ~ group2))

BZML_stats_v2   <- BZML_stats_v2  %>% 
  dplyr::mutate(strain= group1)

dmso_v2_plot <- BZML_total_dsmo_fitness_v2 %>%
  # dplyr::filter(condition == "DMSO") %>%
  dplyr::mutate(strain = factor(strain, levels = c("N2","ECA882", "PHX2647", "PHX2502", "PHX2481", "ECA2687", "ECA2688", "ECA2689", "ECA2690")))%>%
  ggplot() +
  aes(x = strain, y = (fitness), fill = strain) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  geom_hline(yintercept = 0) +
  scale_x_discrete(breaks = c("N2","ECA882", "PHX2647","PHX2502","PHX2481","ECA2687","ECA2688","ECA2689","ECA2690"),
                   labels = c("N2\n \n ", "ben-1\n \n ","avr-14\n \n ","avr-15\n \n ","glc-1\n \n ","avr-14\navr-15\n ","avr-14\nglc-1\n ","avr-15\nglc-1\n ","avr-14\navr-15\nglc-1")) +
  stat_pvalue_manual(BZML_stats_v2, label = "p.adj.signif", y.position = c(2), xmax="group1", remove.bracket = TRUE) +
  scale_y_continuous(sec.axis = dup_axis(name = "DMSO")) +
  scale_fill_manual(values = col_comp) +  # Set fill colors using col_comp
  cowplot::theme_cowplot(12) +
  ylab("Competitive fitness") +
  theme(legend.position = "none",
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

dmso_v2_plot


###################################
#                                 #
#       Fig 1E + F  v2 IVM        #
#                                 #
###################################


# Fig 1 E - IVM allele freq 
allelefreq_IVM <- comp_assays_ivm_mean_sd_df  %>%
  dplyr::filter(condition == "Ivermectin") %>%
  ggplot() +
  aes(x = generation, y = mean, color = strain) +
  geom_point(size = 3) +
  ylab("Relative allele frequency") +
  xlab("Generation") +
  geom_errorbar(
    aes(ymin = (mean - sd), ymax = (mean + sd)),
    width = 0.25,
    size = 0.5
  ) +
  geom_line(size = 1) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), limits = c(0, 1.0)) +
  scale_color_manual(
    name = "Strain",
    labels = c(
      "N2" = "N2",
      "PHX2647" = "avr-14",
      "PHX2502" = "avr-15",
      "PHX2481" = "glc-1",
      "ECA2687" = "avr-14;avr-15",
      "ECA2688" = "avr-14;glc-1",
      "ECA2689" = "avr-15,glc-1",
      "ECA2690" = "avr-14;avr-15,glc-1"
    ),
    values = col_comp
  ) +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 7)) +
  cowplot::theme_cowplot(12) +
  theme(legend.position = "none", plot.background = element_rect(fill = "white"))

allelefreq_IVM



# Fig 1 F - IVM comp fitness
BZML_stats_v2_IVM <- BZML_total_ivm_fitness %>%
  dplyr::filter(condition == "Ivermectin")%>%
  aov(fitness ~ strain, data = .)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(group1 == 'N2' | group2 == 'N2')

# Flip to always have N2 represented in the group 2 column
BZML_stats_v2_IVM   <- BZML_stats_v2_IVM   %>%
  mutate(
    temp = group1,
    group1 = ifelse(group1 == 'N2', group2, group1),
    group2 = ifelse(temp == 'N2', temp, group2)
  ) %>%
  dplyr::select(-temp)

View(BZML_stats_v2_IVM )

BZML_stats_v2_IVM    <- BZML_stats_v2_IVM   %>%
  mutate(
    group2 = case_when(
      group2 %in% c("PHX2481", "PHX2502", "PHX2647") ~ "N2",
      TRUE ~ group2))

BZML_stats_v2_IVM    <- BZML_stats_v2_IVM   %>% 
  dplyr::mutate(strain= group1)

IVM_v2_plot <- BZML_total_ivm_fitness %>%
  # dplyr::filter(condition == "DMSO") %>%
  dplyr::mutate(strain = factor(strain, levels = c("N2","ECA882", "PHX2647", "PHX2502", "PHX2481", "ECA2687", "ECA2688", "ECA2689", "ECA2690")))%>%
  ggplot() +
  aes(x = strain, y = (fitness), fill = strain) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  geom_hline(yintercept = 0) +
  scale_x_discrete(breaks = c("N2","ECA882", "PHX2647","PHX2502","PHX2481","ECA2687","ECA2688","ECA2689","ECA2690"),
                   labels = c("N2\n \n ", "ben-1\n \n ","avr-14\n \n ","avr-15\n \n ","glc-1\n \n ","avr-14\navr-15\n ","avr-14\nglc-1\n ","avr-15\nglc-1\n ","avr-14\navr-15\nglc-1")) +
  stat_pvalue_manual(BZML_stats_v2_IVM, label = "p.adj.signif", y.position = c(7), xmax="group1", remove.bracket = TRUE) +
  scale_y_continuous(sec.axis = dup_axis(name = "Ivermectin")) +
  scale_fill_manual(values = col_comp) +  # Set fill colors using col_comp
  cowplot::theme_cowplot(12) +
  ylab("Competitive fitness") +
  theme(legend.position = "none",
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        # axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 1, margin = unit(c(0, 0, 0, 0), units = "in"), face = "italic"),
        axis.title.x = element_blank())

IVM_v2_plot




# Fig 1 C - ABZ allele freq 
allelefreq_v3_ABZ <- v3_grouped %>%
  dplyr::filter(condition == "albendazole") %>%
  ggplot() +
  aes(x = generation, y = mean_allele_freq, color = strain) +
  geom_point(size = 3) +
  ylab("Relative allele frequency") +
  xlab("Generation") +
  geom_errorbar(
    aes(ymin = (mean_allele_freq - st_dev_allele_freq), ymax = (mean_allele_freq + st_dev_allele_freq)),
    width = 0.25,
    size = 0.5
  ) +
  geom_line(size = 1) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), limits = c(0, 1.0)) +
  scale_color_manual(
    name = "Strain",
    labels = c(
      "N2" = "N2",
      "PHX2647" = "avr-14",
      "PHX2502" = "avr-15",
      "PHX2481" = "glc-1",
      "ECA2687" = "avr-14;avr-15",
      "ECA2688" = "avr-14;glc-1",
      "ECA2689" = "avr-15,glc-1",
      "ECA2690" = "avr-14;avr-15,glc-1"
    ),
    values = col_comp
  ) +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 7)) +
  cowplot::theme_cowplot(12) +
  theme(legend.position = "none", plot.background = element_rect(fill = "white"))


allelefreq_v3_ABZ

# stats for ABZ comp fitness 

# Calculate v3 DMSO stats for comp fitness assay
BZML_stats_ABZ <- BZML_total_abz_fitness %>%
  dplyr::filter(condition == "albendazole")%>%
  aov(fitness ~ strain, data = .)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(group1 == 'N2' | group2 == 'N2')

# Flip to always have N2 represented in the group 2 column
BZML_stats_ABZ  <- BZML_stats_ABZ  %>%
  mutate(
    temp = group1,
    group1 = ifelse(group1 == 'N2', group2, group1),
    group2 = ifelse(temp == 'N2', temp, group2)
  ) %>%
  dplyr::select(-temp)

View(BZML_stats_ABZ)

BZML_stats_ABZ   <- BZML_stats_ABZ   %>%
  mutate(
    group2 = case_when(
      group2 %in% c("PHX2481", "PHX2502", "PHX2647") ~ "N2",
      TRUE ~ group2))

BZML_stats_ABZ   <- BZML_stats_ABZ  %>% 
  dplyr::mutate(strain= group1)




# Fig 1D Plot ABZ competitive fitness 
abz_plot <- BZML_total_abz_fitness %>%
  dplyr::filter(condition == "albendazole") %>%
  dplyr::mutate(strain = factor(strain, levels = c("N2","ECA882", "PHX2647", "PHX2502", "PHX2481", "ECA2687", "ECA2688", "ECA2689", "ECA2690")))%>%
  ggplot() +
  aes(x = strain, y = (fitness), fill = strain) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  geom_hline(yintercept = 0) +
  scale_x_discrete(breaks = c("N2","ECA882", "PHX2647","PHX2502","PHX2481","ECA2687","ECA2688","ECA2689","ECA2690"),
                   labels = c("N2\n \n ", "ben-1\n \n ","avr-14\n \n ","avr-15\n \n ","glc-1\n \n ","avr-14\navr-15\n ","avr-14\nglc-1\n ","avr-15\nglc-1\n ","avr-14\navr-15\nglc-1"))+
  stat_pvalue_manual(BZML_stats_ABZ, label = "p.adj.signif", y.position = c(4), xmax="group1", remove.bracket = TRUE) +
  scale_y_continuous(sec.axis = dup_axis(name = "Albendazole")) +
  scale_fill_manual(values = col_comp) +  # Set fill colors using col_comp
  cowplot::theme_cowplot(12) +
  ylab("Competitive fitness") +
  theme(legend.position = "none",
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

abz_plot



################################
#                              #
#       Fig 1  cowplot         #
#                              #
################################          
#BZML manuscript cowplot 
BZMLcompplot <- cowplot::plot_grid(allelefreq_DMSO, # Fig 1 A - DMSO allele freq
                                   dmso_v2_plot, # Fig 1 B - DMSO comp fit
                                   allelefreq_v3_ABZ, # Fig 1 C - ABZ allele freq
                                   abz_plot, # Fig 1 D - ABZ comp fit 
                                   allelefreq_IVM, # Fig 1 E - IVM allele freq
                                   IVM_v2_plot, # Fig 1 F - IVM comp fit 
                                   ncol = 2,nrow = 3,labels = c("A","B","C","D", "E", "F"),
                                   align = "vh",axis = "lrbt",
                                   label_size = 12,
                                   label_fontfamily = "Arial",
                                   rel_widths = c(1,1,1,1),rel_heights = c(1,1,1,1))

BZMLcompplot

ggsave(filename = "figures/BZML_compfitness_cowplot.eps",plot = BZMLcompplot ,width = 7.5,height = 8,units = "in")



################################
#                              #
#       Fig SX  v3 DMSO        #
#                              #
################################
#v2
v3_grouped <- group_by(v3_total, strain, condition, generation)
v3_mean_sd_df <- summarise(v3_grouped, mean=mean(allele_freq), sd=sd(allele_freq))
# mean_sd_df <- mutate(mean_sd_df, mean_percent = mean * 100)
# mean_sd_df <- mutate(mean_sd_df, sd_percent = sd * 100)

#remove columns
v3_grouped <- v3_grouped %>%
  select(-plate, -well)

#remove empty row with NAs
v3_grouped <-v3_grouped %>%
  filter(!is.na())

write.csv(v3_grouped, "data/competition_assays/v3_grouped.csv", row.names = FALSE)

v3_grouped <- read.csv("data/competition_assays/v3_grouped.csv")

# rm(v3_grouped)

v3_grouped <- v3_grouped %>%
  filter(strain != "PTM229")

allelefreq_v3_DMSO <- v3_grouped %>%
  dplyr::filter(condition == "DMSO") %>%
  ggplot() +
  aes(x = generation, y = mean_allele_freq, color = strain) +
  geom_point(size = 3) +
  ylab("Relative allele frequency") +
  xlab("Generation") +
  geom_errorbar(
    aes(ymin = (mean_allele_freq - st_dev_allele_freq), ymax = (mean_allele_freq + st_dev_allele_freq)),
    width = 0.25,
    size = 0.5
  ) +
  geom_line(size = 1) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), limits = c(0, 1.0)) +
  scale_color_manual(
    name = "Strain",
    labels = c(
      "N2" = "N2",
      "PHX2647" = "avr-14",
      "PHX2502" = "avr-15",
      "PHX2481" = "glc-1",
      "ECA2687" = "avr-14;avr-15",
      "ECA2688" = "avr-14;glc-1",
      "ECA2689" = "avr-15,glc-1",
      "ECA2690" = "avr-14;avr-15,glc-1"
    ),
    values = col_comp
  ) +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 7)) +
  cowplot::theme_cowplot(12) +
  theme(legend.position = "none", plot.background = element_rect(fill = "white"))


allelefreq_v3_DMSO


# Calculate v3 DMSO stats for comp fitness assay
BZML_total_fitness_DMSO_v3 <- BZML_total_dsmo_fitness_v3 %>%
  dplyr::filter(condition == "DMSO")%>%
  aov(fitness ~ strain, data = .)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(group1 == 'N2' | group2 == 'N2')

# Flip to always have N2 represented in the group 2 column
BZML_total_fitness_DMSO_v3 <- BZML_total_fitness_DMSO_v3 %>%
  mutate(
    temp = group1,
    group1 = ifelse(group1 == 'N2', group2, group1),
    group2 = ifelse(temp == 'N2', temp, group2)
  ) %>%
  dplyr::select(-temp)

View(BZML_total_fitness_DMSO_v3)

BZML_total_fitness_DMSO_v3  <- BZML_total_fitness_DMSO_v3  %>%
  mutate(
    group2 = case_when(
      group2 %in% c("PHX2481", "PHX2502", "PHX2647") ~ "N2",
      TRUE ~ group2))

BZML_total_fitness_DMSO_v3  <- BZML_total_fitness_DMSO_v3 %>% 
  dplyr::mutate(strain= group1)

View(BZML_total_fitness_DMSO_v3)


# FIG SXB - DMSO v3 Competitive fitness plot 
dmso_v3_plot <- BZML_total_dsmo_fitness_v3 %>%
  dplyr::filter(condition == "DMSO") %>%
  dplyr::mutate(strain = factor(strain, levels = c("N2","ECA882", "PHX2647", "PHX2502", "PHX2481", "ECA2687", "ECA2688", "ECA2689", "ECA2690")))%>%
  ggplot() +
  aes(x = strain, y = (fitness), fill = strain) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  geom_hline(yintercept = 0) +
  scale_x_discrete(breaks = c("N2","ECA882", "PHX2647","PHX2502","PHX2481","ECA2687","ECA2688","ECA2689","ECA2690"),
                   labels = c("N2\n \n ", "ben-1\n \n ","avr-14\n \n ","avr-15\n \n ","glc-1\n \n ","avr-14\navr-15\n ","avr-14\nglc-1\n ","avr-15\nglc-1\n ","avr-14\navr-15\nglc-1")) +
  stat_pvalue_manual(BZML_total_fitness_DMSO_v3, label = "p.adj.signif", y.position = c(1), xmax="group1", remove.bracket = TRUE) +
  scale_y_continuous(sec.axis = dup_axis(name = "DMSO")) +
  scale_fill_manual(values = col_comp) +  # Set fill colors using col_comp
  cowplot::theme_cowplot(12) +
  ylab("Competitive fitness") +
  theme(legend.position = "none",
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        # axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 0.5, vjust = 1, margin = unit(c(0, 0, 0, 0), units = "in"), face = "italic"),
        axis.title.x = element_blank())

dmso_v3_plot


# Fig SX DMSO cowplot
FigS2 <- cowplot::plot_grid(allelefreq_v3_DMSO , dmso_v3_plot, ncol = 2,nrow = 1,labels = c("A","B"),align = "vh",axis = "lrbt",label_size = 12,label_fontfamily = "Arial",rel_widths = c(1,1),rel_heights = c(1,1))

FigS2

ggsave(filename = "figures/FigS2_v3_DMSO_allele_freq_comp_assaay.eps",plot = FigS2 ,width = 7.5,height = 3,units = "in")


