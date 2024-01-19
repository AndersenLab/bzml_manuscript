# AOS, last edited 2023 Dec 24

#load libraries 
#install.packages("devtools")
devtools::install_github("AndersenLab/easyXpress")
library(tidyverse)
library(dplyr)
library(easyXpress)
library(drc)
library(knitr)
library(RCurl)
library(RColorBrewer)
library(cowplot)
library(ggbeeswarm)
library(ggrepel)
library(ggnewscale)
library(scales)
require(magrittr)
require(grid)
require(gridExtra)
require(sommer)
library(cowplot)
library(gbm)

# set working directory
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

# get the date
today <- format(Sys.Date(), "%Y%m%d")

# figure directory 
figure_dir <- "/projects/b1059/projects/Amanda/projects/BZML/figures"

# assign strain colors
col_comp <- c("N2" = "#FFA500", #orange
              "ECA882" = "#FFE400", #yellow
              "PHX2647" = "#9DDD00", #lime green
              "PHX2502" = "#00B558", #medium green
              "PHX2481" = "#008878", #teal
              "ECA2687" = "#14517F", #lightest dark blue
              "ECA2688" = "#00C3D7", #medium blue
              "ECA2689" = "#2D52E4", #dark blue
              "ECA2690" = "#3A0045") #dark blue

# assign strain labels
strain_labels <- c("N2" = "N2", "ECA882" = "ben-1", "PHX2647" = "avr-14",
                   "PHX2502" = "avr-15", "PHX2481" = "glc-1", "ECA2687" = "avr-14;avr-15",
                   "ECA2688" = "avr-14;glc-1", "ECA2689" = "avr-15,glc-1", "ECA2690" = "avr-14;avr-15,glc-1")

###########################################
#                                         #
#     BZML  Easyxpress - ABZ and DMSO     #
#                                         #
###########################################
# drugs in assay: albendazole (30 uM), DMSO (control)
# all assays combined (1-3)
# DMSO = albendazole (condition), 0 (concentration_um)

# Define experimental directory and file name 
dirsBZML_ABZ <- "/projects/b1059/projects/Amanda/projects/cellprofiler-nf/projects/20231225_BZML_HTA_ABZ/Analysis-20231225" #data and design match
datafileBZML_ABZ <- "20231225_BZML_HTA_ABZ_Analysis-20231225.RData"

BZML_ABZ <- easyXpress::readXpress(
  filedir = dirsBZML_ABZ,
  rdafile = datafileBZML_ABZ,
  design = TRUE,
  doseR = FALSE)

# model selection 
ms_ABZ <- easyXpress::modelSelection(BZML_ABZ$raw_data)

# flag objects
ef_ABZ <- edgeOF(data = ms_ABZ)

# Use the clusterOF function
cf_ABZ <- clusterOF(data = ef_ABZ)

# Plot flagged objects 
c1_ABZ <- checkOF(data = cf_ABZ, strain, drug)

# loop at the plot returned
flagged_plot_ABZ <- c1_ABZ$p

# save flagged object plot 
ggsave(
  filename = glue::glue("{figure_dir}/{today}_flagged_plot.png"),
  plot = flagged_plot_ABZ,
  width = 10,
  height = 10,
  units = "in",
  dpi = 300)

# flagged objects by condition - DMSO 
c1_of_strain_DMSO <- checkOF(
  data = dplyr::filter(cf_ABZ, concentration_um == "0"),
  strain)

c1_of_strain_DMSOplot <- c1_of_strain_DMSO$p

# save plot 
ggsave(
  filename = glue::glue("{figure_dir}/{today}_flagged_plot_by_strain_DMSO.png"),
  plot = c1_of_strain_DMSOplot,
  width = 10,
  height = 10,
  units = "in",
  dpi = 300)

# flagged objects by condition - albendzole 
c1_of_strain_ABZ <- checkOF(
  data = dplyr::filter(cf_ABZ, concentration_um == "30"),
  strain)

c1_of_strain_ABZplot <- c1_of_strain_ABZ$p


# save plot 
ggsave(
  filename = glue::glue("{figure_dir}/{today}_flagged_plot_by_strain_ABZ.png"),
  plot = c1_of_strain_ABZplot,
  width = 10,
  height = 10,
  units = "in",
  dpi = 300)


#####################################
##                                 ##
##  Check the size distribution    ##
##                                 ##
#####################################
c2_ABZ <- checkObjs(data = cf_ABZ, OF = 'filter', strain, drug, concentration_um)

cluster_plot_ABZ <- c2_ABZ
cluster_plot_ABZ

ggsave(
  filename = glue::glue("{figure_dir}/{today}_size_distribution.png"),
  plot = cluster_plot_ABZ,
  width = 10,
  height = 10,
  units = "in",
  dpi = 300)

#####################
##                 ##
##  Check models   ##
##                 ##
#####################
# Create a variable to specify the path to the images
checkmodels_ABZ <- cf_ABZ %>%
  mutate(
    i.dir = "/projects/b1059/projects/Amanda/projects/cellprofiler-nf/projects/20231225_BZML_HTA_ABZ/Analysis-20231225/processed_images",
    w.lab = paste(drug, strain, sep = "_"))

# Now we can run the checkModels function
cm_ABZ_out <- checkModels(data = checkmodels_ABZ ,
                      # the grouping vars (...), make a plot for each.
                      drug,
                      proc.img.dir = "i.dir",
                      well.label = "w.lab",
                      # save in the repo you cloned
                      out.dir = glue::glue("{figure_dir}/{today}_checkModels"))


#########################
##                     ##
##  Get rid of debris  ##
##                     ##
#########################
# add the user variable that will be converted to an object flag
u_ABZ = checkmodels_ABZ %>%
  dplyr::mutate(user = dplyr::case_when(drug == "albendazole" &
                                          model == "MDHD" ~ "junk",
                                        drug == "albendazole" &
                                        worm_length_um < 165 ~ "junk",
                                        TRUE ~ NA_character_))

# Run the userOF function and specify user variable as the flag
uf_ABZ <- easyXpress::userOF(data = u_ABZ, user)

# Check the object data again to see if the bimodal distributions are resolved.
checkObjs(data = uf_ABZ, OF = "filter", drug, concentration_um)

#> 3 ObjectFlags detected in data. They were applied in the following order:
#> edge_ObjectFlag
#> cluster_ObjectFlag
#> user_ObjectFlag
#> The flagged objects will be filtered from the plot.

########################################
##                                    ##
##  Apply the classifierOF function   ##
##                                    ##
########################################
cl_ABZ <- classifierOF(data = uf_ABZ) #7 rows with NAs

####################################
##                                ##
##  Apply the outlierOF function  ##
##                                ##
####################################
# apply the outlierOF function
o_ABZ <- easyXpress::outlierOF(data = cl_ABZ)


#####################
##                 ##
##  checkOBjs      ##
##                 ##
#####################
# check objects again, notice how there are 5 ObjectFlags detected now.
checkObjs(data = o_ABZ, OF = 'filter', drug, concentration_um)

###############
##           ##
##  checkOF  ##
##           ##
###############
co2_ABZ <- checkOF(data = o_ABZ, drug, concentration_um) # check object flags again

co2_ABZ$p # check the plot

################
##            ##
##  filterOF  ##
##            ##
################
# finalize the object data made with outlierOF function above.
proc.objs_ABZ <- easyXpress::filterOF(o_ABZ, rmVars = T)

#####################
##                 ##
##  summarizeWells ##
##                 ##
#####################
# remove all flags, summarize wells, and drop object vars all in one function
raw.wells_ABZ <- easyXpress::summarizeWells(data = o_ABZ)


################
##            ##
## FLAG WELLS ##
##            ##
################
# titerWF 
# Use the titerWF function within a particular bleach and drug for each strain
tf_ABZ <- easyXpress::titerWF(data = raw.wells_ABZ,
                          Metadata_Experiment, bleach, strain, drug, 
                          doseR = F)

# look at the diagnostic plot
tf_ABZ$p

# nwf 
# use the nWF function with the data output from the titerWF function above
n_ABZ <- easyXpress::nWF(data = tf_ABZ$d, drug, concentration_um, max = 30, min = 5)

# Look at the diagnostic plot
n_ABZ$p

# outlierWF
# use the outlierWF function to flag outliers within groups
# then use dplyr::mutate to add a variable to indicate independent bleaches.
ow_ABZ <- easyXpress::outlierWF(data = n_ABZ$d,
                            Metadata_Experiment, bleach, drug,
                            concentration_um, strain) %>%
  dplyr::mutate(assay_bleach = paste(Metadata_Experiment, bleach, sep = "_"))

# checkWF
# Use the checkWF function
cw1_ABZ <- easyXpress::checkWF(data = ow_ABZ, drug, concentration_um)

# look at the plot
cw1_ABZ$p

# filterWF
# Use the filterWF function and drop the flagging variables afterward
fw_ABZ <- filterWF(data = ow_ABZ, rmVars = T)

# checkBalance
# use the checkBalance function and add assay_bleach var to design 
# ISSUE - IT CANT FIND Metadata_Experiment but it does exist 
# cb_ABZ <- checkBalance(data = fw_ABZ, drug, concentration_um,
#                    design = BZML_ABZ$design %>%
#                      dplyr::mutate(assay_bleach =
#                                      paste(bleach)),
#                    x = assay_bleach)

# Look at the plot and add a nicer x-axis with ggplot2 package
# cb_assay1$p +
#   ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5)) +
#   ggplot2::geom_hline(yintercept = 0.75, linetype = 2, color = "red")

# checkEff - Assay1
# filter the unbalanced bleaches - nearly half the data in this example.
# drop_ABZ <- fw_ABZ %>%
#   dplyr::mutate(b.filter =
#                   dplyr::case_when(concentration_um == "0" & strain == "ECA2687"
#                                    & mean_wormlength_um > 800) ~  "drop",
#                                    TRUE ~ "keep") %>%
#   dplyr::filter(b.filter == "keep") %>%
#   dplyr::select(-b.filter)

drop_ABZ <- fw_ABZ %>%
  filter(!(concentration_um == "0" & strain == "ECA2687" & mean_wormlength_um > 800))


# use the checkEff function
ce1_ABZ <- easyXpress::checkEff(data = drop_ABZ,
                            drug, strain,
                            x = concentration_um,
                            y = median_wormlength_um,
                            fill = assay_bleach,
                            scales = "free_x")

# look at the plot
ce1_ABZ

# Get rid of rows with NAs 
# fw_assay1_clean <- fw_assay1 %>%
#   filter(!is.na(drug))  # Remove rows where drug is NA

####################################
##                                ##
## Finalize Results - All Assays  ##
##                                ##
####################################
# regEff
# Regress the effect of independent bleaches for each drug using regEff()
reg_ABZ <- easyXpress::regEff(data = drop_ABZ, #drop
                          drug,
                          d.var = median_wormlength_um,
                          c.var = assay_bleach)

# Look at the regression coefficients in the diagnostic plot p2
reg_ABZ$p2


# delta 
# use the delta() function
del_ABZ <- easyXpress::delta(data = reg_ABZ$d,
                         assay_bleach, drug, strain, # group with ...
                         WF = "filter",
                         doseR = TRUE,
                         vars = "median_wormlength_um_reg")

View(del_assay)

# check the finalized data
checkEff(data = del_ABZ, drug, strain, x = concentration_um,
         y = median_wormlength_um_reg_delta,
         fill = assay_bleach,
         scales = "free_x")

write.csv(del_ABZ, "data/20231227_BZML_del_ABZ_delta_regression.csv", row.names = FALSE)


#############################
#                           #
#   Stats + Final Plots     #
#                           #
#############################
# Make sure there are no NAs in delta reg df
del_ABZ_clean <- del_ABZ %>%
  filter(!is.na(drug))  # Remove rows where drug is NA

View(del_ABZ_clean)



# Box plot and stats for ABZ
# # ANOVA and Tukey HSD for abz
statsabz <- del_ABZ_clean %>%
  dplyr::filter(drug == "albendazole")%>%
  dplyr::filter(concentration_um == "30")%>%
  aov(median_wormlength_um_reg_delta ~ strain, data = .) %>%
  rstatix::tukey_hsd() %>%
  dplyr::filter(group1 == 'N2' | group2 == 'N2')

View(statsabz)
# 
# # Flip to always have N2 represented in the group 2 column 
statsabz <- statsabz %>%
  mutate(
    temp = ifelse(group1 == 'N2', group2, group1),
    group2 = ifelse(group1 == 'N2', 'N2', group2),
    group1 = temp
  ) %>%
  select(-temp)

# Box plot for ABZ
HTA.fig.abz <- del_ABZ %>%  
  dplyr::mutate(strain=factor(strain, levels = c("N2","ECA882","PHX2647","PHX2502","PHX2481","ECA2687","ECA2688","ECA2689","ECA2690")))%>%
  dplyr::filter(drug %in% c("albendazole")) %>% 
  dplyr::filter(concentration_um %in% c("30")) %>% 
  ggplot +
  aes(x=strain, y= median_wormlength_um_reg_delta) +
  geom_jitter(width=0.1, size = 0.3) +
  geom_boxplot(aes(fill = strain),alpha = 0.8,outlier.shape = NA) +
  xlab(" ")+
  ylab("Normalized Animal Length") +
  # ggtitle("Albendazole")+
  scale_fill_manual(values  = c("N2"="#FFA500","ECA882"="#FFE400","PHX2647"="#9DDD00","PHX2502"="#00B558","PHX2481"="#008878","ECA2687"="#14517F","ECA2688"="#00C3D7","ECA2689"="#2D52E4","ECA2690"="#3A0045"))+
  scale_x_discrete(breaks = c("N2","ECA882", "PHX2647","PHX2502","PHX2481","ECA2687","ECA2688","ECA2689","ECA2690"),
                   labels = c("N2\n \n ", "ben-1\n \n ","avr-14\n \n ","avr-15\n \n ","glc-1\n \n ","avr-14\navr-15\n ","avr-14\nglc-1\n ","avr-15\nglc-1\n ","avr-14\navr-15\nglc-1"))+
  stat_pvalue_manual(statsabz, label = "p.adj.signif",y.position = c(10), xmax="group1", remove.bracket = TRUE)+
  cowplot::theme_cowplot(12) +
  theme(axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 1, margin = unit(c(0, 0, 0, 0), units = "in"), face = "italic"),
  # theme(axis.text.y = element_text(size = 12),
        plot.background = element_rect(fill="white"),
        strip.background = element_blank(),
        # axis.title.x =   element_text(size = 12),
        strip.text = element_text(face = "bold", color = "black"),
        plot.title = element_text(face="bold"),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        # axis.text.x = element_text(angle = 0, vjust = 2),  # Rotate x-axis labels
        # plot.title = element_text(vjust = 2)  # Adjust title position using vjust
  )

HTA.fig.abz
ggsave("figures/HTA.abz-20231227.png", plot = HTA.fig.abz, width = 8, height = 4, units = "in")


# Box plot and stats for DMSO
# # ANOVA and Tukey HSD for DMSO
statsDMSO <- del_ABZ_clean %>%
  dplyr::filter(drug == "albendazole")%>%
  dplyr::filter(concentration_um == "0")%>%
  aov(median_wormlength_um ~ strain, data = .) %>%
  rstatix::tukey_hsd() %>%
  dplyr::filter(group1 == 'N2' | group2 == 'N2')

View(statsDMSO)

# # Flip to always have N2 represented in the group 2 column 
statsDMSO  <- statsDMSO  %>%
  mutate(
    temp = ifelse(group1 == 'N2', group2, group1),
    group2 = ifelse(group1 == 'N2', 'N2', group2),
    group1 = temp
  ) %>%
  select(-temp)

# Box plot for DMSO
HTA.fig.DMSO <- del_ABZ %>%  
  dplyr::mutate(strain=factor(strain, levels = c("N2","ECA882","PHX2647","PHX2502","PHX2481","ECA2687","ECA2688","ECA2689","ECA2690")))%>%
  dplyr::filter(drug %in% c("albendazole")) %>% 
  dplyr::filter(concentration_um %in% c("0")) %>% 
  ggplot +
  aes(x=strain, y= median_wormlength_um) +
  geom_jitter(width=0.1, size = 0.3) +
  geom_boxplot(aes(fill = strain),alpha = 0.8,outlier.shape = NA) +
  xlab("Strain")+
  ylab("Median Animal Length") +
  # ggtitle("DMSO")+
  scale_fill_manual(values  = c("N2"="#FFA500","ECA882"="#FFE400","PHX2647"="#9DDD00","PHX2502"="#00B558","PHX2481"="#008878","ECA2687"="#14517F","ECA2688"="#00C3D7","ECA2689"="#2D52E4","ECA2690"="#3A0045"))+
  scale_x_discrete(breaks = c("N2","ECA882", "PHX2647","PHX2502","PHX2481","ECA2687","ECA2688","ECA2689","ECA2690"),
                   labels = c("N2\n \n ", "ben-1\n \n ","avr-14\n \n ","avr-15\n \n ","glc-1\n \n ","avr-14\navr-15\n ","avr-14\nglc-1\n ","avr-15\nglc-1\n ","avr-14\navr-15\nglc-1"))+
  stat_pvalue_manual(statsDMSO, label = "p.adj.signif",y.position = c(1000), xmax="group1", remove.bracket = TRUE)+
  cowplot::theme_cowplot(12) +
  theme(axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 1, margin = unit(c(0, 0, 0, 0), units = "in"), face = "italic"),
        # theme(axis.text.y = element_text(size = 12),
        plot.background = element_rect(fill="white"),
        strip.background = element_blank(),
        # axis.title.x =   element_text(size = 12),
        strip.text = element_text(face = "bold", color = "black"),
        plot.title = element_text(face="bold"),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        # axis.text.x = element_text(angle = 0, vjust = 2),  # Rotate x-axis labels
        # plot.title = element_text(vjust = 2)  # Adjust title position using vjust
  )

HTA.fig.DMSO

ggsave("figures/HTA.fig.DMSO_20230831.png", plot = HTA.fig.DMSO, width = 8, height = 4, units = "in")


HTA.fig.DMSO <- del_ABZ %>%  
  dplyr::mutate(strain = factor(strain, levels = c("N2", "ECA882", "PHX2647", "PHX2502", "PHX2481", "ECA2687", "ECA2688", "ECA2689", "ECA2690"))) %>%
  dplyr::filter(drug %in% c("albendazole")) %>% 
  dplyr::filter(concentration_um %in% c("0")) %>% 
  ggplot() +
  aes(x = strain, y = median_wormlength_um) +
  geom_jitter(width = 0.1, size = 0.3) +
  geom_boxplot(aes(fill = strain), alpha = 0.8, outlier.shape = NA) +
  xlab("Strain") +
  ylab("Normalized Animal Length") +
  # ggtitle("DMSO") +
  scale_fill_manual(values = c("N2" = "#FFA500", "ECA882" = "#FFE400", "PHX2647" = "#9DDD00", "PHX2502" = "#00B558", "PHX2481" = "#008878", "ECA2687" = "#14517F", "ECA2688" = "#00C3D7", "ECA2689" = "#2D52E4", "ECA2690" = "#3A0045")) +
  scale_x_discrete(
    breaks = c("N2", "ECA882", "PHX2647", "PHX2502", "PHX2481", "ECA2687", "ECA2688", "ECA2689", "ECA2690"),
    labels = c("N2", "ben-1", "avr-14", "avr-15", "glc-1", "avr-14\navr-15", "avr-14\nglc-1", "avr-15\nglc-1", "avr-14\navr-15\nglc-1")
  ) +
  stat_pvalue_manual(statsDMSO, label = "p.adj.signif", y.position = c(1000), xmax = "group1", remove.bracket = TRUE) +
  cowplot::theme_cowplot(12) +
  theme(
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      hjust = 0.5,
      vjust = 1,
      margin = unit(c(0, 0, 0, 0), units = "in")
    ),
    axis.title.x = element_text(size = 12),
    plot.background = element_rect(fill = "white"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", color = "black"),
    plot.title = element_text(face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(face = ifelse(c("N2", "ECA882", "PHX2647", "PHX2502", "PHX2481", "ECA2687", "ECA2688", "ECA2689", "ECA2690") == "N2", "plain", "italic"))
  )


HTA.fig.DMSO




































##### OLD 
Fig3_label_positions <- data.frame(label = c("A", "B"),
                                   x = c(0.014, 0.51),  # X-coordinate
                                   y = c(0.98, 0.98)) # Y-coordinate

Fig3_with_labels <- Fig3 +
  lapply(1:nrow(Fig3_label_positions), function(i) {
    draw_label(Fig3_label_positions$label[i], x = Fig3_label_positions$x[i], y = Fig3_label_positions$y[i],
               size = 16, fontface = "bold")
  })


Fig3_with_labels

Fig3 <- cowplot::plot_grid(HTA.fig.abz,HTA.fig.ivm25, labels = c("A","B"))
Fig3

# Supplemental Figure 2 - DMSO HTA 
# ANOVA and Tukey HSD for DMSO control  
statsDMSO <- ass.effs %>% 
  dplyr::filter(drug == "DMSO")%>%
  aov(median_wormlength_um ~ strain, data = .) %>% 
  rstatix::tukey_hsd() %>% 
  dplyr::filter(group1 == 'N2' | group2 == 'N2')

View(statsDMSO)

statsDMSO <- statsDMSO %>%
  mutate(
    temp = ifelse(group1 == 'N2', group2, group1),
    group2 = ifelse(group1 == 'N2', 'N2', group2),
    group1 = temp
  ) %>%
  select(-temp)



# Supplemental Figure 3 - Ivermectin 0.50 uM
# ANOVA and Tukey HSD for second concentration of Ivermectin (0.50 uM) 
statsivm50 <- BZML.cleaned.ivm50 %>% 
  dplyr::filter(drug == "ivermectin_50")%>%
  aov(median_wormlength_um ~ strain, data = .) %>% 
  rstatix::tukey_hsd() %>% 
  dplyr::filter(group1 == 'N2' | group2 == 'N2')

View(statsivm50)

statsivm50<- statsivm50 %>%
  mutate(
    temp = ifelse(group1 == 'N2', group2, group1),
    group2 = ifelse(group1 == 'N2', 'N2', group2),
    group1 = temp
  ) %>%
  select(-temp)

# Box plot for Ivermectin 50 uM
HTA.fig.ivm50 <- BZML.cleaned.ivm50 %>%  
  dplyr::mutate(strain=factor(strain, levels = c("N2","ECA882","PHX2647","PHX2502","PHX2481","ECA2687","ECA2688","ECA2689","ECA2690")))%>%
  dplyr::filter(drug %in% c("ivermectin_50")) %>% 
  ggplot() +
  aes(x=strain, y= median_wormlength_um_delta_reg) +
  geom_jitter(width=0.1, size = 0.3) +
  geom_boxplot(aes(fill = strain),alpha = 0.8,outlier.shape = NA) +
  xlab(" ") +
  ylab("Normalized Animal Length") +
  ggtitle("Ivermectin") +
  scale_fill_manual(values  = c("N2"="#FFA500","ECA882"="#FFE400","PHX2647"="#9DDD00","PHX2502"="#00B558","PHX2481"="#008878","ECA2687"="#14517F","ECA2688"="#00C3D7","ECA2689"="#2D52E4","ECA2690"="#3A0045"))+
  scale_x_discrete(breaks = c("N2","ECA882","PHX2647","PHX2502","PHX2481","ECA2687","ECA2688","ECA2689","ECA2690"),  
                   labels = strain_labels) + 
  stat_pvalue_manual(statsivm50, label = "p.adj.signif",y.position = c(210), xmax="group1", remove.bracket = TRUE) +
  cowplot::theme_cowplot(12) +
  theme(axis.text.y = element_text(size = 12),
        plot.background = element_rect(fill="white"),
        strip.background = element_blank(),
        axis.title.x =   element_text(size = 12),
        strip.text = element_text(face = "bold", color = "black"),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
        plot.title = element_text(face="bold", vjust = 2)  # Adjust title position using vjust
  )

HTA.fig.ivm50

ggsave("figures/HTA.fig.ivm50_20230831.png", plot = HTA.fig.ivm50, width = 7.5, height = 3.5, units = "in")


# Figure for Erik's R01 
BZML.cleaned.ivm50.R01 <- BZML.cleaned.ivm50 %>%
  dplyr::filter(strain != "ECA882")

# ANOVA and Tukey HSD for IVM 0.25 uM 
statsivm50.R01 <- BZML.cleaned.ivm50.R01  %>% 
  dplyr::filter(drug == "ivermectin_50")%>%
  aov(median_wormlength_um_delta_reg ~ strain, data = .) %>% 
  rstatix::tukey_hsd() %>% 
  dplyr::filter(group1 == 'N2' | group2 == 'N2')


# Flip to always have N2 represented in the group 2 column 
statsivm50.R01 <- statsivm50.R01 %>%
  mutate(
    temp = ifelse(group1 == 'N2', group2, group1),
    group2 = ifelse(group1 == 'N2', 'N2', group2),
    group1 = temp
  ) %>%
  select(-temp)

HTA.fig.ivm50_noben1 <- BZML.cleaned.ivm50.R01 %>%  
  dplyr::mutate(strain=factor(strain, levels = c("N2","PHX2647","PHX2502","PHX2481","ECA2687","ECA2688","ECA2689","ECA2690")))%>%
  dplyr::filter(drug %in% c("ivermectin_50")) %>% 
  ggplot+
  aes(x=strain, y= median_wormlength_um_delta_reg)+
  geom_jitter(width=0.1, size = 0.3)+
  geom_boxplot(aes(fill = strain),alpha = 0.8,outlier.shape = NA)+
  xlab(" ")+
  ylab("Normalized Animal Length") +
  ggtitle("Ivermectin")+
  scale_fill_manual(values  = c("N2"="#FFA500","PHX2647"="#9DDD00","PHX2502"="#00B558","PHX2481"="#008878","ECA2687"="#14517F","ECA2688"="#00C3D7","ECA2689"="#2D52E4","ECA2690"="#3A0045"))+
  scale_x_discrete(breaks = c("N2","PHX2647","PHX2502","PHX2481","ECA2687","ECA2688","ECA2689","ECA2690"),  #,labels = c(x1,x2,x3,x4,x5,x6,x7,x8,x9))+
                   labels = strain_labels_2) + 
  stat_pvalue_manual(statsivm50.R01, label = "p.adj.signif",y.position = c(210), xmax="group1", remove.bracket = TRUE)+
  cowplot::theme_cowplot(12)+
  theme(axis.text.y = element_text(size = 12),
        plot.background = element_rect(fill="white"),
        strip.background = element_blank(),
        axis.title.x =   element_text(size = 12),
        strip.text = element_text(face = "bold", color = "black"),
        plot.title = element_text(face="bold", vjust = 2),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

HTA.fig.ivm50_noben1

ggsave("figures/HTA.fig.ivm50_noben1.png", plot = HTA.fig.ivm50_noben1, width = 7.5, height = 3.5, units = "in")


# Box plot for IVM 0.25 Erik's RO1 
# Remove ECA882 from df 
BZML.cleaned.ivm25.R01 <- BZML.cleaned.ivm25 %>%
  dplyr::filter(strain != "ECA882")

# ANOVA and Tukey HSD for IVM 0.25 uM 
statsivm25.R01 <- BZML.cleaned.ivm25.R01  %>% 
  dplyr::filter(drug == "ivermectin_25")%>%
  aov(median_wormlength_um_delta_reg ~ strain, data = .) %>% 
  rstatix::tukey_hsd() %>% 
  dplyr::filter(group1 == 'N2' | group2 == 'N2')


# Flip to always have N2 represented in the group 2 column 
statsivm25.R01 <- statsivm25.R01 %>%
  mutate(
    temp = ifelse(group1 == 'N2', group2, group1),
    group2 = ifelse(group1 == 'N2', 'N2', group2),
    group1 = temp
  ) %>%
  select(-temp)


HTA.fig.ivm25_noben1 <- BZML.cleaned.ivm25.R01 %>%  
  dplyr::mutate(strain=factor(strain, levels = c("N2","PHX2647","PHX2502","PHX2481","ECA2687","ECA2688","ECA2689","ECA2690")))%>%
  dplyr::filter(drug %in% c("ivermectin_25")) %>% 
  ggplot+
  aes(x=strain, y= median_wormlength_um_delta_reg)+
  geom_jitter(width=0.1, size = 0.3)+
  geom_boxplot(aes(fill = strain),alpha = 0.8,outlier.shape = NA)+
  xlab(" ")+
  ylab("Normalized Animal Length") +
  ggtitle("Ivermectin")+
  scale_fill_manual(values  = c("N2"="#FFA500","PHX2647"="#9DDD00","PHX2502"="#00B558","PHX2481"="#008878","ECA2687"="#14517F","ECA2688"="#00C3D7","ECA2689"="#2D52E4","ECA2690"="#3A0045"))+
  scale_x_discrete(breaks = c("N2","ECA882", "PHX2647","PHX2502","PHX2481","ECA2687","ECA2688","ECA2689","ECA2690"),
                   labels = c("N2\n \n ", "ben-1\n \n ","avr-14\n \n ","avr-15\n \n ","glc-1\n \n ","avr-14\navr-15\n ","avr-14\nglc-1\n ","avr-15\nglc-1\n ","avr-14\navr-15\nglc-1"))+
  stat_pvalue_manual(statsivm25.R01, label = "p.adj.signif",y.position = c(210), xmax="group1", remove.bracket = TRUE)+
  cowplot::theme_cowplot(12)+
  theme(axis.text.y = element_text(size = 12),
        plot.background = element_rect(fill="white"),
        strip.background = element_blank(),
        axis.title.x =   element_text(size = 12),
        strip.text = element_text(face = "bold", color = "black"),
        plot.title = element_text(face="bold", vjust = 2),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1))  

HTA.fig.ivm25_noben1

ggsave("figures/HTA.fig.ivm25_noben1.png", plot = HTA.fig.ivm25_noben1, width = 8, height = 4, units = "in")


