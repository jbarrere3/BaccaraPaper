#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name _targets.R  
#' @description R script to launch the target pipeline
#' @author Julien BARRERE
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Options and packages ----------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Load targets
library(targets)
# Load functions
lapply(grep("R$", list.files("R"), value = TRUE), function(x) source(file.path("R", x)))
# install if needed and load packages
packages.in <- c("dplyr", "ggplot2", "RCurl", "httr", "tidyr", "data.table", 
                 "lme4", "cowplot", "archive", "terra", "knitr", "car", "sf", 
                 "rnaturalearth", "elevatr", "raster")
for(i in 1:length(packages.in)) if(!(packages.in[i] %in% rownames(installed.packages()))) install.packages(packages.in[i])
# Targets options
options(tidyverse.quiet = TRUE)
tar_option_set(packages = packages.in)
set.seed(2)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Targets workflow --------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

list(
  
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##' - DATA FORMATTING
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # - Raw files
  tar_target(data_seedlings_file, "data/data_seedlings.csv", format = "file"),
  
  # - Read raw files
  tar_target(data_seedlings, fread(data_seedlings_file, h = T, sep = ";")),
  
  # - Format raw files
  tar_target(df_br, (data_seedlings %>% dplyr::select(-LLS) %>% drop_na())), 
  tar_target(df_gr, (data_seedlings %>% drop_na())), 
  
  # - Get precipitation from chelsa bioclim
  tar_target(plot_coord, extract_coord(df_gr)), 
  tar_target(chelsa_files, download_CHELSA(bioclim = 12, path = "data/CHELSA"), format = "file"), 
  tar_target(map_per_plot, get_chelsa_map_per_plot(chelsa_files, plot_coord)), 
  
  
  
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##' - FIGURES
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
 
  # - Plot growth and browsing with an effect of elevation
  tar_target(fig_gr_elev, plot_gr_elev(
    subset(df_gr, Site != "Gallivare"), dir.in = "fig/elevation"), format = "file"),
  tar_target(fig_br_elev, plot_br_elev(
    subset(df_br, Site != "Gallivare"), dir.in = "fig/elevation"), format = "file"),
  
  # - Plot growth and browsing with elevation and Gallivare
  tar_target(fig_gr_elev_withG, plot_gr_elev(
    df_gr, dir.in = "fig/elevation_withG"), format = "file"),
  tar_target(fig_br_elev_withG, plot_br_elev(
    df_br, dir.in = "fig/elevation_withG"), format = "file"), 
  
  # - Plot the elevation and climatic gradient in each site
  tar_target(plot_Zclimate, plot_climate(
    df_br, map_per_plot, "fig/supplementary/climate.jpg"), format = "file"), 
  
  # - Plot the height distribution of seedlings
  tar_target(fig_height, plot_height(df_br, "fig/supplementary/height.jpg"), 
             format = "file"), 
  
  # - Map of the study sites
  tar_target(fig_map, map_sites(df_br, "fig/map"), format = "file")
  
  
)