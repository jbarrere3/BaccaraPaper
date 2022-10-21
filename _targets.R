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
packages.in <- c("dplyr", "ggplot2", "RCurl", "httr", "tidyr", "data.table", "lme4", "cowplot", "archive", "terra", "knitr")
for(i in 1:length(packages.in)) if(!(packages.in[i] %in% rownames(installed.packages()))) install.packages(packages.in[i])
# Targets options
options(tidyverse.quiet = TRUE)
tar_option_set(packages = packages.in)
set.seed(2)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Targets workflow --------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

list(
  # - Raw files
  tar_target(upi_file, "data/field/UPI.csv", format = "file"), 
  tar_target(bdd_file, "data/field/JdD_MB.csv", format = "file"), 
  
  # - Read raw files
  tar_target(upi, fread(upi_file, h = T, sep = ";", dec = ",")), 
  tar_target(bdd, read.table(bdd_file, h = T, sep = ";")), 
  
  # - Format raw files
  tar_target(df_br, format_browsing_rate(bdd, upi)), 
  tar_target(df_gr, format_growth(df_br)), 
  
  # - Get precipitation from chelsa bioclim
  tar_target(plot_coord, extract_coord(df_gr)), 
  tar_target(chelsa_files, download_CHELSA(bioclim = 12, path = "data/CHELSA"), format = "file"), 
  tar_target(map_per_plot, get_chelsa_map_per_plot(chelsa_files, plot_coord)), 
  
  # Get ungulate mass density index (umdi) for each plot
  tar_target(ungulate_files, list.files(path = "data/ungulates", full.names = TRUE), format = "file"),
  tar_target(umdi_per_plot, get_umdi_per_plot(plot_coord, ungulate_files)),
  
  # - Fit and plot models for browsing and growth with Swedish site
  # --- Model fit
  tar_target(browsing_models, fit_browsing(df_br, map_per_plot)),
  # --- Figure with one model including both map and temperature
  tar_target(fig_browsing, plot_browsingproba(
    df_br, map_per_plot, browsing_models, "fig/browsing.pdf"), format = "file"),
  # --- Figure with separate models for map and temperature
  tar_target(fig_browsing_sep, plot_browsingproba_sep(
    df_br, map_per_plot, browsing_models, "fig/browsing_sep.pdf"), format = "file"), 
  
  # - Fit and plot models for browsing without Swedish site
  # --- Model fit
  tar_target(browsing_models_noG, fit_browsing(subset(df_br, Site != "Gallivare"), map_per_plot)),
  # --- Figure with one model including both map and temperature
  tar_target(fig_browsing_noG, plot_browsingproba(
    subset(df_br, Site != "Gallivare"), map_per_plot, browsing_models_noG, "fig/browsing_noG.pdf"), format = "file"),
  # --- Figure with separate models for map and temperature
  tar_target(fig_browsing_sep_noG, plot_browsingproba_sep(
    subset(df_br, Site != "Gallivare"), map_per_plot, browsing_models_noG, "fig/browsing_sep_noG.pdf"), format = "file"), 
  
  # - Fit and plot models for browsing with Swedish site and winter temperature instead of mean
  # --- Model fit
  tar_target(browsing_models_hiv, fit_browsing_hiv(df_br, map_per_plot)),
  # --- Figure with one model including both map and temperature
  tar_target(fig_browsing_hiv, plot_browsingproba_hiv(
    df_br, map_per_plot, browsing_models_hiv, "fig/browsing_hiv.pdf"), format = "file"), 
  
  # - Fit and plot model for growth without Swedish site
  # --- Model fit
  tar_target(growth_models_noG, fit_growth(subset(df_gr, Site != "Gallivare"), map_per_plot)), 
  # --- Figure with one model including both map and temperature
  tar_target(fig_growth_noG, plot_growth_ms(
    subset(df_gr, Site != "Gallivare"), map_per_plot, growth_models_noG, "fig/growth_noG.pdf"), format = "file"),
  
  # - Plot climatic variables
  tar_target(fig_tmean_vs_map, plot_climatic_var(df_gr, map_per_plot, "Tmean", "map", "fig/tmean_vs_map.pdf")),
  tar_target(fig_tmean_vs_thiv, plot_climatic_var(df_gr, map_per_plot, "Tmean", "Tm_hiv", "fig/tmean_vs_thiv.pdf")),
  
  # - Knit the Rmarkdown document presenting the analyses
  tar_target(rmd_file, "analysis.Rmd", format = "file"), 
  tar_target(fig_files, c(fig_browsing, fig_tmean_vs_map, fig_tmean_vs_thiv, fig_browsing_sep, fig_browsing_hiv), format = "file"), 
  tar_target(knitted_pdf, knit_rmd(rmd_file, fig_files), format = "file")
  
)