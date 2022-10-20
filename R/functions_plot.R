#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_plot.R  
#' @description R script containing all functions relative to data
#               importation and formatting
#' @author Julien Barrere, Marianne Bernard
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





#' Plot the effect of temperature and mapipitation on browsing probability
#' @param df_br Dataset for browsing
#' @param map_per_plot mapippitation per plot
#' @param browsing_models Fitted models for browsing
#' @param file.in Name (including path) of the file to save
plot_browsingproba <- function(df_br, map_per_plot, browsing_models, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Identify the species for which to plot the model outputs
  species.in <- names(browsing_models[[2]])
  
  # Add mapipitation to the dataset
  df_br <- df_br %>%
    left_join(map_per_plot, by = c("Site", "Plot", "Subplot"))
  
  # Initialize the basis of the dataset for plotting
  data.plot.0 <- expand.grid(Tmean = seq(quantile(df_br$Tmean, 0.05), quantile(df_br$Tmean, 0.95), length.out = 100), 
                             map = mean(df_br$map, na.rm = T), 
                             HT = as.numeric(quantile(df_br$HT, probs = c(0.1, 0.5, 0.9)))) %>%
    rbind.data.frame(expand.grid(Tmean = mean(df_br$Tmean, na.rm = T), 
                                 map = seq(quantile(df_br$map, 0.05), quantile(df_br$map, 0.95), length.out = 100), 
                                 HT = as.numeric(quantile(df_br$HT, probs = c(0.1, 0.5, 0.9))))) %>%
    # Add scaled variables
    mutate(Tmean.scaled = predict(lm(scale(Tmean) ~ Tmean, data = df_br), newdata = .), 
           map.scaled = predict(lm(scale(map) ~ map, data = df_br), newdata = .), 
           HT.scaled = predict(lm(scale(HT) ~ HT, data = df_br), newdata = .))
  
  # Loop on all species
  for(i in 1:length(species.in)){
    
    ## - To plot the predicted probabilities
    
    # Data for species i
    data.plot.i <- data.plot.0  %>%
      # Add mean and confidence interval
      left_join(predict.mbr(model = browsing_models[[2]][[i]][[6]], model.type = "all", data = .), 
                by = c("Tmean", "map", "HT")) %>%
      mutate(species = species.in[i])
    
    # Add to the final dataset
    if(i == 1) data.plot <- data.plot.i
    else data.plot <- rbind(data.plot, data.plot.i)
    
    
    ## - To plot the effect of each variable
    
    # Parameters for species i
    param.confint.i <- data.frame(
      param = c("Int.", "T", "Ht", "P", "T:Ht", "P:Ht"), 
      coef = as.numeric(coefficients(summary(browsing_models[[2]][[i]][[6]]))[, 1]), 
      low = as.numeric(confint(browsing_models[[2]][[i]][[6]], method = "Wald")[-c(1, 2), 1]),
      high = as.numeric(confint(browsing_models[[2]][[i]][[6]], method = "Wald")[-c(1, 2), 2]), 
      species = species.in[i]
    )
    
    # Add to the final param dataset
    if(i == 1) param.confint = param.confint.i
    else param.confint = rbind(param.confint, param.confint.i)
  }
  
  # Plot for temperature
  plot.temp <- data.plot %>%
    left_join(data.frame(species = c("ABAL", "ACPS", "FASY", "PIAB"), 
                         sp = c("A. alba", "A. pseudoplatanus", "F. sylvatica", "P. abies")), 
              by = "species") %>%
    filter(map == mean(df_br$map, na.rm = TRUE)) %>%
    mutate(height = paste0("Height = ", HT, " cm")) %>%
    mutate(height = factor(height, levels = paste0("Height = ", unique(.$HT)[order(unique(.$HT))], " cm"))) %>%
    ggplot(aes(x = Tmean, y = fit.mean, group = height, color = height)) + 
    geom_line() + 
    geom_ribbon(aes(ymin = fit.low, ymax = fit.high, fill = height), color = NA, alpha = 0.2) +
    facet_wrap(~ sp, nrow = 1) + 
    scale_color_manual(values = c("#335C67", "#E09F3E", "#9E2A2B")) +
    scale_fill_manual(values = c("#335C67", "#E09F3E", "#9E2A2B")) +
    ylab("Browsing\n probability") + xlab(expression(atop(paste("Mean Annual Temperature (", degree,"C)")))) +
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          legend.position = "none")
  
  # Plot for mapipitation
  plot.map <- data.plot %>%
    left_join(data.frame(species = c("ABAL", "ACPS", "FASY", "PIAB"), 
                         sp = c("A. alba", "A. pseudoplatanus", "F. sylvatica", "P. abies")), 
              by = "species") %>%
    filter(Tmean == mean(df_br$Tmean, na.rm = TRUE)) %>%
    mutate(height = paste0("Height = ", HT, " cm")) %>%
    mutate(height = factor(height, levels = paste0("Height = ", unique(.$HT)[order(unique(.$HT))], " cm"))) %>%
    ggplot(aes(x = map, y = fit.mean, group = height, color = height)) + 
    geom_line() + 
    geom_ribbon(aes(ymin = fit.low, ymax = fit.high, fill = height), color = NA, alpha = 0.2) +
    facet_wrap(~ sp, nrow = 1) + 
    scale_color_manual(values = c("#335C67", "#E09F3E", "#9E2A2B")) +
    scale_fill_manual(values = c("#335C67", "#E09F3E", "#9E2A2B")) +
    ylab("Browsing\n probability") + xlab("Annual precipitation (mm)") +
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          legend.position = "none", 
          axis.text.x = element_text(size = 7))
  
  # Plot the effect of each variable
  plot.param <- param.confint  %>%
    left_join(data.frame(species = c("ABAL", "ACPS", "FASY", "PIAB"), 
                         sp = c("A. alba", "A. pseudoplatanus", "F. sylvatica", "P. abies")), 
              by = "species") %>%
    mutate(param = factor(param, levels = c("P:Ht", "P", "T:Ht", "T", "Ht", "Int."))) %>%
    mutate(signif = ifelse((low > 0 | high < 0), "Yes", "No")) %>%
    ggplot(aes(x = param, y = coef, color = signif)) + 
    geom_point() +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0) +
    scale_color_manual(values = c("black", "red")) +
    facet_wrap(~ sp, nrow = 1) + 
    geom_hline(yintercept = 0, linetype = "dashed") +
    ylab("Effect on browsing probability") + xlab("") +
    coord_flip() +
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          legend.position = "none")
  
  # Final plot
  plot.out <- plot_grid(plot_grid(plot.temp, plot.map, plot.param, ncol = 1, align = "v", labels = c("(a)", "(b)", "(c)")), 
                        get_legend(plot.temp + theme(legend.position = "left", legend.title = element_blank())), 
                        nrow = 1, rel_widths = c(1, 0.3))
  
  # Save the plot
  ggsave(file.in, plot = plot.out, width = 18, height = 15, units = "cm", dpi = 600)
  
  # Return plot
  return(file.in)
}


#' Plot the effect of WINTER temperature and precipitation on browsing probability
#' @param df_br Dataset for browsing
#' @param map_per_plot precipitation per plot
#' @param browsing_models Fitted models for browsing
#' @param file.in Name (including path) of the file to save
plot_browsingproba_hiv <- function(df_br, map_per_plot, browsing_models, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Identify the species for which to plot the model outputs
  species.in <- names(browsing_models[[2]])
  
  # Add mapipitation to the dataset
  df_br <- df_br %>%
    left_join(map_per_plot, by = c("Site", "Plot", "Subplot"))
  
  # Initialize the basis of the dataset for plotting
  data.plot.0 <- expand.grid(Tm_hiv = seq(quantile(df_br$Tm_hiv, 0.05), quantile(df_br$Tm_hiv, 0.95), length.out = 100), 
                             map = mean(df_br$map, na.rm = T), 
                             HT = as.numeric(quantile(df_br$HT, probs = c(0.1, 0.5, 0.9)))) %>%
    rbind.data.frame(expand.grid(Tm_hiv = mean(df_br$Tm_hiv, na.rm = T), 
                                 map = seq(quantile(df_br$map, 0.05), quantile(df_br$map, 0.95), length.out = 100), 
                                 HT = as.numeric(quantile(df_br$HT, probs = c(0.1, 0.5, 0.9))))) %>%
    # Add scaled variables
    mutate(Tm_hiv.scaled = predict(lm(scale(Tm_hiv) ~ Tm_hiv, data = df_br), newdata = .), 
           map.scaled = predict(lm(scale(map) ~ map, data = df_br), newdata = .), 
           HT.scaled = predict(lm(scale(HT) ~ HT, data = df_br), newdata = .))
  
  # Loop on all species
  for(i in 1:length(species.in)){
    
    ## - To plot the predicted probabilities
    
    # Data for species i
    data.plot.i <- data.plot.0  %>%
      # Add mean and confidence interval
      left_join(predict.mbr_hiv(model = browsing_models[[2]][[i]][[6]], model.type = "all", data = .), 
                by = c("Tm_hiv", "map", "HT")) %>%
      mutate(species = species.in[i])
    
    # Add to the final dataset
    if(i == 1) data.plot <- data.plot.i
    else data.plot <- rbind(data.plot, data.plot.i)
    
    
    ## - To plot the effect of each variable
    
    # Parameters for species i
    param.confint.i <- data.frame(
      param = c("Int.", "T", "Ht", "P", "T:Ht", "P:Ht"), 
      coef = as.numeric(coefficients(summary(browsing_models[[2]][[i]][[6]]))[, 1]), 
      low = as.numeric(confint(browsing_models[[2]][[i]][[6]], method = "Wald")[-c(1, 2), 1]),
      high = as.numeric(confint(browsing_models[[2]][[i]][[6]], method = "Wald")[-c(1, 2), 2]), 
      species = species.in[i]
    )
    
    # Add to the final param dataset
    if(i == 1) param.confint = param.confint.i
    else param.confint = rbind(param.confint, param.confint.i)
  }
  
  # Plot for temperature
  plot.temp <- data.plot %>%
    left_join(data.frame(species = c("ABAL", "ACPS", "FASY", "PIAB"), 
                         sp = c("A. alba", "A. pseudoplatanus", "F. sylvatica", "P. abies")), 
              by = "species") %>%
    filter(map == mean(df_br$map, na.rm = TRUE)) %>%
    mutate(height = paste0("Height = ", HT, " cm")) %>%
    mutate(height = factor(height, levels = paste0("Height = ", unique(.$HT)[order(unique(.$HT))], " cm"))) %>%
    ggplot(aes(x = Tm_hiv, y = fit.mean, group = height, color = height)) + 
    geom_line() + 
    geom_ribbon(aes(ymin = fit.low, ymax = fit.high, fill = height), color = NA, alpha = 0.2) +
    facet_wrap(~ sp, nrow = 1) + 
    scale_color_manual(values = c("#335C67", "#E09F3E", "#9E2A2B")) +
    scale_fill_manual(values = c("#335C67", "#E09F3E", "#9E2A2B")) +
    ylab("Browsing\n probability") + xlab(expression(atop(paste("Mean Winter Temperature (", degree,"C)")))) +
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          legend.position = "none")
  
  # Plot for mapipitation
  plot.map <- data.plot %>%
    left_join(data.frame(species = c("ABAL", "ACPS", "FASY", "PIAB"), 
                         sp = c("A. alba", "A. pseudoplatanus", "F. sylvatica", "P. abies")), 
              by = "species") %>%
    filter(Tm_hiv == mean(df_br$Tm_hiv, na.rm = TRUE)) %>%
    mutate(height = paste0("Height = ", HT, " cm")) %>%
    mutate(height = factor(height, levels = paste0("Height = ", unique(.$HT)[order(unique(.$HT))], " cm"))) %>%
    ggplot(aes(x = map, y = fit.mean, group = height, color = height)) + 
    geom_line() + 
    geom_ribbon(aes(ymin = fit.low, ymax = fit.high, fill = height), color = NA, alpha = 0.2) +
    facet_wrap(~ sp, nrow = 1) + 
    scale_color_manual(values = c("#335C67", "#E09F3E", "#9E2A2B")) +
    scale_fill_manual(values = c("#335C67", "#E09F3E", "#9E2A2B")) +
    ylab("Browsing\n probability") + xlab("Annual precipitation (mm)") +
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          legend.position = "none", 
          axis.text.x = element_text(size = 7))
  
  # Plot the effect of each variable
  plot.param <- param.confint  %>%
    left_join(data.frame(species = c("ABAL", "ACPS", "FASY", "PIAB"), 
                         sp = c("A. alba", "A. pseudoplatanus", "F. sylvatica", "P. abies")), 
              by = "species") %>%
    mutate(param = factor(param, levels = c("P:Ht", "P", "T:Ht", "T", "Ht", "Int."))) %>%
    mutate(signif = ifelse((low > 0 | high < 0), "Yes", "No")) %>%
    ggplot(aes(x = param, y = coef, color = signif)) + 
    geom_point() +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0) +
    scale_color_manual(values = c("black", "red")) +
    facet_wrap(~ sp, nrow = 1) + 
    geom_hline(yintercept = 0, linetype = "dashed") +
    ylab("Effect on browsing probability") + xlab("") +
    coord_flip() +
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          legend.position = "none")
  
  # Final plot
  plot.out <- plot_grid(plot_grid(plot.temp, plot.map, plot.param, ncol = 1, align = "v", labels = c("(a)", "(b)", "(c)")), 
                        get_legend(plot.temp + theme(legend.position = "left", legend.title = element_blank())), 
                        nrow = 1, rel_widths = c(1, 0.3))
  
  # Save the plot
  ggsave(file.in, plot = plot.out, width = 18, height = 15, units = "cm", dpi = 600)
  
  # Return plot
  return(file.in)
}



#' Plot the effect of temperature and mapipitation on browsing probability with two separate models
#' @param df_br Dataset for browsing
#' @param map_per_plot mapippitation per plot
#' @param browsing_models Fitted models for browsing
#' @param file.in Name (including path) of the file to save
plot_browsingproba_sep <- function(df_br, map_per_plot, browsing_models, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Identify the species for which to plot the model outputs
  species.in <- names(browsing_models[[2]])
  
  # Add mapipitation to the dataset
  df_br <- df_br %>%
    left_join(map_per_plot, by = c("Site", "Plot", "Subplot"))
  
  # Initialize the basis of the dataset for temperature
  data.temp.0 <- expand.grid(Tmean = seq(quantile(df_br$Tmean, 0.05), quantile(df_br$Tmean, 0.95), length.out = 100), 
                             map = mean(df_br$map, na.rm = T), 
                             HT = as.numeric(quantile(df_br$HT, probs = c(0.1, 0.5, 0.9)))) %>%
    # Add scaled variables
    mutate(Tmean.scaled = predict(lm(scale(Tmean) ~ Tmean, data = df_br), newdata = .), 
           map.scaled = predict(lm(scale(map) ~ map, data = df_br), newdata = .), 
           HT.scaled = predict(lm(scale(HT) ~ HT, data = df_br), newdata = .))
  
  # Initialize the basis of the dataset for precipitation
  data.map.0 <- expand.grid(Tmean = mean(df_br$Tmean, na.rm = T), 
                            map = seq(quantile(df_br$map, 0.05), quantile(df_br$map, 0.95), length.out = 100), 
                            HT = as.numeric(quantile(df_br$HT, probs = c(0.1, 0.5, 0.9)))) %>%
    # Add scaled variables
    mutate(Tmean.scaled = predict(lm(scale(Tmean) ~ Tmean, data = df_br), newdata = .), 
           map.scaled = predict(lm(scale(map) ~ map, data = df_br), newdata = .), 
           HT.scaled = predict(lm(scale(HT) ~ HT, data = df_br), newdata = .))
  
  # Loop on all species
  for(i in 1:length(species.in)){
    
    ## - To plot the predicted probabilities
    
    # Temperature data for species i
    data.temp.i <- data.temp.0  %>%
      # Add mean and confidence interval
      left_join(predict.mbr(model = browsing_models[[2]][[i]][[4]], model.type = "temp", data = .), 
                by = c("Tmean", "map", "HT")) %>%
      mutate(species = species.in[i])
    
    # Precipitation data for species i
    data.map.i <- data.map.0  %>%
      # Add mean and confidence interval
      left_join(predict.mbr(model = browsing_models[[2]][[i]][[5]], model.type = "prec", data = .), 
                by = c("Tmean", "map", "HT")) %>%
      mutate(species = species.in[i])
    
    # Parameters for species i
    param.confint.i <- data.frame(
      param = c("Int.", "T", "Ht", "T:Ht"), 
      coef = as.numeric(coefficients(summary(browsing_models[[2]][[i]][[4]]))[, 1]), 
      low = as.numeric(confint(browsing_models[[2]][[i]][[4]], method = "Wald")[-c(1, 2), 1]),
      high = as.numeric(confint(browsing_models[[2]][[i]][[4]], method = "Wald")[-c(1, 2), 2]), 
      species = species.in[i], 
      model.type = "temp.\n model"
    ) %>%
      rbind.data.frame(data.frame(
        param = c("Int.", "P", "Ht", "P:Ht"), 
        coef = as.numeric(coefficients(summary(browsing_models[[2]][[i]][[5]]))[, 1]), 
        low = as.numeric(confint(browsing_models[[2]][[i]][[5]], method = "Wald")[-c(1, 2), 1]),
        high = as.numeric(confint(browsing_models[[2]][[i]][[5]], method = "Wald")[-c(1, 2), 2]), 
        species = species.in[i], 
        model.type = "prec.\n model"
      ))
    
    # Add to the final dataset
    if(i == 1){
      data.temp <- data.temp.i
      data.map = data.map.i
      param.confint = param.confint.i
    } else {
      data.temp <- rbind(data.temp, data.temp.i)
      data.map <- rbind(data.map, data.map.i)
      param.confint = rbind(param.confint, param.confint.i)
    } 
    
    
  }
  
  # Plot for temperature
  plot.temp <- data.temp %>%
    left_join(data.frame(species = c("ABAL", "ACPS", "FASY", "PIAB"), 
                         sp = c("A. alba", "A. pseudoplatanus", "F. sylvatica", "P. abies")), 
              by = "species") %>%
    mutate(height = paste0("Height = ", HT, " cm")) %>%
    mutate(height = factor(height, levels = paste0("Height = ", unique(.$HT)[order(unique(.$HT))], " cm"))) %>%
    ggplot(aes(x = Tmean, y = fit.mean, group = height, color = height)) + 
    geom_line() + 
    geom_ribbon(aes(ymin = fit.low, ymax = fit.high, fill = height), color = NA, alpha = 0.2) +
    facet_wrap(~ sp, nrow = 1) + 
    scale_color_manual(values = c("#335C67", "#E09F3E", "#9E2A2B")) +
    scale_fill_manual(values = c("#335C67", "#E09F3E", "#9E2A2B")) +
    ylab("Browsing\n probability") + xlab(expression(atop(paste("Mean Annual Temperature (", degree,"C)")))) +
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          legend.position = "none")
  
  # Plot for precipitation
  plot.map <- data.map %>%
    left_join(data.frame(species = c("ABAL", "ACPS", "FASY", "PIAB"), 
                         sp = c("A. alba", "A. pseudoplatanus", "F. sylvatica", "P. abies")), 
              by = "species") %>%
    mutate(height = paste0("Height = ", HT, " cm")) %>%
    mutate(height = factor(height, levels = paste0("Height = ", unique(.$HT)[order(unique(.$HT))], " cm"))) %>%
    ggplot(aes(x = map, y = fit.mean, group = height, color = height)) + 
    geom_line() + 
    geom_ribbon(aes(ymin = fit.low, ymax = fit.high, fill = height), color = NA, alpha = 0.2) +
    facet_wrap(~ sp, nrow = 1) + 
    scale_color_manual(values = c("#335C67", "#E09F3E", "#9E2A2B")) +
    scale_fill_manual(values = c("#335C67", "#E09F3E", "#9E2A2B")) +
    ylab("Browsing\n probability") + xlab("Annual precipitation (mm)") +
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          legend.position = "none", 
          axis.text.x = element_text(size = 7))
  
  # Plot the effect of each variable
  plot.param <- param.confint  %>%
    left_join(data.frame(species = c("ABAL", "ACPS", "FASY", "PIAB"), 
                         sp = c("A. alba", "A. pseudoplatanus", "F. sylvatica", "P. abies")), 
              by = "species") %>%
    filter(param != "Int.") %>%
    mutate(signif = ifelse((low > 0 | high < 0), "Yes", "No")) %>%
    ggplot(aes(x = param, y = coef, color = signif)) + 
    geom_point() +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0) +
    scale_color_manual(values = c("black", "red")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ylab("Effect on browsing probability") + xlab("") +
    facet_grid(model.type ~ sp, scales = "free_y") +
    coord_flip() +
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          legend.position = "none")
  
  # Final plot
  plot.out <- plot_grid(plot_grid(plot.temp, plot.map, plot.param, ncol = 1, align = "v", labels = c("(a)", "(b)", "(c)")), 
                        get_legend(plot.temp + theme(legend.position = "left", legend.title = element_blank())), 
                        nrow = 1, rel_widths = c(1, 0.3))
  
  # Save the plot
  ggsave(file.in, plot = plot.out, width = 18, height = 15, units = "cm", dpi = 600)
  
  # Return plot
  return(file.in)
}





#' Function to plot two climatic variables together
#' @param df_gr data set to fit growth models
#' @param map_per_plot mean annual precipitation per plot
#' @param var.1 climatic variable to plot ("map", "Tmean" or "Tm_hiv)
#' @param var.2 climatic variable to plot ("map", "Tmean" or "Tm_hiv)
#' @param file.in name (including path) of the file to save
plot_climatic_var <- function(df_gr, map_per_plot, var.1, var.2, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Make the plot
  plot.out <- map_per_plot %>%
    left_join((df_gr %>% 
                 dplyr::select("Site", "Plot", "Subplot", "Tmean", "Tm_hiv") %>%
                 distinct()), 
              by = c("Site", "Plot", "Subplot")) %>%
    dplyr::select("Site", "var.1" = var.1, "var.2" = var.2) %>%
    ggplot(aes(x = var.1, y = var.2, color = Site)) + 
    geom_point() + 
    xlab(var.1) + ylab(var.2) +
    theme_bw()
  
  # Save the plot
  ggsave(file.in, plot = plot.out, width = 12, height = 8, units = "cm", dpi = 600)
  
  # Return plot
  return(file.in)
}




#' Plot the effect of temperature and map on growth
#' @param df_gr Data set for growth
#' @param map_per_plot mapippitation per plot
#' @param browsing_models Fitted models for browsing
#' @param file.in Name (including path) of the file to save
plot_growth_ms <- function(df_gr, map_per_plot, growth_models, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Identify the species for which to plot the model outputs
  species.in <- names(growth_models[[2]])
  
  # Add precipitation to the dataset
  df_gr <- df_gr %>% left_join(map_per_plot, by = c("Site", "Plot", "Subplot"))
  
  
  # Loop on all species
  for(i in 1:length(species.in)){
    
    ## - To plot the predicted probabilities
    
    # Dataset for species i
    df_gr.i <- subset(df_gr, Species == species.in[i])
    
    # Data for species i
    data.plot.i <- expand.grid(Tmean = seq(quantile(df_gr.i$Tmean, 0.05), quantile(df_gr.i$Tmean, 0.95), length.out = 100), 
                               map = mean(df_gr.i$map, na.rm = T), 
                               HT = mean(df_gr.i$HT, na.rm = TRUE), 
                               B = c(0, 1)) %>%
      rbind.data.frame(expand.grid(Tmean = mean(df_gr.i$Tmean, na.rm = T), 
                                   map = seq(quantile(df_gr.i$map, 0.05), quantile(df_gr.i$map, 0.95), length.out = 100), 
                                   HT = mean(df_gr.i$HT, na.rm = TRUE), 
                                   B = c(0, 1))) %>%
      # Add scaled variables
      mutate(Tmean.scaled = predict(lm(scale(Tmean) ~ Tmean, data = df_gr.i), newdata = .), 
             map.scaled = predict(lm(scale(map) ~ map, data = df_gr.i), newdata = .), 
             HT.scaled = predict(lm(scale(HT) ~ HT, data = df_gr.i), newdata = .))  %>%
      # Add mean and confidence interval
      left_join(predict.mgr(model = growth_models[[2]][[i]][[3]], model.type = "all", data = .), 
                by = c("Tmean", "map", "HT", "B")) %>%
      mutate(species = species.in[i])
    
    # Add to the final dataset
    if(i == 1) data.plot <- data.plot.i
    else data.plot <- rbind(data.plot, data.plot.i)
    
    
    ## - To plot the effect of each variable
    
    # Parameters for species i
    param.confint.i <- data.frame(
      param = c("Int.", "Ht", "T", "Br", "P", "T:Br", "P:Br"), 
      coef = as.numeric(coefficients(summary(growth_models[[2]][[i]][[3]]))[, 1]), 
      low = as.numeric(confint(growth_models[[2]][[i]][[3]], method = "Wald")[-c(1:3), 1]),
      high = as.numeric(confint(growth_models[[2]][[i]][[3]], method = "Wald")[-c(1:3), 2]), 
      species = species.in[i]
    )
    
    # Add to the final param dataset
    if(i == 1) param.confint = param.confint.i
    else param.confint = rbind(param.confint, param.confint.i)
  }
  
  # Raw data to plot for temperature
  data.real.temp <- df_gr %>%
    mutate(Tmean = round(Tmean, digits = 0)) %>%
    mutate(browsing = ifelse(B == 1, "browsed", "unbrowsed")) %>%
    mutate(Tmean = ifelse(browsing == "browsed", Tmean + 0.15, Tmean - 0.15)) %>%
    group_by(Species) %>%
    mutate(sp.HTlow = quantile(HT, probs = 0.1), 
           sp.HThigh = quantile(HT, probs = 0.9)) %>%
    filter(HT >= sp.HTlow & HT <= sp.HThigh) %>%
    ungroup() %>%
    left_join(data.frame(Species = c("ABAL", "ACPS", "FASY", "PIAB"), 
                         sp = c("A. alba", "A. pseudoplatanus", "F. sylvatica", "P. abies")), 
              by = "Species") %>%
    group_by(sp, Tmean, browsing) %>%
    summarize(fit.mean = mean(LLS, na.rm = TRUE), 
              fit.sd = sd(LLS, na.rm = TRUE), 
              n = n()) %>%
    mutate(label = paste0("(", n, ")")) %>%
    filter(n > 9)
  
  # Raw data to plot for precipitation
  data.real.map <- df_gr %>%
    mutate(map = round(map/100, digits = 0)*100) %>%
    mutate(browsing = ifelse(B == 1, "browsed", "unbrowsed")) %>%
    mutate(map = ifelse(browsing == "browsed", map + 15, map - 15)) %>%
    group_by(Species) %>%
    mutate(sp.HTlow = quantile(HT, probs = 0.1), 
           sp.HThigh = quantile(HT, probs = 0.9)) %>%
    filter(HT >= sp.HTlow & HT <= sp.HThigh) %>%
    ungroup() %>%
    left_join(data.frame(Species = c("ABAL", "ACPS", "FASY", "PIAB"), 
                         sp = c("A. alba", "A. pseudoplatanus", "F. sylvatica", "P. abies")), 
              by = "Species") %>%
    group_by(sp, map, browsing) %>%
    summarize(fit.mean = mean(LLS, na.rm = TRUE), 
              fit.sd = sd(LLS, na.rm = TRUE), 
              n = n()) %>%
    mutate(label = paste0("(", n, ")")) %>%
    filter(n > 9)
  
  # Plot for temperature
  plot.temp <- data.plot %>%
    left_join(data.frame(species = c("ABAL", "ACPS", "FASY", "PIAB"), 
                         sp = c("A. alba", "A. pseudoplatanus", "F. sylvatica", "P. abies")), 
              by = "species") %>%
    filter(map %in% (df_gr %>% group_by(Species) %>% summarize(p = mean(map)))$p) %>%
    mutate(browsing = ifelse(B == 1, "browsed", "unbrowsed")) %>%
    ggplot(aes(x = Tmean, y = fit.mean, group = browsing, color = browsing)) + 
    geom_line() + 
    geom_ribbon(aes(ymin = fit.low, ymax = fit.high, fill = browsing), color = NA, alpha = 0.2) + 
    geom_point(data = data.real.temp, inherit.aes = TRUE) + 
    geom_errorbar(data = data.real.temp, inherit.aes = TRUE, width = 0, 
                  aes(ymin = fit.mean - fit.sd, ymax = fit.mean + fit.sd)) +
    facet_wrap(~ sp, nrow = 1) + 
    scale_color_manual(values = c("#335C67", "#9E2A2B")) +
    scale_fill_manual(values = c("#335C67", "#9E2A2B")) +
    ylab("Last shoot\n length (cm)") + xlab(expression(atop(paste("Mean Annual Temperature (", degree,"C)")))) +
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          legend.position = "none")
  
  # Plot for mapipitation
  plot.map <- data.plot %>%
    left_join(data.frame(species = c("ABAL", "ACPS", "FASY", "PIAB"), 
                         sp = c("A. alba", "A. pseudoplatanus", "F. sylvatica", "P. abies")), 
              by = "species") %>%
    filter(Tmean %in% (df_gr %>% group_by(Species) %>% summarize(p = mean(Tmean)))$p) %>%
    mutate(browsing = ifelse(B == 1, "browsed", "unbrowsed")) %>%
    ggplot(aes(x = map, y = fit.mean, group = browsing, color = browsing)) + 
    geom_line() + 
    geom_ribbon(aes(ymin = fit.low, ymax = fit.high, fill = browsing), color = NA, alpha = 0.2) + 
    geom_point(data = data.real.map, inherit.aes = TRUE) + 
    geom_errorbar(data = data.real.map, inherit.aes = TRUE, width = 0, 
                  aes(ymin = fit.mean - fit.sd, ymax = fit.mean + fit.sd)) +
    facet_wrap(~ sp, nrow = 1) + 
    scale_color_manual(values = c("#335C67", "#9E2A2B")) +
    scale_fill_manual(values = c("#335C67", "#9E2A2B")) +
    ylab("Last shoot\n length (cm)") + xlab("Annual precipitation (mm)") +
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          legend.position = "none", 
          axis.text.x = element_text(size = 7))
  
  # Plot the effect of each variable
  plot.param <- param.confint  %>%
    left_join(data.frame(species = c("ABAL", "ACPS", "FASY", "PIAB"), 
                         sp = c("A. alba", "A. pseudoplatanus", "F. sylvatica", "P. abies")), 
              by = "species") %>%
    filter(param != "Int.") %>%
    mutate(param = factor(param, levels = c("P:Br", "P", "T:Br", "T", "Ht", "Br"))) %>%
    mutate(signif = ifelse((low > 0 | high < 0), "Yes", "No")) %>%
    ggplot(aes(x = param, y = coef, color = signif)) + 
    geom_point() +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0) +
    scale_color_manual(values = c("black", "red")) +
    facet_wrap(~ sp, nrow = 1) + 
    geom_hline(yintercept = 0, linetype = "dashed") +
    ylab("Effect on last shoot length") + xlab("") +
    coord_flip() +
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          legend.position = "none", 
          axis.text.x = element_text(size = 7))
  
  # Final plot
  plot.out <- plot_grid(plot_grid(plot.temp, plot.map, plot.param, ncol = 1, align = "v", labels = c("(a)", "(b)", "(c)")), 
                        get_legend(plot.temp + theme(legend.position = "left", legend.title = element_blank())), 
                        nrow = 1, rel_widths = c(1, 0.3))
  
  # Save the plot
  ggsave(file.in, plot = plot.out, width = 18, height = 15, units = "cm", dpi = 600)
  
  # Return plot
  return(file.in)
}


