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







#' Plot the effect of elevation on growth
#' @param df_gr Dataset for growth
#' @param dir.in Name (including path) of the directory where to save files
plot_gr_elev = function(df_gr, dir.in){
  
  # Name of the files to export
  file.fig.in = paste0(dir.in, "/fig_growth.jpg")
  file.table.in = paste0(dir.in, "/table_growth.csv")
  
  # Create directory if needed
  create_dir_if_needed(file.fig.in)
  
  # Vector of species present in the dataset
  sp <- unique(df_gr$Species)
  
  # Initialize final list containing the models
  list.out <- list()
  
  # Initialize the basis of the dataset for plotting
  data.plot.0 <- expand.grid(Z = seq(quantile(df_gr$Z, 0.05), quantile(df_gr$Z, 0.95), length.out = 100), 
                             HT = as.numeric(mean(df_gr$HT, na.rm = TRUE)), 
                             B = c(0, 1)) %>%
    # Add scaled variables
    mutate(Z.scaled = predict(lm(scale(Z) ~ Z, data = df_gr), newdata = .), 
           HT.scaled = predict(lm(scale(HT) ~ HT, data = df_gr), newdata = .))
  
  # Initialize the table with statistics 
  table.stat = data.frame(col1 = "", col2 = "Estimate (se)", col3 = "p-value")
  
  # Loop on all species
  for(i in 1:length(sp)){
    
    # Initialize model list for species i
    model.i = lmer(log(LLS+1) ~ scale(HT) + scale(Z)*B + (1|Site/Z.Classe),
                   data =subset(df_gr, Species == sp[i]))
    
    
    # Extract lower and upper CI for each param
    param.confint.i <- data.frame(param = letters[c(1:dim(coefficients(summary(model.i)))[1])], 
                                  low = as.numeric(confint(model.i, method = "Wald")[-c(1:3), 1]),
                                  high = as.numeric(confint(model.i, method = "Wald")[-c(1:3), 2]))
    
    # Initialize dataset for parameters
    param.i <- data.frame(sim = c(1:5000))
    
    # Add parameters value for each simulation
    for(j in 1:dim(param.confint.i)[1]){
      eval(parse(text = paste0("param.i$", letters[j], " <- rnorm(5000, mean(c(", 
                               param.confint.i$low[j], ", ", param.confint.i$high[j], 
                               ")), (", param.confint.i$high[j], " - ", param.confint.i$low[j], 
                               ")/(2*1.96))")))
    }
    
    # Calculate confidence interval and mean
    data.plot.i <- full_join(data.plot.0, data.frame(sim = c(1:5000)), by = character()) %>%
      left_join(param.i, by = "sim") %>%
      mutate(fit.raw = a + b*HT.scaled + c*Z.scaled + d*B + e*Z.scaled*B) %>%
      mutate(fit = exp(fit.raw) - 1) %>%
      group_by(Z, HT, B) %>%
      summarize(fit.low = quantile(fit, probs = 0.025), 
                fit.mean = mean(fit, na.rm = TRUE),
                fit.high = quantile(fit, probs = 0.975)) %>%
      mutate(species = sp[i])
    
    # Parameters for species i
    coef.est.i <- data.frame(
      param = c("Int.", "Ht", "Z", "B", "Z:B"), 
      coef = as.numeric(coefficients(summary(model.i))[, 1]), 
      se = as.numeric(coefficients(summary(model.i))[, 2]), 
      low = as.numeric(confint(model.i, method = "Wald")[-c(1:3), 1]),
      high = as.numeric(confint(model.i, method = "Wald")[-c(1:3), 2]), 
      p = as.numeric(c(NA_real_, car::Anova(model.i)[, 3])), 
      species = sp[i]
    )
    
    # Complete the stat table 
    table.stat = table.stat %>%
      rbind(data.frame(
        col1 = c("", "", "Int.", "Ht", "Z", "B", "Z:B"), 
        col2 = c("", sp[i], paste0(round(coef.est.i$coef, 2), " (", 
                                   round(coef.est.i$se, 2), ")")), 
        col3 = c("", "", "", scales::pvalue(coef.est.i$p[-1]))
      ))
    
    # Add to the final dataset
    if(i == 1){
      data.plot <- data.plot.i
      coef.est = coef.est.i
    } else {
      data.plot <- rbind(data.plot, data.plot.i)
      coef.est <- rbind(coef.est, coef.est.i)
    } 
    
  }
  
  # Raw data to plot for temperature
  data.real.Z <- df_gr %>%
    mutate(Z = round(Z/250, digits = 0)*250) %>%
    mutate(browsing = ifelse(B == 1, "browsed", "unbrowsed")) %>%
    mutate(Z = ifelse(browsing == "browsed", Z + 10, Z - 10)) %>%
    group_by(Species) %>%
    mutate(sp.HTlow = quantile(HT, probs = 0.1), 
           sp.HThigh = quantile(HT, probs = 0.9)) %>%
    filter(HT >= sp.HTlow & HT <= sp.HThigh) %>%
    ungroup() %>%
    left_join(data.frame(Species = c("ABAL", "ACPS", "FASY", "PIAB"), 
                         sp = c("A. alba", "A. pseudoplatanus", "F. sylvatica", "P. abies")), 
              by = "Species") %>%
    group_by(sp, Z, browsing) %>%
    summarize(fit.mean = mean(LLS, na.rm = TRUE), 
              fit.sd = sd(LLS, na.rm = TRUE), 
              n = n()) %>%
    mutate(label = paste0("(", n, ")")) %>%
    filter(n > 9)
  
  # Plot the predictions
  plot.prediction <- data.plot %>%
    left_join(data.frame(species = c("ABAL", "ACPS", "FASY", "PIAB"), 
                         sp = c("A. alba", "A. pseudoplatanus", "F. sylvatica", "P. abies")), 
              by = "species") %>%
    mutate(browsing = ifelse(B == 1, "browsed", "unbrowsed")) %>%
    ggplot(aes(x = Z, y = fit.mean, group = browsing, color = browsing)) + 
    geom_line() + 
    geom_ribbon(aes(ymin = fit.low, ymax = fit.high, fill = browsing), 
                color = NA, alpha = 0.2) + 
    geom_point(data = data.real.Z, inherit.aes = TRUE) + 
    geom_errorbar(data = data.real.Z, inherit.aes = TRUE, width = 0, 
                  aes(ymin = fit.mean - fit.sd, ymax = fit.mean + fit.sd)) +
    facet_wrap(~ sp, nrow = 1) + 
    scale_color_manual(values = c("#335C67", "#9E2A2B")) +
    scale_fill_manual(values = c("#335C67", "#9E2A2B")) +
    ylab("Last shoot\n length (cm)") + xlab("elevation (m)") +
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          axis.text.x = element_text(size = 5),
          strip.background = element_blank(), 
          legend.position = "none")
  
  # Plot the effect of each variable
  plot.param <- coef.est  %>%
    left_join(data.frame(species = c("ABAL", "ACPS", "FASY", "PIAB"), 
                         sp = c("A. alba", "A. pseudoplatanus", "F. sylvatica", "P. abies")), 
              by = "species") %>%
    filter(param != "Int.") %>%
    mutate(param = factor(param, levels = c("Ht", "Z", "B", "Z:B"))) %>%
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
          axis.text.x = element_text(size = 5),
          legend.position = "none")
  
  # Final plot
  plot.out <- plot_grid(plot_grid(plot.prediction, plot.param, ncol = 1, align = "v", labels = c("(a)", "(b)")), 
                        get_legend(plot.prediction + theme(legend.position = "left", legend.title = element_blank())), 
                        nrow = 1, rel_widths = c(1, 0.3))
  
  # Save the plot
  ggsave(file.fig.in, plot = plot.out, width = 18, height = 10, units = "cm", dpi = 600)
  
  # Save the stat table
  write.table(table.stat, file.table.in, row.names = F, col.names = F, quote = F)
  
  # Return plot
  return(c(file.fig.in, file.table.in))
}




#' Plot the effect of elevation on browsing probability
#' @param df_br Dataset for browsing
#' @param dir.in Name (including path) of the directory where to save files
plot_br_elev = function(df_br, dir.in){
  
  # Name of the files to export
  file.fig.in = paste0(dir.in, "/fig_browsing.jpg")
  file.table.in = paste0(dir.in, "/table_browsing.csv")
  
  # Create directory if needed
  create_dir_if_needed(file.fig.in)
  
  # Vector of species present in the dataset
  sp <- unique(df_br$Species)
  
  # Initialize final list containing the models
  list.out <- list()
  
  # Initialize the basis of the dataset for plotting
  data.plot.0 <- expand.grid(Z = seq(quantile(df_br$Z, 0.05), quantile(df_br$Z, 0.95), length.out = 100), 
                             HT = as.numeric(quantile(df_br$HT, probs = c(0.1, 0.5, 0.9)))) %>%
    # Add scaled variables
    mutate(Z.scaled = predict(lm(scale(Z) ~ Z, data = df_br), newdata = .), 
           HT.scaled = predict(lm(scale(HT) ~ HT, data = df_br), newdata = .))
  
  # Initialize the table with statistics 
  table.stat = data.frame(col1 = "", col2 = "Estimate (se)", col3 = "p-value")
  
  # Loop on all species
  for(i in 1:length(sp)){
    
    # Initialize model list for species i
    model.i = glmer(B ~ scale(Z)*scale(HT) + (1|Site/Z.Classe), family = 'binomial', 
                    data =subset(df_br, Species == sp[i]))
    
    
    # Extract lower and upper CI for each param
    param.confint.i <- data.frame(param = letters[c(1:dim(coefficients(summary(model.i)))[1])], 
                                  low = as.numeric(confint(model.i, method = "Wald")[-c(1, 2), 1]),
                                  high = as.numeric(confint(model.i, method = "Wald")[-c(1, 2), 2]))
    
    # Initialize dataset for parameters
    param.i <- data.frame(sim = c(1:5000))
    
    # Add parameters value for each simulation
    for(j in 1:dim(param.confint.i)[1]){
      eval(parse(text = paste0("param.i$", letters[j], " <- rnorm(5000, mean(c(", 
                               param.confint.i$low[j], ", ", param.confint.i$high[j], 
                               ")), (", param.confint.i$high[j], " - ", param.confint.i$low[j], 
                               ")/(2*1.96))")))
    }
    
    # Calculate confidence interval and mean
    data.plot.i <- full_join(data.plot.0, data.frame(sim = c(1:5000)), by = character()) %>%
      left_join(param.i, by = "sim") %>%
      mutate(fit = a + b*Z.scaled + c*HT.scaled + d*Z.scaled*HT.scaled) %>%
      mutate(fit = plogis(fit)) %>%
      group_by(Z, HT) %>%
      summarize(fit.low = quantile(fit, probs = 0.025), 
                fit.mean = mean(fit, na.rm = TRUE),
                fit.high = quantile(fit, probs = 0.975)) %>%
      mutate(species = sp[i])
    
    # Parameters for species i
    coef.est.i <- data.frame(
      param = c("Int.", "Z", "Ht", "Z:Ht"), 
      coef = as.numeric(coefficients(summary(model.i))[, 1]), 
      se = as.numeric(coefficients(summary(model.i))[, 2]), 
      low = as.numeric(confint(model.i, method = "Wald")[-c(1, 2), 1]),
      high = as.numeric(confint(model.i, method = "Wald")[-c(1, 2), 2]), 
      p = as.numeric(coefficients(summary(model.i))[, 4]), 
      species = sp[i]
    )
    
    # Identify the name of sp i
    species.i = (data.frame(species = sp[i]) %>%
                   left_join(data.frame(species = c("ABAL", "ACPS", "FASY", "PIAB"), 
                                        sp = c("A. alba", "A. pseudoplatanus", 
                                               "F. sylvatica", "P. abies")), 
                             by = "species"))$sp
    
    # Complete the stat table 
    table.stat = table.stat %>%
      rbind(data.frame(
        col1 = c("", "", "Int.", "Z", "Ht", "Z:Ht"), 
        col2 = c("", species.i, paste0(round(coef.est.i$coef, 2), " (", 
                                       round(coef.est.i$se, 2), ")")), 
        col3 = c("", "", scales::pvalue(coef.est.i$p))
      ))
    # Add to the final dataset
    if(i == 1){
      data.plot <- data.plot.i
      coef.est = coef.est.i
    } else {
      data.plot <- rbind(data.plot, data.plot.i)
      coef.est <- rbind(coef.est, coef.est.i)
    } 
    
  }
  
  
  
  # Plot for temperature
  plot.prediction <- data.plot %>%
    left_join(data.frame(species = c("ABAL", "ACPS", "FASY", "PIAB"), 
                         sp = c("A. alba", "A. pseudoplatanus", "F. sylvatica", "P. abies")), 
              by = "species") %>%
    mutate(height = paste0("Height = ", HT, " cm")) %>%
    mutate(height = factor(height, levels = paste0("Height = ", unique(.$HT)[order(unique(.$HT))], " cm"))) %>%
    ggplot(aes(x = Z, y = fit.mean, group = height, color = height)) + 
    geom_line() + 
    geom_ribbon(aes(ymin = fit.low, ymax = fit.high, fill = height), color = NA, alpha = 0.2) +
    facet_wrap(~ sp, nrow = 1) + 
    scale_color_manual(values = c("#335C67", "#E09F3E", "#9E2A2B")) +
    scale_fill_manual(values = c("#335C67", "#E09F3E", "#9E2A2B")) +
    ylab("Browsing\n probability") + xlab("Altitude (m)") +
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          axis.text.x = element_text(size = 5),
          strip.background = element_blank(), 
          legend.position = "none")
  
  # Plot the effect of each variable
  plot.param <- coef.est  %>%
    left_join(data.frame(species = c("ABAL", "ACPS", "FASY", "PIAB"), 
                         sp = c("A. alba", "A. pseudoplatanus", "F. sylvatica", "P. abies")), 
              by = "species") %>%
    filter(param != "Int.") %>%
    mutate(param = factor(param, levels = c("Z", "Ht", "Z:Ht"))) %>%
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
          axis.text.x = element_text(size = 5),
          strip.background = element_blank(), 
          legend.position = "none")
  
  # Final plot
  plot.out <- plot_grid(plot_grid(plot.prediction, plot.param, ncol = 1, align = "v", labels = c("(a)", "(b)")), 
                        get_legend(plot.prediction + theme(legend.position = "left", legend.title = element_blank())), 
                        nrow = 1, rel_widths = c(1, 0.3))
  
  # Save the plot
  ggsave(file.fig.in, plot = plot.out, width = 18, height = 10, units = "cm", dpi = 600)
  
  # Save the stat table
  write_on_disk(table.stat, file.table.in)
  
  # Return plot
  return(c(file.fig.in, file.table.in))
  
}




#' Function to plot the climatic gradient
#' @param df_br browsing data
#' @param map_per_plot precipitation data
#' @param file.in Name of the file to save, including path
plot_climate = function(df_br, map_per_plot, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Make the plot
  plot.out = df_br %>%
    dplyr::select(Site, Plot, Subplot, Z, Tmean) %>%
    left_join((map_per_plot %>% select(Site, Plot, map)), 
              by = c("Site", "Plot")) %>%
    dplyr::select(-Subplot) %>%
    distinct() %>%
    ggplot(aes(x = Tmean, y = Z, fill = map)) + 
    geom_point(shape = 21, color = "black") +
    facet_wrap(~ Site)
  
  # Save the plot
  ggsave(file.in, plot = plot.out, width = 14, height = 9, units = "cm", dpi = 600)
  
  # Return file
  return(file.in)
}

