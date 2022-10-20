#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_analysis.R  
#' @description R script containing all functions relative to data
#               analysis
#' @author Julien Barrere
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#' Function to apply different models to predict the effect of climate on browsing probability
#' @param df_br Dataset formatted for the browsing model
fit_browsing <- function(df_br, map_per_plot){
  # Vector of species present in the dataset
  sp <- unique(df_br$Species)
  
  # Add precipitation to the dataset
  df_br <- df_br %>%
    left_join(map_per_plot, by = c("Site", "Plot", "Subplot"))
  
  # Initialize table of the best model
  table.out <- data.frame(
    Model = paste0("MBr", c(1:6)), 
    Formulation = c("T + H + H*T + SUA", "Pr + H + H*Pr + SUA", 
                    "Pr + H*Pr + T + H*T + H + SUA", "T + H + H*T", 
                    "Pr + H + H*Pr", "Pr + H*Pr + T + H*T + H")
  )
  
  # Initialize final list containing the models
  list.out <- list()
  
  # Loop on all species
  for(i in 1:length(sp)){
    
    # Initialize model list for species i
    model.list <- list()
    
    # Model with temperature
    model.list$mbr1 <- glmer(B ~ scale(Tmean)*scale(HT) + scale(UPI) + (1|Site/Z.Classe),
                             family = 'binomial', data =subset(df_br, Species == sp[i]))
    
    # Model with precipitation
    model.list$mbr2 <- glmer(B ~ scale(map)*scale(HT) + scale(UPI) + (1|Site/Z.Classe),
                             family = 'binomial', data =subset(df_br, Species == sp[i]))
    
    # Model with precipitation
    model.list$mbr3 <- glmer(B ~ scale(Tmean)*scale(HT) + scale(map)*scale(HT) + scale(UPI) + (1|Site/Z.Classe),
                             family = 'binomial', data =subset(df_br, Species == sp[i]))
    
    # Model with temperature and no ungulate index
    model.list$mbr4 <- glmer(B ~ scale(Tmean)*scale(HT) + (1|Site/Z.Classe),
                             family = 'binomial', data =subset(df_br, Species == sp[i]))
    
    # Model with precipitation and no ungulate index
    model.list$mbr5 <- glmer(B ~ scale(map)*scale(HT) + (1|Site/Z.Classe),
                             family = 'binomial', data =subset(df_br, Species == sp[i]))
    
    # Model with precipitation and no ungulate index
    model.list$mbr6 <- glmer(B ~ scale(Tmean)*scale(HT) + scale(map)*scale(HT) + (1|Site/Z.Classe),
                             family = 'binomial', data =subset(df_br, Species == sp[i]))
    
    # Add to the model list
    eval(parse(text = paste0("list.out$", sp[i], " <- model.list")))
    
    # Add AIC to the final table
    eval(parse(text = paste0("table.out$", sp[i], "_AIC <- as.numeric(unlist(lapply(model.list, AIC)))")))
    eval(parse(text = paste0("table.out$", sp[i], "_deltaAIC <- table.out$", sp[i], "_AIC - min(table.out$", sp[i], "_AIC)")))
  }
  
  # Final output
  out <- list(table.out, list.out)
  return(out)
}



#' Function to apply different models to predict the effect of climate on browsing probability
#' @param df_br Dataset formatted for the browsing model
fit_browsing_hiv <- function(df_br, map_per_plot){
  # Vector of species present in the dataset
  sp <- unique(df_br$Species)
  
  # Add precipitation to the dataset
  df_br <- df_br %>%
    left_join(map_per_plot, by = c("Site", "Plot", "Subplot"))
  
  # Initialize table of the best model
  table.out <- data.frame(
    Model = paste0("MBr", c(1:6)), 
    Formulation = c("T + H + H*T + SUA", "Pr + H + H*Pr + SUA", 
                    "Pr + H*Pr + T + H*T + H + SUA", "T + H + H*T", 
                    "Pr + H + H*Pr", "Pr + H*Pr + T + H*T + H")
  )
  
  # Initialize final list containing the models
  list.out <- list()
  
  # Loop on all species
  for(i in 1:length(sp)){
    
    # Initialize model list for species i
    model.list <- list()
    
    # Model with temperature
    model.list$mbr1 <- glmer(B ~ scale(Tm_hiv)*scale(HT) + scale(UPI) + (1|Site/Z.Classe),
                             family = 'binomial', data =subset(df_br, Species == sp[i]))
    
    # Model with precipitation
    model.list$mbr2 <- glmer(B ~ scale(map)*scale(HT) + scale(UPI) + (1|Site/Z.Classe),
                             family = 'binomial', data =subset(df_br, Species == sp[i]))
    
    # Model with precipitation
    model.list$mbr3 <- glmer(B ~ scale(Tm_hiv)*scale(HT) + scale(map)*scale(HT) + scale(UPI) + (1|Site/Z.Classe),
                             family = 'binomial', data =subset(df_br, Species == sp[i]))
    
    # Model with temperature and no ungulate index
    model.list$mbr4 <- glmer(B ~ scale(Tm_hiv)*scale(HT) + (1|Site/Z.Classe),
                             family = 'binomial', data =subset(df_br, Species == sp[i]))
    
    # Model with precipitation and no ungulate index
    model.list$mbr5 <- glmer(B ~ scale(map)*scale(HT) + (1|Site/Z.Classe),
                             family = 'binomial', data =subset(df_br, Species == sp[i]))
    
    # Model with precipitation and no ungulate index
    model.list$mbr6 <- glmer(B ~ scale(Tm_hiv)*scale(HT) + scale(map)*scale(HT) + (1|Site/Z.Classe),
                             family = 'binomial', data =subset(df_br, Species == sp[i]))
    
    # Add to the model list
    eval(parse(text = paste0("list.out$", sp[i], " <- model.list")))
    
    # Add AIC to the final table
    eval(parse(text = paste0("table.out$", sp[i], "_AIC <- as.numeric(unlist(lapply(model.list, AIC)))")))
    eval(parse(text = paste0("table.out$", sp[i], "_deltaAIC <- table.out$", sp[i], "_AIC - min(table.out$", sp[i], "_AIC)")))
  }
  
  # Final output
  out <- list(table.out, list.out)
  return(out)
}



#' Function to compute prediction of model 6 from range of parameters
#' @param model model fit
#' @param model.type character with three possible choices: "temp", "prec", "all" 
#' @param data dataset containing scaled value of Tmean, Prec and HT
predict.mbr <- function(model, model.type = "all", data){
  
  # Extract lower and upper CI for each param
  param.confint <- data.frame(param = letters[c(1:dim(coefficients(summary(model)))[1])], 
                              low = as.numeric(confint(model, method = "Wald")[-c(1, 2), 1]),
                              high = as.numeric(confint(model, method = "Wald")[-c(1, 2), 2]))
  
  # Initialize dataset for parameters
  param <- data.frame(sim = c(1:5000))
  
  # Add parameters value for each simulation
  for(i in 1:dim(param.confint)[1]){
    eval(parse(text = paste0("param$", letters[i], " <- rnorm(5000, mean(c(", 
                             param.confint$low[i], ", ", param.confint$high[i], 
                             ")), (", param.confint$high[i], " - ", param.confint$low[i], 
                             ")/(2*1.96))")))
  }
  
  # Calculate confidence interval and mean
  out <- full_join(data, data.frame(sim = c(1:5000)), by = character()) %>%
    left_join(param, by = "sim")
  if(model.type == "all") out <- mutate(
    out, fit = a + b*Tmean.scaled + c*HT.scaled + d*map.scaled +  e*Tmean.scaled*HT.scaled + f*map.scaled*HT.scaled)
  if(model.type == "temp") out <- mutate(out, fit = a + b*Tmean.scaled + c*HT.scaled + d*Tmean.scaled*HT.scaled)
  if(model.type == "prec") out <- mutate(out, fit = a + b*map.scaled + c*HT.scaled + d*map.scaled*HT.scaled)
  out <- out %>%
    mutate(fit = plogis(fit)) %>%
    group_by(Tmean, map, HT) %>%
    summarize(fit.low = quantile(fit, probs = 0.025), 
              fit.mean = mean(fit, na.rm = TRUE),
              fit.high = quantile(fit, probs = 0.975))
  
  return(out)
}



#' Function to compute prediction of model 6 from range of parameters
#' @param model model fit
#' @param model.type character with three possible choices: "temp", "prec", "all" 
#' @param data dataset containing scaled value of Tmean, Prec and HT
predict.mbr_hiv <- function(model, model.type = "all", data){
  
  # Extract lower and upper CI for each param
  param.confint <- data.frame(param = letters[c(1:dim(coefficients(summary(model)))[1])], 
                              low = as.numeric(confint(model, method = "Wald")[-c(1, 2), 1]),
                              high = as.numeric(confint(model, method = "Wald")[-c(1, 2), 2]))
  
  # Initialize dataset for parameters
  param <- data.frame(sim = c(1:5000))
  
  # Add parameters value for each simulation
  for(i in 1:dim(param.confint)[1]){
    eval(parse(text = paste0("param$", letters[i], " <- rnorm(5000, mean(c(", 
                             param.confint$low[i], ", ", param.confint$high[i], 
                             ")), (", param.confint$high[i], " - ", param.confint$low[i], 
                             ")/(2*1.96))")))
  }
  
  # Calculate confidence interval and mean
  out <- full_join(data, data.frame(sim = c(1:5000)), by = character()) %>%
    left_join(param, by = "sim")
  if(model.type == "all") out <- mutate(
    out, fit = a + b*Tm_hiv.scaled + c*HT.scaled + d*map.scaled +  e*Tm_hiv.scaled*HT.scaled + f*map.scaled*HT.scaled)
  if(model.type == "temp") out <- mutate(out, fit = a + b*Tm_hiv.scaled + c*HT.scaled + d*Tm_hiv.scaled*HT.scaled)
  if(model.type == "prec") out <- mutate(out, fit = a + b*map.scaled + c*HT.scaled + d*map.scaled*HT.scaled)
  out <- out %>%
    mutate(fit = plogis(fit)) %>%
    group_by(Tm_hiv, map, HT) %>%
    summarize(fit.low = quantile(fit, probs = 0.025), 
              fit.mean = mean(fit, na.rm = TRUE),
              fit.high = quantile(fit, probs = 0.975))
  
  return(out)
}


# Function to calculate VIF per predictor in a mixed model
vif.mer <- function (fit) {
  
  v <- vcov(fit)
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}


#' Function to apply different models to predict the effect of climate and browsing on growth
#' @param df_gr Dataset formatted for the browsing model
#' @param map_per_plot precipitation per plot
fit_growth <- function(df_gr, map_per_plot){
  
  # Vector of species present in the dataset
  sp <- unique(df_gr$Species)
  
  # Add precipitation to the dataset
  df_gr <- df_gr %>%
    left_join(map_per_plot, by = c("Site", "Plot", "Subplot"))
  
  # Initialize table of the best model
  table.out <- data.frame(
    Model = paste0("Mgr", c(1:3)), 
    Formulation = c("H + Br + T + T*Br", "H + Br + P + P*Br", 
                    "H + Br + T + T*Br + P + P*Br")
  )
  
  # Initialize final list containing the models
  list.out <- list()
  
  # Loop on all species
  for(i in 1:length(sp)){
    
    # Initialize model list for species i
    model.list <- list()
    
    # Model with temperature
    model.list$mgr1 <- lmer(log(LLS+1) ~ scale(HT) + scale(Tmean)*B + (1|Site/Z.Classe), 
                            data = subset(df_gr, Species == sp[i]))
    
    # Model with precipitation
    model.list$mgr2 <- lmer(log(LLS+1) ~ scale(HT) + scale(map)*B + (1|Site/Z.Classe), 
                            data = subset(df_gr, Species == sp[i]))
    
    # Model with precipitation
    model.list$mgr3 <- lmer(log(LLS+1) ~ scale(HT) + scale(Tmean)*B + scale(map)*B + (1|Site/Z.Classe), 
                            data = subset(df_gr, Species == sp[i]))
    
    # Add to the model list
    eval(parse(text = paste0("list.out$", sp[i], " <- model.list")))
    
    # Add AIC to the final table
    eval(parse(text = paste0("table.out$", sp[i], "_AIC <- as.numeric(unlist(lapply(model.list, AIC)))")))
    eval(parse(text = paste0("table.out$", sp[i], "_deltaAIC <- table.out$", sp[i], "_AIC - min(table.out$", sp[i], "_AIC)")))
  }
  
  # Final output
  out <- list(table.out, list.out)
  return(out)
}


#' Function to compute prediction of growth model from range of parameters
#' @param model model fit
#' @param model.type character with three possible choices: "temp", "prec", "all" 
#' @param data dataset containing scaled value of Tmean, Prec and HT
predict.mgr <- function(model, model.type = "all", data){
  
  # Extract lower and upper CI for each param
  param.confint <- data.frame(param = letters[c(1:dim(coefficients(summary(model)))[1])], 
                              low = as.numeric(confint(model, method = "Wald")[-c(1:3), 1]),
                              high = as.numeric(confint(model, method = "Wald")[-c(1:3), 2]))
  
  # Initialize dataset for parameters
  param <- data.frame(sim = c(1:5000))
  
  # Add parameters value for each simulation
  for(i in 1:dim(param.confint)[1]){
    eval(parse(text = paste0("param$", letters[i], " <- rnorm(5000, mean(c(", 
                             param.confint$low[i], ", ", param.confint$high[i], 
                             ")), (", param.confint$high[i], " - ", param.confint$low[i], 
                             ")/(2*1.96))")))
  }
  
  # Calculate confidence interval and mean
  out <- full_join(data, data.frame(sim = c(1:5000)), by = character()) %>%
    left_join(param, by = "sim")
  if(model.type == "all") out <- mutate(
    out, fit = a + b*HT.scaled + c*Tmean.scaled + d*B + e*map.scaled + f*Tmean.scaled*B + g*map.scaled*B)
  if(model.type == "temp") out <- mutate(out, fit = a + b*HT.scaled + c*Tmean.scaled + d*B + e*Tmean.scaled*B)
  if(model.type == "prec") out <- mutate(out, fit = a + b*HT.scaled + c*map.scaled + d*B + e*map.scaled*B)
  out <- out %>%
    mutate(fit = exp(fit) - 1) %>%
    group_by(Tmean, map, HT, B) %>%
    summarize(fit.low = quantile(fit, probs = 0.025), 
              fit.mean = mean(fit, na.rm = TRUE),
              fit.high = quantile(fit, probs = 0.975))
  
  return(out)
}
