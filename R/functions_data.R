#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_data.R  
#' @description R script containing all functions relative to data
#               importation and formatting
#' @author Julien Barrere, Marianne Bernard
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## -- GENERIC FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Get file from its url and write it on disk, at a specified location. 
#' @param dir.name.in Directory where the file should be written (ex: "data/BWI")
#' @param url.in URL where to download the file.
get_and_write <- function(dir.name.in, url.in){
  
  # Write directories if they do not exist
  path.in <- strsplit(dir.name.in, "/")[[1]]
  for(i in 1:length(path.in)){
    if(i == 1) path.in_i <- path.in[i]
    else path.in_i <- paste(path.in_i, path.in[i], sep = "/")
    if(!dir.exists(path.in_i)) dir.create(path.in_i)
  }
  
  # Write file on the disk
  url.in_split <- strsplit(url.in, "/")[[1]]
  file.in <- paste(dir.name.in, url.in_split[length(url.in_split)], sep = "/")
  if(!file.exists(file.in)){
    try(GET(url.in, authenticate('guest', ""), write_disk(file.in, overwrite = TRUE)))
    # Specific case of zip file: unzip and delete zip file
    if("zip" %in% strsplit(file.in, split = "\\.")[[1]]){
      unzip(file.in, exdir = dir.name.in, overwrite = T)
      print(paste0("---Getting and unzipping ", file.in))
      unlink(file.in)
    }else{print(paste0("---Getting ", file.in))}
  } 
}



#' Function to get the path of a file, and create directories if they don't exist
#' @param file.in character: path of the file, filename included (ex: "plot/plot.png")
create_dir_if_needed <- function(file.in){
  
  path.in <- strsplit(file.in, "/")[[1]]
  if(length(path.in) > 1){
    for(i in 1:(length(path.in)-1)){
      if(i == 1) path.in_i <- path.in[i]
      else path.in_i <- paste(path.in_i, path.in[i], sep = "/")
      if(!dir.exists(path.in_i)) dir.create(path.in_i)
    }
  }
}




#' Write a table on disk
#' @param table.in dataframe to write on the disk
#' @param file.in Name (and path) of the file on the disk
write_on_disk <- function(table.in, file.in){
  create_dir_if_needed(file.in)
  write.table(table.in, file = file.in, row.names = F)
  return(file.in)
}





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## -- FUNCTIONS FOR FORMATTING FIELD DATA ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Format the data set to fit browsing rates model
#' @param bdd original data set
#' @param upi data set with the ungulate index
#' @author Marianne Bernard
format_browsing_rate <- function(bdd, upi){
  bdd$site_alt <- paste0(bdd$Site, bdd$Z)
  bdd$SpeciesUPI <- "Deciduous"
  bdd$SpeciesUPI[bdd$Species == "ABAL"] <- "Fir"
  bdd$SpeciesUPI[bdd$Species == "PIAB"] <- "Spruce"
  bdd <- left_join(bdd,upi, by =  c("SpeciesUPI" = "Species", "Site" = "Site"))  
  bdd$Tm_hiv<-(bdd$tm_new_dec+bdd$tm_new_jan+bdd$tm_new_fev)/3
  bdd$Tm_hiv_CR <- scale(bdd$Tm_hiv)
  bdd$HT_CR <- scale(bdd$HT)
  bdd <-  bdd[ !is.na(bdd$HT) & !is.na(bdd$Tm_hiv) & !is.na(bdd$B) & !is.na(bdd$Site) & !is.na(bdd$Z.Classe), ]
  bdd$Z_F <- factor(bdd$Z.Classe)
  return(bdd)
}


#' Format the data set to fit growth model
#' @param bdd original data set 
#' @author Marianne Bernard
format_growth <- function(bdd){
  bdd$B<-as.factor(bdd$B)
  bdd$Press.ong<-as.numeric(as.character(bdd$Press.ong))
  bdd$Thiv09_10 <- bdd$Tm_hiv
  return(bdd)
}



#' Function to calculate an index of ungulate mass density (umdi) for each plot coordinate
#' @param plot_coord Coordinates of each plot
#' @param ungulate_files list of the files contained in the data/ungulates folder
get_umdi_per_plot <- function(plot_coord, ungulate_files){
  
  # Ungulate body mass from https://naturalis.github.io/trait-organismal-ungulates/data/
  bodymass <- fread(ungulate_files[grep("csv", ungulate_files)]) %>%
    mutate(species = paste(Genus, Species, sep = " ")) %>%
    filter(species %in% c("Cervus elaphus", "Alces alces", "Rupicapra rupicapra", "Capreolus capreolus", 
                          "Capra ibex", "Dama dama", "Ovis aries")) %>%
    left_join(data.frame(species = c("Cervus elaphus", "Alces alces", "Rupicapra rupicapra", "Capreolus capreolus", 
                                     "Capra ibex", "Dama dama", "Ovis aries"), 
                         name = c("RedDeer", "Moose", "NorthernChamois", "RoeDeer", "AlpineIbex", "FallowDeer", "Mouflon")), 
              by = "species") %>%
    mutate(bodymass_kg = as.numeric(X5.1_AdultBodyMass_g)/1000) %>%
    dplyr::select(name, bodymass_kg)
  
  # Restrict ungulate_files to raster files
  ungulate_rast_files <- ungulate_files[grep("tif", ungulate_files)]
  
  # Vector containing the name of each ungulate species
  ungulate_species <- gsub("\\_.+", "", gsub(".+10km\\_", "", ungulate_rast_files))
  ungulate_species <- gsub("Red", "RedDeer", ungulate_species)
  ungulate_species <- gsub("Roe", "RoeDeer", ungulate_species)
  ungulate_species <- gsub("Fd", "FallowDeer", ungulate_species)
  
  # Initialize output data set
  out <- plot_coord
  
  # Loop on all ungulate species
  for(i in 1:length(ungulate_species)){
    
    # Raster for ungulate species i
    rast.i <- project(rast(ungulate_rast_files[i]), y = "epsg:4326")
    
    # Extract raster value for each plot
    out$ungulate.i <- as.numeric(terra::extract(rast.i, cbind(out$X, out$Y))[, 1])
    
    # Replace NA by 0
    out <- out %>% mutate(ungulate.i = ifelse(is.na(ungulate.i), 0, ungulate.i))
    
    # Rename newly created column with ungulate species name
    colnames(out)[dim(out)[2]] <- ungulate_species[i]
    
  }
  
  # calculate ungulate bodymass density index
  out <- out %>%
    gather(key = "name", value = "density", bodymass$name) %>%
    left_join(bodymass, by = "name") %>%
    mutate(bodymass_kg.km2 = bodymass_kg*density) %>%
    group_by(Site, Plot, Subplot) %>%
    summarize(umdi = sum(bodymass_kg.km2, na.rm = TRUE))
  
  # Return output
  return(out)
  
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## -- FUNCTIONS FOR CLIMATE DATA ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Extract only the plot coordinates from the dataset
#' @param df_gr formatted dataset for growth
extract_coord <- function(df_gr){
  df_gr %>%
    dplyr::select(Site, Plot, Subplot, X, Y, Z) %>%
    distinct()
}


#' Download chelsa files
#' @param bioclim numeric vector of the bioclim codes to extract
#' @param path character: directory where the data should be stored
#' @return 
#' @author G. Kunstler, J. Barrere
download_CHELSA <- function(bioclim, path){
  
  # Loop on all bioclim codes
  for(b in bioclim){
    if(! b %in% c(1:19)) stop("Argument bioclim needs to be in range 1:19")
    # Name of the file to extract
    filename.ymv <- paste("CHELSA_bio10", b, "land.7z", sep = "_")
    # Name of the url to extract
    url.ymv <- paste("https://zenodo.org/record/4996318/files", 
                     filename.ymv, sep = "/")
    # Download the file
    get_and_write(path, url.ymv)
  }
  
  # Identify archived files
  chelsa_files_archived <- paste(path, list.files(path), sep = "/")
  
  # Loop on all chelsa files archived to unzip them
  for(i in 1:length(chelsa_files_archived)){
    print(paste0("--- Unarchiving file ", chelsa_files_archived[i]))
    # Extract the file(s) from the archive
    archive_extract(chelsa_files_archived[i], dir = path)
    # Remove the archive
    unlink(chelsa_files_archived[i])
  }
  
  
  return(paste(path, list.files(path), sep = "/"))
}



#' Extract climate variables for each location 
#' @param chelsa_files character vector containing all chelsa files
#' @param plot_coord Cooordinates per plot
get_chelsa_map_per_plot <- function(chelsa_files, plot_coord){
  
  # Initialize output
  out <- plot_coord
  
  # Extract mean annual precipitation
  chelsa_file_map <- grep("bio10_12.tif", chelsa_files, value = TRUE)
  raster_map <- terra::rast(chelsa_file_map)
  out$map <- as.numeric(terra::extract(raster_map, cbind(out$X, out$Y))[, 1])
  
  # Return output
  return(out)
  
}




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## -- FUNCTIONS FOR RMARKDOWN ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




#' Function to knit the rmarkdown doocument
#' @param rmd_file Name including path of the rmarkdown document to knit
#' @param fig_files List of the figures that should be updated in the document
#' @return the name of the pdf file generated
knit_rmd <- function(rmd_file, fig_files){
  # Convert rmd file content in a character vector
  rmd_character <- readLines(con = rmd_file)
  
  # Count the presence of the figures in the rmd file
  k = 0
  for(i in 1:length(fig_files)) k = k + length(grep(fig_files[i], rmd_character))
  
  # Name of the output file
  file.out <- gsub("Rmd", "pdf", rmd_file)
  
  # If some of the figures are used in the rmd file, then knit it
  if(k > 0) rmarkdown::render(rmd_file, output_file = file.out, quiet = TRUE)
  
  # Return the pdf file 
  return(file.out)
  
}


#' Function to format the model outputs for the browsing model
#' @param browsing_models Fitted models for browsing
#' @param index.name Name of the index used for ungulate pressure
format_browsing_models <- function(browsing_models, index.name){
  
  # Species included in the model
  species.in <- names(browsing_models[[2]])
  
  # Convert species abbreviation in correct name
  species.name.in <- (data.frame(sp = species.in) %>%
                        left_join(data.frame(sp = c("ABAL", "ACPS", "FASY", "PIAB"), 
                                             species = c("A. alba", "A. pseudoplatanus", "F. sylvatica", "P. abies")), 
                                  by = "sp"))$species
  
  # Loop on all species
  for(i in 1:length(species.in)){
    # Results with ungulate index
    data.addvar.i = data.frame(var = c("Int", "T", "Ht", "P", index.name, "T:Ht", "P:Ht"), 
                               Est.1 = round(summary(browsing_models[[2]][[i]][[3]])$coefficients[, 1], digits = 2), 
                               Est.sd.1 = round(summary(browsing_models[[2]][[i]][[3]])$coefficients[, 2], digits = 2), 
                               p.1 = scales::pvalue(summary(browsing_models[[2]][[i]][[3]])$coefficients[, 4], accuracy = 0.01), 
                               VIF.1 = c("", as.character(round(vif.mer(browsing_models[[2]][[i]][[3]]), digits = 1))))
    rownames(data.addvar.i) <- NULL
    # Results without ungulate index
    data.novar.i = data.frame(var = c("Int", "T", "Ht", "P", "T:Ht", "P:Ht"), 
                              Est.2 = round(summary(browsing_models[[2]][[i]][[6]])$coefficients[, 1], digits = 2), 
                              Est.sd.2 = round(summary(browsing_models[[2]][[i]][[6]])$coefficients[, 2], digits = 2), 
                              p.2 = scales::pvalue(summary(browsing_models[[2]][[i]][[6]])$coefficients[, 4], accuracy = 0.01), 
                              VIF.2 = c("", as.character(round(vif.mer(browsing_models[[2]][[i]][[6]]), digits = 1))))
    rownames(data.novar.i) <- NULL
    
    # Merge the two dataset in a matrix
    matrix.i <- as.matrix(left_join(data.addvar.i, data.novar.i, by = "var"))
    matrix.i[is.na(matrix.i)] <- ""
    
    # Arrange the final table for species i
    data.i <- data.frame(
      col1 = c("Var", "",  "", matrix.i[, 1]), 
      col2 = c("Est. (sd)", species.name.in[i], paste0("Model with ", index.name), paste0(matrix.i[, 2], " (", matrix.i[, 3], ")")), 
      col3 = c("p value", "", paste0("(AIC=", round(AIC(browsing_models[[2]][[i]][[3]]), digits = 0), ")"), matrix.i[, 4]), 
      col4 = c("VIF", "", "", matrix.i[, 5]), 
      col5 = "",
      col6 = c("Est. (sd)", "", paste0("Model without ", index.name), paste0(matrix.i[, 6], " (", matrix.i[, 7], ")")), 
      col7= c("p value", "", paste0("(AIC=", round(AIC(browsing_models[[2]][[i]][[6]]), digits = 0), ")"), matrix.i[, 8]), 
      col8 = c("VIF", "", "", matrix.i[, 9])
    )
    
    # Add to the final table
    if(i == 1) out <- data.i
    else out <- rbind(out, matrix("", nrow = 1, ncol = 8, dimnames = list(NULL, colnames(data.i))), data.i[-1, ])
    
  }
  
  # Remove colnames of the output
  colnames(out) = NULL
  
  # Return output
  return(out)
}


#' Function to format the model outputs for the browsing model
#' @param growth_models Fitted models for growth
format_growth_models <- function(growth_models){
  
  # Species included in the model
  species.in <- names(growth_models[[2]])
  
  # Convert species abbreviation in correct name
  species.name.in <- (data.frame(sp = species.in) %>%
                        left_join(data.frame(sp = c("ABAL", "ACPS", "FASY", "PIAB"), 
                                             species = c("A. alba", "A. pseudoplatanus", "F. sylvatica", "P. abies")), 
                                  by = "sp"))$species
  
  # Loop on all species
  for(i in 1:length(species.in)){
    # Results with ungulate index
    data.i = data.frame(var = c("Int", "Ht", "T", "Br", "P", "T:Br", "P:Br"), 
                        Est. = round(summary(growth_models[[2]][[i]][[3]])$coefficients[, 1], digits = 2), 
                        Est.sd = round(summary(growth_models[[2]][[i]][[3]])$coefficients[, 2], digits = 2), 
                        p. = c("", scales::pvalue(car::Anova(growth_models[[2]][[i]][[3]])[, 3], accuracy = 0.01)), 
                        VIF = c("", as.character(round(vif.mer(growth_models[[2]][[i]][[3]]), digits = 1))))
    rownames(data.i) <- NULL
    
    # Merge the dataset in a matrix
    matrix.i <- as.matrix(data.i)
    matrix.i[is.na(matrix.i)] <- ""
    
    # Arrange the final table for species i
    data.i <- data.frame(
      col1 = c("Var", "",  matrix.i[, 1]), 
      col2 = c("Est. (sd)", species.name.in[i], paste0(matrix.i[, 2], " (", matrix.i[, 3], ")")), 
      col3 = c("p value", paste0("(AIC=", round(AIC(growth_models[[2]][[i]][[3]]), digits = 0), ")"), matrix.i[, 4]), 
      col4 = c("VIF", "", matrix.i[, 5]))
    
    # Add to the final table
    if(i == 1) out <- data.i
    else out <- rbind(out, matrix("", nrow = 1, ncol = 4, dimnames = list(NULL, colnames(data.i))), data.i[-1, ])
    
  }
  
  # Remove colnames of the output
  colnames(out) = NULL
  
  # Return output
  return(out)
}


