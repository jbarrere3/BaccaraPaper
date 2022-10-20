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
