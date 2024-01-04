# BaccaraPaper

R script to reproduce the analyses of the paper *Elevation affects both the occurrence of ungulate browsing and its effect on tree seedling growth for four major tree species in European mountain forests* by Marianne Bernard, Julien Barrere, Xavier Morin, Sonia Saïd, Vincent Boulanger, Elena Granda, Raquel Benavides, Hervé Jactel, Marco Heurich, Sonia G. Rabasa, Fernando Valladares, and Georges Kunstler. This script is based on initial analyses by Marianne Bernard and Georges Kunstler, later modified by Julien Barrere. 

Before running the script, some data are required. The list of files needed and the structure that the "data" folder must have is visible in lines 43 to 44  of the ```_targets.R``` script. These files can be downloaded freely on Zenodo, at https://zenodo.org/records/10370097

Package ```targets``` is needed to run the script. Once the package is installed, just run ```targets::tar_make()``` from R and the script should install automatically the packages missing, and run the analyses. The figures will be stored in the directory "fig". 

For any question about the script or analyses, contact Julien BARRERE (julien.barrere@inrae.fr). 
