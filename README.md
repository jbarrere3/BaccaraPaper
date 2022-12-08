# DisturbancePaper

R script to reproduce the analyses of the paper *Climate is a driver of ungulate browsing intensity and of its effect on tree seedling growth* by JMarianne Bernard, Xavier Morin, Sonia Saïd, Vincent Boulanger, Elena Granda, Julien Barrere, Raquel Benavides, Hervé Jactel, Marco Heurich, Sonia G. Rabasa, Fernando Valladares, and Georges Kunstler. This script is based on initial analyses by Marianne Bernard and Georges Kunstler, later modified by Julien Barrere. 

Before running the script, some data are required. The list of files needed and the structure that the "data" folder must have is visible in lines 43 to 44  of the ```_targets.R``` script. These files can be obtained upon request to Julien BARRERE (julien.barrere@inrae.fr) or Georges Kunstler (georges.kunstler@inrae.fr). 

Package ```targets``` is needed to run the script. Once the package is indeed, just run ```targets::tar_make()``` from R and the script should install automatically the packages missing, and run the analyses. The figures will be stored in the directory "fig". 

For any question about the script or analyses, contact Julien BARRERE (julien.barrere@inrae.fr). 
