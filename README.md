# EpiVisR: A Tool for Exploratory Data Analysis and Visualisation within Epigenome Wide Methylation Analyses (EWAS)

EpiVisR is a tool for visualisation of differentially methylated probes (DMP).

The EpiVisR package is designed to visualize results from EWAS analyses. It supports the user by providing
- Manhattan and volcano plots for CpG selection;
- Trait Methylation plots;
- Methylation profile plots and
- Correlation plots.

## Installation

library('remotes')  
install_github('SteRoe/EpiVisR')

## Usage

Point your working directory to the location with your files to visualize. This folder should also contain a valid yml-file with all necessary settings for

 betaFileName  
 traitFileName  
 genderFileName  
 dataDir  
 EWAScatalogFileName  
 baseURL_EWASDataHub  
 baseURL_MRCEWASCatalog  
 probeAttribut  
 mergeAttribut  
 genderAttribut  
 genderFemaleValue  
 genderMaleValue  

An example is given in the projects folder.  
setwd("your.working.directory")  
Start the App using EpiVisR::EpiVisRApp()

## License
[MIT](https://choosealicense.com/licenses/mit/)
