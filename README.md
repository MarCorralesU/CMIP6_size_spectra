### <div style="clear: both;">  <font size= “1”>Code repository for the manuscript 'Assessment of plankton size structure from CMIP6 Earth System Models with a novel pelagic size structure database'</font></div>  
###### Marco Corrales-Ugalde, Jessica Y. Luo, Colleen M. Petrik, Charles A. Stock, Mathilde Dugenne, Rainer Kiko, Fabien Lombard, Todd O'Brien and Lars Stemmann
<a href="mailto:mcugalde88@gmail.com">Contact</a>


This repository documents the code used to generate the results for Corrales-Ugalde et al., usign the CMIP6 Earth System Model outputs from [ESGF-CoG](https://esgfnode.llnl.gov/search/cmip6/ "") and [ISIMIP DKRZ](https://www.isimip.org/dashboard/accessing-isimip-data-dkrz-server/ ""), and the  latest version of the [bulk PSSdb data products](https://zenodo.org/records/11050013 ""). 


1. Scripts 1.1-1.7 use the ESM biomass outputs to calculate the normalized biovolume, the size spectra coefficients, and group these results by biome following [Petrik et al. (2022)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2022GB007367 "").
there is a separate script for each ESM given the model-specific size ranges of their plankton functional group (PFT) and the approach to use the allometric equations to convert biomass to biovolume.

2. Script 2 uses the intermediate, weekly averaged PSSdb products with the biovolume information (not included in PSSdb) and the PSSdb-bulk data products to integrate the data by depth and sort the data by biomes and seasons to perform the skill assessment of the ESMs size spectra.
PSSdb weekly data products can be made available upon request.

3. Script 3 processes the .nc datafiles generated in scripts 1.1-1.7 and performs the data agreggation by biome and season to perform the model assessment skill with PSSdb

4. Scripts 4.1-4.x generate the figures in the main text and supplementary material. 

5. Script 5.1 contains a modified version of the taylor.diagram function of the [plotrix package by Lemon et al. (2006)](https://cran.r-project.org/web/packages/plotrix/index.html ""), and script 5.2 generated the coefficients included in tables S3-S7 and plots in figures 8 and S7




Acknowledgment: This work is funded by NOAA (Award #NA21OAR4310254)
