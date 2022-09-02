# oriT-Strast

### Multiple plasmid origin-of-transfer regions might aid the spread of antimicrobial resistance to human pathogens
Link to paper: [10.1002/mbo3.1129](http://dx.doi.org/10.1002/mbo3.1129)

---------------

<img src=https://github.com/JanZrimec/oriT-Strast/blob/master/docs/Fig_2b.png alt="drawing" width="400"> <img src=https://github.com/JanZrimec/oriT-Strast/blob/master/docs/strast_logo.png alt="drawing" width="400">

Figure 1. We develop a structural alignment procedure for typing plasmid-borne origin-of-transfer substrates that enables finding these regions in large plasmid datasets. Thousands of putative DNA transfer substrates are identified, showing that plasmid mobility can be 2-fold higher and span almost 2-fold more host species than is currently known. Over half of all the putative mobile plasmids carry conjugation systems belonging to different mobility groups, which might facilitate the transfer of antimicrobial resistance from environmental genetic reservoirs to human pathogens.

---------------

This repository contains scripts to reproduce the analysis and figures. 
The data for running the scripts is available at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3990610.svg)](https://doi.org/10.5281/zenodo.3990610), extract the archive to the folder named 'data'.
Dependencies are provided in the conda environment.yml file in the 'scripts' folder.

The algorithm scripts and readme are in the folder 'strast'.
Datasets S1, S2 and S3 are located in the folder 'data'.

See also the repositories [DNA_structural_variables predictor](https://github.com/JanZrimec/DNA_structural_variables) for prediction of DNA structural properties and [smer_acm_bcb_20](https://github.com/JanZrimec/smer_acm_bcb_20) for computing structural representations.
