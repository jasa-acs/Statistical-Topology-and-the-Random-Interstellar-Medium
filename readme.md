# Statistical Topology and the Random Interstellar Medium

# Author Contributions Checklist Form

## Data 

### Abstract 

Emission of neutral atomic hydrogen in a region of the southern sky from the Galactic All Sky Survey.  

### Availability

Available.

### Description

The source data are publicly available from 

https://www.astro.uni-bonn.de/hisurvey/gass/

The data originators  request that users  of GASS data products reference the following:

Kalberla, P.M.W. and Haud, U. (2015) A&A 578, A78 (Highlight: Gass: The Parkes Galactic All-Sky Survey. III)

Kalberla, P.M.W., McClure-Griffiths, N.M., Pisano, D. J. Calabretta, M. R., Ford, H. Alyson, Lockman, Felix J., Staveley-Smith, L., Kerp, J., Winkel, B., Murphy, T., Newton-McGee, K. (2010) A&A, 512, A14 (Highlight: Gass: The Parkes Galactic All-Sky Survey. II)

McClure-Griffiths, N. M.; Pisano, D. J.; Calabretta, M. R.; Ford, H. Alyson; Lockman, Felix J.; Staveley-Smith, L.; Kalberla, P. M. W.; Bailin, J.; Dedes, L.; Janowiecki, S , Gibson, B.K., Murpy, T., Newton-McGee, K. (2009) ApJS, 181, 398 (http://adsabs.harvard.edu/abs/2009ApJS..181..398M)

We have given detailed instructions on how to download the original data and matlab code to create the integrated version used for analysis.


## Code

### Abstract 

Code to read the data, produce figures, tables and simulation results for the submission are provided as supplementary material.

### Description 

R code for all analyses within the paper, together with instructions.

### Optional Information 

Requires 

library(TDA)
library(RandomFields)
library(mvtnorm) 
library(mgcv)
library(plot3D

The R version we used was 3.5.3.  The package versions  are   TDA1.6.5,    RandomFields 3.3.6,   mvtnorm 1.0-10,  mgcv 1.8-28, plot3D 1.1.1


## Instructions for Use

### Reproducibility

All tables and figures from the paper can be reproduced by following Instructions.pdf  

### Replication 

The code is not written for replication but will work for any square matrix of observations.

