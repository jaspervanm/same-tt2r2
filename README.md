# Introduction
This folder contains everything that is needed to replicate the findings from the paper "Clinical benefit of the SAMe-TT2R2 score: a systematic review and simulation meta-analysis".

# How to use this
If you just want to take a look at the R scripts, without actually running them, you can find them in  `.scripts/`.  
If you want to run the simulation yourself, follow these steps:

1. Extract this archive or clone this repository
2. In R, set the working directory to the directory where you extracted the files, e.g. `setwd("C:/SAMe-TT2R2")`
3. Optional: if you want to install all dependencies to an isolated library: in R, run `source(".Rprofile")`. This will install all packages needed in the library in `packrat/`. If you want to use your own library, navigate to `scripts/init.R` and install all packages mentioned there.
4. To run the scripts, `source("scripts/manuscript.R")` to execute all scripts, or just `source` the one you're interested in