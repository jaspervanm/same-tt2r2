# Introduction
This folder contains everything that is needed to replicate the findings from the paper "Van Miert JHA, Bos S, Veeger NJGM, Meijer K (2018). Clinical usefulness of the SAMe-TT2R2 score: a systematic review and simulation meta-analysis. PLOS ONE.".

# How to use this
The folder `scripts/` contains all the scripts needed to replicate our study.  
The data that are used by the analyses, can be found in `search/` and `scripts/studies.R`.

To run the analysis, set the working directory to the directory of this file (e.g. `setwd("C:/same-tt2r2")`) and `source("scripts/manuscript.R")`.

One can install the packages needed manually by examining `scripts/init.R`, or install everything using *packrat*: `source(".Rprofile")`. The latter will download the same version of the dependencies as used when generating the article.