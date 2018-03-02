# SAMe-TT2R2 Evaluation
#
# Purpose of this file:
# - load all required libraries
# - define custom functions

# misc:
library(plyr)
library(tidyverse)
library(purrrlyr)

# calculations:
library(caret)

# plots: 
library(ggplot2)
library(ggthemes)
library(gridExtra)

# read input:
library(readxl)

# meta-analysis:
library(metafor)

# format output:
library(pander)

#------------------------------------------------------------------------------
# FORMATTING
#------------------------------------------------------------------------------
Variable_names <- c(
	prev      = "Low TTR prevalence",
	sens      = "Sensitivity",
	spec      = "Specificity",
	PSEP      = "PSEP",
	LRp       = "LR+",
	LRn       = "LR-",
	prev.pred = "Prev. of high score",
	npv       = "NPV",
	ppv       = "PPV"
)
Reference_lines <- c(
	LRp  = 0,
	LRn  = 0,
	PSEP = 0
)
DPI <- 100

colours <- c( # PLOS
	"#f8af2d", # yellow
	"#3c63af", # blue
	"#16a127", # green
	"#333333", # black
	"#891fb1", # purple
	"#667C7E", # gray
	neg = "#333333", pos = "#16a127",
	LRn = "#f8af2d", LRp = "#3c63af",
	PSEP = "#f8af2d", sens = "#3c63af", spec = "#16a127",
	pooled = "#333333"
)
LaTeX_escape_math <- function(x) {
	gsub(">=", "\\\\geq", x) %>%
		gsub("<=", "\\\\leq", . ) %>%
		paste0("$", . , "$")
}

#------------------------------------------------------------------------------
# CUSTOM FUNCTIONS
# generating reference intervals
#------------------------------------------------------------------------------
RIlow <- function(x, na.rm = FALSE) {
	quantile(x, probs = 0.025, na.rm = na.rm)
}
RIhigh <- function(x, na.rm = FALSE) {
	quantile(x, probs = 0.975, na.rm = na.rm)
}

#------------------------------------------------------------------------------
# CUSTOM FUNCTIONS
# formatting
#------------------------------------------------------------------------------
Pct <- function(x) {
	round(x * 100) %>%
		paste0( collapse = "-") %>%
		paste0( . , "%")
}
Collapse <- function(x, y = "and") {
	if(length(x) > 2) {
		x[ length(x)] <- paste0(y, " ", x[length(x)])
		paste0(x, collapse = ", ")
	} else if( length(x) == 2) {
		paste(x[1], y, x[2])
	} else {
		x[1]
	}
}
meanCI <- function(x, ci = 0.95, digits = 2) {
	if(ci > 1) { ci <- ci / 100 }
	
	probs  <- c( (1 - ci)/2, 1 - (1-ci)/2 )
	cis <- format( round(quantile(x, probs, na.rm = TRUE), digits), nsmall = digits )
	paste0(
		format( round(mean(x, na.rm = TRUE), digits), nsmall = digits),
		" (", cis[1], " - ", cis[2], ")"
	)
}
Paste <- function(x, digits = 2, ...) {
	if(is.numeric(x)) {
		x <- format(x, digits = digits)
	}
	paste(x, ...)
}


sumDF <- function(df) {
	df %>%
		select( which( sapply(., class) %in% c("numeric", "integer") ), study, contains('cutoff') ) %>%
		gather( "key", "value", -study, -contains('cutoff') ) %>%
		group_by( study, score_cutoff, qoac_cutoff, key ) %>%
		summarise( Mean = mean( value, na.rm = TRUE), CI = meanCI( value ) ) %>%
		ungroup()
}
