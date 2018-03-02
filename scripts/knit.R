# The file below is used to automatically render the manuscript.
# Not to be included in final version
knitr::knit("clinical_usefulness_of_the_samett2r2_score.Rmd",
			"manuscript/generated/body.md")
knitr::knit("abstract.Rmd",
			"manuscript/generated/abstract.md")

# Write information about the current session
sink("scripts/session_info.txt")
sessionInfo()
sink()
