source("scripts/simulation.R")

# fix random number generator to get consistent results
set.seed(123)

# generate a sample of TTRs (scripts/simulation.R)
# based on study population (scripts/studies.R)
ttr <- tibble(
	study = names( study_populations ),
	data  = map( study_populations, generate_sample)
) %>%
	unnest() %>%
	select( study, samett2r2_score, value ) %>%
	# see if the test is positive
	# i.e. SAMe-TT2R2 >= cutoff
	mutate(
		study     = factor( study ),
		two_pos   = samett2r2_score >= 2,
		three_pos = samett2r2_score >= 3
	)

# Compare the effect of SAMe-TT2R2 group on continuous TTR
compare_ttr <- function( data, grouping_variable ) {
	data <- as.data.frame( data )
	data$group <- factor( data[, grouping_variable] )
	
	data %>%
		group_by( study ) %>%
		filter( length( unique( group ) ) > 1 ) %>% # to exclude studies that don't have both groups
		group_by( study, group ) %>%
		summarise(
			Mean = mean( value ),
			Sd   = sd( value ),
			N    = n(),
			se   = Sd / sqrt( N )
		) %>%
		do({
			# Calculate Hedges' g
			escalc( "SMDH", Mean/Sd ~ group | study, weights = N, data = ., slab = study) 
		}) %>%
		rma( yi, vi, data = ., slab = study ) %>%
		(function(x) {
			forest( x, xlab = paste0("Hedges' g of TTR based on ", Grouping_variable[grouping_variable]) )
			x
		})
}

Grouping_variable <- c(
	"two_pos" = "SAMe-TT2R2 >= 2",
	"three_pos" = "SAMe-TT2R2 >= 3"
)

# Create a forest plot for standardised mean difference (Hedges' g)
# of TTR based on SAMe-TT2R2-score
pdf("plots/forest_plot_ttr.pdf", width = 29.7 / 2.54, height = 21 / 2.54 )
par(mfcol=c(1,2))

ttr_smd <- list(
	cutoff.two   = compare_ttr( ttr, "two_pos"),
	cutoff.three = compare_ttr( ttr, "three_pos")
)

par(mfcol=c(1,1))
dev.off()

ttr_comp <- map( ttr_smd, predict) %>%
	do.call( rbind.data.frame, . ) %>%
	mutate(
		cutoff = c(2, 3),
		report = paste0(
			round(pred, 2),
			" (", round(ci.lb, 2), "-", round(ci.ub, 2), ")")
	)
