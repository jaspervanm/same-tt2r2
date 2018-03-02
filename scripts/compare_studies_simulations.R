source("scripts/simulation.R")

# Combining the results from
# original data (study_findings), and 
# simulated data (simulation_summary)
sim_stud_comparison <- bind_rows(simulation_summary, study_findings) %>%
	ungroup() %>%
	filter(
		variable     %in% unique(study_findings$variable),
		qoac_cutoff  %in% study_findings$qoac_cutoff,
		score_cutoff %in% study_findings$score_cutoff,
		variable     %in% c("npv", "ppv", "sens", "spec", "prev")
	) %>%
	mutate(
		variable = Variable_names[variable]
	) %>%
	select(study, score_cutoff, qoac_cutoff, origin, variable, pe, ll, ul) %>%
	gather("key", "value", -study, -contains("cutoff"), -origin, -variable ) %>%
	spread(origin, value) %>%
	filter(
		# only keep cutoffs evaluated in both the original study and the simulation
		!is.na(original_data),
		!is.na(simulation)
	)

# Calculating parameters for regression line
sim_stud_lm <- sim_stud_comparison %>%
	group_by(variable) %>%
	nest() %>%
	mutate(
		lm = map(data, ~lm(simulation ~ original_data, . )),
		lm_coef = map(lm, broom::tidy ),
		corr = map_dbl(data, ~cor(.$simulation, .$original_data))
	) %>%
	unnest(lm_coef) %>%
	select(variable, corr, term, estimate) %>%
	spread(term, estimate) %>%
	as.data.frame() %>%
	rename( intercept = `(Intercept)`, slope = original_data )

# Calibration plot
sim_stud_comparison %>%
	filter(key == "pe") %>% # only plot point estimates
	ggplot(aes(x     = original_data,
	   	       y     = simulation,
	   	       color = variable) ) +
	geom_point() + # < points in graph
	geom_abline(aes( # < regression line
		colour = variable,
		intercept = intercept,
		slope = slope,
		linetype = "regression line"),
		data = sim_stud_lm) +
	geom_abline(aes( # < line y=x 
		linetype = "y = x", 
		slope = 1, 
		intercept = 0) ) + 
	geom_text( aes(
		x = 1, y = 0,
		label = paste("r =", format(corr, digits = 2) ),
		hjust = "right", vjust = "bottom"),
		data = sim_stud_lm, size = 2.5) +
	# cosmetic settings
	theme_few(base_size = 8) +
	scale_colour_manual(values = unname(colours), guide = FALSE ) +
	lims( x = c(0, 1),
		  y = c(0, 1) ) +
	labs(
		x     = "Value in the original data",
		y     = "Value in the simulation") +
	facet_grid( ~ variable )

ggsave("correlation_simulation_original.png", path = "plots/", width = 7.5, height = 1.7, dpi = 300, units = "in")
ggsave("correlation_simulation_original.eps", path = "plots/", width = 19.05, height = 4.3, units = "cm")

# Bland Altman plot
sim_stud_comp_summary <- sim_stud_comparison %>%
	group_by(variable) %>%
	nest() %>%
	mutate(
		meanVal = map(data, ~with( . , (simulation + original_data) / 2) ),
		diffVal = map(data, ~with( . , simulation - original_data)),
		Mean = map_dbl(diffVal, mean),
		Sd   = map_dbl( diffVal, sd),
		ul   = Mean + 1.96 * Sd,
		ll   = Mean - 1.96 * Sd,
		Cor = map2( diffVal, meanVal,
					~cor.test( .x, .y) %>%
						broom::tidy( . ) %>%
						with( data.frame(cor = estimate, conf.low, conf.high, p.value))	),
		Lm   = map2( diffVal, meanVal, function(y, x) {
			lin_reg <- data.frame(x, y) %>%
				lm( y ~ x, . )
			
			cbind(
				broom::tidy( lin_reg ),
				confint( lin_reg )
			) %>%
				mutate(
					term = c("(Intercept)" = "intercept", "x" = "slope")[term]
				) %>%
				select(term, estimate, p.value, ll = `2.5 %`, ul = `97.5 %`) %>%
				gather(key, value, -term) %>%
				unite(key, term, key) %>%
				spread(key, value)
		} )
	) %>%
	unnest(Lm, Cor) %>%
	mutate(
		lin_reg_label = paste0("Intercept: ", round(intercept_estimate, 2), " (", round(intercept_ll, 2), " - ", round(intercept_ul, 2), ")\n",
							   "Slope: ", round(slope_estimate, 2), " (", round(slope_ll, 2)," - ", round(slope_ul, 2), ")")
	)

sim_stud_comparison %>%
	filter(key == "pe") %>% # only plot point estimates
	ggplot(aes( x      = ( simulation + original_data ) / 2,
			    y      = simulation - original_data,
			    shape = factor(qoac_cutoff),
			    colour  = factor(score_cutoff) ) ) +
	geom_point() + # < add points to the plot
	# add reference lines
	geom_hline(aes( yintercept = Mean, linetype = "Mean"),
			   data = sim_stud_comp_summary) +
	geom_hline(aes( yintercept = ll, linetype = "Mean ± 1.96 × sd"),
			   data = sim_stud_comp_summary) +
	geom_hline(aes( yintercept = ul, linetype = "Mean ± 1.96 × sd"),
			   data = sim_stud_comp_summary) +
	# add text with results from linear regression
	geom_text( aes( x = 0, y = -0.11, vjust = "bottom", hjust = "left",
					label = lin_reg_label ),
					data = sim_stud_comp_summary,
					inherit.aes = FALSE,
					size = 2.5) +
	# cosmetic settings
	theme_few() +
	theme(legend.position = "bottom",
		  legend.box = "horizontal",
		  legend.margin = margin(0,0,0,0),
		  legend.spacing.y = unit(0, "cm"),
		  legend.spacing.x = unit(2.5, "cm")) +
	scale_colour_manual(values = unname(colours), name = "SAMe-TT2R2 cutoff" ) +
	scale_shape_discrete(name = "TTR cutoff") +
	scale_linetype_discrete(guide = FALSE) +
	facet_grid( ~ variable ) +
	ylim(c(-0.11, 0.1)) +
	labs(
		x     = "Mean of simulated and original data",
		y     = "Difference between simulated and original data" )
ggsave("bland_altman_plot.png", path = "plots/", width = 1300/DPI, height = 450/DPI, dpi = DPI)
ggsave("bland_altman_plot.pdf", path = "plots/", width = 1300/DPI, height = 450/DPI, dpi = DPI)


studies_simulation_stats <- with(
	sim_stud_comparison, # < this one includes 95% conf/ref interval
	{ 
		c(
			cor  = cor(  simulation, original_data ),
			mean = mean( original_data - simulation ),
			sd   = sd(   original_data - simulation )
		)
	}
)


# how often is the point estimate of the simulation outside the original CI?
outside_CI <- sim_stud_comparison %>%
	group_by( study, qoac_cutoff, score_cutoff, variable ) %>%
	rename( original = original_data) %>%
	gather( source, value, original, simulation ) %>%
	unite( temp, source, key ) %>%
	spread( temp, value ) %>%
	mutate(
		outside_CI = ( simulation_pe > original_ul | simulation_pe < original_ll )
	)
outside_CI_mean <- mean( outside_CI$outside_CI, na.rm = TRUE )
