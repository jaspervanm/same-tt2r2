source("scripts/simulation.R")
source("scripts/meta_analysis.R")

# Function to apply LR to prior probability
do_LR <- function( LR, p ) {
	pre_odds <- p / (1 - p)
	post_odds <- pre_odds * LR
	
	output <- (post_odds / (post_odds + 1))
	names(output) <- p
	output
}

# Determine plot labels
plot_labels <- c(
	"2" = "SAMe-TT2R2 cutoff = 2",
	"3" = "SAMe-TT2R2 cutoff = 3",
	neg = "score < cutoff",
	pos = "score >= cutoff"
)

x <- 0:100 / 100 # the range of prior-probabilites (0-1) we want in the graph

# Construct prior-posterior probabilities data
# Each prior probability X is multiplied with the LR found
LR_plot_data <- MetaAnalysis %>%
	filter( variable %in% c("LRn", "LRp") ) %>%
	select( -meta_analysis ) %>%
	gather( "est", "value", pe, ll, ul ) %>%
	group_by( score_cutoff, variable, qoac_cutoff ) %>%
	mutate(
		y    = map( value, do_LR, x ),
		x    = map( y, ~ as.numeric( names( .x ) ) ),
		test = ifelse( variable == "LRp", "pos", "neg")
	) %>%
	unnest( x, y ) %>%
	select( -value ) %>%
	group_by( score_cutoff, variable, qoac_cutoff, x ) %>%
	spread( est, y ) %>%
	as.data.frame() %>%
	filter( !is.na(pe) )

# Construct dataframe with information about the original studies
plot_pre_post_data_studies <- simulation_and_studies %>%
	filter(
		qoac_cutoff == 70, score_cutoff %in% 2:3,
		variable %in% c("ppv", "npv", "prev")
	) %>%
	group_by( study, score_cutoff, qoac_cutoff ) %>%
	select( study, variable, pe, ll, ul, contains("cutoff") ) %>%
	mutate(
		prev_pe = pe[ variable == "prev"],
		prev_ll = ll[ variable == "prev"],
		prev_ul = ul[ variable == "prev"],
		pe   = ifelse( variable == "npv", 1 - pe, pe ),
		ll   = ifelse( variable == "npv", 1 - ll, ll ),
		ul   = ifelse( variable == "npv", 1 - ul, ul ),
		test = ifelse( variable == "ppv", "pos", "neg")
	) %>%
	filter( variable != "prev") 

# Create the plot
ggplot( plot_pre_post_data_studies, # < data for last part of the plot
		aes( x = prev_pe, y = pe, ymin = ll, ymax = ul) ) +
	geom_ribbon( # < Prior probability * confidence interval of LR
		aes( x = x, ymin = ll, ymax = ul, colour = test),
		LR_plot_data,
		linetype = "dotted", size = 0.8,
		fill = "white"
	) +
	geom_line( # < Prior probability * point estimate of LR
		aes( x = x, y = pe, colour = test),
		LR_plot_data,
		size = 1.2
	) +
	# Unity line (y=x, no additional information by test)
	geom_abline( aes( intercept = 0, slope = 1), linetype = "dashed") +
	# Information from the studies
	geom_errorbar( aes(colour = test ) ) + # < vertical error bars based on 95% CI of LRs
	geom_errorbarh( aes(xmin = prev_ll, xmax = prev_ul, colour = test) ) + # < horizontal error bars based on 95% CI of prevalence
	geom_point( data = plot_pre_post_data_studies) +
	# Make separate plot for each test outcome and SAMe-TT2R2 cutoff
	facet_grid( test ~ score_cutoff,
				labeller = labeller(test = plot_labels, score_cutoff = plot_labels) ) +
	# Cosmetic options
	guides( fill = FALSE, colour = FALSE ) +
	theme_few(base_size = 10) +
	scale_colour_manual( values = colours) +
	scale_fill_manual( values = colours) +
	labs(
		x = "Prior probability (prevalence) of a TTR < 70%",
		y = "Post test probability of a TTR < 70%"
	)

ggsave( "pretest_posttest.png", path = "plots/", width = 5.2, height = 4.5, dpi = 300)
ggsave( "pretest_posttest.eps", path = "plots/", width = 5.2, height = 4.5, dpi = 300)
