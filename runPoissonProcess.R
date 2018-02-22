require(rstan)
require(bayesplot)

# Get the path of this script - the stan file should be here too
# https://github.com/molgenis/molgenis-pipelines.wiki.git
LocationOfThisScript = function() # Function LocationOfThisScript returns the location of this .R script (may be needed to source other files in same dir)
{
	this.file = NULL
	# This file may be 'sourced'
	for (i in -(1:sys.nframe())) {
		if (identical(sys.function(i), base::source)) this.file = (normalizePath(sys.frame(i)$ofile))
	}

	if (!is.null(this.file)) return(dirname(this.file))

	# But it may also be called from the command line
	cmd.args = commandArgs(trailingOnly = FALSE)
	cmd.args.trailing = commandArgs(trailingOnly = TRUE)
	cmd.args = cmd.args[seq.int(from=1, length.out=length(cmd.args) - length(cmd.args.trailing))]
	res = gsub("^(?:--file=(.*)|.*)$", "\\1", cmd.args)

	# If multiple --file arguments are given, R uses the last one
	res = tail(res[res != ""], 1)
	if (0 < length(res)) return(dirname(res))

	# Both are not the case. Maybe we are in an R GUI?
	return(NULL)
}
current.dir <- normalizePath(LocationOfThisScript())
setwd(current.dir)

# Compile the model
my_model <- rstan::stan_model(file = "poissonProcess.stan",
                              model_name = "CTVT poisson process")

# List of model parameters + input data
dat <- list(
    n = 2,  # number of post-MRCA mutation estimates
    m = c(23, 27),  # post-MRCA mutation estimates
    M = 221370,  # pre-MRCA mutation estimates (use array to pass a single value)
    min_permitted_lower_bound = 200,  # lower bound on uniform prior on 'T'
    max_permitted_upper_bound = 25000,  # upper bound on uniform prior on 'T'
    expon_rate = 0.025  # parameter of exponential prior on mean mutation rate (corresponds to ~= 40)
)

# Get MCMC samples from the model
trace <- rstan::sampling(my_model, data = dat, iter = 50000, chains = 4)

# List parameter estimates
summary(trace, prob = c(0.025, 0.975))$summary

# Plot some stuff
bayesplot::mcmc_areas(as.array(trace), pars = "t_mrca", prob = .95) + ggplot2::labs(title = "t_mrca")
bayesplot::mcmc_areas(as.array(trace), pars = "t[1]", prob = .95) + ggplot2::labs(title = "t_608")
bayesplot::mcmc_areas(as.array(trace), pars = "t[2]", prob = .95) + ggplot2::labs(title = "t_609")
bayesplot::mcmc_areas(as.array(trace), pars = "t_origin", prob = .95) + ggplot2::labs(title = "Time to CTVT origin (years)"); #ggsave("timeToOrigin.pdf")
bayesplot::mcmc_areas(as.array(trace), pars = "rate", prob = .95) + ggplot2::labs(title = "Mutation rate (total substitutions per year)"); #ggsave("mutationRate.pdf")
bayesplot::mcmc_areas(as.array(trace), pars = "i[1]", prob = .95) + ggplot2::labs(title = "i_608")
bayesplot::mcmc_areas(as.array(trace), pars = "i[2]", prob = .95) + ggplot2::labs(title = "i_609")
