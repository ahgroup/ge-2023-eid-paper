# Replicating Yang's Models ####################################################
# Project: Analyzing the association between inoculum dose and Norovirus
#   infection outcomes
# Author: Zane Billings
# Since: 2022-03-11
# Desc: For this paper, Yang fit a bunch of bayesian regression models using
#   the brms package. I will implement his models using rethinking::ulam()
#   and check to see if the model results are the same.

# Dependencies =================================================================

# install.packages(c("coda","mvtnorm","devtools","loo","dagitty","shape"))
# devtools::install_github("rmcelreath/rethinking")

# some stuff doesn't work when I don't load the entire rethinking namespace
library(rethinking)

options(mc.cores = parallel::detectCores() - 1)

# Load data ####################################################################

# Yang's cleaned data -- only needs minimal preprocessing now
# TODO update this with correct data path
dat_orig <- readRDS(here::here("2_supplementary_materials/data/norodata.rds"))

# Apply preprocessing steps to all list elements except data_tab1
norovirus_data <-
  dat_orig[-1] |>
  # Apply same preprocessing steps to every list element
  purrr::map(
    ~ .x |>
      # Filter out the person in the low dose group
      dplyr::filter(dose > 0.48) |>
      dplyr::mutate(
        # drop the very low dose factor level
        dose_cat = forcats::fct_drop(dose_cat),
        # convert ID to numeric
        ID_old = factor(ID, levels = paste0("ID", 1:19)),
        ID = as.integer(droplevels(ID_old))
      ) |>
      dplyr::arrange(ID)
  )

# Create the output folders
suppressWarnings({
  dir.create(here::here("2_supplementary_materials", "analyses", "rethinking-results"))
  dir.create(here::here("2_supplementary_materials", "analyses", "rethinking-results", "data"))
  dir.create(here::here("2_supplementary_materials", "analyses", "rethinking-results", "res"))
  dir.create(here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds"))
  dir.create(here::here("2_supplementary_materials", "analyses", "rethinking-results", "figs"))
})

# Global constants #############################################################
# These are set the same as in the brms code

# Target acceptance rate (adaptive delta)
DELTA <- 0.999999
# maximum tree depth
MAX_TD <- 20
# number of chains to sample from
N_CHAINS <- 4
# number of cores to distribute chains over
N_CORES <- 4
# Number of warmup sampling iterations
N_WU <- 25000
# Number of sampling iterations
N_ITER <- 30000
# Should log likelihood be calculated for each iteration?
RTN_LL <- TRUE
# Should cmdstan be used instead of rstan?
USE_CMDSTAN <- TRUE
# Seed for Stan
SEED <- 123456789

# Total fecal shedding <= 96 hours - Continuous Dose ###########################

# Select only the necessary columns to pass to Stan
tfs_96_dat <-
  norovirus_data$feces96 |>
  dplyr::select(ID, dose_log, outcome_log)

# Save processed data to file
saveRDS(
  object = tfs_96_dat,
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "data", "tfs_96_con_data.Rds")
)

# Fitting the model as described in Yang's supplement.
tfs_96_model <- rethinking::ulam(
  data = tfs_96_dat,
  # This part defines the model
  flist = alist(
    # Distribution of outcome
    outcome_log ~ dnorm(mu, sigma),
    # Mu_i is a linear function of the dose, intercepts are random
    mu <- alpha[ID] + beta * dose_log,
    # Use an adaptive prior for random intercepts
    alpha[ID] ~ dnorm(mu_alpha, gamma),
    # Fixed priors for other parameters
    mu_alpha ~ dnorm(25, 5),
    gamma ~ dcauchy(0, 2),
    beta ~ dnorm(0, 1),
    sigma ~ dcauchy(0, 2)
  ),
  # Variance parameters must be constrained to be > 0. This is equivalent
  # to specifying them with a half Cauchy distribution.
  constraints = list(
    sigma = "lower=0",
    gamma = "lower=0"
  ),
  # These are Stan parameters that I defined as global constants.
  log_lik = RTN_LL,
  cmdstan = USE_CMDSTAN,
  control = list(
    adapt_delta = DELTA,
    max_treedepth = MAX_TD
  ),
  chains = N_CHAINS,
  cores = N_CORES,
  warmup = N_WU,
  iter = N_ITER,
  seed = SEED
)

# Save model results to disk
saveRDS(
  object = tfs_96_model,
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "res", "tfs_96_model.Rds")
)

# Total fecal shedding full time - Continuous Dose #############################

# Read in processed data
tfs_all_dat <-
  norovirus_data$fecesall |>
  dplyr::select(ID, dose_log, outcome_log)

# Save processed data to file
saveRDS(
  object = tfs_all_dat,
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "data", "tfs_all_con_data.Rds")
)

# Next fitting the model
tfs_model_all <- rethinking::ulam(
  data = tfs_all_dat,
  # This part defines the model
  flist = alist(
    # Distribution of outcome
    outcome_log ~ dnorm(mu, sigma),
    # Mu_i is a linear function of the dose, intercepts are random
    mu <- alpha[ID] + beta * dose_log,
    # Use an adaptive prior for random in tercepts
    alpha[ID] ~ dnorm(mu_alpha, gamma),
    # Fixed priors for other parameters
    mu_alpha ~ dnorm(25, 5),
    gamma ~ dcauchy(0, 2),
    beta ~ dnorm(0, 1),
    sigma ~ dcauchy(0, 2)
  ),
  # Variance parameters must be constrained to be > 0. This is equivalent
  # to specifying them with a half Cauchy distribution.
  constraints = list(
    sigma = "lower=0",
    gamma = "lower=0"
  ),
  # These are Stan parameters that I defined as global constants.
  log_lik = RTN_LL,
  cmdstan = USE_CMDSTAN,
  control = list(
    adapt_delta = DELTA,
    max_treedepth = MAX_TD
  ),
  chains = N_CHAINS,
  cores = N_CORES,
  warmup = N_WU,
  iter = N_ITER,
  seed = SEED
)

# Save model results to disk
saveRDS(
  object = tfs_model_all,
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "res", "tfs_all_model.Rds")
)

# Total vomit shedding - Continuous dose #######################################

# Read in processed data
vvshed_dat <-
  norovirus_data$vomit |>
  dplyr::select(ID, dose_log, outcome_log)

# Save processed data to file
saveRDS(
  object = vvshed_dat,
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "data", "vvshed_con_data.Rds")
)

# Next fitting the model
vvshed_model_all <- rethinking::ulam(
  data = vvshed_dat,
  # This part defines the model
  flist = alist(
    # Distribution of outcome
    outcome_log ~ dnorm(mu, sigma),
    # Mu_i is a linear function of the dose, intercepts are random
    mu <- alpha[ID] + beta * dose_log,
    # Use an adaptive prior for random in tercepts
    alpha[ID] ~ dnorm(mu_alpha, gamma),
    # Fixed priors for other parameters
    mu_alpha ~ dnorm(20, 5),
    gamma ~ dcauchy(0, 2),
    beta ~ dnorm(0, 1),
    sigma ~ dcauchy(0, 2)
  ),
  # Variance parameters must be constrained to be > 0. This is equivalent
  # to specifying them with a half Cauchy distribution.
  constraints = list(
    sigma = "lower=0",
    gamma = "lower=0"
  ),
  # These are Stan parameters that I defined as global constants.
  log_lik = RTN_LL,
  cmdstan = USE_CMDSTAN,
  control = list(
    adapt_delta = DELTA,
    max_treedepth = MAX_TD
  ),
  chains = N_CHAINS,
  cores = N_CORES,
  warmup = N_WU,
  iter = N_ITER,
  seed = SEED
)

# Save model results to disk
saveRDS(
  object = vvshed_model_all,
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "res", "vvshed_model_all.Rds")
)

# Fecal shedding <= 96 hours - Categorical Dose ################################

tfs_96_cat_dat <-
  norovirus_data$feces96 |>
  dplyr::select(ID, dose_cat, outcome_log)

# Save processed data to file
saveRDS(
  object = tfs_96_cat_dat,
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "data", "tfs_96_cat_data.Rds")
)

# Fitting the model as described in Yang's supplement.
tfs_96_cat_model <- rethinking::ulam(
  data = tfs_96_cat_dat,
  # This part defines the model
  flist = alist(
    # Distribution of outcome
    outcome_log ~ dnorm(mu, sigma),
    # Mu_i is a linear function of the dose, intercepts are random
    mu <- alpha[ID] + beta[dose_cat],
    # Use an adaptive prior for random in tercepts
    alpha[ID] ~ dnorm(0, gamma),
    # Fixed priors for other parameters
    gamma ~ dcauchy(0, 2),
    beta[dose_cat] ~ dnorm(25, 5),
    sigma ~ dcauchy(0, 2)
  ),
  # Variance parameters must be constrained to be > 0. This is equivalent
  # to specifying them with a half Cauchy distribution.
  constraints = list(
    sigma = "lower=0",
    gamma = "lower=0"
  ),
  # These are Stan parameters that I defined as global constants.
  log_lik = RTN_LL,
  cmdstan = USE_CMDSTAN,
  control = list(
    adapt_delta = DELTA,
    max_treedepth = MAX_TD
  ),
  chains = N_CHAINS,
  cores = N_CORES,
  warmup = N_WU,
  iter = N_ITER,
  seed = SEED
)

# Save model results to disk
saveRDS(
  object = tfs_96_cat_model,
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "res", "tfs_96_cat_model.Rds")
)

# Fecal shedding full time - Categorical dose ##################################

# Read in processed data
tfs_all_cat_dat <-
  norovirus_data$fecesall |>
  dplyr::select(ID, dose_cat, outcome_log)

# Save processed data to file
saveRDS(
  object = tfs_all_cat_dat,
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "data", "tfs_all_cat_data.Rds")
)

# Next fitting the model
tfs_model_all_cat <- rethinking::ulam(
  data = tfs_all_cat_dat,
  # This part defines the model
  flist = alist(
    # Distribution of outcome
    outcome_log ~ dnorm(mu, sigma),
    # Mu_i is a linear function of the dose, intercepts are random
    mu <- alpha[ID] + beta[dose_cat],
    # Use an adaptive prior for random in tercepts
    alpha[ID] ~ dnorm(0, gamma),
    # Fixed priors for other parameters
    gamma ~ dcauchy(0, 2),
    beta[dose_cat] ~ dnorm(25, 5),
    sigma ~ dcauchy(0, 2)
  ),
  # Variance parameters must be constrained to be > 0. This is equivalent
  # to specifying them with a half Cauchy distribution.
  constraints = list(
    sigma = "lower=0",
    gamma = "lower=0"
  ),
  # These are Stan parameters that I defined as global constants.
  log_lik = RTN_LL,
  cmdstan = USE_CMDSTAN,
  control = list(
    adapt_delta = DELTA,
    max_treedepth = MAX_TD
  ),
  chains = N_CHAINS,
  cores = N_CORES,
  warmup = N_WU,
  iter = N_ITER,
  seed = SEED
)

# Save model results to disk
saveRDS(
  object = tfs_model_all_cat,
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "res", "tfs_model_all_cat.Rds")
)

# Total vomit shedding - Categorical dose ######################################

vvshed_dat_cat <-
  norovirus_data$vomit |>
  dplyr::select(ID, dose_cat, outcome_log)

# Save processed data to file
saveRDS(
  object = vvshed_dat_cat,
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "data", "vvshed_cat_data.Rds")
)

# Next fitting the model
vvshed_model_cat <- rethinking::ulam(
  data = vvshed_dat_cat,
  # This part defines the model
  flist = alist(
    # Distribution of outcome
    outcome_log ~ dnorm(mu, sigma),
    # Mu_i is a linear function of the dose, intercepts are random
    mu <- alpha[ID] + beta[dose_cat],
    # Use an adaptive prior for random in tercepts
    alpha[ID] ~ dnorm(0, gamma),
    # Fixed priors for other parameters
    gamma ~ dcauchy(0, 2),
    beta[dose_cat] ~ dnorm(20, 5),
    sigma ~ dcauchy(0, 2)
  ),
  # Variance parameters must be constrained to be > 0. This is equivalent
  # to specifying them with a half Cauchy distribution.
  constraints = list(
    sigma = "lower=0",
    gamma = "lower=0"
  ),
  # These are Stan parameters that I defined as global constants.
  log_lik = RTN_LL,
  cmdstan = USE_CMDSTAN,
  control = list(
    adapt_delta = DELTA,
    max_treedepth = MAX_TD
  ),
  chains = N_CHAINS,
  cores = N_CORES,
  warmup = N_WU,
  iter = N_ITER,
  seed = SEED
)

# Save model results to disk
saveRDS(
  object = vvshed_model_cat,
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "res", "vvshed_model_cat.Rds")
)

# Viral kinetics (longitudinal analysis) - Continuous Dose #####################

timeseries_con_dat <-
  norovirus_data$ftsshedding |>
  dplyr::select(ID, dose_log, days, logyc) |>
  dplyr::filter(days > 0)

# Save processed data to file
saveRDS(
  object = timeseries_con_dat,
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "data", "timeseries_con_data.Rds")
)

# Stupid priors explanation for gamma prior
# I think rethinking is mangling my gamma parameters, and I'm not sure how to
# make it stop. ulam() forces all gamma priors to go through
# rethinking::dgamma2(), which internally does a division and passes the results
# to stats::dgamma() . It is apparently impossible to use dgamma() directly in
# ulam() instead, and when you specify a gamma prior using Stan's gamma(alpha,
# beta) syntax, it gets passed to dgamma2, divided (even though this step
# doesn't need to happen), and passed to dgamma().
#
# I didn't notice this before because I had the parameters misspecified and they
# were both decimals. Now I'm specifying the prior as dgamma2(20, 10), which is
# equivalent to the prior Yang is using with brms. But I get a warning from the
# Stan code that tells me it's automatically rounding the integer division 20 /
# 10 and if I want to avoid rounding, I have to use 20.0 / 10 instead. But
# passing in dgamma2(20.0, 10.0) doesn't fix it, I still get the same warning.
# (That makes sense because in R, every number is a float by default.)
#
# But anyways the 20 / 10 division isn't the problem, the problem is that it has
# to compute the inverse scale parameter, which is 1 / 10. It rounds this to
# zero and then the model fails because the parameter has to be positive.
#
# So if I set the prior to dgamma2(20.01, 10.01), hopefully this distribution
# will be close enough to Yang's distribution not to make a difference even
# if the models are sensitive to the priors. But adding decimals forces ulam
# to write stan code that declares the parameters as type real instead of as
# type int, forcing them to actually be divided instead of rounding.

# Next fit the model. Note that this one takes forever to run and has a lot of
# divergent transitions during warmup.
timeseries_con_model <- rethinking::ulam(
  data = timeseries_con_dat,
  # This part defines the model
  flist = alist(
    # Distribution of outcome
    logyc ~ dstudent(k, mu, sigma),
    # Overall time series equation (I let u = T since capital T is a special
    # character in R and I don't want to mess it up).
    mu <- log((2 * exp(p)) / bottom),
    bottom <- exp(-exp(g) * (days - exp(u))) + exp(exp(d) * (days - exp(u))),
    # Parameters are linearly dependent on dose
    p <- p0[ID] + p1 * dose_log,
    g <- g0[ID] + g1 * dose_log,
    u <- u0[ID] + u1 * dose_log,
    d <- d0[ID] + d1 * dose_log,
    # Adaptive priors for random intercepts
    p0[ID] ~ dnorm(p0mu, p0sd),
    g0[ID] ~ dnorm(g0mu, g0sd),
    u0[ID] ~ dnorm(u0mu, u0sd),
    d0[ID] ~ dnorm(d0mu, d0sd),
    # Fixed priors for remaining parameters
    p0mu ~ dnorm(25, 5),
    p0sd ~ dcauchy(0, 1),
    g0mu ~ dnorm(3, 1),
    g0sd ~ dcauchy(0, 1),
    u0mu ~ dnorm(0, 1),
    u0sd ~ dcauchy(0, 1),
    d0mu ~ dnorm(-1, 0.5),
    d0sd ~ dcauchy(0, 1),
    p1 ~ dnorm(0, 1),
    g1 ~ dnorm(0, 0.5),
    u1 ~ dnorm(0, 0.5),
    d1 ~ dnorm(0, 0.5),
    k ~ dgamma2(20.01, 10.01),
    sigma ~ dcauchy(0, 1)
  ),
  constraints = list(
    # Variance parameters must be constrained to be > 0.
    p0sd  = "lower=0",
    g0sd  = "lower=0",
    u0sd  = "lower=0",
    d0sd  = "lower=0",
    sigma = "lower=0",
    # brms and Andrew Gelman recommend a >1 constraint for k to ensure
    # the first moment exists.
    k     = "lower=1"
  ),
  # These are Stan parameters that I defined as global constants.
  log_lik = RTN_LL,
  cmdstan = USE_CMDSTAN,
  control = list(
    adapt_delta = DELTA,
    max_treedepth = MAX_TD
  ),
  chains = N_CHAINS,
  cores = N_CORES,
  warmup = N_WU,
  iter = N_ITER,
  seed = SEED
)

# Save model results to disk
saveRDS(
  object = timeseries_con_model,
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "res", "timeseries_con_model.Rds")
)

# Viral kinetics (longitudinal analysis) - Categorical Dose ####################

timeseries_cat_dat <-
  norovirus_data$ftsshedding |>
  dplyr::select(ID, dose_cat, days, logyc) |>
  dplyr::filter(days > 0)

# Save processed data to file
saveRDS(
  object = timeseries_cat_dat,
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "data", "timeseries_cat_data.Rds")
)

# Fit the model
timeseries_cat_model <- rethinking::ulam(
  data = timeseries_cat_dat,
  # This part defines the model
  flist = alist(
    # Distribution of outcome
    logyc ~ dstudent(k, mu, sigma),
    # Overall time series equation (I let u = T since I don't want to
    # overwrite T)
    mu <- log(top / bottom),
    top <- 2 * exp( p ),
    bottom <- L + R,
    L <- exp( -exp( g ) * ( days - exp( u ) ) ),
    R <- exp(  exp( d ) * ( days - exp( u ) ) ),
    # Parameters are linearly dependent on dose
    p <- p0[ID] + p1[dose_cat],
    g <- g0[ID] + g1[dose_cat],
    u <- u0[ID] + u1[dose_cat],
    d <- d0[ID] + d1[dose_cat],
    # Adaptive priors for random intercepts
    p0[ID] ~ dcauchy(0, gamma_p),
    g0[ID] ~ dcauchy(0, gamma_g),
    u0[ID] ~ dcauchy(0, gamma_u),
    d0[ID] ~ dcauchy(0, gamma_d),
    # Fixed priors for remaining parameters
    gamma_p ~ dcauchy(0, 1),
    gamma_g ~ dcauchy(0, 1),
    gamma_u ~ dcauchy(0, 1),
    gamma_d ~ dcauchy(0, 1),
    p1[dose_cat] ~ dnorm(25,   5),
    g1[dose_cat] ~ dnorm( 3,   1),
    u1[dose_cat] ~ dnorm( 0,   1),
    d1[dose_cat] ~ dnorm(-1, 0.5),
    # Rethinking uses mu, scale parametrization of gamma distribution vs.
    # shape, scale in brms. transformation is mu = shape/scale.
    k ~ dgamma2(mu = 20.01, scale = 10.01),
    sigma ~ dcauchy(0, 1)
  ),
  # Variance parameters must be constrained to be > 0.
  constraints = list(
    gamma_p = "lower=0",
    gamma_g = "lower=0",
    gamma_u = "lower=0",
    gamma_d = "lower=0",
    sigma = "lower=0"
  ),
  # These are Stan parameters that I defined as global constants.
  log_lik = RTN_LL,
  cmdstan = USE_CMDSTAN,
  control = list(
    adapt_delta = DELTA,
    max_treedepth = MAX_TD
  ),
  chains = N_CHAINS,
  cores = N_CORES,
  warmup = N_WU,
  iter = N_ITER,
  seed = SEED
)

# Save model results to disk
saveRDS(
  object = timeseries_cat_model,
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "res", "timeseries_cat_model.Rds")
)

# Incubation period - Continuous Dose ##########################################

# Read in processed data
incubation_data <-
  norovirus_data$incubation |>
  dplyr::select(ID, dose_log, incubation)

# Save processed data to file
saveRDS(
  object = incubation_data,
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "data", "incubation_con_data.Rds")
)

# Fit the model
incubation_model <- rethinking::ulam(
  data = incubation_data,
  # This part defines the model
  flist = alist(
    # Distribution of outcome
    incubation ~ dlnorm(mu, sigma),
    mu <- alpha[ID] + beta * dose_log,
    # Use an adaptive prior for random intercepts
    alpha[ID] ~ dnorm(delta, gamma),
    # Fixed priors for other parameters
    delta ~ dnorm(1, 5),
    gamma ~ dcauchy(0, 2),
    beta ~ dnorm(0, 1),
    sigma ~ dcauchy(0, 2)
  ),
  # Variance parameters must be constrained to be > 0. This is equivalent
  # to specifying them with a half Cauchy distribution.
  constraints = list(
    sigma = "lower=0",
    gamma = "lower=0"
  ),
  # These are Stan parameters that I defined as global constants.
  log_lik = RTN_LL,
  cmdstan = USE_CMDSTAN,
  control = list(
    adapt_delta = DELTA,
    max_treedepth = MAX_TD
  ),
  chains = N_CHAINS,
  cores = N_CORES,
  warmup = N_WU,
  iter = N_ITER,
  seed = SEED
)

# Save model results to disk
saveRDS(
  object = incubation_model,
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "res", "incubation_model.Rds")
)

# Modified Vesikari Score (MVS) - Continuous dose ##############################

mvs_dat <-
  norovirus_data$mvs |>
  dplyr::select(ID, dose_log, mvs)

# Save processed data to file
saveRDS(
  object = mvs_dat,
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "data", "mvs_con_data.Rds")
)

# Fit the model
mvs_model <- rethinking::ulam(
  data = mvs_dat,
  # This part defines the model
  flist = alist(
    # Distribution of outcome
    mvs ~ dgampois(exp(lambda), phi),
    # Lambda is a linear function of the dose, intercepts are random
    lambda <- alpha[ID] + beta * dose_log,
    # Use an adaptive prior for random intercepts
    alpha[ID] ~ dnorm(mu, gamma),
    # Fixed priors for other parameters
    mu ~ dnorm(1, 5),
    gamma ~ dcauchy(0, 2),
    beta ~ dnorm(0, 1),
    phi ~ dcauchy(0, 2)
  ),
  # Variance parameters must be constrained to be > 0. This is equivalent
  # to specifying them with a half Cauchy distribution.
  constraints = list(
    phi = "lower=0",
    gamma = "lower=0"
  ),
  # These are Stan parameters that I defined as global constants.
  log_lik = RTN_LL,
  cmdstan = USE_CMDSTAN,
  control = list(
    adapt_delta = DELTA,
    max_treedepth = MAX_TD
  ),
  chains = N_CHAINS,
  cores = N_CORES,
  warmup = N_WU,
  iter = N_ITER,
  seed = SEED
)

# Save model results to disk
saveRDS(
  object = mvs_model,
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "res", "mvs_model.Rds")
)

# Comprehensive Symptom Score (CSS) - Continuous dose###########################

css_dat <-
  norovirus_data$css |>
  dplyr::group_by(ID, dose_log) |>
  dplyr::summarize(outcome = sum(css), .groups = "drop")

# Save processed data to file
saveRDS(
  object = css_dat,
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "data", "css_con_data.Rds")
)

# Fit the model
css_model <- rethinking::ulam(
  data = css_dat,
  # This part defines the model
  flist = alist(
    # Distribution of outcome
    outcome ~ dgampois(exp(lambda), phi),
    # Lambda is a linear function of the dose, intercepts are random
    lambda <- alpha[ID] + beta * dose_log,
    # Use an adaptive prior for random intercepts
    alpha[ID] ~ dnorm(mu, gamma),
    # Fixed priors for other parameters
    mu ~ dnorm(1, 5),
    gamma ~ dcauchy(0, 2),
    beta ~ dnorm(0, 1),
    phi ~ dcauchy(0, 2)
  ),
  # Variance parameters must be constrained to be > 0. This is equivalent
  # to specifying them with a half Cauchy distribution.
  constraints = list(
    phi = "lower=0",
    gamma = "lower=0"
  ),
  # These are Stan parameters that I defined as global constants.
  log_lik = RTN_LL,
  cmdstan = USE_CMDSTAN,
  control = list(
    adapt_delta = DELTA,
    max_treedepth = MAX_TD
  ),
  chains = N_CHAINS,
  cores = N_CORES,
  warmup = N_WU,
  iter = N_ITER,
  seed = SEED
)

# Save model results to disk
saveRDS(
  object = css_model,
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "res", "css_model.Rds")
)

# Incubation period  - Categorical Dose ########################################

incubation_cat_dat <-
  norovirus_data$incubation |>
  dplyr::select(ID, dose_cat, incubation)

# Save processed data to file
saveRDS(
  object = incubation_cat_dat,
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "data", "incubation_cat_data.Rds")
)

# Fit the model
incubation_cat_model <- rethinking::ulam(
  data = incubation_cat_dat,
  # This part defines the model
  flist = alist(
    # Distribution of outcome
    incubation ~ dlnorm(mu, sigma),
    # Mu is a linear function of the dose, intercepts are random
    mu <- alpha[ID] + beta[dose_cat],
    # Use an adaptive prior for random intercepts
    alpha[ID] ~ dnorm(0, gamma),
    # Fixed priors for other parameters
    gamma ~ dcauchy(0, 2),
    beta[dose_cat] ~ dnorm(1, 5),
    sigma ~ dcauchy(0, 2)
  ),
  # Variance parameters must be constrained to be > 0. This is equivalent
  # to specifying them with a half Cauchy distribution.
  constraints = list(
    sigma = "lower=0",
    gamma = "lower=0"
  ),
  # These are Stan parameters that I defined as global constants.
  log_lik = RTN_LL,
  cmdstan = USE_CMDSTAN,
  control = list(
    adapt_delta = DELTA,
    max_treedepth = MAX_TD
  ),
  chains = N_CHAINS,
  cores = N_CORES,
  warmup = N_WU,
  iter = N_ITER,
  seed = SEED
)

# Save model results to disk
saveRDS(
  object = incubation_cat_model,
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "res", "incubation_cat_model.Rds")
)

# Modified Vesikari score (MVS) - Categorical dose #############################

# Read in processed data
mvs_cat_dat <-
  norovirus_data$mvs |>
  dplyr::select(ID, dose_cat, mvs)

# Save processed data to file
saveRDS(
  object = mvs_cat_dat,
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "data", "mvs_cat_data.Rds")
)

# Fit the model
mvs_cat_model <- rethinking::ulam(
  data = mvs_cat_dat,
  # This part defines the model
  flist = alist(
    # Distribution of outcome
    mvs ~ dgampois(exp(lambda), phi),
    # Lambda is a linear function of the dose, intercepts are random
    lambda <- alpha[ID] + beta[dose_cat],
    # Use an adaptive prior for random intercepts
    alpha[ID] ~ dnorm(0, gamma),
    # Fixed priors for other parameters
    gamma ~ dcauchy(0, 2),
    beta[dose_cat] ~ dnorm(1, 5),
    phi ~ dcauchy(0, 2)
  ),
  # Variance parameters must be constrained to be > 0. This is equivalent
  # to specifying them with a half Cauchy distribution.
  constraints = list(
    phi = "lower=0",
    gamma = "lower=0"
  ),
  # These are Stan parameters that I defined as global constants.
  log_lik = RTN_LL,
  cmdstan = USE_CMDSTAN,
  control = list(
    adapt_delta = DELTA,
    max_treedepth = MAX_TD
  ),
  chains = N_CHAINS,
  cores = N_CORES,
  warmup = N_WU,
  iter = N_ITER,
  seed = SEED
)

# Save model results to disk
saveRDS(
  object = mvs_cat_model,
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "res", "mvs_cat_model.Rds")
)

# Comprehensive Symptom Score (CSS) - Categorical data #########################

css_cat_dat <-
  norovirus_data$css |>
  dplyr::group_by(ID, dose_cat) |>
  dplyr::summarize(outcome = sum(css), .groups = "drop")

# Save processed data to file
saveRDS(
  object = css_cat_dat,
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "data", "css_cat_data.Rds")
)

# Fit the model
css_cat_model <- rethinking::ulam(
  data = css_cat_dat,
  # This part defines the model
  flist = alist(
    # Distribution of outcome
    outcome ~ dgampois(exp(lambda), phi),
    # Lambda is a linear function of the dose, intercepts are random
    lambda <- alpha[ID] + beta[dose_cat],
    # Use an adaptive prior for random intercepts
    alpha[ID] ~ dnorm(0, gamma),
    # Fixed priors for other parameters
    gamma ~ dcauchy(0, 2),
    beta[dose_cat] ~ dnorm(0, 1),
    phi ~ dcauchy(0, 2)
  ),
  # Variance parameters must be constrained to be > 0. This is equivalent
  # to specifying them with a half Cauchy distribution.
  constraints = list(
    phi = "lower=0",
    gamma = "lower=0"
  ),
  # These are Stan parameters that I defined as global constants.
  log_lik = RTN_LL,
  cmdstan = USE_CMDSTAN,
  control = list(
    adapt_delta = DELTA,
    max_treedepth = MAX_TD
  ),
  chains = N_CHAINS,
  cores = N_CORES,
  warmup = N_WU,
  iter = N_ITER,
  seed = SEED
)

# Save model results to disk
saveRDS(
  object = css_cat_model,
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "res", "css_cat_model.Rds")
)

###
# End of File ##################################################################
###
