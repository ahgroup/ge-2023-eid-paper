# Getting predictions ##########################################################
# Project: Analyzing the association between inoculum dose and Norovirus
#   infection outcomes
# Author: Zane Billings
# Since: 2022-04-11
# Desc: For this paper, Yang fit a bunch of bayesian regression models using
#   the brms package. I fitted the models and extracted predictions in the
#   models.R script. In this script, I will use the model results to get
#   predicted fits.

# Setting up dependencies ######################################################

# If I don't load the entire rethinking namespace, it will mess up some stan
# functionality. But declare all other dependencies without loading the
# namespace
library(rethinking)
box::use(
  arrow,
  dplyr,
  here,
  tibble,
  tidyr
)

# Seed variable
S <- 109189

# Function definitions #########################################################

# OLD FUNCTION FROM PREVIOUS CI METHOD -- NEED FOR TIME SERIES MODEL
# This function is necessary for getting predictions for each individual for
# the time series model. So it has not been deleted.
get_preds <- function(model, dat = NULL, ...) {
  
  # If .data is not specified, use the data stored in the Stan model object.
  if (is.null(dat)) {dat <- model@data}
  
  # Use rethinking::link to compute the predictions. Then get the
  # mean and credible interval of the estimate over all the samples?
  # I guess. I am not really sure what it is doing.
  # TODO make this a function instead of loose
  linkmod <- rethinking::link(model, data = dat, ...)
  if (mode(linkmod) == "list") {linkmod <- linkmod$mu}
  estimate_mean <- apply(linkmod, 2, mean)
  estimate_ci   <- apply(linkmod, 2, rethinking::PI, prob = 0.95)
  
  fit <- data.frame(
    estimate = estimate_mean,
    lwr = estimate_ci[1, ],
    upr = estimate_ci[2, ]
  )
  
  # Use the estimates to make a data frame of predictions. This is dependent
  # on the order that rethinking::link outputs things, so it has to be in the
  # same order as preddat since nothing in the output is labeled.
  out <- cbind(dat, fit) |> tibble::tibble()
  
  return(out)
  
}

# NEW FUNCTION USING BRMS CI METHOD

# According to brms documentation, "If "uncertainty" (default), each posterior
# sample for a new level is drawn from the posterior draws of a randomly chosen
# existing level. Each posterior sample for a new level may be drawn from a
# different existing level such that the resulting set of new posterior draws
# represents the variation across existing levels."
# So to calculate the mean and credible interval, I will create an overall
# alpha vector by randomly sampling across the 19 individual random intercepts
# for each posterior draw.
# Apparently this approximates integrating out the random effect.

get_t3_ci <- function(model, newdata, re, seed) {
  # ARGS
  #   Model is the result of running rethinking::ulam().
  #   Newdata is the data to be used for predictions.
  #   re is a character vector giving the names of the random effects
  #     components of the posterior sample to modify
  #   seed is the seed for the PRNG that will sample across levels.
  # VALUE a tibble containing posterior mean estimates with CI's of
  #   type 3.
  
  # Extract the posterior samples from the model
  post <- extract.samples(model, n = Inf)
  
  # Randomly sample the random effect
  set.seed(seed)
  for (i in 1:length(re)) {
    # Get next level
    this_re <- re[[i]]
    
    # Number of posterior samples
    N <- nrow(post[[this_re]])
    
    # Number of random effects levels
    k <- ncol(post[[this_re]])
    
    # Randomly sample the current re
    sampled_vec <-
      purrr::map_dbl(
        .x = 1:N,
        .f = ~post[[this_re]][.x, sample(1:k, size = 1)]
      )
    
    
    # Modify posterior samples to trick rethinking::link
    post[[this_re]] <- array(sampled_vec, dim = dim(post[[this_re]]))
  }
  
  # Compute the linear model term from the model
  linkmod <-
    rethinking::link(
      model,
      data = newdata,
      post = post
    )
  
  return(linkmod)
}

# Estimate the mean and CI of "Type 1"
# based on Andrew Heiss' blog post (link)
# We calculate the width of the credible interval for the "grand mean"
# (AKA estimate of hypothetical population mean) by using the estimated
# mean of the random effects, rather than accounting for variation in the
# individual random effects estimate.
get_t1_ci <- function(model, newdata, re, re_mean) {
  # ARGS
  #   Model is the result of running rethinking::ulam().
  #   Newdata is the data to be used for predictions.
  #   re is a character vector giving the names of the random effects
  #     components of the posterior sample to modify
  #   re_mean is a character vector of equal length giving the name of the
  #     component that should be used to replace the corresponding component
  #     of the re vector in the modified posterior
  # VALUE a tibble containing posterior mean estimates with CI's of
  #   type 1.
  
  # Extract the posterior samples from the model
  post <- extract.samples(model, n = Inf)
  

  
  # Modify the posterior to replace the random effects with their mean
  for (i in 1:length(re)) {
    # Get next level
    this_re <- re[[i]]
    
    # We need a special case where re_mean is allowed to be a number,
    # specifically the number 0, to handle some of the categorical models.
    if (mode(re_mean) == "numeric" && length(re_mean) == 1) {
      post[["spam"]] <- re_mean
      this_re_mean <- "spam"
    } else {
      this_re_mean <- re_mean[[i]]
    }
    
    # Number of posterior samples
    N <- nrow(post[[this_re]])
    
    # Number of random effects levels
    k <- ncol(post[[this_re]])
    
    # Modify posterior samples to trick rethinking::link
    post[[this_re]] <- array(post[[this_re_mean]], dim = dim(post[[this_re]]))
  }
  
  # Compute the linear model term from the model
  linkmod <-
    rethinking::link(
      model,
      data = newdata,
      post = post
    )
  
  return(linkmod)
}

get_mean_ci <- function(type, model, newdata, re, re_mean = NULL, seed = NULL) {
  
  if(type == 1) {
    linkmod <- get_t1_ci(model, newdata, re, re_mean)
  } else if (type == 3) {
    linkmod <- get_t3_ci(model, newdata, re, seed)
  } else {
    stop("type must be 1 or 3. see blog post for details.")
  }
  
  # If there is more than one deterministic expression in the ulam model code,
  # rethinking::link() returns a list with all of them calculated. But
  # we only want the one for mu.
  if (mode(linkmod) == "list") {
    linkmod <- linkmod$mu
  }
  
  # Compute mu and 89% equal-tail credible interval
  mean_est <- colMeans(linkmod)
  CI_est   <- apply(linkmod, 2, rethinking::PI, prob = 0.95)
  
  # Bind results into a data frame
  out <-
    dplyr::bind_cols(
      tibble::as_tibble(newdata),
      tibble::tibble(
        mu = mean_est,
        lwr = CI_est[1, ],
        upr = CI_est[2, ]
      )
    ) |>
    dplyr::select(-ID)
  
  return(out)
}

# Function for processing time series outcomes data
get_ts_outcomes <- function(model, data, type, re, re_mean = NULL, S = NULL) {
  
  # First use the designated CI type to get the estimates -- we need the
  # raw link output though.
  if (type == 1) {
    ts_mu_samples <-
      get_t1_ci(
        model = model,
        newdata = data,
        re = re,
        re_mean = re_mean
      )$mu
  } else if (type == 3) {
    ts_mu_samples <-
      get_t3_ci(
        model = model,
        newdata = data,
        re = re,
        seed = S
      )$mu
  } else {
    stop("type must be 1 or 3")
  }
  
  # Process the inputted data -- this does require the data to have columns
  # with specific names, but since I only need this function twice I am not
  # going to refactor it
  vec_days <- data$days |> unique()
  cat_chk <- model@model |> grepl(pattern = "dose_cat", fixed = TRUE)
  if (isTRUE(cat_chk)) {
    vec_dose <- data$dose_cat
  } else {
    vec_dose <- data$dose_log
  }
  
  # Nest the samples
  nested_samples <-
    ts_mu_samples |>
    exp() |>
    t() |>
    tibble::as_tibble(.name_repair = "none") |>
    rlang::set_names(paste0("sample_", 1:nrow(ts_mu_samples))) |>
    dplyr::mutate(
      dose = vec_dose,
      .before = tidyselect::everything()
    ) |>
    tidyr::nest(data = -dose)
  
  # Calculate the outcomes of interest using the nested samples
  auc_vec <-
    nested_samples$data |>
    purrr::map(~apply(.x, 2, function(k) log(pracma::trapz(vec_days, k))))
  
  peak_vec <-
    nested_samples$data |>
    purrr::map(~apply(.x, 2, function(k) vec_days[which.max(k)]))
  
  onset_vec <-
    nested_samples$data |>
    purrr::map(~apply(.x, 2, function(k) vec_days[which(k >= 225000)][1]))
  
  end_vec <-
    nested_samples$data |>
    purrr::map(~apply(.x, 2, function(k) rev(vec_days[which(k >= 225000)])[1]))
  
  duration_vec <- purrr::map2(.x = onset_vec, .y = end_vec, ~.y - .x)
  
  # Process the results
  out <-
    tibble::tibble(
      dose = nested_samples$dose,
      auc = auc_vec,
      peak = peak_vec,
      onset = onset_vec,
      duration = duration_vec
    ) |>
    dplyr::rowwise(dose) |>
    dplyr::summarise(
      dplyr::across(
        c(auc, peak, onset, duration),
        ~brms::posterior_summary(.x, probs = c(0.025, 0.975)) |>
          tibble::as_tibble() |>
          list()
      ),
      .groups = "drop"
    ) |>
    tidyr::pivot_longer(
      cols = -dose,
      names_to = "stat",
      values_to = "value"
    ) |>
    tidyr::unnest(value)
  
  if (isTRUE(cat_chk)) {
    out <- rlang::set_names(out,
                            c("dose_cat", "stat", "mean", "se", "lwr", "upr"))
  } else {
    out <- rlang::set_names(out,
                            c("dose_log", "stat", "mean", "se", "lwr", "upr"))
  }
  
  return(out)
}

# Create prediction data frames ################################################

# We need to create "fake" data frames that contain the levels of categories
# that we want to predict on. We'll need to do this twice, once for the
# categorical models and once for the continuous models.
# Note that we can use the same predicted data for all models except for the
# time series models.

# Continuous dose prediction data
con_dose_pred_data <-
  tidyr::expand_grid(
    ID = 1,
    dose_log = seq(from = log(4.8), to=log(4800), length = 4000) - log(48)
  )

# Categorical dose prediction data
cat_dose_pred_data <-
  tidyr::expand_grid(
    ID = 1,
    dose_cat = 1:3
  )

# Dose vector to predict on
dose_vec <- seq(from = log(4.8), to=log(4800), length.out = 4000) - log(48)

# Continuous dose models ######################################################

## Total fecal shedding <= 96 hours   ###########################

# Read in model object
tfs_96_con_model <- readRDS(here::here("2_supplementary_materials", "analyses", "rethinking-results", "res", "tfs_96_model.Rds"))

# Calculate the mean and credible interval.
tfs_96_con_fit <-
  get_mean_ci(
    model = tfs_96_con_model,
    newdata = con_dose_pred_data,
    type = 1,
    re = "alpha",
    re_mean = "mu_alpha",
    seed = S
  )

# Write to disk
arrow::write_parquet(
  tfs_96_con_fit,
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "tfs_96_con_fit.parquet")
)

## Total fecal shedding full time   #############################

tfs_all_con_model <- readRDS(here::here("2_supplementary_materials", "analyses", "rethinking-results", "res", "tfs_all_model.Rds"))

# Get marginal predictions from model -- see first model for explanation
tfs_all_con_fit <-
  get_mean_ci(
    model = tfs_all_con_model,
    newdata = con_dose_pred_data,
    type = 1,
    re = "alpha",
    re_mean = "mu_alpha",
    seed = S
  )

# Write to disk
arrow::write_parquet(
  tfs_all_con_fit,
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "tfs_all_con_fit.parquet")
)

## Total vomit shedding   #######################################

vvshed_model_all <- readRDS(here::here("2_supplementary_materials", "analyses", "rethinking-results", "res", "vvshed_model_all.Rds"))

# Get marginal predictions from model -- see first model for explanation
vvshed_con_fit <-
  get_mean_ci(
    model = vvshed_model_all,
    newdata = con_dose_pred_data,
    type = 1,
    re = "alpha",
    re_mean = "mu_alpha",
    seed = S
  )

arrow::write_parquet(
  vvshed_con_fit,
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "vvshed_con_fit.parquet")
)

## Viral kinetics (longitudinal analysis)   #####################

rm(list = ls.str(mode = "S4"))
gc()

timeseries_con_model <- readRDS(
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "res", "timeseries_con_model.Rds")
)

### Individual level predictions ################################################
# We use the posterior samples to predict the outcomes
# I based this on what Andreas did in his blog posts
# First create the "fake data" (analogous to newdata in stats::predict) that
# we will feed into the fitted model to get predictions.
ts_con_ind_pred_data <-
  expand.grid(
    ID = 1:19,
    # I think starting time at 0 will not work, so add a small number and then
    # create a regular time series.
    days = c(seq(from = 0.1, to = 7, by = 0.1),
             seq(from = 7, to = 90, by = 1))
  ) |>
  # Get the dose for each person
  dplyr::left_join(
    timeseries_con_model@data |>
      as.data.frame() |>
      dplyr::select(ID, dose_log) |>
      dplyr::distinct(),
    by = "ID"
  )

timeseries_con_fit_ind <- get_preds(timeseries_con_model, ts_con_ind_pred_data)

# Save predictions to disk
arrow::write_parquet(
  timeseries_con_fit_ind,
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "timeseries_con_ind_fit.parquet")
)

### Population level (marginal) predictions #####################################
# Get population level curve fits

ts_con_pop_pred_data <- ts_con_ind_pred_data
ts_con_pop_pred_data$ID <- rep(1, length = nrow(ts_con_pop_pred_data))

timeseries_con_fit_pop <-
  get_mean_ci(
    model = timeseries_con_model,
    newdata = ts_con_pop_pred_data,
    type = 1,
    re = c("p0", "g0", "u0", "d0"),
    re_mean = c("p0mu", "g0mu", "u0mu", "d0mu"),
    seed = S
  )

# Save population level curve fits to disk
arrow::write_parquet(
  timeseries_con_fit_pop,
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "timeseries_con_pop_fit.parquet")
)

### Outcomes derived from pop. level curves #####################################

# Get value of days to use
day_range <- range(timeseries_con_model@data$days)

# Create new predicted data with all the dose values of interest
gc()
outcomes_pred_data <-
  tidyr::expand_grid(
    ID = 1,
    dose_log = seq(from = log(4.8), to = log(4800), length.out = 50) - log(48),
    days = c(seq(from = 1/24, to = (7 - 0.01), by = 0.01),
             seq(from = 7, to = (30 - 0.25), by = 0.25),
             seq(from = 30, to = 90, by = 1))
  )

ts_outcomes_fitted <-
  get_ts_outcomes(
    model = timeseries_con_model,
    data = outcomes_pred_data,
    type = 1,
    re = c("p0", "g0", "u0", "d0"),
    re_mean = c("p0mu", "g0mu", "u0mu", "d0mu"),
    S = S
  )

# Save the results
arrow::write_parquet(
  ts_outcomes_fitted,
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "timeseries_con_stats.parquet")
)

gc()

## Incubation period   ##########################################

# Read in model object
incubation_model <- readRDS(here::here("2_supplementary_materials", "analyses", "rethinking-results", "res", "incubation_model.Rds"))

# We can use the shedding fake data for getting these predictions

# Calculate the mean and credible interval.
incubation_con_fit <-
  get_mean_ci(
    model = incubation_model,
    newdata = con_dose_pred_data,
    type = 1,
    re = "alpha",
    re_mean = "delta",
    seed = S
  )

arrow::write_parquet(
  incubation_con_fit,
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "incubation_con_fit.parquet")
)

## Modified Vesikari Score (MVS)   ##############################

mvs_model <- readRDS(here::here("2_supplementary_materials", "analyses", "rethinking-results", "res", "mvs_model.Rds"))

# Extract model samples, get predictions, and save to disk
# Calculate the mean and credible interval.
mvs_con_fit <-
  get_mean_ci(
    model = mvs_model,
    newdata = con_dose_pred_data,
    type = 1,
    re = "alpha",
    re_mean = "mu",
    seed = S
  )

arrow::write_parquet(
  mvs_con_fit,
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "mvs_con_fit.parquet")
)

## Comprehensive Symptom Score (CSS)  ###########################

css_model <- readRDS(here::here("2_supplementary_materials", "analyses", "rethinking-results", "res", "css_model.Rds"))

# Calculate the mean and credible interval.
css_con_fit <-
  get_mean_ci(
    model = css_model,
    newdata = con_dose_pred_data,
    type = 1,
    re = "alpha",
    re_mean = "mu",
    seed = S
  )


arrow::write_parquet(
  css_con_fit,
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "css_con_fit.parquet")
)


# Categorical dose models #####################################################

## Total fecal shedding <= 96 hours   ###########################

# Read in model object
tfs_96_cat_model <- readRDS(here::here("2_supplementary_materials", "analyses", "rethinking-results", "res", "tfs_96_cat_model.Rds"))

# Calculate the mean and credible interval.
tfs_96_cat_fit <-
  get_mean_ci(
    model = tfs_96_cat_model,
    newdata = cat_dose_pred_data,
    type = 1,
    re = "alpha",
    re_mean = 0, # The RIs have an assume mean of 0 in this model
    seed = S
  )

# Write to disk
arrow::write_parquet(
  tfs_96_cat_fit,
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "tfs_96_cat_fit.parquet")
)

## Total fecal shedding full time   #############################

tfs_all_cat_model <- readRDS(here::here("2_supplementary_materials", "analyses", "rethinking-results", "res", "tfs_model_all_cat.Rds"))

# Get marginal predictions from model -- see first model for explanation
tfs_all_cat_fit <-
  get_mean_ci(
    model = tfs_all_cat_model,
    newdata = cat_dose_pred_data,
    type = 1,
    re = "alpha",
    re_mean = 0, # The RIs have an assume mean of 0 in this model
    seed = S
  )

# Write to disk
arrow::write_parquet(
  tfs_all_cat_fit,
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "tfs_all_cat_fit.parquet")
)

## Total vomit shedding   #######################################

vvshed_model_all_cat <- readRDS(here::here("2_supplementary_materials", "analyses", "rethinking-results", "res", "vvshed_model_cat.Rds"))

# Get marginal predictions from model -- see first model for explanation
vvshed_cat_fit <-
  get_mean_ci(
    model = vvshed_model_all_cat,
    newdata = cat_dose_pred_data,
    type = 1,
    re = "alpha",
    re_mean = 0, # The RIs have an assume mean of 0 in this model
    seed = S
  )

arrow::write_parquet(
  vvshed_cat_fit,
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "vvshed_cat_fit.parquet")
)

## Viral kinetics (longitudinal analysis)   #####################

rm(list = ls.str(mode = "S4"))
gc()

timeseries_cat_model <- readRDS(
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "res", "timeseries_cat_model.Rds")
)

### Individual level predictions ################################################
# We use the posterior samples to predict the outcomes
# I based this on what Andreas did in his blog posts
# First create the "fake data" (analogous to newdata in stats::predict) that
# we will feed into the fitted model to get predictions.
ts_cat_ind_pred_data <-
  expand.grid(
    ID = 1:19,
    # I think starting time at 0 will not work, so add a small number and then
    # create a regular time series.
    days = c(seq(from = 0.1, to = 7, by = 0.1),
             seq(from = 7, to = 90, by = 1))
  ) |>
  # Get the dose for each person
  dplyr::left_join(
    timeseries_cat_model@data |>
      as.data.frame() |>
      dplyr::select(ID, dose_cat) |>
      dplyr::distinct(),
    by = "ID"
  )

timeseries_cat_fit_ind <- get_preds(timeseries_cat_model, ts_cat_ind_pred_data)

# Save predictions to disk
arrow::write_parquet(
  timeseries_cat_fit_ind,
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "timeseries_cat_ind_fit.parquet")
)

### Population level (marginal) predictions #####################################
# Get population level curve fits

ts_cat_pop_pred_data <- ts_cat_ind_pred_data
ts_cat_pop_pred_data$ID <- rep(1, length = nrow(ts_cat_pop_pred_data))

timeseries_cat_fit_pop <-
  get_mean_ci(
    model = timeseries_cat_model,
    newdata = ts_cat_pop_pred_data,
    type = 1,
    re = c("p0", "g0", "u0", "d0"),
    re_mean = 0,
    seed = S
  )

# Save population level curve fits to disk
arrow::write_parquet(
  timeseries_cat_fit_pop,
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "timeseries_cat_pop_fit.parquet")
)

### Outcomes derived from pop. level curves #####################################

# gc is used intermittently because my computer only had enough free RAM
# after I manually called the garbage collector. probably unnecessary but
# I can evaluate a longer time sequence this way

# Create new predicted data with all the dose values of interest
gc()
outcomes_cat_pred_data <-
  tidyr::expand_grid(
    ID = 1,
    dose_cat = 1:3,
    days = c(seq(from = 1/24, to = (7 - 0.01), by = 0.01),
             seq(from = 7, to = (30 - 0.25), by = 0.25),
             seq(from = 30, to = 90, by = 1))
  )

# Get the samples for the mean of each dose category
ts_cat_outcomes_fitted <-
  get_ts_outcomes(
    model = timeseries_cat_model,
    data = outcomes_cat_pred_data,
    type = 1,
    re = c("p0", "g0", "u0", "d0"),
    re_mean = 0,
    S = S
  )

# Save the results
arrow::write_parquet(
  ts_cat_outcomes_fitted,
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "timeseries_cat_stats.parquet")
)

rm(list = ls.str(mode = "S4"))
gc()

## Incubation period   ##########################################

# Read in model object
incubation_cat_model <- readRDS(here::here("2_supplementary_materials", "analyses", "rethinking-results", "res", "incubation_cat_model.Rds"))

# We can use the shedding fake data for getting these predictions

# Calculate the mean and credible interval.
incubation_cat_fit <-
  get_mean_ci(
    model = incubation_cat_model,
    newdata = cat_dose_pred_data,
    type = 1,
    re = "alpha",
    re_mean = 0, # The RIs have an assume mean of 0 in this model
    seed = S
  )

arrow::write_parquet(
  incubation_cat_fit,
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "incubation_cat_fit.parquet")
)

## Modified Vesikari Score (MVS)   ##############################

mvs_cat_model <- readRDS(here::here("2_supplementary_materials", "analyses", "rethinking-results", "res", "mvs_cat_model.Rds"))

# Extract model samples, get predictions, and save to disk
# Calculate the mean and credible interval.
mvs_cat_fit <-
  get_mean_ci(
    model = mvs_cat_model,
    newdata = cat_dose_pred_data,
    type = 1,
    re = "alpha",
    re_mean = 0, # The RIs have an assume mean of 0 in this model
    seed = S
  )

arrow::write_parquet(
  mvs_cat_fit,
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "mvs_cat_fit.parquet")
)

## Comprehensive Symptom Score (CSS)  ###########################

css_cat_model <- readRDS(here::here("2_supplementary_materials", "analyses", "rethinking-results", "res", "css_cat_model.Rds"))

# Calculate the mean and credible interval.
css_cat_fit <-
  get_mean_ci(
    model = css_cat_model,
    newdata = cat_dose_pred_data,
    type = 1,
    re = "alpha",
    re_mean = 0, # The RIs have an assume mean of 0 in this model
    seed = S
  )


arrow::write_parquet(
  css_cat_fit,
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "css_cat_fit.parquet")
)

# End of File ##################################################################
