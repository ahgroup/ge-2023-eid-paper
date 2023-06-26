# Packages
# ===============================================
library(dplyr)
library(tidyr)
library(here)
library(brms)
library(foreach)
library(doParallel)
library(pracma)

# ===============================================
rm(list=ls())

source(here("2_supplementary_materials/analyses/common_functions/brms_setting.R"))
source(here("2_supplementary_materials/analyses/common_functions/newdata4mods.R"))
source(here("2_supplementary_materials/analyses/common_functions/holder_brms_fit_con.R"))
source(here("2_supplementary_materials/analyses/common_functions/holder_brms_fit_op_con.R"))
source(here("2_supplementary_materials/analyses/common_functions/fun_virus_growth_rate.R"))

norovirus_data <- readRDS(here("2_supplementary_materials/data/norodata.rds"))

prior_holder <- c(prior(cauchy(0,1),   class = "sigma"),
                  prior(gamma(2, 0.1), class = "nu"), # https://mc-stan.org/docs/2_24/functions-reference/gamma-distribution.html
                  # peak v
                  prior(normal(25,5), class = "b",  nlpar = "Vp", coef = "Intercept"),
                  prior(normal(0,1),   class = "b",  nlpar = "Vp", coef = "dose_log"),
                  prior(cauchy(0,1),   class = "sd", nlpar = "Vp"),
                  # lambdag
                  prior(normal(3,1), class = "b",  nlpar = "lambdag", coef = "Intercept"),
                  prior(normal(0,0.5),  class = "b",  nlpar = "lambdag", coef = "dose_log"),
                  prior(cauchy(0,1),  class = "sd", nlpar = "lambdag"),
                  # peak time
                  prior(normal(0,1), class = "b",  nlpar = "tpeak", coef = "Intercept"),
                  prior(normal(0,0.5), class = "b",  nlpar = "tpeak", coef = "dose_log"),
                  prior(cauchy(0,1), class = "sd", nlpar = "tpeak"),
                  # lambdad
                  prior(normal(-1,0.5), class = "b",  nlpar = "lambdad", coef = "Intercept"),
                  prior(normal(0,0.5),  class = "b",  nlpar = "lambdad", coef = "dose_log"),
                  prior(cauchy(0,1),  class = "sd", nlpar = "lambdad"))

my_days <- c(seq(1/24,7,1/24), seq(8, 90, 1))

# ===============================================
# Without 0.48
# ===============================================
feces_ts19 <- norovirus_data$ftsshedding %>% filter(dose != 0.48, days > 0) %>% droplevels()

# start to fit the model
brms_holder19 <- brm(bf(logyc|cens(y_censor,y2=logy) ~ log(2*exp(Vp) / (exp(-exp(lambdag)*(days-exp(tpeak))) + exp(exp(lambdad)*(days-exp(tpeak))))),
                        lambdad ~ dose_log + (1|ID),
                        Vp ~      dose_log + (1|ID),
                        lambdag ~ dose_log + (1|ID),
                        tpeak ~   dose_log + (1|ID),
                        nl = TRUE),
                     data = feces_ts19,
                     family = student(),  
                     prior = prior_holder,
                     warmup = mywarmup, iter = myiter,
                     chains = mychains, cores = mycores, 
                     seed = myseed,
                     backend = "cmdstanr", 
                     sample_prior = TRUE,
                     silent = 1, refresh = 2000,
                     control = list(adapt_delta = myadapt_delta, 
                                    max_treedepth = mymax_treedepth))

saveRDS(brms_holder19, here("2_supplementary_materials/analyses/rds/brms_holder_con19.rds"))
# brms_holder19 <- readRDS(here("2_supplementary_materials/analyses/rds/brms_holder_con19.rds"))


# to take the individual randomness into account, for individual goodness of
# fit, we used IDs in the raw data, for the model prediction, we used a
# hypothetical individual (ID0).

# prediction
data4prediction19 <- crossing(days = my_days, newdata_cnt19)
holder_curves_mu_con19 <- holder_brms_fit_con(data = data4prediction19, 
                                              brms_fit = brms_holder19, 
                                              lod1 = 15e3, 
                                              lod2 = 40e6, 
                                              mod = "Posterior_grandmean")
saveRDS(holder_curves_mu_con19, here("2_supplementary_materials/analyses/rds/holder_curves_mu_con19.rds"))

# goodness of fit
data4gof19 <- crossing(days = my_days, feces_ts19 %>% distinct(dose_log, dose, dose_cat, ID))
holder_curves_mu_id_con19 <- holder_brms_fit_con(data = data4gof19, 
                                                 brms_fit = brms_holder19, 
                                                 lod1 = 15e3, 
                                                 lod2 = 40e6, 
                                                 mod = "GoodnessOfFit")

saveRDS(holder_curves_mu_id_con19, here("2_supplementary_materials/analyses/rds/holder_curves_mu_id_con19.rds"))
# holder_curves_mu_id_con19 <- readRDS(here("2_supplementary_materials/analyses/rds/holder_curves_mu_id_con19.rds"))

# ---
# Calculate auc for id ----
# ---

fun_auc_id4days <- function(x){
  x4days   <- x %>% filter(days <= 4) # need to be same to figure 1: <= 4 days
  x4days_y <- x4days %>% filter(Estimate >= 0) # auc curve for the part virus >= 0
  auc4days <- trapz(x4days_y$days, x4days_y$Estimate)
  op <- data.frame(ID = unique(x4days_y$ID), AUC = auc4days)
  return(op)
} 

# ---
# Calculate auc for id ----
# ---
auc_id4days <- lapply(holder_curves_mu_id_con19, fun_auc_id4days) %>% # unit = GEC
  bind_rows()

# ---
# total virus shed in feces ----
# ---
feces96 <- norovirus_data$feces96 %>% filter(dose != 0.48) %>% droplevels()


# ---
# corrleation between auc and total shed ----
# ---
feces96_2shed <- left_join(feces96, auc_id4days)

saveRDS(feces96_2shed, here("2_supplementary_materials/analyses/rds/feces96_2shed.rds"))

myprior = c(prior(cauchy(0,2), class = "sigma"),
            prior(normal(0,5), class = "b"),
            prior(normal(35,10), class = "Intercept"),
            prior(cauchy(0,2),  class = "sd"))

feces96_2shed_brms <- brm(AUC ~ outcome_log + (1 | ID), data = feces96_2shed, 
                          family = gaussian(), prior = myprior,
                          backend = "cmdstanr", silent = 1, refresh = 5000,
                          warmup = 25000, iter = 30000, 
                          chains = mychains, cores = mycores, seed = 123456789,
                          control=list(adapt_delta=0.999999, max_treedepth = 20))

saveRDS(feces96_2shed_brms, here("2_supplementary_materials/analyses/rds/feces96_2shed_cor_brms.rds"))


# rate of growth for each ID. log(peak value within 4 days)/(time of peak within 4 days)
virus_growth_rate_19_list <- lapply(X = holder_curves_mu_id_con19, FUN = fun_virus_growth_rate)
virus_growth_rate_19 <- bind_rows(virus_growth_rate_19_list)

saveRDS(virus_growth_rate_19, here("2_supplementary_materials/analyses/rds/virus_growth_rate_con19.rds"))


# model output 2nd results, similar to holder_curves_mu19, but with time points
# (days) in more fine grids.
holder_brms_fit_op_con19 <- holder_brms_fit_op_con(data = feces_ts19, 
                                                   brms_fit = brms_holder19, 
                                                   dose_log = newdata_cnt19)

saveRDS(holder_brms_fit_op_con19, here("2_supplementary_materials/analyses/rds/holder_brms_fit_op_con19.rds"))


# Prior v.s. Posteriors
mod_prior19 <- prior_draws(brms_holder19) %>% as_tibble() %>% select(-sigma,-nu,-sd_ID,-sd_ID__1,-sd_ID__2,-sd_ID__3) %>% mutate(type = "Prior")
mod_post_prior19 <- as_draws_df(brms_holder19) %>% as_tibble() %>% mutate(type = "Posteriors") %>% select(colnames(mod_prior19))


model_con_priors19_compare <- bind_rows(mod_prior19, mod_post_prior19) %>%
  pivot_longer(cols = b_lambdad_Intercept:b_tpeak_dose_log,
               names_to = "Parameters",
               values_to = "value")

saveRDS(model_con_priors19_compare, here("2_supplementary_materials/analyses/rds/model_con_priors19_compare.rds"))


# ===============================================
# With 0.48
# ===============================================
feces_ts20 <- norovirus_data$ftsshedding %>% filter(days > 0) %>% droplevels()

# start to fit the model
brms_holder20 <- brm(bf(logyc|cens(y_censor,y2=logy) ~ log(2*exp(Vp) / (exp(-exp(lambdag)*(days-exp(tpeak))) + exp(exp(lambdad)*(days-exp(tpeak))))),
                        lambdad ~ dose_log + (1|ID),
                        Vp ~      dose_log + (1|ID),
                        lambdag ~ dose_log + (1|ID),
                        tpeak ~   dose_log + (1|ID),
                        nl = TRUE),
                     data = feces_ts20,
                     family = student(),  
                     prior = prior_holder,
                     warmup = mywarmup, iter = myiter,
                     chains = mychains, cores = mycores, 
                     seed = myseed,
                     sample_prior = TRUE,
                     backend = "cmdstanr", 
                     silent = 1, refresh = 2000,
                     control = list(adapt_delta = myadapt_delta, 
                                    max_treedepth = mymax_treedepth))

saveRDS(brms_holder20, here("2_supplementary_materials/analyses/rds/brms_holder_con20.rds"))
# brms_holder20 <- readRDS(here("2_supplementary_materials/analyses/rds/brms_holder_con20.rds"))


# prediction
data4prediction20 <- crossing(days = my_days, newdata_cnt20)
holder_curves_mu_con20 <- holder_brms_fit_con(data = data4prediction20, 
                                              brms_fit = brms_holder20, 
                                              lod1 = 15e3, 
                                              lod2 = 40e6, 
                                              mod = "Posterior_grandmean")

saveRDS(holder_curves_mu_con20, here("2_supplementary_materials/analyses/rds/holder_curves_mu_con20.rds"))

# goodness of fit
data4gof20 <- crossing(days = my_days, feces_ts20 %>% distinct(dose_log, dose, dose_cat, ID))
holder_curves_mu_id_con20 <- holder_brms_fit_con(data = data4gof20, 
                                                 brms_fit = brms_holder20, 
                                                 lod1 = 15e3, 
                                                 lod2 = 40e6, 
                                                 mod = "GoodnessOfFit")

saveRDS(holder_curves_mu_id_con20, here("2_supplementary_materials/analyses/rds/holder_curves_mu_id_con20.rds"))
# holder_curves_mu_id_con20 <- readRDS(here("2_supplementary_materials/analyses/rds/holder_curves_mu_id_con20.rds"))

# rate of growth for each ID. log(peak value within 4 days)/(time of peak within 4 days)
virus_growth_rate_20_list <- lapply(X = holder_curves_mu_id_con20, FUN = fun_virus_growth_rate)
virus_growth_rate_20 <- bind_rows(virus_growth_rate_20_list)

saveRDS(virus_growth_rate_20, here("2_supplementary_materials/analyses/rds/virus_growth_rate_con20.rds"))


# model output 2nd results, similar to holder_curves_mu19, but with time points
# (days) in more fine grids.
holder_brms_fit_op_con20 <- holder_brms_fit_op_con(data = feces_ts20, 
                                                   brms_fit = brms_holder20, 
                                                   dose_log = newdata_cnt20)

saveRDS(holder_brms_fit_op_con20, here("2_supplementary_materials/analyses/rds/holder_brms_fit_op_con20.rds"))

# Prior v.s. Posteriors

mod_prior20 <- prior_draws(brms_holder20) %>% as_tibble() %>% select(-sigma,-nu,-sd_ID,-sd_ID__1,-sd_ID__2,-sd_ID__3) %>% mutate(type = "Prior")
mod_post_prior20 <- as_draws_df(brms_holder20) %>% as_tibble() %>% mutate(type = "Posteriors") %>% select(colnames(mod_prior20))

model_con_priors20_compare <- bind_rows(mod_prior20, mod_post_prior20) %>%
  pivot_longer(cols = b_lambdad_Intercept:b_tpeak_dose_log,
               names_to = "Parameters",
               values_to = "value")

saveRDS(model_con_priors20_compare, here("2_supplementary_materials/analyses/rds/model_con_priors20_compare.rds"))
