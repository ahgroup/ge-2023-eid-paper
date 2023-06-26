# Packages
# ===============================================
library(dplyr)
library(tidyr)
library(here)
library(brms)
library(foreach)
library(doParallel)
library(pracma)
library(stringr)

# ===============================================
rm(list=ls())

source(here("2_supplementary_materials/analyses/common_functions/brms_setting.R"))
source(here("2_supplementary_materials/analyses/common_functions/holder_brms_fit_cat.R"))
source(here("2_supplementary_materials/analyses/common_functions/holder_brms_fit_op_cat.R"))
source(here("2_supplementary_materials/analyses/common_functions/fun_virus_growth_rate.R"))

norovirus_data <- readRDS(here("2_supplementary_materials/data/norodata.rds"))

prior_holder_cat <- c(prior(cauchy(0,1),   class = "sigma"),
                      prior(gamma(2, 0.1), class = "nu"),
                      # peak v
                      prior(normal(25,5), class = "b",  nlpar = "Vp"),
                      prior(cauchy(0,1),   class = "sd", nlpar = "Vp"),
                      # lambdag
                      prior(normal(3,1), class = "b",  nlpar = "lambdag"),
                      prior(cauchy(0,1),  class = "sd", nlpar = "lambdag"),
                      # peak time
                      prior(normal(0,1), class = "b",  nlpar = "tpeak"),
                      prior(cauchy(0,1), class = "sd", nlpar = "tpeak"),
                      # lambdad
                      prior(normal(-1,0.5), class = "b",  nlpar = "lambdad"),
                      prior(cauchy(0,1),  class = "sd", nlpar = "lambdad"))

my_dose_cat <-  c("Low","Medium","High")

my_days <- c(seq(1/24,7,1/24), seq(8, 90, 1))

# ===============================================
# Without 0.48
# ===============================================

feces_ts19 <- norovirus_data$ftsshedding %>% filter(dose != 0.48, days > 0) %>% droplevels()

# start to fit the model
brms_holder_cat <- brm(bf(logyc|cens(y_censor,y2=logy) ~ log(2*exp(Vp) / (exp(-exp(lambdag)*(days-exp(tpeak))) + exp(exp(lambdad)*(days-exp(tpeak))))),
                          lambdad ~ 0 + dose_cat + (1|ID),
                          Vp ~      0 + dose_cat + (1|ID),
                          lambdag ~ 0 + dose_cat + (1|ID),
                          tpeak ~   0 + dose_cat + (1|ID),
                          nl = TRUE),
                       data = feces_ts19,
                       family = student(),  
                       prior = prior_holder_cat,
                       warmup = mywarmup, iter = myiter, 
                       chains = mychains, cores = mycores, 
                       seed = myseed,
                       sample_prior = TRUE,
                       backend = "cmdstanr", silent = 1, refresh = 2000,
                       control = list(adapt_delta = myadapt_delta, 
                                      max_treedepth = mymax_treedepth))

saveRDS(brms_holder_cat, here("2_supplementary_materials/analyses/rds/brms_holder_cat19.rds"))
# brms_holder_cat <- readRDS(here("2_supplementary_materials/analyses/rds/brms_holder_cat19.rds"))

# to take the individual randomness into account, for individual goodness of
# fit, we used IDs in the raw data, for the model prediction, we used a
# hypothetical individual (ID0).

# prediction
data4prediction <- crossing(days = my_days, dose_cat = my_dose_cat)
holder_curves_mu19 <- holder_brms_fit_cat(data = data4prediction, 
                                          brms_fit = brms_holder_cat, 
                                          lod1 = 15e3, 
                                          lod2 = 40e6, 
                                          mod = "Posterior_grandmean")
saveRDS(holder_curves_mu19, here("2_supplementary_materials/analyses/rds/holder_curves_mu_cat19.rds"))

# goodness of fit
data4gof19 <- crossing(days = my_days, feces_ts19 %>% distinct(dose_cat, ID))
holder_curves_mu_id_cat19 <- holder_brms_fit_cat(data = data4gof19, 
                                             brms_fit = brms_holder_cat, 
                                             lod1 = 15e3, 
                                             lod2 = 40e6, 
                                             mod = "GoodnessOfFit")
saveRDS(holder_curves_mu_id_cat19, here("2_supplementary_materials/analyses/rds/holder_curves_mu_id_cat19.rds"))
# holder_curves_mu_id_cat19 <- readRDS(here("2_supplementary_materials/analyses/rds/holder_curves_mu_id_cat19.rds"))

# rate of growth for each ID. log(peak value within 4 days)/(time of peak within 4 days)
virus_growth_rate_19_list <- lapply(X = holder_curves_mu_id_cat19, FUN = fun_virus_growth_rate)
virus_growth_rate_19 <- bind_rows(virus_growth_rate_19_list)

saveRDS(virus_growth_rate_19, here("2_supplementary_materials/analyses/rds/virus_growth_rate_cat19.rds"))


# model output 2nd results, similar to holder_curves_mu19, but with time points
# (days) in more fine grids.
holder_brms_fit_op_cat19 <- holder_brms_fit_op_cat(data = feces_ts19, 
                                                   brms_fit = brms_holder_cat, 
                                                   dose_cat = my_dose_cat)
saveRDS(holder_brms_fit_op_cat19, here("2_supplementary_materials/analyses/rds/holder_brms_fit_op_cat19.rds"))


# Prior v.s. Posteriors

brms_holder_cat_prior <- prior_draws(brms_holder_cat) %>% 
  as_tibble() %>% 
  select(-sigma,-nu,-sd_ID,-sd_ID__1,-sd_ID__2,-sd_ID__3) %>% 
  mutate(type = "Prior") %>% 
  pivot_longer(cols = b_lambdad:b_tpeak,
               names_to = "Parameters",
               values_to = "value")

brms_holder_cat_post_prior <- as_draws_df(brms_holder_cat) %>% 
  as_tibble() %>% 
  select(b_lambdad_dose_catLow:b_tpeak_dose_catHigh) %>% 
  pivot_longer(cols = b_lambdad_dose_catLow:b_tpeak_dose_catHigh,
               names_to = "Outputs",
               values_to = "value") %>% 
  mutate(Parameters = str_remove(Outputs,"_dose_catLow|_dose_catMedium|_dose_catHigh"),
         Group = str_remove(Outputs,"b_lambdad_dose_cat|b_Vp_dose_cat|b_lambdag_dose_cat|b_tpeak_dose_cat"),
         type = paste0("Posteriors-",Group)) %>% 
  select(colnames(brms_holder_cat_prior)) %>% 
  mutate(type = factor(type,levels = c("Posteriors-Low","Posteriors-Medium","Posteriors-High")))

model_con_priors_cat_compare <- bind_rows(brms_holder_cat_prior, brms_holder_cat_post_prior)

saveRDS(model_con_priors_cat_compare, here("2_supplementary_materials/analyses/rds/model_con_priors_cat_compare.rds"))
