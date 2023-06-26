# Packages
# ===============================================
library(dplyr)
library(tidyr)
library(here)
library(brms)

# ===============================================
rm(list=ls())

source(here("2_supplementary_materials/analyses/common_functions/brms_setting.R"))
norodata <- readRDS(here("2_supplementary_materials/data/norodata.rds"))

# ===============================================
# Without 0.48
# ===============================================

# incubation
incubation <- norodata$incubation %>% filter(dose != 0.48) %>% droplevels()

icb_myprior = c(prior(cauchy(0,2), class = "sigma"),
                prior(normal(0,1), class = "b"),
                prior(normal(1,5), class = "Intercept"),
                prior(cauchy(0,2),  class = "sd"))

# outcome time unit is day
icb_brms <- brm(incubation ~ dose_log + (1 | ID), data = incubation, 
                family = lognormal(), prior = icb_myprior,
                backend = "cmdstanr", silent = 1, refresh = 5000,
                warmup = 25000, iter = 30000, 
                chains = mychains, cores = mycores, seed = 123456789,
                control=list(adapt_delta=0.999999, max_treedepth = 20))

saveRDS(icb_brms, here("2_supplementary_materials/analyses/rds/icb_brms_con19.rds"))

# css
css <- norodata$css %>% filter(dose != 0.48) %>% droplevels()

ttcss <- css %>% 
  group_by(ID,dose,dose_cat,dose_log) %>% 
  dplyr::summarize(ttcss = sum(css, na.rm = T)) %>% 
  ungroup()

ttcss_myprior = c(prior(gamma(0.01, 0.01), class = "shape"),
                  prior(normal(0,1), class = "b"),
                  prior(normal(1,5), class = "Intercept"),
                  prior(cauchy(0,2),  class = "sd"))

ttcss_brms <- brm(ttcss ~ dose_log + (1 | ID), data = ttcss, 
                  family = negbinomial(), prior = ttcss_myprior,
                  backend = "cmdstanr", silent = 1, refresh = 5000,
                  warmup = 25000, iter = 30000, 
                  chains = mychains, cores = mycores, seed = 123456789,
                  control=list(adapt_delta=0.999999, max_treedepth = 20))

saveRDS(ttcss_brms, here("2_supplementary_materials/analyses/rds/ttcss_brms_con19.rds"))

# ttcss ~ virus growth_rate, need to first run `ttshedding_con.R` to have virus_growth_rate_con19.rds
virus_growth_rate_con19 <- readRDS(here("2_supplementary_materials/analyses/rds/virus_growth_rate_con19.rds"))

virus_growth_rate_con19_ttcss <- left_join(virus_growth_rate_con19, ttcss) 

ttcss_growth_rate_brms19 <- brm(ttcss ~ growth_rate + (1 | ID), 
                              data = virus_growth_rate_con19_ttcss, 
                              family = negbinomial(), prior = ttcss_myprior,
                              backend = "cmdstanr", silent = 1, refresh = 5000,
                              warmup = 25000, iter = 30000, 
                              chains = mychains, cores = mycores, seed = 123456789,
                              control=list(adapt_delta=0.999999, max_treedepth = 20))

saveRDS(ttcss_growth_rate_brms19, here("2_supplementary_materials/analyses/rds/virus_growth_rate_ttcss_brms_con19.rds"))

# css css_10symptoms
css_10symptoms <- norodata$css_10symptoms %>% filter(dose != 0.48) %>% droplevels()

tt_MaxOT_score_brms <- brm(tt_MaxOT_score ~ dose_log + (1 | ID), data = css_10symptoms, 
                           family = negbinomial(), prior = ttcss_myprior,
                           backend = "cmdstanr", silent = 1, refresh = 5000,
                           warmup = 25000, iter = 30000, 
                           chains = mychains, cores = mycores, seed = 123456789,
                           control=list(adapt_delta=0.999999, max_treedepth = 20))

# css tt_MaxMA_score
tt_MaxMA_score_brms <- brm(tt_MaxMA_score ~ dose_log + (1 | ID), data = css_10symptoms, 
                           family = negbinomial(), prior = ttcss_myprior,
                           backend = "cmdstanr", silent = 1, refresh = 5000,
                           warmup = 25000, iter = 30000, 
                           chains = mychains, cores = mycores, seed = 123456789,
                           control=list(adapt_delta=0.999999, max_treedepth = 20))

# css tt_MaxMuA_score
tt_MaxMuA_score_brms <- brm(tt_MaxMuA_score ~ dose_log + (1 | ID), data = css_10symptoms, 
                            family = negbinomial(), prior = ttcss_myprior,
                            backend = "cmdstanr", silent = 1, refresh = 5000,
                            warmup = 25000, iter = 30000, 
                            chains = mychains, cores = mycores, seed = 123456789,
                            control=list(adapt_delta=0.999999, max_treedepth = 20))

# css tt_MaxHE_score
tt_MaxHE_score_brms <- brm(tt_MaxHE_score ~ dose_log + (1 | ID), data = css_10symptoms, 
                           family = negbinomial(), prior = ttcss_myprior,
                           backend = "cmdstanr", silent = 1, refresh = 5000,
                           warmup = 25000, iter = 30000, 
                           chains = mychains, cores = mycores, seed = 123456789,
                           control=list(adapt_delta=0.999999, max_treedepth = 20))

# css tt_MaxNA_score
tt_MaxNA_score_brms <- brm(tt_MaxNA_score ~ dose_log + (1 | ID), data = css_10symptoms, 
                           family = negbinomial(), prior = ttcss_myprior,
                           backend = "cmdstanr", silent = 1, refresh = 5000,
                           warmup = 25000, iter = 30000, 
                           chains = mychains, cores = mycores, seed = 123456789,
                           control=list(adapt_delta=0.999999, max_treedepth = 20))

# css tt_MaxCH_score
tt_MaxCH_score_brms <- brm(tt_MaxCH_score ~ dose_log + (1 | ID), data = css_10symptoms, 
                           family = negbinomial(), prior = ttcss_myprior,
                           backend = "cmdstanr", silent = 1, refresh = 5000,
                           warmup = 25000, iter = 30000, 
                           chains = mychains, cores = mycores, seed = 123456789,
                           control=list(adapt_delta=0.999999, max_treedepth = 20))

# css tt_MaxAN_score
tt_MaxAN_score_brms <- brm(tt_MaxAN_score ~ dose_log + (1 | ID), data = css_10symptoms, 
                           family = negbinomial(), prior = ttcss_myprior,
                           backend = "cmdstanr", silent = 1, refresh = 5000,
                           warmup = 25000, iter = 30000, 
                           chains = mychains, cores = mycores, seed = 123456789,
                           control=list(adapt_delta=0.999999, max_treedepth = 20))

# css tt_MaxCR_score
tt_MaxCR_score_brms <- brm(tt_MaxCR_score ~ dose_log + (1 | ID), data = css_10symptoms, 
                           family = negbinomial(), prior = ttcss_myprior,
                           backend = "cmdstanr", silent = 1, refresh = 5000,
                           warmup = 25000, iter = 30000, 
                           chains = mychains, cores = mycores, seed = 123456789,
                           control=list(adapt_delta=0.999999, max_treedepth = 20))

# css tt_MaxSF_score
tt_MaxSF_score_brms <- brm(tt_MaxSF_score ~ dose_log + (1 | ID), data = css_10symptoms, 
                           family = negbinomial(), prior = ttcss_myprior,
                           backend = "cmdstanr", silent = 1, refresh = 5000,
                           warmup = 25000, iter = 30000, 
                           chains = mychains, cores = mycores, seed = 123456789,
                           control=list(adapt_delta=0.999999, max_treedepth = 20))

# css tt_MaxVT_score
tt_MaxVT_score_brms <- brm(tt_MaxVT_score ~ dose_log + (1 | ID), data = css_10symptoms, 
                           family = negbinomial(), prior = ttcss_myprior,
                           backend = "cmdstanr", silent = 1, refresh = 5000,
                           warmup = 25000, iter = 30000, 
                           chains = mychains, cores = mycores, seed = 123456789,
                           control=list(adapt_delta=0.999999, max_treedepth = 20))

css_10symptoms_brms_list <- list("Body temperature" = tt_MaxOT_score_brms,
                                 Malaise = tt_MaxMA_score_brms,
                                 "Muscle aches" = tt_MaxMuA_score_brms,
                                 Headache = tt_MaxHE_score_brms,
                                 Nausea = tt_MaxNA_score_brms,
                                 Chills = tt_MaxCH_score_brms,
                                 Anorexia = tt_MaxAN_score_brms,
                                 Cramps = tt_MaxCR_score_brms,
                                 "Feces form" = tt_MaxSF_score_brms,
                                 Vomit = tt_MaxVT_score_brms)

saveRDS(css_10symptoms_brms_list, here("2_supplementary_materials/analyses/rds/css_10symptoms_brms_list_con19.rds"))

# mvs
mvs <- norodata$mvs %>% filter(dose != 0.48) %>% droplevels() 

mvs_myprior = c(prior(gamma(0.01, 0.01), class = "shape"),
                prior(normal(0,1), class = "b"),
                prior(normal(1,5), class = "Intercept"),
                prior(cauchy(0,2),  class = "sd"))

mvs_brms <- brm(mvs ~ dose_log + (1 | ID), data = mvs, 
                family = negbinomial(), prior = mvs_myprior,
                backend = "cmdstanr", silent = 1, refresh = 5000,
                warmup = 25000, iter = 30000, 
                chains = mychains, cores = mycores, seed = 123456789,
                control=list(adapt_delta=0.999999, max_treedepth = 20))

saveRDS(mvs_brms, here("2_supplementary_materials/analyses/rds/mvs_brms_con19.rds"))

# need to first run `ttshedding_con.R` to have virus_growth_rate_con19.rds
virus_growth_rate_con19_mvs <- left_join(virus_growth_rate_con19, mvs) 

mvs_growth_rate_brms19 <- brm(mvs ~ growth_rate + (1 | ID), 
                              data = virus_growth_rate_con19_mvs, 
                              family = negbinomial(), prior = ttcss_myprior,
                              backend = "cmdstanr", silent = 1, refresh = 5000,
                              warmup = 25000, iter = 30000, 
                              chains = mychains, cores = mycores, seed = 123456789,
                              control=list(adapt_delta=0.999999, max_treedepth = 20))

saveRDS(mvs_growth_rate_brms19, here("2_supplementary_materials/analyses/rds/virus_growth_rate_mvs_brms_con19.rds"))

# ===============================================
# With 0.48
# ===============================================

# incubation
incubation <- norodata$incubation

# outcome time unit is day
icb_brms <- brm(incubation ~ dose_log + (1 | ID), data = incubation, 
                family = lognormal(), prior = icb_myprior,
                backend = "cmdstanr", silent = 1, refresh = 5000,
                warmup = 25000, iter = 30000, 
                chains = mychains, cores = mycores, seed = 123456789,
                control=list(adapt_delta=0.999999, max_treedepth = 20))

saveRDS(icb_brms, here("2_supplementary_materials/analyses/rds/icb_brms_con20.rds"))

# css
css <- norodata$css

ttcss <- css %>% 
  group_by(ID,dose,dose_cat,dose_log) %>% 
  dplyr::summarize(ttcss = sum(css, na.rm = T)) %>% 
  ungroup()

ttcss_brms <- brm(ttcss ~ dose_log + (1 | ID), data = ttcss, 
                  family = negbinomial(), prior = ttcss_myprior,
                  backend = "cmdstanr", silent = 1, refresh = 5000,
                  warmup = 25000, iter = 30000, 
                  chains = mychains, cores = mycores, seed = 123456789,
                  control=list(adapt_delta=0.999999, max_treedepth = 20))

saveRDS(ttcss_brms, here("2_supplementary_materials/analyses/rds/ttcss_brms_con20.rds"))

# css css_10symptoms
css_10symptoms <- norodata$css_10symptoms

tt_MaxOT_score_brms <- brm(tt_MaxOT_score ~ dose_log + (1 | ID), data = css_10symptoms, 
                           family = negbinomial(), prior = ttcss_myprior,
                           backend = "cmdstanr", silent = 1, refresh = 5000,
                           warmup = 25000, iter = 30000, 
                           chains = mychains, cores = mycores, seed = 123456789,
                           control=list(adapt_delta=0.999999, max_treedepth = 20))

# css tt_MaxMA_score
tt_MaxMA_score_brms <- brm(tt_MaxMA_score ~ dose_log + (1 | ID), data = css_10symptoms, 
                           family = negbinomial(), prior = ttcss_myprior,
                           backend = "cmdstanr", silent = 1, refresh = 5000,
                           warmup = 25000, iter = 30000, 
                           chains = mychains, cores = mycores, seed = 123456789,
                           control=list(adapt_delta=0.999999, max_treedepth = 20))

# css tt_MaxMuA_score
tt_MaxMuA_score_brms <- brm(tt_MaxMuA_score ~ dose_log + (1 | ID), data = css_10symptoms, 
                            family = negbinomial(), prior = ttcss_myprior,
                            backend = "cmdstanr", silent = 1, refresh = 5000,
                            warmup = 25000, iter = 30000, 
                            chains = mychains, cores = mycores, seed = 123456789,
                            control=list(adapt_delta=0.999999, max_treedepth = 20))

# css tt_MaxHE_score
tt_MaxHE_score_brms <- brm(tt_MaxHE_score ~ dose_log + (1 | ID), data = css_10symptoms, 
                           family = negbinomial(), prior = ttcss_myprior,
                           backend = "cmdstanr", silent = 1, refresh = 5000,
                           warmup = 25000, iter = 30000, 
                           chains = mychains, cores = mycores, seed = 123456789,
                           control=list(adapt_delta=0.999999, max_treedepth = 20))

# css tt_MaxNA_score
tt_MaxNA_score_brms <- brm(tt_MaxNA_score ~ dose_log + (1 | ID), data = css_10symptoms, 
                           family = negbinomial(), prior = ttcss_myprior,
                           backend = "cmdstanr", silent = 1, refresh = 5000,
                           warmup = 25000, iter = 30000, 
                           chains = mychains, cores = mycores, seed = 123456789,
                           control=list(adapt_delta=0.999999, max_treedepth = 20))

# css tt_MaxCH_score
tt_MaxCH_score_brms <- brm(tt_MaxCH_score ~ dose_log + (1 | ID), data = css_10symptoms, 
                           family = negbinomial(), prior = ttcss_myprior,
                           backend = "cmdstanr", silent = 1, refresh = 5000,
                           warmup = 25000, iter = 30000, 
                           chains = mychains, cores = mycores, seed = 123456789,
                           control=list(adapt_delta=0.999999, max_treedepth = 20))

# css tt_MaxAN_score
tt_MaxAN_score_brms <- brm(tt_MaxAN_score ~ dose_log + (1 | ID), data = css_10symptoms, 
                           family = negbinomial(), prior = ttcss_myprior,
                           backend = "cmdstanr", silent = 1, refresh = 5000,
                           warmup = 25000, iter = 30000, 
                           chains = mychains, cores = mycores, seed = 123456789,
                           control=list(adapt_delta=0.999999, max_treedepth = 20))

# css tt_MaxCR_score
tt_MaxCR_score_brms <- brm(tt_MaxCR_score ~ dose_log + (1 | ID), data = css_10symptoms, 
                           family = negbinomial(), prior = ttcss_myprior,
                           backend = "cmdstanr", silent = 1, refresh = 5000,
                           warmup = 25000, iter = 30000, 
                           chains = mychains, cores = mycores, seed = 123456789,
                           control=list(adapt_delta=0.999999, max_treedepth = 20))

# css tt_MaxSF_score
tt_MaxSF_score_brms <- brm(tt_MaxSF_score ~ dose_log + (1 | ID), data = css_10symptoms, 
                           family = negbinomial(), prior = ttcss_myprior,
                           backend = "cmdstanr", silent = 1, refresh = 5000,
                           warmup = 25000, iter = 30000, 
                           chains = mychains, cores = mycores, seed = 123456789,
                           control=list(adapt_delta=0.999999, max_treedepth = 20))

# css tt_MaxVT_score
tt_MaxVT_score_brms <- brm(tt_MaxVT_score ~ dose_log + (1 | ID), data = css_10symptoms, 
                           family = negbinomial(), prior = ttcss_myprior,
                           backend = "cmdstanr", silent = 1, refresh = 5000,
                           warmup = 25000, iter = 30000, 
                           chains = mychains, cores = mycores, seed = 123456789,
                           control=list(adapt_delta=0.999999, max_treedepth = 20))

css_10symptoms_brms_list <- list("Body temperature" = tt_MaxOT_score_brms,
                                 Malaise = tt_MaxMA_score_brms,
                                 "Muscle aches" = tt_MaxMuA_score_brms,
                                 Headache = tt_MaxHE_score_brms,
                                 Nausea = tt_MaxNA_score_brms,
                                 Chills = tt_MaxCH_score_brms,
                                 Anorexia = tt_MaxAN_score_brms,
                                 Cramps = tt_MaxCR_score_brms,
                                 "Feces form" = tt_MaxSF_score_brms,
                                 Vomit = tt_MaxVT_score_brms)

saveRDS(css_10symptoms_brms_list, here("2_supplementary_materials/analyses/rds/css_10symptoms_brms_list_con20.rds"))

# mvs
mvs <- norodata$mvs

mvs_brms <- brm(mvs ~ dose_log + (1 | ID), data = mvs, 
                family = negbinomial(), prior = mvs_myprior,
                backend = "cmdstanr", silent = 1, refresh = 5000,
                warmup = 25000, iter = 30000, 
                chains = mychains, cores = mycores, seed = 123456789,
                control=list(adapt_delta=0.999999, max_treedepth = 20))

saveRDS(mvs_brms, here("2_supplementary_materials/analyses/rds/mvs_brms_con20.rds"))