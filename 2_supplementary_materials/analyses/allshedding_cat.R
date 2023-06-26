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
# Without 0.48 ----
# ===============================================
## data in 96 hours
feces96 <- norodata$feces96 %>% filter(dose != 0.48) %>% droplevels()

myprior = c(prior(cauchy(0,2),  class = "sigma"),
            prior(normal(25,5), class = "b"),
            prior(cauchy(0,2),  class = "sd"))

feces96_brms <- brm(outcome_log ~ 0 + dose_cat + (1 | ID), data = feces96, 
                    family = gaussian(), prior = myprior,
                    backend = "cmdstanr", silent = 1, refresh = 5000,
                    warmup = 55000, iter = 60000, 
                    chains = mychains, cores = mycores, seed = 123456789,
                    control=list(adapt_delta=0.999999, max_treedepth = 20))

saveRDS(feces96_brms, here("2_supplementary_materials/analyses/rds/feces96_brms_cat19.rds"))

## all feces data
fecesall <- norodata$fecesall %>% filter(dose != 0.48) %>% droplevels()

fecesall_brms <- brm(outcome_log ~ 0 + dose_cat + (1 | ID), data = fecesall, 
                     family = gaussian(), prior = myprior,
                     backend = "cmdstanr", silent = 1, refresh = 5000,
                     warmup = 25000, iter = 30000, 
                     chains = mychains, cores = mycores, seed = 123456789,
                     control=list(adapt_delta=0.999999, max_treedepth = 20))

saveRDS(fecesall_brms, here("2_supplementary_materials/analyses/rds/fecesall_brms_cat19.rds"))

## all data lod1 = 0, lod2 = lod1 value
fvshed_gmu_all_LOD12_0_15000 <- norodata$fvshed_gmu_all_LOD12_0_15000 %>% filter(dose != 0.48) %>% droplevels()

fecesall_brms_LOD12_0_15000 <- brm(outcome_log ~ 0 + dose_cat + (1 | ID), data = fvshed_gmu_all_LOD12_0_15000, 
                                   family = gaussian(), prior = myprior,
                                   backend = "cmdstanr", silent = 1, refresh = 5000,
                                   warmup = 25000, iter = 30000, 
                                   chains = mychains, cores = mycores, seed = 123456789,
                                   control=list(adapt_delta=0.9999999, max_treedepth = 20))

saveRDS(fecesall_brms_LOD12_0_15000, here("2_supplementary_materials/analyses/rds/fecesall_brms_LOD12_0_15000_cat19.rds"))

## all data lod1 = 0, lod2 = lod2 value
fvshed_gmu_all_LOD12_0_4e7 <- norodata$fvshed_gmu_all_LOD12_0_4e7 %>% filter(dose != 0.48) %>% droplevels()

fecesall_brms_LOD12_0_4e7 <- brm(outcome_log ~ 0 + dose_cat + (1 | ID), data = fvshed_gmu_all_LOD12_0_4e7, 
                                 family = gaussian(), prior = myprior,
                                 backend = "cmdstanr", silent = 1, refresh = 5000,
                                 warmup = 25000, iter = 30000, 
                                 chains = mychains, cores = mycores, seed = 123456789,
                                 control=list(adapt_delta=0.9999999, max_treedepth = 20))

saveRDS(fecesall_brms_LOD12_0_4e7, here("2_supplementary_materials/analyses/rds/fecesall_brms_LOD12_0_4e7_cat19.rds"))

## all vomit data
vomitall <- norodata$vomit %>% filter(dose != 0.48) %>% droplevels()

myprior = c(prior(cauchy(0,2),  class = "sigma"),
            prior(normal(20,5), class = "b"),
            prior(cauchy(0,2),  class = "sd"))

vshedall_brms <- brm(outcome_log ~ 0 + dose_cat + (1 | ID), data = vomitall, 
                     family = gaussian(), prior = myprior,
                     backend = "cmdstanr", silent = 1, refresh = 5000,
                     warmup = 25000, iter = 30000, 
                     chains = mychains, cores = mychains, seed = 123456789,
                     control=list(adapt_delta=0.99999999, max_treedepth = 20))

saveRDS(vshedall_brms, here("2_supplementary_materials/analyses/rds/vshedall_brms_cat19.rds"))
