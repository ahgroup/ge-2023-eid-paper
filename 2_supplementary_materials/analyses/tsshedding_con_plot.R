# Packages
# ===============================================
library(dplyr)
library(tidyr)
library(here)
library(ggplot2)
library(ggthemes)
library(patchwork)
library(ggpubr)

# ===============================================
rm(list=ls())
options(warn = -1)
# as_tibble or visualizations have variable name warnings which can be ignored

source(here("2_supplementary_materials/analyses/common_functions/plot_theme.R"))
source(here("2_supplementary_materials/analyses/common_functions/plot_holder_curves_mu_con.R"))
source(here("2_supplementary_materials/analyses/common_functions/plot_holder_curves_mu_id.R"))
source(here("2_supplementary_materials/analyses/common_functions/plots_ggpubr_id19.R"))
source(here("2_supplementary_materials/analyses/common_functions/plots_ggpubr_id20.R"))
source(here("2_supplementary_materials/analyses/common_functions/plot_holder_brms_fit_op_con.R"))

norovirus_data <- readRDS(here("2_supplementary_materials/data/norodata.rds"))

# ===============================================
# Without 0.48
# ===============================================
feces_ts19 <- norovirus_data$ftsshedding %>% 
  filter(dose != 0.48, days > 0) %>% 
  droplevels()

# prediction
holder_curves_mu19 <- readRDS(here("2_supplementary_materials/analyses/rds/holder_curves_mu_con19.rds"))

holder_curves_mu19 <- holder_curves_mu19 %>% 
  bind_rows() %>%
  mutate(Dose = factor(dose, level = c("4.8","48","4800"))) %>% 
  drop_na(Dose) 

holder_curves_mu_con19_plot <- plot_holder_curves_mu_con(data = holder_curves_mu19, lod1 = 15e3, lod2 = 40e6, idnum = "19")

ggsave(file = here("2_supplementary_materials/analyses/plots/holder_curves_mu_con19_plot.png"), 
       holder_curves_mu_con19_plot, dpi=dpiset, units="in", width=6, height=4) 

# Goodness of fit on raw data
holder_curves_mu_id_con19 <- readRDS(here("2_supplementary_materials/analyses/rds/holder_curves_mu_id_con19.rds"))

holder_curves_mu_id_con19 <- bind_rows(holder_curves_mu_id_con19)

holder_curves_mu_id_con19_plot <- plot_holder_curves_mu_id(data = feces_ts19, 
                                                           brms_fit = holder_curves_mu_id_con19, 
                                                           lod1 = 15e3, 
                                                           lod2 = 40e6)

holder_curves_mu_id_con19_plot_combine <- plots_ggpubr_id19(plot_list = holder_curves_mu_id_con19_plot, width = 1.52)

ggsave(file = here("2_supplementary_materials/analyses/plots/holder_curves_mu_id_con19_plot_combine.png"), 
       holder_curves_mu_id_con19_plot_combine, dpi=dpiset, units="in", width=9, height=5)

# model output 2nd results
holder_brms_fit_op_con19 <- readRDS(here("2_supplementary_materials/analyses/rds/holder_brms_fit_op_con19.rds"))

holder_brms_fit_op_con19_plot <- plot_holder_brms_fit_op_con(brms_fit_op = holder_brms_fit_op_con19, idnum = "19")

ggsave(file = here("2_supplementary_materials/analyses/plots/holder_brms_fit_op_con19_plot.png"), 
       holder_brms_fit_op_con19_plot, dpi = dpiset, units = "in", width = 7, height = 5)

# ===============================================
# With 0.48
# ===============================================
feces_ts20 <- norovirus_data$ftsshedding %>% 
  filter(days > 0) %>% 
  droplevels()

# fixed level predict
holder_curves_mu20 <- readRDS(here("2_supplementary_materials/analyses/rds/holder_curves_mu_con20.rds"))

holder_curves_mu20 <- holder_curves_mu20 %>% 
  bind_rows() %>%
  mutate(Dose = factor(dose, level = c("0.48","4.8","48","4800"))) %>% 
  drop_na(Dose) 

holder_curves_mu_con20_plot <- plot_holder_curves_mu_con(data = holder_curves_mu20, 
                                                         lod1 = 15e3, 
                                                         lod2 = 40e6, 
                                                         idnum = "20")

ggsave(file = here("2_supplementary_materials/analyses/plots/holder_curves_mu_con20_plot.png"), 
       holder_curves_mu_con20_plot, dpi=dpiset, units="in", width=6, height=4) 

# random level predict
holder_curves_mu_id_con20 <- readRDS(here("2_supplementary_materials/analyses/rds/holder_curves_mu_id_con20.rds"))

holder_curves_mu_id_con20 <- bind_rows(holder_curves_mu_id_con20)

holder_curves_mu_id_con20_plot <- plot_holder_curves_mu_id(data = feces_ts20, 
                                                           brms_fit = holder_curves_mu_id_con20, 
                                                           lod1 = 15e3, 
                                                           lod2 = 40e6)

holder_curves_mu_id_con20_plot_combine <- plots_ggpubr_id20(plot_list = holder_curves_mu_id_con20_plot, width = 1.5)

ggsave(file = here("2_supplementary_materials/analyses/plots/holder_curves_mu_id_con20_plot_combine.png"), 
       holder_curves_mu_id_con20_plot_combine, dpi=dpiset, units="in", width=9, height=5)

# model output 2nd results
holder_brms_fit_op_con20 <- readRDS(here("2_supplementary_materials/analyses/rds/holder_brms_fit_op_con20.rds"))
holder_brms_fit_op_con20_plot <- plot_holder_brms_fit_op_con(brms_fit_op = holder_brms_fit_op_con20, idnum = "20")

ggsave(file = here("2_supplementary_materials/analyses/plots/holder_brms_fit_op_con20_plot.png"), 
       holder_brms_fit_op_con20_plot, dpi = dpiset, units = "in", width = 7, height = 5)



# Prior v.s. Posteriors
model_con_priors19_compare <- readRDS(here("2_supplementary_materials/analyses/rds/model_con_priors19_compare.rds"))
model_con_priors20_compare <- readRDS(here("2_supplementary_materials/analyses/rds/model_con_priors20_compare.rds"))

# make plot
model_con_priors19_compare_plot <- model_con_priors19_compare %>%
  ggplot() +
  geom_density(aes(x = value, color = type), size = 1) +
  theme_minimal() +
  facet_wrap(~Parameters, scales = "free", nrow = 4)

ggsave(file = here("2_supplementary_materials/analyses/plots/model_con_priors19_compare_plot.png"), 
       model_con_priors19_compare_plot, dpi = dpiset, units = "in", width = 5, height = 5)

# make plot
model_con_priors20_compare_plot <- model_con_priors20_compare %>%
  ggplot() +
  geom_density(aes(x = value, color = type), size = 1) +
  theme_minimal() +
  facet_wrap(~Parameters, scales = "free", nrow = 4)

ggsave(file = here("2_supplementary_materials/analyses/plots/model_con_priors20_compare_plot.png"), 
       model_con_priors20_compare_plot, dpi = dpiset, units = "in", width = 5, height = 5)
