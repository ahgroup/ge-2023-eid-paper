# Packages
# ===============================================
library(dplyr)
library(tidyr)
library(here)
library(brms)
library(ggplot2)
library(patchwork)

# ===============================================
rm(list=ls())
options(warn = -1)
# as_tibble or visualizations have variable name warnings which can be ignored

norodata <- readRDS(here("2_supplementary_materials/data/norodata.rds"))
source(here("2_supplementary_materials/analyses/common_functions/plot_theme.R"))
source(here("2_supplementary_materials/analyses/common_functions/newdata4mods.R"))
source(here("2_supplementary_materials/analyses/common_functions/mycoef.R"))

# ===============================================
# Without 0.48 ----
# ===============================================
## data in 96 hours
feces96 <- norodata$feces96 %>% filter(dose != 0.48) %>% droplevels()
feces96_brms <- readRDS(here("2_supplementary_materials/analyses/rds/feces96_brms_cat19.rds"))

slope <- mycoef(feces96_brms) %>% pull(op)

feces96_fit <- posterior_epred(object = feces96_brms, 
                               newdata = newdata_cat %>% filter(dose_cat != "Low_L"),
                               re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range) %>% 
  as.data.frame() %>% 
  mutate(dose_cat = newdata_cat %>% filter(dose_cat != "Low_L") %>% pull(dose_cat)) %>% 
  mutate(outcome = "feces96",
         dose_cat = factor(dose_cat, levels = c("Low","Medium","High")))

feces96_fit_plot <- ggplot(data = feces96_fit, aes(x = dose_cat, y = exp(Estimate))) +
  geom_pointrange(aes(ymin=exp(Q2.5), ymax=exp(Q97.5)), size = 0.8) +
  geom_point(data = feces96, aes(x = dose_cat, y = outcome), shape = 1, size = 3, col = "blue") +
  scale_y_log10(label=scientific_10) +
  scale_x_discrete(breaks = c("Low","Medium","High"), labels = c("Low","Medium","High")) +
  labs(y = "Virus shedding (GEC)",
       x = "Dose") +
  theme_minimal() +
  annotate("text", x= c("Low","Medium","High"), y = exp(feces96_fit$Estimate), label = paste0(slope), 
           size = 3, angle = 90, vjust = 2) +
  plot_theme +
  theme(axis.text.y.left = element_text(size=9, colour = "black"),
        axis.text.x = element_text(size=9, colour = "black")) 

## all feces data
fecesall <- norodata$fecesall %>% filter(dose != 0.48) %>% droplevels()
fecesall_brms <- readRDS(here("2_supplementary_materials/analyses/rds/fecesall_brms_cat19.rds"))

slope <- mycoef(fecesall_brms) %>% pull(op)

fecesall_fit <- posterior_epred(object = fecesall_brms, 
                                newdata = newdata_cat %>% filter(dose_cat != "Low_L"),
                                re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range) %>% 
  as.data.frame() %>% 
  mutate(dose_cat = newdata_cat %>% filter(dose_cat != "Low_L") %>% pull(dose_cat)) %>% 
  mutate(outcome = "fecesall",
         dose_cat = factor(dose_cat, levels = c("Low","Medium","High")))

fecesall_fit_plot <- ggplot(data = fecesall_fit, aes(x = dose_cat, y = exp(Estimate))) +
  geom_pointrange(aes(ymin=exp(Q2.5), ymax=exp(Q97.5)), size = 0.8) +
  geom_point(data = fecesall, aes(x = dose_cat, y = outcome), shape = 1, size = 3, col = "blue") +
  scale_y_log10(label=scientific_10) +
  scale_x_discrete(breaks = c("Low","Medium","High"), labels = c("Low","Medium","High")) +
  labs(y = "Virus shedding (GEC)",
       x = "Dose") +
  theme_minimal() +
  annotate("text", x= c("Low","Medium","High"), y = exp(fecesall_fit$Estimate), label = paste0(slope), 
           size = 3, angle = 90, vjust = 2) +
  plot_theme +
  theme(axis.text.y.left = element_text(size=9, colour = "black"),
        axis.text.x = element_text(size=9, colour = "black")) 

## all vomit data
vomitall <- norodata$vomit %>% filter(dose != 0.48) %>% droplevels()
vshedall_brms <- readRDS(here("2_supplementary_materials/analyses/rds/vshedall_brms_cat19.rds"))

slope <- mycoef(vshedall_brms) %>% pull(op)

vshedall_fit <- posterior_epred(object = vshedall_brms, 
                                newdata = newdata_cat %>% filter(dose_cat != "Low_L"),
                                re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range) %>% 
  as.data.frame() %>% 
  mutate(dose_cat = newdata_cat %>% filter(dose_cat != "Low_L") %>% pull(dose_cat)) %>%  
  mutate(outcome = "vshedall",
         dose_cat = factor(dose_cat, levels = c("Low","Medium","High")))

vshedall_fit_plot <- ggplot(data = vshedall_fit, aes(x = dose_cat, y = exp(Estimate))) +
  geom_pointrange(aes(ymin=exp(Q2.5), ymax=exp(Q97.5)), size = 0.8) +
  geom_point(data = vomitall, aes(x = dose_cat, y = outcome), shape = 1, size = 3, col = "blue") +
  scale_y_log10(label=scientific_10) +
  scale_x_discrete(breaks = c("Low","Medium","High"), labels = c("Low","Medium","High")) +
  labs(y = "Virus shedding (copies)",
       x = "Dose") +
  theme_minimal() +
  annotate("text", x= c("Low","Medium","High"), y = exp(vshedall_fit$Estimate)*40, label = paste0(slope), 
           size = 3, angle = 90, vjust = 2) +
  plot_theme +
  theme(axis.text.y.left = element_text(size=9, colour = "black"),
        axis.text.x = element_text(size=9, colour = "black")) 

# ===============================================
feces96_fit_plot <- feces96_fit_plot + 
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5))

fecesall_fit_plot <- fecesall_fit_plot + 
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5), axis.title.y.left = element_blank(), axis.text.y.left = element_blank())

vshedall_fit_plot <- vshedall_fit_plot + 
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5), axis.title.y.left = element_blank(), axis.text.y.left = element_blank())

allshedding_con19 <- ((feces96_fit_plot | fecesall_fit_plot | vshedall_fit_plot) & 
                        scale_y_log10(limits = c(1e3,1e16),label=scientific_10)) + 
  plot_annotation(tag_levels = "A")

ggsave(file = here("2_supplementary_materials/analyses/plots/allshedding_cat19.png"), 
       allshedding_con19, dpi = dpiset, units = "in", width = 7, height = 5)