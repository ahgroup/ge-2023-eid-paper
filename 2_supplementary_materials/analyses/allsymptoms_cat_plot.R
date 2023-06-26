# Packages
# ===============================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(here)
library(brms)
library(ggthemes)

# ===============================================
rm(list=ls())
options(warn = -1)
# as_tibble or visualizations have variable name warnings which can be ignored

source(here("2_supplementary_materials/analyses/common_functions/brms_setting.R"))
source(here("2_supplementary_materials/analyses/common_functions/plot_theme.R"))
source(here("2_supplementary_materials/analyses/common_functions/newdata4mods.R"))
norodata <- readRDS(here("2_supplementary_materials/data/norodata.rds"))

# ===============================================
# Without 0.48
# ===============================================

# incubation
incubation <- norodata$incubation %>% filter(dose != 0.48) %>% droplevels()

icb_brms <- readRDS(here("2_supplementary_materials/analyses/rds/icb_brms_cat19.rds"))

incubation_fit <- posterior_epred(object = icb_brms, 
                                  newdata = newdata_cat %>% filter(dose_cat != "Low_L"),
                                  re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range) %>% 
  as.data.frame() %>% 
  mutate(dose_cat = newdata_cat %>% filter(dose_cat != "Low_L") %>% pull(dose_cat)) %>% 
  mutate(outcome = "A: Incubation period (day)")

# css
css <- norodata$css %>% filter(dose != 0.48) %>% droplevels()

ttcss <- css %>% 
  group_by(ID,dose,dose_cat,dose_log) %>% 
  dplyr::summarize(ttcss = sum(css, na.rm = T)) %>% 
  ungroup()

ttcss_brms <- readRDS(here("2_supplementary_materials/analyses/rds/ttcss_brms_cat19.rds"))

css_fit <- posterior_epred(object = ttcss_brms, 
                           newdata = newdata_cat %>% filter(dose_cat != "Low_L"),
                           re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range) %>% 
  as.data.frame() %>% 
  mutate(dose_cat = newdata_cat %>% filter(dose_cat != "Low_L") %>% pull(dose_cat)) %>% 
  mutate(outcome = "C: Comprehensive symptom score")

# mvs
mvs <- norodata$mvs %>% filter(dose != 0.48) %>% droplevels() 

mvs_brms <- readRDS(here("2_supplementary_materials/analyses/rds/mvs_brms_cat19.rds"))

mvs_fit <- posterior_epred(object = mvs_brms, 
                           newdata = newdata_cat %>% filter(dose_cat != "Low_L"),
                           re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range) %>% 
  as.data.frame() %>% 
  mutate(dose_cat = newdata_cat %>% filter(dose_cat != "Low_L") %>% pull(dose_cat)) %>% 
  mutate(outcome = "B: Modified Vesikari score")

# =========
all_fit <- bind_rows(incubation_fit, css_fit, mvs_fit) %>% 
  mutate(dose_cat = factor(dose_cat, levels = c("Low","Medium","High")))

incubation_plot <- all_fit %>% 
  filter(outcome == "A: Incubation period (day)") %>% 
  ggplot(aes(x = dose_cat, y = Estimate)) +
  geom_pointrange(aes(ymin=Q2.5, ymax=Q97.5), size = 1) +
  geom_point(data = incubation, aes(x = dose_cat, y = incubation), 
             shape = 1, size = 3, col = "blue",
             position = position_jitter(w = 0.2, h = 0)) +
  scale_x_discrete(breaks = c("Low","Medium","High"), labels = c("Low","Medium","High")) +
  labs(x = "Dose", y = "Incubation period (day)") +
  theme_minimal() +
  plot_theme

mvs_plot <- all_fit %>% 
  filter(outcome == "B: Modified Vesikari score") %>% 
  ggplot(aes(x = dose_cat, y = Estimate)) +
  geom_pointrange(aes(ymin=Q2.5, ymax=Q97.5), size = 1) +
  geom_point(data = mvs, aes(x = dose_cat, y = mvs), 
             shape = 1, size = 3, col = "blue",
             position = position_jitter(w = 0.2, h = 0)) +
  scale_x_discrete(breaks = c("Low","Medium","High"), labels = c("Low","Medium","High")) +
  labs(x = "Dose", y = "Modified Vesikari score") +
  theme_minimal() +
  plot_theme

css_plot <- all_fit %>% 
  filter(outcome == "C: Comprehensive symptom score") %>% 
  ggplot(aes(x = dose_cat, y = Estimate)) +
  geom_pointrange(aes(ymin=Q2.5, ymax=Q97.5), size = 1) +
  geom_point(data = ttcss, aes(x = dose_cat, y = ttcss), 
             shape = 1, size = 3, col = "blue",
             position = position_jitter(w = 0.2, h = 0)) +
  scale_x_discrete(breaks = c("Low","Medium","High"), labels = c("Low","Medium","High")) +
  labs(x = "Dose", y = "Comprehensive symptom score") +
  theme_minimal() +
  plot_theme

allsymptoms_cat19 <- (incubation_plot | mvs_plot | css_plot) + plot_annotation(tag_levels = "A")

ggsave(file = here("2_supplementary_materials/analyses/plots/allsymptoms_cat19.png"), 
       allsymptoms_cat19, dpi=dpiset, units="in", width=7, height=3) 

# =========
css_10symptoms_brms_list <- readRDS(here("2_supplementary_materials/analyses/rds/css_10symptoms_brms_list_cat19.rds"))
css_10symptoms_brms_posterior_epred_list <- list()

# 1
css_10symptoms_brms_posterior_epred_list$body_temperature <- posterior_epred(object = css_10symptoms_brms_list$`Body temperature`, 
                                        newdata = newdata_cat %>% filter(dose_cat != "Low_L"),
                                        re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range, robust = TRUE) %>% 
  as.data.frame() %>% 
  mutate(dose_cat = newdata_cat %>% filter(dose_cat != "Low_L") %>% pull(dose_cat)) %>% 
  mutate(outcome = "Body temperature")

# 2
css_10symptoms_brms_posterior_epred_list$malaise <- posterior_epred(object = css_10symptoms_brms_list$Malaise, 
                                                                             newdata = newdata_cat %>% filter(dose_cat != "Low_L"),
                                                                             re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range, robust = TRUE) %>% 
  as.data.frame() %>% 
  mutate(dose_cat = newdata_cat %>% filter(dose_cat != "Low_L") %>% pull(dose_cat)) %>% 
  mutate(outcome = "Malaise")

# 3
css_10symptoms_brms_posterior_epred_list$muscle_aches <- posterior_epred(object = css_10symptoms_brms_list$`Muscle aches`, 
                                                                    newdata = newdata_cat %>% filter(dose_cat != "Low_L"),
                                                                    re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range, robust = TRUE) %>% 
  as.data.frame() %>% 
  mutate(dose_cat = newdata_cat %>% filter(dose_cat != "Low_L") %>% pull(dose_cat)) %>% 
  mutate(outcome = "Muscle aches")

# 4
css_10symptoms_brms_posterior_epred_list$headache <- posterior_epred(object = css_10symptoms_brms_list$Headache, 
                                                                         newdata = newdata_cat %>% filter(dose_cat != "Low_L"),
                                                                         re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range, robust = TRUE) %>% 
  as.data.frame() %>% 
  mutate(dose_cat = newdata_cat %>% filter(dose_cat != "Low_L") %>% pull(dose_cat)) %>% 
  mutate(outcome = "Headache")

# 5
css_10symptoms_brms_posterior_epred_list$nausea <- posterior_epred(object = css_10symptoms_brms_list$Nausea, 
                                                                     newdata = newdata_cat %>% filter(dose_cat != "Low_L"),
                                                                     re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range, robust = TRUE) %>% 
  as.data.frame() %>% 
  mutate(dose_cat = newdata_cat %>% filter(dose_cat != "Low_L") %>% pull(dose_cat)) %>% 
  mutate(outcome = "Nausea")

# 6
css_10symptoms_brms_posterior_epred_list$chills <- posterior_epred(object = css_10symptoms_brms_list$Chills, 
                                                                   newdata = newdata_cat %>% filter(dose_cat != "Low_L"),
                                                                   re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range, robust = TRUE) %>% 
  as.data.frame() %>% 
  mutate(dose_cat = newdata_cat %>% filter(dose_cat != "Low_L") %>% pull(dose_cat)) %>% 
  mutate(outcome = "Chills")

# 7
css_10symptoms_brms_posterior_epred_list$anorexia <- posterior_epred(object = css_10symptoms_brms_list$Anorexia, 
                                                                   newdata = newdata_cat %>% filter(dose_cat != "Low_L"),
                                                                   re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range, robust = TRUE) %>% 
  as.data.frame() %>% 
  mutate(dose_cat = newdata_cat %>% filter(dose_cat != "Low_L") %>% pull(dose_cat)) %>% 
  mutate(outcome = "Anorexia")

# 8
css_10symptoms_brms_posterior_epred_list$cramps <- posterior_epred(object = css_10symptoms_brms_list$Cramps, 
                                                                     newdata = newdata_cat %>% filter(dose_cat != "Low_L"),
                                                                     re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range, robust = TRUE) %>% 
  as.data.frame() %>% 
  mutate(dose_cat = newdata_cat %>% filter(dose_cat != "Low_L") %>% pull(dose_cat)) %>% 
  mutate(outcome = "Cramps")

# 9
css_10symptoms_brms_posterior_epred_list$feces_form <- posterior_epred(object = css_10symptoms_brms_list$`Feces form`, 
                                                                   newdata = newdata_cat %>% filter(dose_cat != "Low_L"),
                                                                   re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range, robust = TRUE) %>% 
  as.data.frame() %>% 
  mutate(dose_cat = newdata_cat %>% filter(dose_cat != "Low_L") %>% pull(dose_cat)) %>% 
  mutate(outcome = "Feces form")

# 10
css_10symptoms_brms_posterior_epred_list$vomit <- posterior_epred(object = css_10symptoms_brms_list$Vomit, 
                                                                       newdata = newdata_cat %>% filter(dose_cat != "Low_L"),
                                                                       re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range, robust = TRUE) %>% 
  as.data.frame() %>% 
  mutate(dose_cat = newdata_cat %>% filter(dose_cat != "Low_L") %>% pull(dose_cat)) %>% 
  mutate(outcome = "Vomit")

# all ten combine
css_10symptoms_brms_posterior_epred <- bind_rows(css_10symptoms_brms_posterior_epred_list)

css_10symptoms_plot <- css_10symptoms_brms_posterior_epred %>% 
  mutate(Dose = factor(dose_cat, levels = c("Low", "Medium", "High"))) %>% 
  ggplot(aes(x = Estimate, y = "Dose", col = factor(Dose))) +
  geom_linerange(aes(xmin = Q2.5, xmax = Q97.5), position = position_dodge(width = 1), size = 1) +
  labs(x = "Score", y = "Symptoms", col = "Dose", linetype = "Dose") +
  scale_color_manual(values = c("Low" = colorblind_pal()(8)[1],
                                "Medium" = colorblind_pal()(8)[2],
                                "High" = colorblind_pal()(8)[3])) +
  # geom_vline(xintercept = 0, linetype = 2) + 
  facet_grid(outcome ~ ., switch = "y") +
  theme_bw() +
  plot_theme +
  theme(axis.text.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        strip.text.y.left = element_text(angle = 0), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave(file = here("2_supplementary_materials/analyses/plots/css_10symptoms_plot_cat19.png"), 
       css_10symptoms_plot, dpi=dpiset, units="in", width=6, height=6) 
