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

source(here("2_supplementary_materials/analyses/common_functions/brms_setting.R"))
source(here("2_supplementary_materials/analyses/common_functions/plot_theme.R"))
source(here("2_supplementary_materials/analyses/common_functions/newdata4mods.R"))
norodata <- readRDS(here("2_supplementary_materials/data/norodata.rds"))

virus_growth_rate_con19 <- readRDS(here("2_supplementary_materials/analyses/rds/virus_growth_rate_con19.rds"))

growth_rate_range <- virus_growth_rate_con19 %>% 
  pull(growth_rate) %>% 
  range()

growth_rate_mod_range <- seq(floor(growth_rate_range[[1]]), ceiling(growth_rate_range[[2]]))
# ===============================================
# Without 0.48
# ===============================================

# incubation
incubation <- norodata$incubation %>% filter(dose != 0.48) %>% droplevels()

icb_brms <- readRDS(here("2_supplementary_materials/analyses/rds/icb_brms_con19.rds"))

incubation_fit <- posterior_epred(object = icb_brms, 
                                  newdata = newdata_cnt19,
                                  re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range) %>% 
  as.data.frame() %>% 
  bind_cols(newdata_cnt19) %>% 
  mutate(outcome = "A: Incubation period (day)")

# css
css <- norodata$css %>% filter(dose != 0.48) %>% droplevels()

ttcss <- css %>% 
  group_by(ID,dose,dose_log,dose_cat) %>% 
  dplyr::summarize(ttcss = sum(css, na.rm = T)) %>% 
  ungroup()

plot_css_ts <- css %>% 
  ggplot(aes(which_day, css, shape = factor(dose), col = factor(dose))) +
  geom_point(size = 3, stroke = 0.5, show.legend = F) +
  geom_line(aes(group = ID), linetype=2, lwd=0.3, show.legend = F) +
  scale_shape_manual(values=c(2,1,0,5)) +
  labs(y = "Symptom score", x = "Days postinoculation") +
  scale_y_continuous(limits = c(-2,17)) +
  theme_bw() +
  facet_wrap(dose~.,ncol = 1,
             labeller = labeller(dose = c("4.8" = "Dose: 4.8 RT-PCR units (N = 6)",
                                          "48" = "Dose: 48 RT-PCR units (N = 7)",
                                          "4800" = "Dose: 4800 RT-PCR units (N = 6)"))) +
  plot_theme

ggsave(file = here("2_supplementary_materials/analyses/plots/plot_pkscore_ts_id19.png"), 
       plot_css_ts, dpi=dpiset, units="in", width=4, height=5)

ttcss_brms <- readRDS(here("2_supplementary_materials/analyses/rds/ttcss_brms_con19.rds"))

css_fit <- posterior_epred(object = ttcss_brms, 
                           newdata = newdata_cnt19,
                           re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range) %>% 
  as.data.frame() %>% 
  bind_cols(newdata_cnt19)  %>% 
  mutate(outcome = "C: Comprehensive symptom score")

# ---
virus_growth_rate_ttcss_brms_con19 <- readRDS(here("2_supplementary_materials/analyses/rds/virus_growth_rate_ttcss_brms_con19.rds"))

virus_growth_rate_con19 <- readRDS(here("2_supplementary_materials/analyses/rds/virus_growth_rate_con19.rds"))

virus_growth_rate_con19_ttcss <- left_join(virus_growth_rate_con19, ttcss) 

virus_growth_rate_css <- posterior_epred(object = virus_growth_rate_ttcss_brms_con19, 
                                         newdata = data.frame(growth_rate = growth_rate_mod_range),
                                         re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range) %>% 
  as.data.frame() %>% 
  mutate(growth_rate = growth_rate_mod_range) %>% 
  mutate(outcome = "Comprehensive symptom score")

virus_growth_rate_css_plot <- virus_growth_rate_css %>% 
  ggplot(aes(x = growth_rate, y = Estimate)) +
  geom_line() +
  geom_ribbon(aes(x = growth_rate, ymin=Q2.5, ymax=Q97.5),alpha=0.2, show.legend = F) +
  geom_point(data = virus_growth_rate_con19_ttcss, aes(x = growth_rate, y = ttcss), 
             shape = 1, size = 3, col = "blue") +
  labs(x = "Virus Growth Rate", y = "Comprehensive symptom score") +
  theme_minimal() +
  plot_theme

ggsave(file = here("2_supplementary_materials/analyses/plots/virus_growth_rate_css_con19_plot.png"), 
       virus_growth_rate_css_plot, dpi=dpiset, units="in", width=4, height=3) 


# =========
css_10symptoms_brms_list <- readRDS(here("2_supplementary_materials/analyses/rds/css_10symptoms_brms_list_con19.rds"))
css_10symptoms_brms_posterior_epred_list <- list()

# 1
css_10symptoms_brms_posterior_epred_list$body_temperature <- posterior_epred(object = css_10symptoms_brms_list$`Body temperature`, 
                                                                             newdata = newdata_cnt19,
                                                                             re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range, robust = TRUE) %>% 
  as.data.frame() %>% 
  bind_cols(newdata_cnt19)  %>% 
  mutate(outcome = "Body temperature")

# 2
css_10symptoms_brms_posterior_epred_list$malaise <- posterior_epred(object = css_10symptoms_brms_list$Malaise, 
                                                                    newdata = newdata_cnt19,
                                                                    re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range, robust = TRUE) %>% 
  as.data.frame() %>% 
  bind_cols(newdata_cnt19)  %>% 
  mutate(outcome = "Malaise")

# 3
css_10symptoms_brms_posterior_epred_list$muscle_aches <- posterior_epred(object = css_10symptoms_brms_list$`Muscle aches`, 
                                                                         newdata = newdata_cnt19,
                                                                         re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range, robust = TRUE) %>% 
  as.data.frame() %>% 
  bind_cols(newdata_cnt19)  %>% 
  mutate(outcome = "Muscle aches")

# 4
css_10symptoms_brms_posterior_epred_list$headache <- posterior_epred(object = css_10symptoms_brms_list$Headache, 
                                                                     newdata = newdata_cnt19,
                                                                     re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range, robust = TRUE) %>% 
  as.data.frame() %>% 
  bind_cols(newdata_cnt19)  %>% 
  mutate(outcome = "Headache")

# 5
css_10symptoms_brms_posterior_epred_list$nausea <- posterior_epred(object = css_10symptoms_brms_list$Nausea, 
                                                                   newdata = newdata_cnt19,
                                                                   re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range, robust = TRUE) %>% 
  as.data.frame() %>% 
  bind_cols(newdata_cnt19)  %>% 
  mutate(outcome = "Nausea")

# 6
css_10symptoms_brms_posterior_epred_list$chills <- posterior_epred(object = css_10symptoms_brms_list$Chills, 
                                                                   newdata = newdata_cnt19,
                                                                   re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range, robust = TRUE) %>% 
  as.data.frame() %>% 
  bind_cols(newdata_cnt19)  %>% 
  mutate(outcome = "Chills")

# 7
css_10symptoms_brms_posterior_epred_list$anorexia <- posterior_epred(object = css_10symptoms_brms_list$Anorexia, 
                                                                     newdata = newdata_cnt19,
                                                                     re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range, robust = TRUE) %>% 
  as.data.frame() %>% 
  bind_cols(newdata_cnt19)  %>% 
  mutate(outcome = "Anorexia")

# 8
css_10symptoms_brms_posterior_epred_list$cramps <- posterior_epred(object = css_10symptoms_brms_list$Cramps, 
                                                                   newdata = newdata_cnt19,
                                                                   re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range, robust = TRUE) %>% 
  as.data.frame() %>% 
  bind_cols(newdata_cnt19)  %>% 
  mutate(outcome = "Cramps")

# 9
css_10symptoms_brms_posterior_epred_list$feces_form <- posterior_epred(object = css_10symptoms_brms_list$`Feces form`, 
                                                                       newdata = newdata_cnt19,
                                                                       re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range, robust = TRUE) %>% 
  as.data.frame() %>% 
  bind_cols(newdata_cnt19)  %>% 
  mutate(outcome = "Feces form")

# 10
css_10symptoms_brms_posterior_epred_list$vomit <- posterior_epred(object = css_10symptoms_brms_list$Vomit, 
                                                                  newdata = newdata_cnt19,
                                                                  re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range, robust = TRUE) %>% 
  as.data.frame() %>% 
  bind_cols(newdata_cnt19)  %>% 
  mutate(outcome = "Vomit")

# all ten combine
css_10symptoms_brms_posterior_epred <- bind_rows(css_10symptoms_brms_posterior_epred_list)

css_10symptoms_plot <- css_10symptoms_brms_posterior_epred %>% 
  ggplot(aes(x = dose_log, y = Estimate)) +
  geom_line() +
  geom_ribbon(aes(x=dose_log, ymin=Q2.5, ymax=Q97.5),alpha=0.2, show.legend = F) +
  scale_x_continuous(breaks = c(log(c(4.8,48,4800)) - log(48)), labels = c("Log(4.8)","Log(48)","Log(4800)")) +
  labs(x = "Dose", y = "Symptom score") +
  facet_grid(outcome ~ ., switch = "y") +
  theme_bw() +
  plot_theme +
  theme(strip.text.y.left = element_text(angle = 0), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


ggsave(file = here("2_supplementary_materials/analyses/plots/css_10symptoms_plot_con19.png"), 
       css_10symptoms_plot, dpi=dpiset, units="in", width=6, height=8) 


# mvs
mvs <- norodata$mvs %>% filter(dose != 0.48) %>% droplevels() 

mvs_brms <- readRDS(here("2_supplementary_materials/analyses/rds/mvs_brms_con19.rds"))

mvs_fit <- posterior_epred(object = mvs_brms, 
                           newdata = newdata_cnt19,
                           re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range) %>% 
  as.data.frame() %>% 
  bind_cols(newdata_cnt19)  %>% 
  mutate(outcome = "B: Modified Vesikari score")


# ---
virus_growth_rate_mvs_brms_con19 <- readRDS(here("2_supplementary_materials/analyses/rds/virus_growth_rate_mvs_brms_con19.rds"))

virus_growth_rate_con19_mvs <- left_join(virus_growth_rate_con19, mvs) 

virus_growth_rate_mvs <- posterior_epred(object = virus_growth_rate_mvs_brms_con19, 
                                         newdata = data.frame(growth_rate = growth_rate_mod_range),
                                         re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range) %>% 
  as.data.frame() %>% 
  mutate(growth_rate = growth_rate_mod_range) %>% 
  mutate(outcome = "Modified Vesikari score")

virus_growth_rate_mvs_plot <- virus_growth_rate_mvs %>% 
  ggplot(aes(x = growth_rate, y = Estimate)) +
  geom_line() +
  geom_ribbon(aes(x = growth_rate, ymin=Q2.5, ymax=Q97.5),alpha=0.2, show.legend = F) +
  geom_point(data = virus_growth_rate_con19_mvs, aes(x = growth_rate, y = mvs), 
             shape = 1, size = 3, col = "blue") +
  labs(x = "Virus Growth Rate", y = "Modified Vesikari score") +
  theme_minimal() +
  plot_theme

ggsave(file = here("2_supplementary_materials/analyses/plots/virus_growth_rate_mvs_con19_plot.png"), 
       virus_growth_rate_mvs_plot, dpi=dpiset, units="in", width=4, height=3)

# =========
all_fit <- bind_rows(incubation_fit, css_fit, mvs_fit)

saveRDS(all_fit, here("2_supplementary_materials/analyses/rds/allsymptoms_fit_con19.rds"))

incubation_plot <- all_fit %>% 
  filter(outcome == "A: Incubation period (day)") %>% 
  ggplot(aes(x = dose_log, y = Estimate)) +
  geom_line() +
  geom_ribbon(aes(x=dose_log, ymin=Q2.5, ymax=Q97.5),alpha=0.2, show.legend = F) +
  geom_point(data = incubation, aes(x = dose_log, y = incubation), 
             shape = 1, size = 3, 
             position = position_jitter(w = 0.2, h = 0)) +
  scale_x_continuous(breaks = c(log(c(4.8,48,4800)) - log(48)), labels = c("Log(4.8)","Log(48)","Log(4800)")) +
  labs(x = "Dose", y = "Incubation period (day)") +
  theme_minimal() +
  plot_theme

mvs_plot <- all_fit %>% 
  filter(outcome == "B: Modified Vesikari score") %>% 
  ggplot(aes(x = dose_log, y = Estimate)) +
  geom_line() +
  geom_ribbon(aes(x=dose_log, ymin=Q2.5, ymax=Q97.5),alpha=0.2, show.legend = F) +
  geom_point(data = mvs, aes(x = dose_log, y = mvs), 
             shape = 1, size = 3, 
             position = position_jitter(w = 0.2, h = 0)) +
  scale_x_continuous(breaks = c(log(c(4.8,48,4800)) - log(48)), labels = c("Log(4.8)","Log(48)","Log(4800)")) +
  labs(x = "Dose", y = "Modified Vesikari score") +
  theme_minimal() +
  plot_theme

css_plot <- all_fit %>% 
  filter(outcome == "C: Comprehensive symptom score") %>% 
  ggplot(aes(x = dose_log, y = Estimate)) +
  geom_line() +
  geom_ribbon(aes(x=dose_log, ymin=Q2.5, ymax=Q97.5),alpha=0.2, show.legend = F) +
  geom_point(data = ttcss, aes(x = dose_log, y = ttcss), 
             shape = 1, size = 3, 
             position = position_jitter(w = 0.2, h = 0)) +
  scale_x_continuous(breaks = c(log(c(4.8,48,4800)) - log(48)), labels = c("Log(4.8)","Log(48)","Log(4800)")) +
  labs(x = "Dose", y = "Comprehensive symptom score") +
  theme_minimal() +
  plot_theme

allsymptoms_con19 <- (incubation_plot | mvs_plot | css_plot) + plot_annotation(tag_levels = "A")

ggsave(file = here("2_supplementary_materials/analyses/plots/allsymptoms_con19.png"), 
       allsymptoms_con19, dpi=dpiset, units="in", width=7, height=3) 


# ===============================================
# With 0.48
# ===============================================
# incubation
incubation <- norodata$incubation 

icb_brms <- readRDS(here("2_supplementary_materials/analyses/rds/icb_brms_con20.rds"))

incubation_fit <- posterior_epred(object = icb_brms, 
                                  newdata = newdata_cnt20,
                                  re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range) %>% 
  as.data.frame() %>% 
  bind_cols(newdata_cnt20)  %>% 
  mutate(outcome = "A: Incubation period (day)")

# css
css <- norodata$css 

ttcss <- css %>% 
  group_by(ID,dose,dose_log,dose_cat) %>% 
  dplyr::summarize(ttcss = sum(css, na.rm = T)) %>% 
  ungroup()

plot_css_ts <- css %>% 
  ggplot(aes(which_day, css, shape = factor(dose), col = factor(dose))) +
  geom_point(size = 3, stroke = 0.5, show.legend = F) +
  geom_line(aes(group = ID), linetype=2, lwd=0.3, show.legend = F) +
  scale_shape_manual(values=c(2,1,0,5)) +
  labs(y = "Symptom score", x = "Days postinoculation") +
  scale_y_continuous(limits = c(-2,17)) +
  theme_bw() +
  facet_wrap(dose~.,ncol = 1,
             labeller = labeller(dose = c("0.48" = "Dose: 0.48 RT-PCR units (N = 1)",
                                          "4.8" = "Dose: 4.8 RT-PCR units (N = 6)",
                                          "48" = "Dose: 48 RT-PCR units (N = 7)",
                                          "4800" = "Dose: 4800 RT-PCR units (N = 6)"))) +
  plot_theme

ggsave(file = here("2_supplementary_materials/analyses/plots/plot_pkscore_ts_id20.png"), 
       plot_css_ts, dpi=dpiset, units="in", width=4, height=6)

ttcss_brms <- readRDS(here("2_supplementary_materials/analyses/rds/ttcss_brms_con20.rds"))

css_fit <- posterior_epred(object = ttcss_brms, 
                           newdata = newdata_cnt20,
                           re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range) %>% 
  as.data.frame() %>% 
  bind_cols(newdata_cnt20)  %>% 
  mutate(outcome = "C: Comprehensive symptom score")

# ---
css_10symptoms_brms_list <- readRDS(here("2_supplementary_materials/analyses/rds/css_10symptoms_brms_list_con20.rds"))
css_10symptoms_brms_posterior_epred_list <- list()

# 1
css_10symptoms_brms_posterior_epred_list$body_temperature <- posterior_epred(object = css_10symptoms_brms_list$`Body temperature`, 
                                                                             newdata = newdata_cnt20,
                                                                             re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range, robust = TRUE) %>% 
  as.data.frame() %>% 
  bind_cols(newdata_cnt20)  %>% 
  mutate(outcome = "Body temperature")

# 2
css_10symptoms_brms_posterior_epred_list$malaise <- posterior_epred(object = css_10symptoms_brms_list$Malaise, 
                                                                    newdata = newdata_cnt20,
                                                                    re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range, robust = TRUE) %>% 
  as.data.frame() %>% 
  bind_cols(newdata_cnt20)  %>% 
  mutate(outcome = "Malaise")

# 3
css_10symptoms_brms_posterior_epred_list$muscle_aches <- posterior_epred(object = css_10symptoms_brms_list$`Muscle aches`, 
                                                                         newdata = newdata_cnt20,
                                                                         re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range, robust = TRUE) %>% 
  as.data.frame() %>% 
  bind_cols(newdata_cnt20)  %>% 
  mutate(outcome = "Muscle aches")

# 4
css_10symptoms_brms_posterior_epred_list$headache <- posterior_epred(object = css_10symptoms_brms_list$Headache, 
                                                                     newdata = newdata_cnt20,
                                                                     re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range, robust = TRUE) %>% 
  as.data.frame() %>% 
  bind_cols(newdata_cnt20)  %>% 
  mutate(outcome = "Headache")

# 5
css_10symptoms_brms_posterior_epred_list$nausea <- posterior_epred(object = css_10symptoms_brms_list$Nausea, 
                                                                   newdata = newdata_cnt20,
                                                                   re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range, robust = TRUE) %>% 
  as.data.frame() %>% 
  bind_cols(newdata_cnt20)  %>% 
  mutate(outcome = "Nausea")

# 6
css_10symptoms_brms_posterior_epred_list$chills <- posterior_epred(object = css_10symptoms_brms_list$Chills, 
                                                                   newdata = newdata_cnt20,
                                                                   re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range, robust = TRUE) %>% 
  as.data.frame() %>% 
  bind_cols(newdata_cnt20)  %>% 
  mutate(outcome = "Chills")

# 7
css_10symptoms_brms_posterior_epred_list$anorexia <- posterior_epred(object = css_10symptoms_brms_list$Anorexia, 
                                                                     newdata = newdata_cnt20,
                                                                     re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range, robust = TRUE) %>% 
  as.data.frame() %>% 
  bind_cols(newdata_cnt20)  %>% 
  mutate(outcome = "Anorexia")

# 8
css_10symptoms_brms_posterior_epred_list$cramps <- posterior_epred(object = css_10symptoms_brms_list$Cramps, 
                                                                   newdata = newdata_cnt20,
                                                                   re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range, robust = TRUE) %>% 
  as.data.frame() %>% 
  bind_cols(newdata_cnt20)  %>% 
  mutate(outcome = "Cramps")

# 9
css_10symptoms_brms_posterior_epred_list$feces_form <- posterior_epred(object = css_10symptoms_brms_list$`Feces form`, 
                                                                       newdata = newdata_cnt20,
                                                                       re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range, robust = TRUE) %>% 
  as.data.frame() %>% 
  bind_cols(newdata_cnt20)  %>% 
  mutate(outcome = "Feces form")

# 10
css_10symptoms_brms_posterior_epred_list$vomit <- posterior_epred(object = css_10symptoms_brms_list$Vomit, 
                                                                  newdata = newdata_cnt20,
                                                                  re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range, robust = TRUE) %>% 
  as.data.frame() %>% 
  bind_cols(newdata_cnt20)  %>% 
  mutate(outcome = "Vomit")

# all ten combine
css_10symptoms_brms_posterior_epred <- bind_rows(css_10symptoms_brms_posterior_epred_list)

css_10symptoms_plot <- css_10symptoms_brms_posterior_epred %>% 
  ggplot(aes(x = dose_log, y = Estimate)) +
  geom_line() +
  geom_ribbon(aes(x=dose_log, ymin=Q2.5, ymax=Q97.5),alpha=0.2, show.legend = F) +
  scale_x_continuous(breaks = c(log(c(0.48,4.8,48,4800)) - log(48)), labels = c("Log(0.48)","Log(4.8)","Log(48)","Log(4800)")) +
  labs(x = "Dose", y = "Symptom score") +
  facet_grid(outcome ~ ., switch = "y") +
  theme_bw() +
  plot_theme +
  theme(strip.text.y.left = element_text(angle = 0), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


ggsave(file = here("2_supplementary_materials/analyses/plots/css_10symptoms_plot_con20.png"), 
       css_10symptoms_plot, dpi=dpiset, units="in", width=6, height=8) 
# mvs
mvs <- norodata$mvs

mvs_brms <- readRDS(here("2_supplementary_materials/analyses/rds/mvs_brms_con20.rds"))

mvs_fit <- posterior_epred(object = mvs_brms, 
                           newdata = newdata_cnt20,
                           re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range) %>% 
  as.data.frame() %>% 
  bind_cols(newdata_cnt20)  %>% 
  mutate(outcome = "B: Modified Vesikari score")

# =========
all_fit <- bind_rows(incubation_fit, css_fit, mvs_fit)

incubation_plot <- all_fit %>% 
  filter(outcome == "A: Incubation period (day)") %>% 
  ggplot(aes(x = dose_log, y = Estimate)) +
  geom_line() +
  geom_ribbon(aes(x=dose_log, ymin=Q2.5, ymax=Q97.5),alpha=0.2, show.legend = F) +
  geom_point(data = incubation, aes(x = dose_log, y = incubation), 
             shape = 1, size = 3, 
             position = position_jitter(w = 0.2, h = 0)) +
  scale_x_continuous(breaks = c(log(c(0.48,4.8,48,4800)) - log(48)), labels = c("Log(0.48)","Log(4.8)","Log(48)","Log(4800)")) +
  labs(x = "Dose", y = "Incubation period (day)") +
  theme_minimal() +
  plot_theme

mvs_plot <- all_fit %>% 
  filter(outcome == "B: Modified Vesikari score") %>% 
  ggplot(aes(x = dose_log, y = Estimate)) +
  geom_line() +
  geom_ribbon(aes(x=dose_log, ymin=Q2.5, ymax=Q97.5),alpha=0.2, show.legend = F) +
  geom_point(data = mvs, aes(x = dose_log, y = mvs), 
             shape = 1, size = 3, 
             position = position_jitter(w = 0.2, h = 0)) +
  scale_x_continuous(breaks = c(log(c(0.48,4.8,48,4800)) - log(48)), labels = c("Log(0.48)","Log(4.8)","Log(48)","Log(4800)")) +
  labs(x = "Dose", y = "Modified Vesikari score") +
  theme_minimal() +
  plot_theme

css_plot <- all_fit %>% 
  filter(outcome == "C: Comprehensive symptom score") %>% 
  ggplot(aes(x = dose_log, y = Estimate)) +
  geom_line() +
  geom_ribbon(aes(x=dose_log, ymin=Q2.5, ymax=Q97.5),alpha=0.2, show.legend = F) +
  geom_point(data = ttcss, aes(x = dose_log, y = ttcss), 
             shape = 1, size = 3, 
             position = position_jitter(w = 0.2, h = 0)) +
  scale_x_continuous(breaks = c(log(c(0.48,4.8,48,4800)) - log(48)), labels = c("Log(0.48)","Log(4.8)","Log(48)","Log(4800)")) +
  labs(x = "Dose", y = "Comprehensive symptom score") +
  theme_minimal() +
  plot_theme

allsymptoms_con20 <- (incubation_plot | mvs_plot | css_plot) + plot_annotation(tag_levels = "A")

ggsave(file = here("2_supplementary_materials/analyses/plots/allsymptoms_con20.png"), 
       allsymptoms_con20, dpi=dpiset, units="in", width=10, height=3)
