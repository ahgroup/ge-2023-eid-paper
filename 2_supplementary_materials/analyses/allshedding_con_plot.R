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
# Without 0.48 dose ----
# ===============================================
## data in 96 hours
feces96 <- norodata$feces96 %>% filter(dose != 0.48) %>% droplevels()
feces96_brms <- readRDS(here("2_supplementary_materials/analyses/rds/feces96_brms_con19.rds"))

slope <- mycoef(feces96_brms) %>% pull(op)

feces96_fit <- posterior_epred(object = feces96_brms, 
                               newdata = newdata_cnt19,
                               re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range) %>% 
  as.data.frame() %>% 
  mutate(dose_log = newdata_cnt19$dose_log) %>% 
  mutate(outcome = "feces96")

feces96_fit_plot <- ggplot(data = feces96_fit, aes(x = dose_log, y = exp(Estimate))) +
  geom_line() +
  geom_ribbon(aes(x=dose_log, ymin=exp(Q2.5), ymax=exp(Q97.5)),alpha=0.2, show.legend = F) +
  geom_point(data = feces96, aes(x = dose_log, y = outcome), shape = 1, size = 2) +
  scale_y_log10(label=scientific_10) +
  scale_x_continuous(breaks = c(log(c(4.8,48,4800)) - log(48)), labels = c("log(4.8)","log(48)","log(4800)")) +
  labs(y = "Virus shedding (GEC)",
       x = "Dose") +
  annotate("text", x= log(80) - log(48), y = 1e15, label = paste0("\U03B2: ",slope), size = 2.5) +
  theme_minimal() +
  plot_theme +
  theme(axis.text.y.left = element_text(size=9, colour = "black"),
        axis.text.x = element_text(size=9, colour = "black")) 


## data in 96 hours with data lod1 = 0, lod2 = lod1 value
fvshed_gmu_96h_LOD12_0_15000 <- norodata$fvshed_gmu_96h_LOD12_0_15000 %>% filter(dose != 0.48) %>% droplevels()
feces96_brms_LOD12_0_15000_con19 <- readRDS(here("2_supplementary_materials/analyses/rds/feces96_brms_LOD12_0_15000_con19.rds"))

slope <- mycoef(feces96_brms_LOD12_0_15000_con19) %>% pull(op)

feces96_fit_LOD12_0_15000 <- posterior_epred(object = feces96_brms_LOD12_0_15000_con19, 
                                     newdata = newdata_cnt19,
                                     re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range) %>% 
  as.data.frame() %>% 
  mutate(dose_log = newdata_cnt19$dose_log) %>% 
  mutate(outcome = "LOD12_0_15000")

feces96_fit_LOD12_0_15000_plot <- ggplot(data = feces96_fit_LOD12_0_15000, aes(x = dose_log, y = exp(Estimate))) +
  geom_line() +
  geom_ribbon(aes(x=dose_log, ymin=exp(Q2.5), ymax=exp(Q97.5)),alpha=0.2, show.legend = F) +
  geom_point(data = fvshed_gmu_96h_LOD12_0_15000, aes(x = dose_log, y = outcome), shape = 1, size = 2) +
  scale_y_log10(label=scientific_10, limits = c(1, 1e15)) +
  scale_x_continuous(breaks = c(log(c(4.8,48,4800)) - log(48)), labels = c("log(4.8)","log(48)","log(4800)")) +
  labs(y = "Virus shedding (GEC)",
       x = "Dose") +
  annotate("text", x= log(80) - log(48), y = 1e15, label = paste0("\U03B2: ",slope), size = 3) +
  theme_minimal() +
  plot_theme +
  theme(axis.text.y.left = element_text(size=9, colour = "black"),
        axis.text.x = element_text(size=9, colour = "black")) 

## data in 96 hours with data lod1 = 0, lod2 = lod2 value
fvshed_gmu_96h_LOD12_0_4e7 <- norodata$fvshed_gmu_96h_LOD12_0_4e7 %>% filter(dose != 0.48) %>% droplevels()
feces96_brms_LOD12_0_4e7_con19 <- readRDS(here("2_supplementary_materials/analyses/rds/feces96_brms_LOD12_0_4e7_con19.rds"))

slope <- mycoef(feces96_brms_LOD12_0_4e7_con19) %>% pull(op)

feces96_fit_LOD12_0_4e7 <- posterior_epred(object = feces96_brms_LOD12_0_4e7_con19, 
                                           newdata = newdata_cnt19,
                                           re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range) %>% 
  as.data.frame() %>% 
  mutate(dose_log = newdata_cnt19$dose_log) %>% 
  mutate(outcome = "LOD12_0_4e7")

feces96_fit_LOD12_0_4e7_plot <- ggplot(data = feces96_fit_LOD12_0_4e7, aes(x = dose_log, y = exp(Estimate))) +
  geom_line() +
  geom_ribbon(aes(x=dose_log, ymin=exp(Q2.5), ymax=exp(Q97.5)),alpha=0.2, show.legend = F) +
  geom_point(data = fvshed_gmu_96h_LOD12_0_4e7, aes(x = dose_log, y = outcome), shape = 1, size = 2) +
  scale_y_log10(label=scientific_10, limits = c(1, 1e15)) +
  scale_x_continuous(breaks = c(log(c(4.8,48,4800)) - log(48)), labels = c("log(4.8)","log(48)","log(4800)")) +
  labs(y = "Virus shedding (GEC)",
       x = "Dose") +
  annotate("text", x= log(80) - log(48), y = 1e15, label = paste0("\U03B2: ",slope), size = 3) +
  theme_minimal() +
  plot_theme +
  theme(axis.text.y.left = element_text(size=9, colour = "black"),
        axis.text.x = element_text(size=9, colour = "black")) 


## all feces data
fecesall <- norodata$fecesall %>% filter(dose != 0.48) %>% droplevels()
fecesall_brms <- readRDS(here("2_supplementary_materials/analyses/rds/fecesall_brms_con19.rds"))

slope <- mycoef(fecesall_brms) %>% pull(op)

fecesall_fit <- posterior_epred(object = fecesall_brms, 
                               newdata = newdata_cnt19,
                               re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range) %>% 
  as.data.frame() %>% 
  mutate(dose_log = newdata_cnt19$dose_log) %>% 
  mutate(outcome = "fecesall")

fecesall_fit_plot <- ggplot(data = fecesall_fit, aes(x = dose_log, y = exp(Estimate))) +
  geom_line() +
  geom_ribbon(aes(x=dose_log, ymin=exp(Q2.5), ymax=exp(Q97.5)),alpha=0.2, show.legend = F) +
  geom_point(data = fecesall, aes(x = dose_log, y = outcome), shape = 1, size = 2) +
  scale_y_log10(label=scientific_10) +
  scale_x_continuous(breaks = c(log(c(4.8,48,4800)) - log(48)), labels = c("log(4.8)","log(48)","log(4800)")) +
  labs(y = "Virus shedding (GEC)",
       x = "Dose") +
  theme_minimal() +
  annotate("text", x= log(80) - log(48), y = 1e15, label = paste0("\U03B2: ",slope), size = 2.5) +
  plot_theme +
  theme(axis.text.y.left = element_text(size=9, colour = "black"),
        axis.text.x = element_text(size=9, colour = "black")) 

## all vomit data
lods <- data.frame(x = c(0.5, 0.5), 
                   y = c(1*0.25,2200*0.25), 
                   label = c("LOD1","LOD2"))

plot_vomit <- norodata$vtsshedding %>% 
  filter(dose != 0.48) %>% droplevels() %>% 
  ggplot(aes(days, as.numeric(y))) +
  geom_point(aes(shape = ID), size = 3, stroke = 1,show.legend = F, alpha = 0.8) +
  scale_shape_manual(values=c(0,1,2,5,  15,16,17,   3,4,8)) +
  scale_y_log10(limit = c(0.1,1e10),
                breaks = c(1,2200,1e10), 
                label = scientific_10) +
  scale_x_continuous(breaks = c(0,0.5,1,1.5,2,2.5,3), limits = c(0.5,2.2)) +
  labs(y = "Virus shedding concentration (copy/mL)",
       x = "Days postinoculation") +
  geom_hline(yintercept=c(1,2200), linetype="dashed", lwd=0.5) +
  theme_bw() +
  geom_text(data = lods, mapping = aes(x = x, y = y, label=label), size = 2.5) +
  facet_wrap(dose~.,ncol = 1,
             labeller = labeller(dose = c("4.8" = "Dose: 4.8 RT-PCR units (N = 4)",
                                          "48" = "Dose: 48 RT-PCR units (N = 3)",
                                          "4800" = "Dose: 4800 RT-PCR units (N = 3)"))) +
  plot_theme +
  theme(axis.text.y.left = element_text(size=8, colour = "black"),
        axis.text.x = element_text(size=8, colour = "black")) 

ggsave(file = here("2_supplementary_materials/analyses/plots/plot_vomit_id19.png"), 
       plot_vomit, dpi=dpiset, units="in", width=5, height=5)

# 
vomitall <- norodata$vomit %>% filter(dose != 0.48) %>% droplevels()

vshedall_brms <- readRDS(here("2_supplementary_materials/analyses/rds/vshedall_brms_con19.rds"))

slope <- mycoef(vshedall_brms) %>% pull(op)

vshedall_fit <- posterior_epred(object = vshedall_brms, 
                                newdata = newdata_cnt19,
                                re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range) %>% 
  as.data.frame() %>% 
  mutate(dose_log = newdata_cnt19$dose_log) %>% 
  mutate(outcome = "vshedall")

vshedall_fit_plot <- ggplot(data = vshedall_fit, aes(x = dose_log, y = exp(Estimate))) +
  geom_line() +
  geom_ribbon(aes(x=dose_log, ymin=exp(Q2.5), ymax=exp(Q97.5)),alpha=0.2, show.legend = F) +
  geom_point(data = vomitall, aes(x = dose_log, y = outcome), shape = 1, size = 2) +
  scale_y_log10(label=scientific_10) +
  scale_x_continuous(breaks = c(log(c(0.48,4.8,48,4800)) - log(48)), labels = c("log(0.48)", "log(4.8)","log(48)","log(4800)")) +
  labs(y = "Virus shedding (copies)",
       x = "Dose") +
  theme_minimal() +
  annotate("text", x= log(80) - log(48), y = 1e15, label = paste0("\U03B2: ",slope), size = 2.5) +
  plot_theme +
  theme(axis.text.y.left = element_text(size=9, colour = "black"),
        axis.text.x = element_text(size=9, colour = "black")) 


# ---
# corrleation between auc and total shed ----
# ---
feces96_2shed_brms <- readRDS(here("2_supplementary_materials/analyses/rds/feces96_2shed_cor_brms.rds"))
feces96_2shed <- readRDS(here("2_supplementary_materials/analyses/rds/feces96_2shed.rds"))

slope <- mycoef(feces96_2shed_brms) %>% pull(op)

feces96_2shed_brms_fit <- posterior_epred(object = feces96_2shed_brms, 
                                          newdata = data.frame(outcome_log = seq(min(feces96_2shed$outcome_log), 
                                                                                 max(feces96_2shed$outcome_log))),
                                          re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range) %>% 
  as.data.frame() %>% 
  mutate(outcome_log = seq(min(feces96_2shed$outcome_log), 
                           max(feces96_2shed$outcome_log))) %>% 
  mutate(outcome = "feces96_2shed")

feces96_2shed_brms_fit_plot <- ggplot(data = feces96_2shed_brms_fit, aes(x = outcome_log, y = Estimate)) +
  geom_line() +
  geom_ribbon(aes(x = outcome_log, ymin = Q2.5, ymax = Q97.5), alpha=0.2, show.legend = F) +
  geom_point(data = feces96_2shed, aes(x = outcome_log, y = AUC), shape = 1, size = 2) +
  labs(x = "Total virus shed in feces", y = "Total virus load (AUC)") +
  annotate("text", x= 18, y = 66, label = paste0("\U03B2: ",slope), size = 2.5) +
  theme_minimal() +
  plot_theme +
  theme(axis.text.y.left = element_text(size=9, colour = "black"),
        axis.text.x = element_text(size=9, colour = "black")) 



# ===============================================
feces96_fit_plot <- feces96_fit_plot + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.7))

fecesall_fit_plot <- fecesall_fit_plot + 
  theme(axis.text.x = element_text(angle = 15, vjust = 0.5, hjust=0.7), axis.title.y.left = element_blank())

vshedall_fit_plot <- vshedall_fit_plot + 
  theme(axis.text.x = element_text(angle = 15, vjust = 0.5, hjust=0.7), axis.title.y.left = element_blank(), axis.text.y.left = element_blank())

twoshedding_con19 <- ((feces96_fit_plot | vshedall_fit_plot) & 
                        scale_y_log10(limits = c(1e3,1e16),label=scientific_10)) + 
  plot_annotation(tag_levels = "A")

ggsave(file = here("2_supplementary_materials/analyses/plots/twoshedding_con19.png"), 
       twoshedding_con19, dpi = dpiset, units = "in", width = 5, height = 3)

ggsave(file = here("2_supplementary_materials/analyses/plots/fecesall_fit_plot_con19.png"), 
       fecesall_fit_plot, dpi = dpiset, units = "in", width = 3, height = 3)

allshedding_fit_con19 <- bind_rows(feces96_fit,vshedall_fit,fecesall_fit)

saveRDS(allshedding_fit_con19, file = here("2_supplementary_materials/analyses/rds/allshedding_fit_con19.rds"))

ggsave(file = here("2_supplementary_materials/analyses/plots/feces96_fit_LOD12_0_15000_plot_con19.png"), 
       feces96_fit_LOD12_0_15000_plot + theme(axis.text.x = element_text(angle = 15, vjust = 0.5, hjust=0.7)), 
       dpi = dpiset, units = "in", width = 3, height = 3)

ggsave(file = here("2_supplementary_materials/analyses/plots/feces96_fit_LOD12_0_4e7_plot_con19.png"), 
       feces96_fit_LOD12_0_4e7_plot + theme(axis.text.x = element_text(angle = 15, vjust = 0.5, hjust=0.7)), 
       dpi = dpiset, units = "in", width = 3, height = 3)

ggsave(file = here("2_supplementary_materials/analyses/plots/feces96_2shed_cor_plot.png"), 
       feces96_2shed_brms_fit_plot, 
       dpi = dpiset, units = "in", width = 3, height = 3)

# ===============================================
# With 0.48 dose ----
# ===============================================
feces96 <- norodata$feces96

## data in 96 hours
feces96_brms <- readRDS(here("2_supplementary_materials/analyses/rds/feces96_brms_con20.rds"))

slope <- mycoef(feces96_brms) %>% pull(op)

feces96_fit <- posterior_epred(object = feces96_brms, 
                               newdata = newdata_cnt20,
                               re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range) %>% 
  as.data.frame() %>% 
  mutate(dose_log = newdata_cnt20$dose_log) %>% 
  mutate(outcome = "feces96")

feces96_fit_plot <- ggplot(data = feces96_fit, aes(x = dose_log, y = exp(Estimate))) +
  geom_line() +
  geom_ribbon(aes(x=dose_log, ymin=exp(Q2.5), ymax=exp(Q97.5)),alpha=0.2, show.legend = F) +
  geom_point(data = feces96, aes(x = dose_log, y = outcome), shape = 1, size = 2) +
  scale_y_log10(label=scientific_10) +
  scale_x_continuous(breaks = c(log(c(0.48,4.8,48,4800)) - log(48)), labels = c("log(0.48)", "log(4.8)","log(48)","log(4800)")) +
  labs(y = "Virus shedding (GEC)",
       x = "Dose") +
  theme_minimal() +
  annotate("text", x= log(20) - log(48), y = 1e15, label = paste0("\U03B2: ",slope), size = 2.5) +
  plot_theme +
  theme(axis.text.y.left = element_text(size=9, colour = "black"),
        axis.text.x = element_text(size=9, colour = "black")) 

## all feces data
fecesall <- norodata$fecesall

fecesall_brms <- readRDS(here("2_supplementary_materials/analyses/rds/fecesall_brms_con20.rds"))

slope <- mycoef(fecesall_brms) %>% pull(op)

fecesall_fit <- posterior_epred(object = fecesall_brms, 
                               newdata = newdata_cnt20,
                               re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range) %>% 
  as.data.frame() %>% 
  mutate(dose_log = newdata_cnt20$dose_log) %>% 
  mutate(outcome = "fecesall")

fecesall_fit_plot <- ggplot(data = fecesall_fit, aes(x = dose_log, y = exp(Estimate))) +
  geom_line() +
  geom_ribbon(aes(x=dose_log, ymin=exp(Q2.5), ymax=exp(Q97.5)),alpha=0.2, show.legend = F) +
  geom_point(data = fecesall, aes(x = dose_log, y = outcome), shape = 1, size = 2) +
  scale_y_log10(label=scientific_10) +
  scale_x_continuous(breaks = c(log(c(0.48,4.8,48,4800)) - log(48)), labels = c("log(0.48)", "log(4.8)","log(48)","log(4800)")) +
  labs(y = "Virus shedding (GEC)",
       x = "Dose") +
  theme_minimal() +
  annotate("text", x= log(20) - log(48), y = 1e15, label = paste0("\U03B2: ",slope), size = 2.5) +
  plot_theme +
  theme(axis.text.y.left = element_text(size=9, colour = "black"),
        axis.text.x = element_text(size=9, colour = "black")) 

## all vomit data
lods <- data.frame(x = c(0.5, 0.5), 
                   y = c(1*0.25,2200*0.25), 
                   label = c("LOD1","LOD2"))

plot_vomit <- norodata$vtsshedding %>% 
  ggplot(aes(days, as.numeric(y))) +
  geom_point(aes(shape = ID), size = 3, stroke = 1,show.legend = F, alpha = 0.8) +
  scale_shape_manual(values=c(12,   0,1,2,5,  15,16,17,   3,4,8)) +
  scale_y_log10(limit = c(0.1,1e10),
                breaks = c(1,2200,1e10), 
                label = scientific_10) +
  scale_x_continuous(breaks = c(0,0.5,1,1.5,2,2.5,3), limits = c(0.5,2.2)) +
  labs(y = "Virus shedding concentration (copy/mL)",
       x = "Days postinoculation") +
  geom_hline(yintercept=c(1,2200), linetype="dashed", lwd=0.5) +
  theme_bw() +
  geom_text(data = lods, mapping = aes(x = x, y = y, label=label), size = 2.5) +
  facet_wrap(dose~.,ncol = 1,
             labeller = labeller(dose = c("0.48" = "Dose: 0.48 RT-PCR units (N = 1)",
                                          "4.8" = "Dose: 4.8 RT-PCR units (N = 4)",
                                          "48" = "Dose: 48 RT-PCR units (N = 3)",
                                          "4800" = "Dose: 4800 RT-PCR units (N = 3)"))) +
  plot_theme +
  theme(axis.text.y.left = element_text(size=8, colour = "black"),
        axis.text.x = element_text(size=8, colour = "black")) 

ggsave(file = here("2_supplementary_materials/analyses/plots/plot_vomit_id20.png"), 
       plot_vomit, dpi=dpiset, units="in", width=5, height=6)

#
vomitall <- norodata$vomit

vshedall_brms <- readRDS(here("2_supplementary_materials/analyses/rds/vshedall_brms_con20.rds"))

slope <- mycoef(vshedall_brms) %>% pull(op)

vshedall_fit <- posterior_epred(object = vshedall_brms, 
                                newdata = newdata_cnt20,
                                re_formula = NA) %>% 
  posterior_summary(probs = my_ci_range) %>% 
  as.data.frame() %>% 
  mutate(dose_log = newdata_cnt20$dose_log) %>% 
  mutate(outcome = "vshedall")

vshedall_fit_plot <- ggplot(data = vshedall_fit, aes(x = dose_log, y = exp(Estimate))) +
  geom_line() +
  geom_ribbon(aes(x=dose_log, ymin=exp(Q2.5), ymax=exp(Q97.5)),alpha=0.2, show.legend = F) +
  geom_point(data = vomitall, aes(x = dose_log, y = outcome), shape = 1, size = 2) +
  scale_y_log10(label=scientific_10) +
  scale_x_continuous(breaks = c(log(c(0.48,4.8,48,4800)) - log(48)), labels = c("log(0.48)", "log(4.8)","log(48)","log(4800)")) +
  labs(y = "Virus shedding (copies)",
       x = "Dose") +
  theme_minimal() +
  annotate("text", x= log(20) - log(48), y = 1e15, label = paste0("\U03B2: ",slope), size = 2.5) +
  plot_theme +
  theme(axis.text.y.left = element_text(size=9, colour = "black"),
        axis.text.x = element_text(size=9, colour = "black")) 

# ===============================================
feces96_fit_plot <- feces96_fit_plot + 
  theme(axis.text.x = element_text(angle = 15, vjust = 0.5, hjust=0.7))

fecesall_fit_plot <- fecesall_fit_plot + 
  theme(axis.text.x = element_text(angle = 15, vjust = 0.5, hjust=0.7), axis.title.y.left = element_blank(), axis.text.y.left = element_blank())

vshedall_fit_plot <- vshedall_fit_plot + 
  theme(axis.text.x = element_text(angle = 15, vjust = 0.5, hjust=0.7), axis.title.y.left = element_blank(), axis.text.y.left = element_blank())

allshedding_con20 <- ((feces96_fit_plot | fecesall_fit_plot | vshedall_fit_plot) & 
                        scale_y_log10(limits = c(1e3,1e16),label=scientific_10)) + 
  plot_annotation(tag_levels = "A")

ggsave(file = here("2_supplementary_materials/analyses/plots/allshedding_con20.png"), 
       allshedding_con20, dpi = dpiset, units = "in", width = 7, height = 3)