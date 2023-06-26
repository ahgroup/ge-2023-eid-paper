plot_holder_brms_fit_op_cat <- function(brms_fit_op = NA){
  
  # browser()
  
  # Onset time
  p_onset <- brms_fit_op$ts_onset %>% 
    mutate(dose = factor(dose, levels = c("Low","Medium","High"))) %>% 
    ggplot() +
    geom_linerange(aes(x = dose, y = Estimate, ymin = (Q2.5), ymax= (Q97.5)), size = 0.5) +
    geom_point(aes(x = dose, y = Estimate), show.legend = F, size = 2) +
    scale_y_continuous(limits = c(0,4)) +
    labs(x = "Dose (log scale)",
         title = "Shedding onset",
         y = "Days") +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          # axis.text.x.bottom = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.y.left = element_text(size=8, colour = "black"),
          axis.text.x = element_text(size=8, colour = "black"))
  
  # Shedding duration
  p_duration <- brms_fit_op$ts_duration %>% 
    mutate(dose = factor(dose, levels = c("Low","Medium","High"))) %>% 
    ggplot() +
    geom_linerange(aes(x = dose, y = Estimate, ymin = (Q2.5), ymax= (Q97.5)), size = 0.5) +
    geom_point(aes(x = dose, y = Estimate), show.legend = F, size = 2) +
    scale_y_continuous(limits = c(1,50)) +
    labs(x = "Dose (log scale)",
         title = "Shedding duration",
         y = "Days") +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.y.left = element_text(size=8, colour = "black"),
          axis.text.x = element_text(size=8, colour = "black"))
  
  # Peak shedding time
  p_pkv <- brms_fit_op$ts_tpkv %>% 
    mutate(dose = factor(dose, levels = c("Low","Medium","High"))) %>% 
    ggplot() +
    geom_linerange(aes(x = dose, y = Estimate, ymin = (Q2.5), ymax= (Q97.5)), size = 0.5) +
    geom_point(aes(x = dose, y = Estimate), show.legend = F, size = 2) +
    scale_y_continuous(limits = c(0,4)) +
    labs(x = "Dose (log scale)",
         title = "Time to virus peak shedding",
         y = "Days") +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          # axis.text.x.bottom = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.y.left = element_text(size=8, colour = "black"),
          axis.text.x = element_text(size=8, colour = "black"))
  
  # Area Under The Curve
  p_auc <-  brms_fit_op$ts_auc %>% 
    mutate(dose = factor(dose, levels = c("Low","Medium","High"))) %>% 
    ggplot() +
    geom_linerange(aes(x = dose, y = Estimate, ymin = (Q2.5), ymax= (Q97.5)), size = 0.5) +
    geom_point(aes(x = dose, y = Estimate), show.legend = F, size = 2) +
    scale_y_log10(expand = c(0, 0),
                  limits = c(1e7,1e13),
                  label=scientific_10) +
    labs(x = "Dose (log scale)",
         title = "Total virus load",
         y = "GEC*days/g ") +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.y.left = element_text(size=8, colour = "black"),
          axis.text.x = element_text(size=8, colour = "black"))
  
  
  plots <- (p_onset + p_pkv) / (p_duration + p_auc) + plot_annotation(tag_levels = "A")
  
  return(plots)
}