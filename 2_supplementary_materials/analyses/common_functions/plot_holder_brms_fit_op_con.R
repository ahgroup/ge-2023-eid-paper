plot_holder_brms_fit_op_con <- function(brms_fit_op = NA, idnum = c("19","20")){
  
  # browser()
  
  if (idnum == "19") {
    # Onset time
    p_onset <- ggplot(data = brms_fit_op$ts_onset)+
      geom_line(aes(x = dose_log, y = (Estimate)), size = 0.5) +
      geom_ribbon(aes(x = dose_log, ymin = (Q2.5), ymax= (Q97.5)),alpha = 0.2, show.legend = F) +
      scale_y_continuous(limits = c(0,4)) +
      scale_x_continuous(breaks = c(log(c(4.8,48,4800))-log(48)), labels = c("log(4.8)","log(48)","log(4800)")) +
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
    p_duration <- ggplot(data = brms_fit_op$ts_duration)+
      geom_line(aes(x = dose_log, y = (Estimate)), size = 0.5) +
      geom_ribbon(aes(x = dose_log, ymin = (Q2.5), ymax= (Q97.5)),alpha = 0.2, show.legend = F) +
      scale_y_continuous(limits = c(1,50)) +
      scale_x_continuous(breaks = c(log(c(4.8,48,4800))-log(48)), labels = c("log(4.8)","log(48)","log(4800)")) +
      labs(x = "Dose (log scale)",
           title = "Shedding duration",
           y = "Days") +
      theme_minimal() +
      theme(axis.title.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.text.y.left = element_text(size=8, colour = "black"),
            axis.text.x = element_text(size=8, colour = "black"))
    
    # Peak shedding time
    p_pkv <- ggplot(data = brms_fit_op$ts_tpkv)+
      geom_line(aes(x = dose_log, y = (Estimate)), size = 0.5) +
      geom_ribbon(aes(x = dose_log, ymin = (Q2.5), ymax= (Q97.5)),alpha = 0.2, show.legend = F) +
      scale_y_continuous(limits = c(0,4)) +
      scale_x_continuous(breaks = c(log(c(4.8,48,4800))-log(48)), labels = c("log(4.8)","log(48)","log(4800)")) +
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
    p_auc <- ggplot(data = brms_fit_op$ts_auc)+
      geom_line(aes(x = dose_log, y = (Estimate)), size = 0.5) +
      geom_ribbon(aes(x = dose_log, ymin = (Q2.5), ymax= (Q97.5)),alpha = 0.2, show.legend = F) +
      scale_y_log10(expand = c(0, 0),
                    limits = c(1e7,1e14),
                    label=scientific_10) +
      scale_x_continuous(breaks = c(log(c(4.8,48,4800))-log(48)), labels = c("log(4.8)","log(48)","log(4800)")) +
      labs(x = "Dose (log scale)",
           title = "Total virus load",
           y = "GEC*days/g ") +
      theme_minimal() +
      theme(axis.title.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.text.y.left = element_text(size=8, colour = "black"),
            axis.text.x = element_text(size=8, colour = "black"))
  }
  
  if (idnum == "20") {
    # Onset time
    p_onset <- ggplot(data = brms_fit_op$ts_onset)+
      geom_line(aes(x = dose_log, y = (Estimate)), size = 0.5) +
      geom_ribbon(aes(x = dose_log, ymin = (Q2.5), ymax= (Q97.5)),alpha = 0.2, show.legend = F) +
      scale_y_continuous(limits = c(0,4)) +
      scale_x_continuous(breaks = c(log(c(0.48,4.8,48,4800))-log(48)), labels = c("log(0.48)","log(4.8)","log(48)","log(4800)")) +
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
    p_duration <- ggplot(data = brms_fit_op$ts_duration)+
      geom_line(aes(x = dose_log, y = (Estimate)), size = 0.5) +
      geom_ribbon(aes(x = dose_log, ymin = (Q2.5), ymax= (Q97.5)),alpha = 0.2, show.legend = F) +
      scale_y_continuous(limits = c(1,50)) +
      scale_x_continuous(breaks = c(log(c(0.48,4.8,48,4800))-log(48)), labels = c("log(0.48)","log(4.8)","log(48)","log(4800)")) +
      labs(x = "Dose (log scale)",
           title = "Shedding duration",
           y = "Days") +
      theme_minimal() +
      theme(axis.title.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.text.y.left = element_text(size=8, colour = "black"),
            axis.text.x = element_text(size=8, colour = "black"))
    
    # Peak shedding time
    p_pkv <- ggplot(data = brms_fit_op$ts_tpkv)+
      geom_line(aes(x = dose_log, y = (Estimate)), size = 0.5) +
      geom_ribbon(aes(x = dose_log, ymin = (Q2.5), ymax= (Q97.5)),alpha = 0.2, show.legend = F) +
      scale_y_continuous(limits = c(0,4)) +
      scale_x_continuous(breaks = c(log(c(0.48,4.8,48,4800))-log(48)), labels = c("log(0.48)","log(4.8)","log(48)","log(4800)")) +
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
    p_auc <- ggplot(data = brms_fit_op$ts_auc)+
      geom_line(aes(x = dose_log, y = (Estimate)), size = 0.5) +
      geom_ribbon(aes(x = dose_log, ymin = (Q2.5), ymax= (Q97.5)),alpha = 0.2, show.legend = F) +
      scale_y_log10(expand = c(0, 0),
                    limits = c(1e7,1e14),
                    label=scientific_10) +
      scale_x_continuous(breaks = c(log(c(0.48,4.8,48,4800))-log(48)), labels = c("log(0.48)","log(4.8)","log(48)","log(4800)")) +
      labs(x = "Dose (log scale)",
           title = "Total virus load",
           y = "GEC*days/g ") +
      theme_minimal() +
      theme(axis.title.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.text.y.left = element_text(size=8, colour = "black"),
            axis.text.x = element_text(size=8, colour = "black"))
  }
  
  plots <- (p_onset + p_pkv) / (p_duration + p_auc) + plot_annotation(tag_levels = "A")
  
  return(plots)
}