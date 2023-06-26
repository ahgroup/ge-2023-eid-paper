plot_holder_curves_mu_id <- function(data = NA, brms_fit = NA, lod1 = 15e3, lod2 = 40e6){
  
  
  id_names <- unique(data$ID)
  
  holder_curves_mu_id_con19_list <- list()
  
  # plot per participant
  for (i in 1:length(id_names)) {
    # raw observations
    feces_ts_sub <- data %>% 
      filter(ID == id_names[i]) %>% 
      filter(days > 0) %>% 
      mutate(PCR2 = factor(case_when(y_virus == 15000 ~ "Below LOD1",
                                     y_virus == 40000000 ~ "Between LOD1 and LOD2",
                                     TRUE ~ "Above LOD2"),
                           levels = c("Below LOD1","Between LOD1 and LOD2","Above LOD2")))
    
    # predicted values
    holder_curves_id_sub <- brms_fit %>% 
      filter(ID == id_names[i])
    
    # make plot
    holder_curves_mu_id_con19_list[[i]] <- ggplot()+
      geom_point(data = feces_ts_sub, aes(x=days, y=y_virus, shape = PCR2), 
                 size = 2, stroke = 0.5, alpha = 1, show.legend = T) + 
      geom_line(data = holder_curves_id_sub, aes(x=days, y=exp(Estimate)), size = 0.5, col = "red") +
      geom_ribbon(data = holder_curves_id_sub, 
                  aes(x=days, ymin= exp(Q2.5), ymax= exp(Q97.5)),alpha=0.5, show.legend = F) +
      scale_y_log10(breaks=c(1, lod1, lod2, 1e12),
                    expand = c(0, 0),
                    label=scientific_10) +
      coord_cartesian(ylim = c(1,1e13)) +
      scale_shape_manual(name = "Result of PCR",values=c(0,1,2)) +
      scale_x_continuous(breaks = seq(0,max(data$days)*1.1,30), limits = c(0,max(data$days)*1.1)) +
      geom_hline(yintercept=lod1, linetype="dashed", size = 0.5, alpha = 0.5) +
      geom_hline(yintercept=lod2, linetype="dashed", size = 0.5, alpha = 0.5) +
      labs(x = "Days postinoculation",
           y = "Virus concentration (GEC/g)") +
      theme_classic() +
      guides(shape = guide_legend(override.aes = list(size = 3))) +
      plot_theme +
      theme(legend.position = "bottom") +
      annotate(geom = 'text', label = paste0("Dose: ", unique(feces_ts_sub$dose)), 
               x = 35, y = 1e12, hjust = 0, vjust = 1, fontface = 1, size = 3) +
      annotate(geom = 'text', label = c("LOD1","LOD2"), x = c(50, 50), y = c(lod1*0.4,lod2*0.35), hjust = 0, vjust = 1, size = 3)
    
  }
  
  return(holder_curves_mu_id_con19_list)
}