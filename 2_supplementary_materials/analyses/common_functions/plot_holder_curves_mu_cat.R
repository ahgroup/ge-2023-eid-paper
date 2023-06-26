plot_holder_curves_mu_cat <- function(data = NA, lod1 = 15e3, lod2 = 40e6){
  
  mycols19 = c("Low" = colorblind_pal()(8)[1], 
               "Medium" = colorblind_pal()(8)[2], 
               "High" = colorblind_pal()(8)[3])
  
  annotation <- data.frame(x = c(11,11,15,20),
                           y = c(1000,1e11,8e4,20),
                           label = c("Onset","Peak","Duration","AUC"))
  
  
  holder_curves_mu_plot1 <- ggplot(data = data)+
    geom_ribbon(aes(x=days, ymin= exp(Q2.5), ymax= exp(Q97.5), fill = factor(Dose)),alpha = 0.3, show.legend = F) +
    geom_line(aes(x=days, y= exp(Estimate), col = factor(Dose), linetype = factor(Dose)), size = 1) +
    scale_y_log10(breaks=c(1, lod1, lod2, 1e10,1e13),
                  expand = c(0, 0),
                  label=scientific_10) +
    scale_x_continuous(breaks=c(0, 1, 4, 7, 14), expand = c(0,0)) +
    coord_cartesian(ylim = c(1,1e13),xlim = c(0,7)) +
    labs(x = "Days postinoculation",
         y = "Virus concentration (GEC/g)",
         fill = "Dose",
         col = "Dose",
         linetype = "Dose") +
    geom_hline(yintercept=lod1, linetype="dashed", size = 0.5, alpha = 0.5) +
    geom_hline(yintercept=lod2, linetype="dashed", size = 0.5, alpha = 0.5) +
    annotate(geom = 'text', label = c("LOD1","LOD2"), x = c(5, 5), y = c(lod1*0.5,lod2*0.4), hjust = 0, vjust = 1, size = 3) +
    scale_colour_manual(values = mycols19) +
    scale_fill_manual(values = mycols19) +
    theme_minimal() +
    plot_theme +
    theme(legend.position = "right",
          axis.text.y.left = element_text(size=7, colour = "black"),
          axis.text.x = element_text(size=7, colour = "black"))  
  
  
  holder_curves_mu_plot2 <- ggplot(data = data)+
    geom_ribbon(aes(x=days, ymin= exp(Q2.5), ymax= exp(Q97.5), fill = factor(Dose)),alpha = 0.3, show.legend = F) +
    geom_line(aes(x=days, y= exp(Estimate), col = factor(Dose), linetype = factor(Dose)), size = 1) +
    scale_y_log10(breaks=c(1, lod1, lod2, 1e10,1e13),
                  expand = c(0, 0),
                  label=scientific_10) +
    scale_x_continuous(expand = c(0,0)) +
    coord_cartesian(ylim = c(1,1e13)) +
    labs(x = "Days postinoculation",
         y = "Virus concentration (GEC/g)",
         fill = "Dose",
         col = "Dose",
         linetype = "Dose") +
    geom_hline(yintercept=lod1, linetype="dashed", size = 0.5, alpha = 0.5) +
    geom_hline(yintercept=lod2, linetype="dashed", size = 0.5, alpha = 0.5) +
    annotate(geom = 'text', label = c("LOD1","LOD2"), x = c(70, 70), y = c(lod1*0.5,lod2*0.4), hjust = 0, vjust = 1, size = 3) +
    scale_colour_manual(values = mycols19) +
    scale_fill_manual(values = mycols19) +
    theme_minimal() +
    plot_theme +
    theme(legend.position = "right",
          axis.text.y.left = element_text(size=7, colour = "black"),
          axis.text.x = element_text(size=7, colour = "black")) +
    geom_text(data=annotation, 
              aes(x=x, y=y, label=label),   
              size=2 , angle=0, fontface="bold" ) +
    # "Onset"
    annotate("segment", x = 11, xend = 5, y = 2000, yend = 8000, 
             size=0.4, alpha=1, arrow=arrow(ends = "last", angle = 30, length = unit(.1,"cm"))) +
    # "Peak"
    annotate("segment", x = 11, xend = 5, y = 4e10, yend = 1e10, 
             size=0.4, alpha=1, arrow=arrow(ends = "last", angle = 30, length = unit(.1,"cm"))) +
    # "Duration"
    annotate("segment", x = 3, xend = 28, y = 3e4, yend = 3e4, 
             size=0.4, alpha=1, arrow=arrow(ends = "both", angle = 30, length = unit(.1,"cm"))) +
    # "auc"
    annotate("polygon", x = c(48,0,2.3,48), y = c(0,0,1e9,0), alpha = .08, fill = "red")
  
  holder_curves_mu_plot <- holder_curves_mu_plot2 + holder_curves_mu_plot1 + 
    plot_annotation(tag_levels = 'A') + plot_layout(guides = "collect") & theme(legend.position = 'bottom') 
  
  return(holder_curves_mu_plot)
}


