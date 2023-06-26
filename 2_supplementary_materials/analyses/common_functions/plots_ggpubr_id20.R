plots_ggpubr_id20 <- function(plot_list = NA, width = NA){
  
  # make panels
  width1 = width
  
  source(here("2_supplementary_materials/analyses/common_functions/plot_theme.R"))
  
  dose0.48 <- ggarrange(plot_list[[20]] + clean_theme1, 
                        NA,
                        NA,
                        NA,
                        NA,
                        NA,
                        NA,
                        nrow = 1, widths = c(width1,rep(1,6)))
  
  dose4.8 <- ggarrange(plot_list[[19]] + clean_theme1, 
                       plot_list[[18]] + clean_theme2, 
                       plot_list[[17]] + clean_theme2, 
                       plot_list[[16]] + clean_theme2, 
                       plot_list[[15]] + clean_theme2, 
                       plot_list[[14]] + clean_theme2,
                       NA,
                       nrow = 1, widths = c(width1,rep(1,6)))
  
  dose48 <- ggarrange(plot_list[[13]] + clean_theme1, 
                      plot_list[[12]] + clean_theme2, 
                      plot_list[[11]] + clean_theme2, 
                      plot_list[[10]] + clean_theme2, 
                      plot_list[[9]] + clean_theme2, 
                      plot_list[[8]] + clean_theme2,
                      plot_list[[7]] + clean_theme2,
                      nrow = 1, widths = c(width1,rep(1,6)))
  
  dose4800 <- ggarrange(plot_list[[6]] + clean_theme4, 
                        plot_list[[5]] + clean_theme3, 
                        plot_list[[4]] + clean_theme3, 
                        plot_list[[3]] + clean_theme3, 
                        plot_list[[2]] + clean_theme3, 
                        plot_list[[1]] + clean_theme3,
                        NA,
                        nrow = 1, widths = c(width1,rep(1,6)), common.legend = T, legend = "bottom")
  
  plot_id <- ggarrange(dose0.48, dose4.8, dose48, dose4800, nrow = 4, heights = c(1,1,1,1.55))
  
  plot_id <- annotate_figure(plot_id,
                             bottom = text_grob("Days postinoculation", size = 16, face = "bold"), 
                             left = text_grob("Virus concentration (GEC/g)", size = 16, face = "bold", rot = 90))
  
  return(plot_id)
}