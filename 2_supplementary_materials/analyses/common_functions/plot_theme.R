dpiset=600

# ggplot theme
plot_theme <- theme(axis.title.x = element_text(size = 12, face = "bold"),
                    axis.title.y = element_text(size = 12, face = "bold"),
                    axis.text = element_text(size = 10),
                    axis.text.x = element_text(colour="black"),
                    axis.text.y = element_text(colour="black"),
                    panel.grid.minor = element_blank())

# make panels
clean_theme1 <- theme(legend.position = "none", 
                      axis.text.x = element_blank(),
                      axis.title.x=element_blank(),
                      axis.title.y=element_blank())

clean_theme2 <- theme(legend.position = "none", 
                      axis.text.y.left = element_blank(),
                      axis.text.x.bottom = element_blank(),
                      axis.title.x=element_blank(),
                      axis.title.y=element_blank())

clean_theme3 <- theme(legend.position = "none", 
                      axis.text.y = element_blank(), 
                      axis.title.x=element_blank(),
                      axis.title.y=element_blank())

clean_theme4 <- theme(legend.position = "none", 
                      axis.title.x=element_blank(),
                      axis.title.y=element_blank())

scientific_10 <- function(x) {
  parse(text = sub("\\+","", gsub("e", " %*% 10^", scales::scientific_format()(x))))
}

my_ci_range = c(0.025, 0.975)

