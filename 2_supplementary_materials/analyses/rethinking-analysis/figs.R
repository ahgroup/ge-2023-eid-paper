# Replicating Yang's figures ###################################################
# Project: Analyzing the association between inoculum dose and Norovirus
#   infection outcomes
# Author: Zane Billings
# Since: 2022-03-30
# Desc: For this paper, Yang fit a bunch of bayesian regression models using
#   the brms package. I fitted the models and extracted predictions in the
#   models.R script. In this script, I will recreate the figures.

# Setting up dependencies ######################################################

# Loading entire package namespace for patchwork to combine plots
# It seems like there is not an easy way to only import the plot arithmetic
# operators into the namespace.
library(patchwork)

box::use(
  arrow,
  dplyr,
  ggplot2[...],
  ggtext,
  here,
  colorblindr
)

# colorblindr has specific installation instructions, see the GitHub page
# for more info https://github.com/clauswilke/colorblindr. Note you may also
# need to run R/RStudio as admin on Windows.

# Function for formatting ggplot axes with 10^ format
exp_labs <- function(x) { return( paste0( "10^", log10(x) ) ) }

# Main manuscript Figure 1 - Shedding models ###################################

## Panel A: fecal shedding first 96 hours ######################################

# Load the predictions data file
tfs_96_con_fit <- arrow::read_parquet(
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "tfs_96_con_fit.parquet")
)

# Load the raw data file
tfs_96_con_raw <- readRDS(
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "data", "tfs_96_con_data.Rds")
)

# Make the figure
tfs_96_con_plot <-
  ggplot2::ggplot(
    tfs_96_con_fit,
    ggplot2::aes(x = dose_log, y = exp(mu), ymin = exp(lwr), ymax = exp(upr))
  ) +
  ggplot2::geom_line() +
  ggplot2::geom_ribbon(alpha = 0.2, show.legend = F) +
  ggplot2::geom_point(
    data = tfs_96_con_raw,
    ggplot2::aes(x = dose_log, y = exp(outcome_log)),
    shape = 1, size = 2, inherit.aes = F
  ) +
  ggplot2::scale_y_continuous(
    trans = "log10",
    limits = c(10 ^ 3, 10 ^ 16),
    breaks = 10 ^ seq(to = 16, from = 4, by = 4),
    labels = exp_labs
  ) +
  ggplot2::scale_x_continuous(
    breaks = c(log(c(4.8,48,4800)) - log(48)),
    labels = c("Log(4.8)","Log(48)","Log(4800)")
  ) +
  ggplot2::labs(
    y = "Total virus shed (GEC)\n",
    x = "\nDose"
  ) +
  ggplot2::theme_minimal(base_size = 20)

## Panel B: fecal shedding full time ###########################################

# Load the predictions data file
tfs_all_con_fit <- arrow::read_parquet(
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "tfs_all_con_fit.parquet")
)

# Load the raw data file
tfs_all_con_raw <- readRDS(
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "data", "tfs_all_con_data.Rds")
)

# Make the figure
tfs_all_con_plot <-
  ggplot2::ggplot(
    tfs_all_con_fit,
    ggplot2::aes(x = dose_log, y = exp(mu), ymin = exp(lwr), ymax = exp(upr))
  ) +
  ggplot2::geom_line() +
  ggplot2::geom_ribbon(alpha = 0.2, show.legend = F) +
  ggplot2::geom_point(
    data = tfs_all_con_raw,
    ggplot2::aes(x = dose_log, y = exp(outcome_log)),
    shape = 1, size = 2, inherit.aes = F
  ) +
  ggplot2::scale_y_continuous(
    trans = "log10",
    limits = c(10 ^ 3, 10 ^ 16),
    breaks = 10 ^ seq(to = 16, from = 4, by = 4),
    labels = exp_labs
  ) +
  ggplot2::scale_x_continuous(
    breaks = c(log(c(4.8,48,4800)) - log(48)),
    labels = c("Log(4.8)","Log(48)","Log(4800)")
  ) +
  ggplot2::labs(
    y = "Total virus shed (GEC)\n",
    x = "\nDose"
  ) +
  ggplot2::theme_minimal(base_size = 20)

## Panel C: vomit shedding #####################################################

# Load the predictions data file
vvshed_con_fit <- arrow::read_parquet(
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "vvshed_con_fit.parquet")
)

# Load the raw data file
vvshed_con_raw <- readRDS(
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "data", "vvshed_con_data.Rds")
)

# Make the figure
vvshed_con_plot <-
  ggplot2::ggplot(
    vvshed_con_fit,
    ggplot2::aes(x = dose_log, y = exp(mu), ymin = exp(lwr), ymax = exp(upr))
  ) +
  ggplot2::geom_line() +
  ggplot2::geom_ribbon(alpha = 0.2, show.legend = F) +
  ggplot2::geom_point(
    data = vvshed_con_raw,
    ggplot2::aes(x = dose_log, y = exp(outcome_log)),
    shape = 1, size = 2, inherit.aes = F
  ) +
  ggplot2::scale_y_continuous(
    trans = "log10",
    limits = c(10 ^ 3, 10 ^ 16),
    breaks = 10 ^ seq(to = 16, from = 4, by = 4),
    labels = exp_labs
  ) +
  ggplot2::scale_x_continuous(
    breaks = c(log(c(4.8,48,4800)) - log(48)),
    labels = c("Log(4.8)","Log(48)","Log(4800)")
  ) +
  ggplot2::labs(
    y = "Total virus shed (GEC)\n",
    x = "\nDose"
  ) +
  ggplot2::theme_minimal(base_size = 20)

## Recreating manuscript figure ################################################

fig1 <- (tfs_96_con_plot | tfs_all_con_plot | vvshed_con_plot)

ggplot2::ggsave(
  filename = here::here("2_supplementary_materials", "analyses", "rethinking-results", "figs", "fig-1.png"),
  plot = fig1,
  width = 13, height = 8
)

# Main manuscript Figure 2 - Time series population curves #####################

# Load the predictions data file
ts_con_pop_fit <- arrow::read_parquet(
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "timeseries_con_pop_fit.parquet")
)

## Full data ####
ts_con_pop_plot_all <-
  ts_con_pop_fit |>
  dplyr::mutate(across(c(mu, lwr, upr), exp)) |>
  #dplyr::filter(days <= 7) |>
  ggplot2::ggplot(
    mapping = ggplot2::aes(
      x = days,
      y = mu,
      ymin = lwr,
      ymax = upr,
      color = factor(dose_log)
    )
  ) +
  ggplot2::geom_ribbon(
    ggplot2::aes(fill = factor(dose_log), group = factor(dose_log)),
    color = NA,
    alpha = 0.15
  ) +
  ggplot2::geom_line() +
  ggplot2::scale_y_continuous(
    trans = "log10",
    labels = exp_labs
  ) +
  ggplot2::scale_x_continuous(
    #breaks = seq(0, 7, 1)
  ) +
  ggplot2::coord_cartesian(
    ylim = c(10 ^ 0, 10 ^ 14),
    expand = FALSE
  ) +
  ggplot2::labs(
    y = "Viral shedding (GEC)\n",
    x = "\nTime (day)"
  ) +
  ggplot2::theme_bw(base_size = 20) +
  colorblindr::scale_fill_OkabeIto() +
  colorblindr::scale_color_OkabeIto()

## Zoomed in to first 7 days ####
ts_con_pop_plot_zoom <-
  ts_con_pop_fit |>
  dplyr::mutate(across(c(mu, lwr, upr), exp)) |>
  dplyr::filter(days <= 7) |>
  ggplot2::ggplot(
    mapping = ggplot2::aes(
      x = days,
      y = mu,
      ymin = lwr,
      ymax = upr,
      color = factor(dose_log)
    )
  ) +
  ggplot2::geom_ribbon(
    ggplot2::aes(fill = factor(dose_log), group = factor(dose_log)),
    color = NA,
    alpha = 0.15
  ) +
  ggplot2::geom_line() +
  ggplot2::scale_y_continuous(
    trans = "log10",
    labels = exp_labs
  ) +
  ggplot2::scale_x_continuous(
    #breaks = seq(0, 7, 1)
  ) +
  ggplot2::coord_cartesian(
    ylim = c(10 ^ 0, 10 ^ 14),
    expand = FALSE
  ) +
  ggplot2::labs(
    y = "Viral shedding (GEC)\n",
    x = "\nTime (day)"
  ) +
  ggplot2::theme_bw(base_size = 20) +
  colorblindr::scale_fill_OkabeIto() +
  colorblindr::scale_color_OkabeIto()

## Compose figure ####

ts_con_pop_plot <-
  ts_con_pop_plot_all +
  ts_con_pop_plot_zoom +
  patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = 'bottom')

ggplot2::ggsave(
  filename = here::here("2_supplementary_materials", "analyses", "rethinking-results", "figs", "fig-2.png"),
  plot = ts_con_pop_plot,
  width = 13, height = 8
)

# Main manuscript Figure 3 - Time series outcome plots #########################

outcomes_dat <- arrow::read_parquet(
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "timeseries_con_stats.parquet")
)

# Have to make four independent panel plots since the axes need to be scaled
# differently

## Panel A: Shedding onset #####################################################
onset <- outcomes_dat |>
  dplyr::filter(stat == "onset") |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = dose_log,
      y = mean,
      ymin = lwr,
      ymax = upr
    )
  ) +
  geom_ribbon(alpha = 0.5) +
  geom_line() +
  ggplot2::scale_x_continuous(
    breaks = c(log(c(4.8,48,4800)) - log(48)),
    labels = c("log(4.8)","log(48)","log(4800)")
  ) +
  scale_y_continuous(
    limits = c(0, 3),
    breaks = seq(0, 3, by = 0.5)
  ) +
  labs(
    x = NULL,
    y = "days"
  ) +
  ggplot2::theme_minimal(base_size = 20)

## Panel B: Time to peak #######################################################
peak <- outcomes_dat |>
  dplyr::filter(stat == "peak") |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = dose_log,
      y = mean,
      ymin = lwr,
      ymax = upr
    )
  ) +
  geom_ribbon(alpha = 0.5) +
  geom_line() +
  ggplot2::scale_x_continuous(
    breaks = c(log(c(4.8,48,4800)) - log(48)),
    labels = c("Log(4.8)","Log(48)","Log(4800)")
  ) +
  scale_y_continuous(
    limits = c(0, 3),
    breaks = seq(0, 3, by = 0.5)
  ) +
  labs(
    x = NULL,
    y = "days"
  ) +
  ggplot2::theme_minimal(base_size = 20)

## Panel C: Shedding duration ##################################################
dur <- outcomes_dat |>
  dplyr::filter(stat == "duration") |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = dose_log,
      y = mean,
      ymin = lwr,
      ymax = upr
    )
  ) +
  geom_ribbon(alpha = 0.5) +
  geom_line() +
  ggplot2::scale_x_continuous(
    breaks = c(log(c(4.8,48,4800)) - log(48)),
    labels = c("Log(4.8)","Log(48)","Log(4800)")
  ) +
  scale_y_continuous(
    limits = c(0, 35),
    breaks = seq(0, 35, by = 5)
  ) +
  labs(
    x = NULL,
    y = "days"
  ) +
  ggplot2::theme_minimal(base_size = 20)

## Panel D: total virus load ###################################################
auc <- outcomes_dat |>
  dplyr::filter(stat == "auc") |>
  dplyr::mutate(dplyr::across(c(mean, lwr, upr), exp)) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = dose_log,
      y = mean,
      ymin = lwr,
      ymax = upr
    )
  ) +
  geom_ribbon(alpha = 0.5) +
  geom_line() +
  ggplot2::scale_x_continuous(
    breaks = c(log(c(4.8,48,4800)) - log(48)),
    labels = c("Log(4.8)","Log(48)","Log(4800)")
  ) +
  scale_y_continuous(
    trans = "log10",
    limits = 10 ^ c(0, 12),
    breaks = 10 ^ seq(0, 12, by = 4),
    labels = exp_labs
  ) +
  labs(
    x = NULL,
    y = "days"
  ) +
  ggplot2::theme_minimal(base_size = 20)

## Figure recreation ###########################################################

fig3 <- onset + peak + dur + auc

ggplot2::ggsave(
  filename = here::here("2_supplementary_materials", "analyses", "rethinking-results", "figs", "fig-3.png"),
  plot = fig3,
  width = 13, height = 8
)

# Main manuscript Figure 4 - Symptom outcomes ##################################

## Panel A: Incubation time ####################################################

# Load the predictions data file
incub_con_fit <- arrow::read_parquet(
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "incubation_con_fit.parquet")
)

# Load the raw data file
incub_con_raw <- readRDS(
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "data", "incubation_con_data.Rds")
)

# Make the figure
incubation_con_plot <-
  incub_con_fit |>
  dplyr::mutate(dplyr::across(!dose_log, exp)) |>
  ggplot2::ggplot(
    mapping = ggplot2::aes(x = dose_log, y = mu, ymin = lwr, ymax = upr)
  ) +
  ggplot2::geom_ribbon(alpha = 0.2) +
  ggplot2::geom_line() +
  ggplot2::geom_point(
    data = incub_con_raw,
    ggplot2::aes(x = dose_log, y = incubation),
    inherit.aes = FALSE,
    shape = 21
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 10),
    breaks = seq(0, 3, 1)
  ) +
  ggplot2::scale_x_continuous(
    breaks = c(log(c(4.8,48,4800)) - log(48)),
    labels = c("Log(4.8)","Log(48)","Log(4800)")
  ) +
  ggplot2::labs(
    x = "\nDose",
    y = "Incubation time (days)\n"
  ) +
  ggplot2::theme_minimal(base_size = 20)

## Panel B: Vesikari score ######################################################

# Load the predictions data
mvs_con_fit <- arrow::read_parquet(
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "mvs_con_fit.parquet")
)

# Load the raw data file
mvs_con_raw <- readRDS(
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "data", "mvs_con_data.Rds")
)

mvs_con_plot <-
  mvs_con_fit |>
  dplyr::mutate(dplyr::across(!dose_log, exp)) |>
  ggplot2::ggplot(
    mapping = ggplot2::aes(x = dose_log, y = mu, ymin = lwr, ymax = upr)
  ) +
  ggplot2::geom_ribbon(alpha = 0.2) +
  ggplot2::geom_line() +
  ggplot2::geom_point(
    data = mvs_con_raw,
    mapping = ggplot2::aes(x = dose_log, y = mvs),
    inherit.aes = FALSE,
    shape = 21
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 10),
    breaks = c(0, 5, 10),
    minor_breaks = seq(0, 10, 1)
  ) +
  ggplot2::scale_x_continuous(
    breaks = log(c(4.8,48,4800)) - log(48),
    labels = c("log(4.8)", "log(48)", "log(4800)")
  ) +
  ggplot2::labs(
    x = "\nDose",
    y = "Modified Vesikari score\n"
  ) +
  ggplot2::theme_minimal(base_size = 20)

## Panel C: comprehensive symptom score ########################################

# Load the predictions data
css_con_fit <- arrow::read_parquet(
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "css_con_fit.parquet")
)

# Load the raw data file
css_con_raw <- readRDS(
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "data", "css_con_data.Rds")
)

css_con_plot <-
  css_con_fit |>
  dplyr::mutate(dplyr::across(!dose_log, exp)) |>
  ggplot2::ggplot(
    mapping = ggplot2::aes(x = dose_log, y = mu, ymin = lwr, ymax = upr)
  ) +
  ggplot2::geom_ribbon(alpha = 0.2) +
  ggplot2::geom_line() +
  ggplot2::geom_point(
    data = css_con_raw,
    mapping = ggplot2::aes(x = dose_log, y = outcome),
    inherit.aes = FALSE,
    shape = 21
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 35),
    breaks = seq(0, 35, 5)
  ) +
  ggplot2::scale_x_continuous(
    breaks = log(c(4.8,48,4800)) - log(48),
    labels = c("log(4.8)", "log(48)", "log(4800)")
  ) +
  ggplot2::labs(
    x = "\nDose",
    y = "Comprehensive symptom score\n"
  ) +
  ggplot2::theme_minimal(base_size = 20)

## Recreate the figure #########################################################
fig4 <- (incubation_con_plot | mvs_con_plot | css_con_plot)
ggplot2::ggsave(filename = here::here("2_supplementary_materials", "analyses", "rethinking-results", "figs", "fig-4.png"), plot = fig4,
                width = 13, height = 8)

# Supplement Figure 3: Shedding Models - Categorical Dose ######################

## Panel A: fecal shedding first 96 hours ######################################

# Load the predictions data file
tfs_96_cat_fit <- arrow::read_parquet(
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "tfs_96_cat_fit.parquet")
) |>
  dplyr::mutate(dose_cat = factor(dose_cat,
                                  levels = c(1, 2, 3),
                                  labels = c("Low", "Medium", "High")),
                across(c(mu, lwr, upr), exp))

# Load the raw data file
tfs_96_cat_raw <- readRDS(
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "data", "tfs_96_cat_data.Rds")
)

# Make the figure
tfs_96_cat_plot <-
  ggplot2::ggplot(
    tfs_96_cat_fit,
    ggplot2::aes(
      x = dose_cat,
      y = mu, ymin = lwr, ymax = upr
    )
  ) +
  ggplot2::geom_pointrange(size = 1) +
  ggplot2::geom_point(
    data = tfs_96_cat_raw,
    mapping = ggplot2::aes(
      x = dose_cat,
      y = exp(outcome_log)
    ),
    shape = 1, size = 2, stroke = 1, color = "red",
    inherit.aes = F
  ) +
  ggplot2::scale_x_discrete(
    labels = c("low", "medium", "high")
  ) +
  ggplot2::scale_y_continuous(
    trans = "log10",
    limits = c(10 ^ 3, 10 ^ 15),
    breaks = 10 ^ c(5, 10, 15),
    labels = exp_labs
  ) +
  ggplot2::labs(
    y = "Virus shedding (GEC)\n",
    x = "\nDose"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.title.x     = ggplot2::element_text(size = 12, face = "bold"),
    axis.title.y     = ggplot2::element_text(size = 12, face = "bold"),
    axis.text        = ggplot2::element_text(size = 10),
    axis.text.x      = ggplot2::element_text(colour="black"),
    axis.text.y      = ggplot2::element_text(colour="black"),
    panel.grid.minor = ggplot2::element_blank()
  )

## Panel B: fecal shedding full time ###########################################

# Load the predictions data file
tfs_all_cat_fit <- arrow::read_parquet(
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "tfs_all_cat_fit.parquet")
) |>
  dplyr::mutate(dose_cat = factor(dose_cat,
                                  levels = c(1, 2, 3),
                                  labels = c("Low", "Medium", "High")),
                across(c(mu, lwr, upr), exp))

# Load the raw data file
tfs_all_cat_raw <- readRDS(
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "data", "tfs_all_cat_data.Rds")
)

# Make the figure
tfs_all_cat_plot <-
  ggplot2::ggplot(
    tfs_all_cat_fit,
    ggplot2::aes(
      x = factor(dose_cat),
      y = mu, ymin = lwr, ymax = upr
    )
  ) +
  ggplot2::geom_pointrange(size = 1) +
  ggplot2::geom_point(
    data = tfs_all_cat_raw,
    mapping = ggplot2::aes(
      x = dose_cat,
      y = exp(outcome_log)
    ),
    shape = 1, size = 2, stroke = 1, color = "red",
    inherit.aes = F
  ) +
  ggplot2::scale_x_discrete(
    labels = c("low", "medium", "high")
  ) +
  ggplot2::scale_y_continuous(
    trans = "log10",
    limits = c(10 ^ 3, 10 ^ 15),
    breaks = 10 ^ c(5, 10, 15),
    labels = exp_labs
  ) +
  ggplot2::labs(
    y = "Virus shedding (GEC)\n",
    x = "\nDose"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.title.x     = ggplot2::element_text(size = 12, face = "bold"),
    axis.title.y     = ggplot2::element_text(size = 12, face = "bold"),
    axis.text        = ggplot2::element_text(size = 10),
    axis.text.x      = ggplot2::element_text(colour="black"),
    axis.text.y      = ggplot2::element_text(colour="black"),
    panel.grid.minor = ggplot2::element_blank()
  )

## Panel C: vomit shedding #####################################################

# Load the predictions data file
vvshed_cat_fit <- arrow::read_parquet(
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "vvshed_cat_fit.parquet")
)  |>
  dplyr::mutate(dose_cat = factor(dose_cat,
                                  levels = c(1, 2, 3),
                                  labels = c("Low", "Medium", "High")),
                across(c(mu, lwr, upr), exp))

# Load the raw data file
vvshed_cat_raw <- readRDS(
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "data", "vvshed_cat_data.Rds")
)

# Make the figure
vvshed_cat_plot <-
  ggplot2::ggplot(
    vvshed_cat_fit,
    ggplot2::aes(
      x = factor(dose_cat),
      y = mu, ymin = lwr, ymax = upr
    )
  ) +
  ggplot2::geom_pointrange(size = 1) +
  ggplot2::geom_point(
    data = vvshed_cat_raw,
    mapping = ggplot2::aes(
      x = dose_cat,
      y = exp(outcome_log)
    ),
    shape = 1, size = 2, stroke = 1, color = "red",
    inherit.aes = F
  ) +
  ggplot2::scale_x_discrete(
    labels = c("low", "medium", "high")
  ) +
  ggplot2::scale_y_continuous(
    trans = "log10",
    limits = c(10 ^ 3, 10 ^ 15),
    breaks = 10 ^ c(5, 10, 15),
    labels = exp_labs
  ) +
  ggplot2::labs(
    y = "Virus shedding (GEC)\n",
    x = "\nDose"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.title.x     = ggplot2::element_text(size = 12, face = "bold"),
    axis.title.y     = ggplot2::element_text(size = 12, face = "bold"),
    axis.text        = ggplot2::element_text(size = 10),
    axis.text.x      = ggplot2::element_text(colour="black"),
    axis.text.y      = ggplot2::element_text(colour="black"),
    panel.grid.minor = ggplot2::element_blank()
  )

## Recreating manuscript figure ################################################

fig_s3 <- (tfs_96_cat_plot | tfs_all_cat_plot | vvshed_cat_plot) +
  patchwork::plot_annotation(tag_levels = "A")

ggplot2::ggsave(
  filename = here::here("2_supplementary_materials", "analyses", "rethinking-results", "figs", "fig-s3.png"),
  plot = fig_s3,
  width = 13, height = 8
)


# Supplement Figure 4: Time series individual plots - Categorical Dose #########

# Load the predictions data file
ts_cat_fit <- arrow::read_parquet(
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "timeseries_cat_ind_fit.parquet")
)  |>
  dplyr::mutate(dose_cat = factor(dose_cat,
                                  levels = c(1, 2, 3),
                                  labels = c("Low", "Medium", "High")),
                across(c(estimate, lwr, upr), exp))

# Load the raw data file
ts_cat_raw <- readRDS(
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "data", "timeseries_cat_data.Rds")
)

ts_cat_plot_all <-
  ggplot2::ggplot(
    data = ts_cat_fit,
    mapping = ggplot2::aes(
      x = days,
      y = estimate,
      ymin = lwr,
      ymax = upr,
      color = dose_cat
    )
  ) +
  ggplot2::geom_ribbon(
    color = NA,
    fill = "gray",
    alpha = 0.5
  ) +
  ggplot2::geom_line() +
  ggplot2::facet_wrap(ggplot2::vars(factor(ID, levels = 1:19))) +
  ggplot2::geom_point(
    data = ts_cat_raw,
    mapping = ggplot2::aes(
      x = days,
      y = exp(logyc)
    ),
    shape = 21, size = 1,
    inherit.aes = FALSE
  ) +
  ggplot2::scale_y_continuous(
    trans = "log10",
    labels = exp_labs,
    expand = c(0,0)
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq(0, 100, 20)
  ) +
  ggplot2::coord_cartesian(
    ylim = c(10 ^ 0, 10 ^ 12),
    xlim = c(0, 90)
  ) +
  ggplot2::labs(
    y = "Viral shedding (GEC)\n",
    x = "\nTime (day)",
    color = "log(dose) - log(48)"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x     = ggplot2::element_text(size = 12, face = "bold"),
    axis.title.y     = ggplot2::element_text(size = 12, face = "bold"),
    axis.text        = ggplot2::element_text(size = 10),
    axis.text.x      = ggplot2::element_text(colour="black"),
    axis.text.y      = ggtext::element_markdown(colour="black"),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  filename = here::here("2_supplementary_materials", "analyses", "rethinking-results", "figs", "fig-s4.png"),
  plot = ts_cat_plot_all,
  width = 13, height = 8
)

# Supplement Figure 5: Time series population plots - Categorical Dose #########

# Load the predictions data file
ts_cat_pop_fit <- arrow::read_parquet(
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "timeseries_cat_pop_fit.parquet")
) |>
  dplyr::mutate(dose_cat = factor(dose_cat,
                                  levels = c(1, 2, 3),
                                  labels = c("Low", "Medium", "High")),
                across(c(mu, lwr, upr), exp))

ts_cat_pop_plot_all <-
  ts_cat_pop_fit |>
  ggplot2::ggplot(
    mapping = ggplot2::aes(
      x = days,
      y = mu,
      ymin = lwr,
      ymax = upr,
      color = dose_cat
    )
  ) +
  ggplot2::geom_ribbon(
    ggplot2::aes(fill = dose_cat, group = dose_cat),
    color = NA,
    alpha = 0.15
  ) +
  ggplot2::geom_line() +
  ggplot2::scale_y_continuous(
    trans = "log10",
    labels = exp_labs
  ) +
  ggplot2::scale_x_continuous(
    #breaks = seq(0, 7, 1)
  ) +
  ggplot2::coord_cartesian(
    ylim = c(10 ^ 0, 10 ^ 14),
    expand = FALSE
  ) +
  ggplot2::labs(
    y = "Viral shedding (GEC)\n",
    x = "\nTime (day)"
  ) +
  ggplot2::theme_bw(base_size = 20) +
  colorblindr::scale_fill_OkabeIto() +
  colorblindr::scale_color_OkabeIto()

## Zoomed in to first 7 days ####
ts_cat_pop_plot_zoom <-
  ts_cat_pop_fit |>
  dplyr::filter(days <= 7) |>
  ggplot2::ggplot(
    mapping = ggplot2::aes(
      x = days,
      y = mu,
      ymin = lwr,
      ymax = upr,
      color = dose_cat
    )
  ) +
  ggplot2::geom_ribbon(
    ggplot2::aes(fill = dose_cat, group = dose_cat),
    color = NA,
    alpha = 0.15
  ) +
  ggplot2::geom_line() +
  ggplot2::scale_y_continuous(
    trans = "log10",
    labels = exp_labs
  ) +
  ggplot2::scale_x_continuous(
    #breaks = seq(0, 7, 1)
  ) +
  ggplot2::coord_cartesian(
    ylim = c(10 ^ 0, 10 ^ 14),
    expand = FALSE
  ) +
  ggplot2::labs(
    y = "Viral shedding (GEC)\n",
    x = "\nTime (day)"
  ) +
  ggplot2::theme_bw(base_size = 20) +
  colorblindr::scale_fill_OkabeIto() +
  colorblindr::scale_color_OkabeIto()

## Compose figure ####

ts_con_pop_plot <-
  ts_cat_pop_plot_all +
  ts_cat_pop_plot_zoom +
  patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = 'bottom')

ggplot2::ggsave(
  filename = here::here("2_supplementary_materials", "analyses", "rethinking-results", "figs", "fig-s5.png"),
  plot = ts_con_pop_plot,
  width = 13, height = 8
)

# Supplement Figure 6: Time series outcome plots - categorical dose ############

outcomes_cat <- arrow::read_parquet(
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "timeseries_cat_stats.parquet")
)

# Have to make four independent panel plots since the axes need to be scaled
# differently

## Panel A: Shedding onset #####################################################
onset2 <- outcomes_cat |>
  dplyr::filter(stat == "onset") |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = factor(dose_cat),
      y = mean,
      ymin = lwr,
      ymax = upr
    )
  ) +
  geom_pointrange() +
  scale_x_discrete(
    labels = c("low", "medium", "high")
  ) +
  scale_y_continuous(
    limits = c(0, 3.2),
    breaks = seq(0, 3, by = 1)
  ) +
  labs(
    x = "dose group",
    y = "days"
  ) +
  ggplot2::theme_minimal(base_size = 20)

## Panel B: Time to peak #######################################################
peak2 <- outcomes_cat |>
  dplyr::filter(stat == "peak") |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = factor(dose_cat),
      y = mean,
      ymin = lwr,
      ymax = upr
    )
  ) +
  geom_pointrange() +
  scale_x_discrete(
    labels = c("low", "medium", "high")
  ) +
  scale_y_continuous(
    limits = c(0, 4.0),
    breaks = seq(0, 4, by = 1)
  ) +
  labs(
    x = "dose group",
    y = "days"
  ) +
  ggplot2::theme_minimal(base_size = 20)

## Panel C: Shedding duration ##################################################
dur2 <- outcomes_cat |>
  dplyr::filter(stat == "duration") |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = factor(dose_cat),
      y = mean,
      ymin = lwr,
      ymax = upr
    )
  ) +
  geom_pointrange() +
  scale_x_discrete(
    labels = c("low", "medium", "high")
  ) +
  scale_y_continuous(
    limits = c(0, 45),
    breaks = seq(0, 45, by = 10)
  ) +
  labs(
    x = "dose group",
    y = "days"
  ) +
  ggplot2::theme_minimal(base_size = 20)

## Panel D: total virus load ###################################################
auc2 <- outcomes_cat |>
  dplyr::filter(stat == "auc") |>
  dplyr::mutate(dplyr::across(c(mean, lwr, upr), exp)) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = factor(dose_cat),
      y = mean,
      ymin = lwr,
      ymax = upr
    )
  ) +
  geom_pointrange() +
  scale_x_discrete(
    labels = c("low", "medium", "high")
  ) +
  scale_y_continuous(
    trans = "log10",
    limits = 10 ^ c(0, 12),
    breaks = 10 ^ seq(0, 12, by = 4),
    labels = exp_labs
  ) +
  labs(
    x = "dose group",
    y = "days"
  ) +
  ggplot2::theme_minimal(base_size = 20)

## Figure recreation ###########################################################

fig_s6 <- onset2 + peak2 + dur2 + auc2

ggplot2::ggsave(
  filename = here::here("2_supplementary_materials", "analyses", "rethinking-results", "figs", "fig-s6.png"),
  plot = fig_s6,
  width = 13, height = 8
)

# Supplement Figure 7: Symptom outcomes - Categorical Dose #####################

## Panel A: Incubation time ####################################################

# Load the predictions data file
incub_cat_fit <- arrow::read_parquet(
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "incubation_cat_fit.parquet")
) |>
  dplyr::mutate(dose_cat = factor(dose_cat,
                                  levels = c(1, 2, 3),
                                  labels = c("Low", "Medium", "High")),
                across(c(mu, lwr, upr), exp))

# Load the raw data file
incub_cat_raw <- readRDS(
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "data", "incubation_cat_data.Rds")
)

# Make the plot
incub_cat_plot <-
  ggplot2::ggplot(
    data = incub_cat_fit,
    mapping = ggplot2::aes(
      x = dose_cat,
      y = mu, ymin = lwr, ymax = upr
    )
  ) +
  ggplot2::geom_pointrange(size = 1) +
  ggplot2::geom_point(
    data = incub_cat_raw,
    mapping = ggplot2::aes(
      x = dose_cat,
      y = incubation
    ),
    shape = 1, size = 2, stroke = 1, color = "red",
    inherit.aes = F
  ) +
  ggplot2::scale_x_discrete(
    labels = c("low", "medium", "high")
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 4),
    breaks = seq(0, 4, 1)
  ) +
  ggplot2::labs(
    y = "Incubation time (days)\n",
    x = "\nDose"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.title.x     = ggplot2::element_text(size = 12, face = "bold"),
    axis.title.y     = ggplot2::element_text(size = 12, face = "bold"),
    axis.text        = ggplot2::element_text(size = 10),
    axis.text.x      = ggplot2::element_text(colour="black"),
    axis.text.y      = ggplot2::element_text(colour="black"),
    panel.grid.minor = ggplot2::element_blank()
  )

## Panel B: MVS categorical dose ###############################################

# Load the predictions data file
mvs_cat_fit <- arrow::read_parquet(
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "mvs_cat_fit.parquet")
) |>
  dplyr::mutate(dose_cat = factor(dose_cat,
                                  levels = c(1, 2, 3),
                                  labels = c("Low", "Medium", "High")),
                across(c(mu, lwr, upr), exp))

# Load the raw data file
mvs_cat_raw <- readRDS(
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "data", "mvs_cat_data.Rds")
)

# Make the plot
mvs_cat_plot <-
  ggplot2::ggplot(
    data = mvs_cat_fit,
    mapping = ggplot2::aes(
      x = dose_cat,
      y = mu, ymin = lwr, ymax = upr
    )
  ) +
  ggplot2::geom_pointrange(size = 1) +
  ggplot2::geom_point(
    data = mvs_cat_raw,
    mapping = ggplot2::aes(
      x = dose_cat,
      y = mvs
    ),
    shape = 1, size = 2, stroke = 1, color = "red",
    inherit.aes = F
  ) +
  ggplot2::scale_x_discrete(
    labels = c("low", "medium", "high")
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 10),
    breaks = c(0, 5, 10),
    minor_breaks = seq(0, 10, 1)
  ) +
  ggplot2::labs(
    y = "Modified Vesikari score\n",
    x = "\nDose"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.title.x     = ggplot2::element_text(size = 12, face = "bold"),
    axis.title.y     = ggplot2::element_text(size = 12, face = "bold"),
    axis.text        = ggplot2::element_text(size = 10),
    axis.text.x      = ggplot2::element_text(colour="black"),
    axis.text.y      = ggplot2::element_text(colour="black")
  )

## Panel C: CSS categorical dose ###############################################

# Load the predictions data file
css_cat_fit <- arrow::read_parquet(
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "css_cat_fit.parquet")
) |>
  dplyr::mutate(dose_cat = factor(dose_cat,
                                  levels = c(1, 2, 3),
                                  labels = c("Low", "Medium", "High")),
                across(c(mu, lwr, upr), exp))

# Load the raw data file
css_cat_raw <- readRDS(
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "data", "css_cat_data.Rds")
)

# Make the plot
css_cat_plot <-
  ggplot2::ggplot(
    data = css_cat_fit,
    mapping = ggplot2::aes(
      x = dose_cat,
      y = mu, ymin = lwr, ymax = upr
    )
  ) +
  ggplot2::geom_pointrange(size = 1) +
  ggplot2::geom_point(
    data = css_cat_raw,
    mapping = ggplot2::aes(
      x = dose_cat,
      y = outcome
    ),
    shape = 1, size = 2, stroke = 1, color = "red",
    inherit.aes = F
  ) +
  ggplot2::scale_x_discrete(
    labels = c("low", "medium", "high")
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 35),
    breaks = seq(0, 35, 5),
    minor_breaks = seq(0, 35, 1)
  ) +
  ggplot2::labs(
    y = "Comprehensive symptom score\n",
    x = "\nDose"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.title.x     = ggplot2::element_text(size = 12, face = "bold"),
    axis.title.y     = ggplot2::element_text(size = 12, face = "bold"),
    axis.text        = ggplot2::element_text(size = 10),
    axis.text.x      = ggplot2::element_text(colour="black"),
    axis.text.y      = ggplot2::element_text(colour="black")
  )

## Recreate the figure #########################################################

fig_s7 <- (incub_cat_plot | mvs_cat_plot | css_cat_plot) +
  patchwork::plot_annotation(tag_levels = "A")
ggplot2::ggsave(filename = here::here("2_supplementary_materials", "analyses", "rethinking-results", "figs", "fig-s7.png"), plot = fig_s7,
                width = 13, height = 8)

# former manuscript Figure 2 - Time series individual curves #####################

# Load the predictions data file
ts_con_fit <- arrow::read_parquet(
  file = here::here("2_supplementary_materials", "analyses", "rethinking-results", "preds", "timeseries_con_ind_fit.parquet")
)

# Load the raw data file
ts_con_raw <- readRDS(
  here::here("2_supplementary_materials", "analyses", "rethinking-results", "data", "timeseries_con_data.Rds")
)

ts_con_plot_all <-
  ggplot2::ggplot(
    data = ts_con_fit,
    mapping = ggplot2::aes(
      x = days,
      y = exp(estimate),
      ymin = exp(lwr),
      ymax = exp(upr),
      color = factor(dose_log)
    )
  ) +
  ggplot2::geom_ribbon(
    color = NA,
    fill = "gray",
    alpha = 0.5
  ) +
  ggplot2::geom_line() +
  ggplot2::facet_wrap(ggplot2::vars(factor(ID, levels = 19:1))) +
  ggplot2::geom_point(
    data = ts_con_raw,
    mapping = ggplot2::aes(
      x = days,
      y = exp(logyc)
    ),
    shape = 21, size = 1,
    inherit.aes = FALSE
  ) +
  ggplot2::scale_y_continuous(
    trans = "log10",
    labels = exp_labs,
    expand = c(0,0)
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq(0, 100, 20)
  ) +
  ggplot2::coord_cartesian(
    ylim = c(10 ^ 0, 10 ^ 12),
    xlim = c(0, 90)
  ) +
  ggplot2::labs(
    y = "Viral shedding (GEC)\n",
    x = "\nTime (day)",
    color = "log(dose) - log(48)"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x     = ggplot2::element_text(size = 12, face = "bold"),
    axis.title.y     = ggplot2::element_text(size = 12, face = "bold"),
    axis.text        = ggplot2::element_text(size = 10),
    axis.text.x      = ggplot2::element_text(colour="black"),
    axis.text.y      = ggtext::element_markdown(colour="black"),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  filename = here::here("2_supplementary_materials", "analyses", "rethinking-results", "figs", "fig-f2.png"),
  plot = ts_con_plot_all,
  width = 13, height = 8
)

ts_con_plot_all_zoom <-
  ts_con_fit |>
  dplyr::filter(days <= 7) |>
  ggplot2::ggplot(
    mapping = ggplot2::aes(
      x = days,
      y = exp(estimate),
      ymin = exp(lwr),
      ymax = exp(upr),
      color = factor(dose_log)
    )
  ) +
  ggplot2::geom_ribbon(
    color = NA,
    fill = "gray",
    alpha = 0.5
  ) +
  ggplot2::geom_line() +
  ggplot2::facet_wrap(ggplot2::vars(factor(ID, levels = 19:1))) +
  ggplot2::geom_point(
    data = ts_con_raw |> dplyr::filter(days <= 7),
    mapping = ggplot2::aes(
      x = days,
      y = exp(logyc)
    ),
    shape = 21, size = 1,
    inherit.aes = FALSE
  ) +
  ggplot2::scale_y_continuous(
    trans = "log10",
    labels = exp_labs
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq(0, 7, 1)
  ) +
  ggplot2::coord_cartesian(
    ylim = c(10 ^ 0, 10 ^ 12),
    xlim = c(0, 7),
    expand = FALSE
  ) +
  ggplot2::labs(
    y = "Viral shedding (GEC)\n",
    x = "\nTime (day)",
    color = "log(dose) - log(48)"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x     = ggplot2::element_text(size = 12, face = "bold"),
    axis.title.y     = ggplot2::element_text(size = 12, face = "bold"),
    axis.text        = ggplot2::element_text(size = 10),
    axis.text.x      = ggplot2::element_text(colour="black"),
    axis.text.y      = ggtext::element_markdown(colour="black"),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  filename = here::here("2_supplementary_materials", "analyses", "rethinking-results", "figs", "fig-f2-zoom.png"),
  plot = ts_con_plot_all_zoom,
  width = 13, height = 8
)

# END OF FILE ##################################################################
