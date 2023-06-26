# WRAPPER SCRIPT
# THIS RUNS THE RETHINKING ANALYSIS IN THE CORRECT ORDER
# AND CLEANS THE ENVIRONMENT BETWEEN SCRIPTS

source(
  here::here("2_supplementary_materials/analyses/rethinking-analysis/models.R"),
  local = new.env(),
  echo = TRUE
)

source(
  here::here("2_supplementary_materials/analyses/rethinking-analysis/preds.R"),
  local = new.env(),
  echo = TRUE
)

source(
  here::here("2_supplementary_materials/analyses/rethinking-analysis/figs.R"),
  local = new.env(),
  echo = TRUE
)