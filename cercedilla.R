# libraries
library(sapfluxnetr)
library(tidyverse)
library(taxonlookup)
library(ggridges)
library(viridis)

# data
# filter_by_var(!is.null(si_code), folder = 'plant', .use_cache = TRUE) %>%
#   read_sfn_data(folder = 'plant') %>%
#   sfn_metrics('daily', solar = TRUE,
#               .funs = sapfluxnetr:::.fixed_metrics_funs(c(0.95, 0.99), FALSE),
#               general = TRUE, predawn = FALSE, midday = FALSE, nighttime = FALSE) %>%
#   metrics_tidyfier(metadata = sfn_metadata) -> tidy_metrics_plant

# save(tidy_metrics_plant, file = 'tidy_metrics_plant.RData')

load('tidy_metrics_plant.RData')

#### escalado maximo ####
tidy_metrics_plant %>%
  pull(pl_species) %>%
  unique() %>%
  lookup_table(missing_action = 'NA', by_species = TRUE) %>%
  rownames_to_column('pl_species') %>%
  left_join(tidy_metrics_plant, ., by = 'pl_species') -> tidy_metrics_plant_phylo

tidy_metrics_plant_phylo %>%
  group_by(pl_code) %>%
  summarise(
    sapflow_q_95_plant = quantile(sapflow_q_95, 0.95, na.rm = TRUE),
    pl_dbh = mean(pl_dbh, na.rm = TRUE),
    si_biome = unique(si_biome),
    pl_species = unique(pl_species),
    group = unique(group),
    st_basal_area = mean(st_basal_area, na.rm = TRUE)
  ) %>%
  ggplot(aes(x = pl_dbh, y = sapflow_q_95_plant, colour = group)) +
  geom_point(aes(size = st_basal_area), alpha = 0.4) +
  stat_smooth(method = 'lm', se = FALSE) +
  scale_colour_viridis(discrete = TRUE, option = 'D') +
  # facet_wrap(~ si_biome, scales = 'free') +
  theme(legend.position = 'bottom') -> escalado_max_plot_dbh
escalado_max_plot_dbh

# ggsave('escalado_max_plot_dbh.pdf', escalado_max_plot_dbh, 'pdf')

#### q_95 vs env_vars ####
tidy_metrics_plant_phylo %>%
  ggplot(aes(x = vpd_q_95, y = sapflow_q_95, colour = group)) +
  geom_point(alpha = 0.1) +
  geom_quantile(quantiles = 0.95) +
  scale_colour_viridis(discrete = TRUE, option = 'C') +
  labs
  theme(legend.position = 'bottom') +
  scale_x_continuous(limits = c(0, 10)) + scale_y_continuous(limits = c(0, 130000))

tidy_metrics_plant_phylo %>%
  ggplot(aes(x = swc_shallow_q_95, y = sapflow_q_95, colour = group)) +
  geom_point(alpha = 0.1) +
  geom_quantile(quantiles = 0.95) +
  scale_colour_viridis(discrete = TRUE, option = 'C') +
  theme(legend.position = 'bottom') +
  scale_y_continuous(limits = c(0, 130000))

#### ggjoys per genera ####
