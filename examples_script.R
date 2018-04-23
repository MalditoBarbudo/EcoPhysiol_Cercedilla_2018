### Script análisis de ejemplos

# libraries
library(sapfluxnetr)
library(tidyverse)
library(ggridges)
library(taxonlookup)
library(lubridate)
library(viridis)
library(quantreg)

# data
# Cojo los 75 sitios que ha hecho Rafa
load('tidy_metrics_plant.RData')

#### ggridges by genera vs sp ####

tidy_metrics_plant %>%
  # mutate(pl_functional_group = lookup_table(pl_species)[['group']]) %>%
  filter(
    str_detect(pl_species, 'Pinus') | 
      str_detect(pl_species, 'Quercus') | 
      str_detect(pl_species, 'Eucalyptus')
  ) -> genera_tidy

genera_tidy %>%
  filter(str_detect(pl_species, 'Pinus')) %>%
  mutate(month = month(TIMESTAMP, label = TRUE)) %>%
  ggplot(aes(x = sapflow_q_95, y = month, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  # geom_density_ridges_gradient(scale = 3) +
  scale_fill_viridis(name = 'Sap Flow [cm3/h]', option = 'C') +
  labs(
    x = 'Maximum Sap Flow as Q95 [cm3/h]',
    y = 'Month',
    title = 'Month variation of daily maximum flows',
    subtitle = 'Pinus species'
  ) +
  facet_wrap(~ pl_species, ncol = 2, scales = 'free_x') -> pinus_joy
pinus_joy

genera_tidy %>%
  filter(str_detect(pl_species, 'Quercus')) %>%
  mutate(month = month(TIMESTAMP, label = TRUE)) %>%
  ggplot(aes(x = sapflow_q_95, y = month, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  # geom_density_ridges_gradient(scale = 3) +
  scale_fill_viridis(name = 'Sap Flow [cm3/h]', option = 'C') +
  labs(
    x = 'Maximum Sap Flow as Q95 [cm3/h]',
    y = 'Month',
    title = 'Month variation of daily maximum flows',
    subtitle = 'Quercus species'
  ) +
  facet_wrap(~ pl_species, ncol = 2, scales = 'free_x') -> quercus_joy
quercus_joy

genera_tidy %>%
  filter(str_detect(pl_species, 'Eucalyptus')) %>%
  mutate(month = month(TIMESTAMP, label = TRUE)) %>%
  ggplot(aes(x = sapflow_q_95, y = month, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  # geom_density_ridges_gradient(scale = 3) +
  scale_fill_viridis(name = 'Sap Flow [cm3/h]', option = 'C') +
  labs(
    x = 'Maximum Sap Flow as Q95 [cm3/h]',
    y = 'Month',
    title = 'Month variation of daily maximum flows',
    subtitle = 'Eucalyptus species'
  ) +
  facet_wrap(~ pl_species, ncol = 2, scales = 'free_x') -> eucalyptus_joy
eucalyptus_joy

#### vpd & swc vs max Jt ####

tidy_metrics_taxonized <- lookup_table(
  unique(tidy_metrics_plant[['pl_species']]), by_species = TRUE, missing_action = 'NA'
) %>%
  rownames_to_column('pl_species') %>%
  left_join(tidy_metrics_plant, ., by = 'pl_species') %>%
  filter(si_country != 'RUS',
         !si_code %in% c('USA_HIL_HF2'),
         sapflow_q_95 > 0,
         vpd_q_95 > 0,
         genus %in% c('Eucalyptus', 'Quercus', 'Pinus'),
         !pl_species %in% c('Quercus rotundifolia'))

quantreg::nlrq(
  # sapflow_q_95 ~ SSasympOrig(vpd_q_95, Asym, lrc),
  sapflow_q_95 ~ Asym * (1 - exp(-exp(lrc) * vpd_q_95)),
  start = list(Asym = 100000, lrc = 0.9),
  tidy_metrics_taxonized %>% filter(genus == 'Eucalyptus'),
  tau = 0.99
) -> euc_model

quantreg::nlrq(
  sapflow_q_95 ~ Asym * (1 - exp(-exp(lrc) * vpd_q_95)),
  tidy_metrics_taxonized %>% filter(genus == 'Quercus'),
  start = list(Asym = 25000, lrc = 0.9),
  tau = 0.99
) -> que_model

quantreg::nlrq(
  sapflow_q_95 ~ Asym * (1 - exp(-exp(lrc) * vpd_q_95)),
  start = list(Asym = 25000, lrc = 0.9),
  tidy_metrics_taxonized %>% filter(genus == 'Pinus'),
  tau = 0.99
) -> pin_model

seqvpd <- data.frame(vpd = seq(0, 10, by = 0.05))
euc_line <- data_frame(
  vpd_q_95 = tidy_metrics_taxonized %>%
    filter(genus == 'Eucalyptus') %>%
    pull(vpd_q_95),
  genus = 'Eucalyptus',
  fit = predict(euc_model, newdata = seqvpd)
)
pin_line <- data_frame(
  vpd_q_95 = tidy_metrics_taxonized %>%
    filter(genus == 'Pinus') %>%
    pull(vpd_q_95),
  genus = 'Pinus',
  fit = predict(pin_model, newdata = seqvpd)
)
que_line <- data_frame(
  vpd_q_95 = tidy_metrics_taxonized %>%
    filter(genus == 'Quercus') %>%
    pull(vpd_q_95),
  genus = 'Quercus',
  fit = predict(que_model, newdata = seqvpd)
)

vpd_response <- tidy_metrics_taxonized %>%
  ggplot(aes(x = vpd_q_95, y = sapflow_q_95, colour = genus)) +
  geom_point(alpha = 0.1) +
  geom_line(aes(x = vpd_q_95, y = fit), data = euc_line, size = 2) +
  geom_line(aes(x = vpd_q_95, y = fit), data = pin_line, size = 2) +
  geom_line(aes(x = vpd_q_95, y = fit), data = que_line, size = 2) +
  # stat_smooth(se = FALSE, size = 2) +
  # facet_wrap(~ pl_species, ncol = 3, scales = 'fixed') +
  scale_colour_viridis(discrete = TRUE, option = 'C') +
  labs(
    colour = '',
    x = 'VPD max [kPA]',
    y = 'Jt max [cm3 h]',
    title = 'Sap flow responses to VPD',
    subtitle = 'by genera'
  ) +
  theme(legend.position = 'bottom')
vpd_response

# swc_shallow_response <- tidy_metrics_taxonized %>%
#   ggplot(aes(x = swc_shallow_q_95, y = sapflow_q_95, colour = group)) +
#   geom_point(alpha = 0.1) +
#   # stat_smooth(se = FALSE, size = 2) +
#   # facet_wrap(~ si_biome, ncol = 2, scales = 'free') +
#   scale_colour_viridis(discrete = TRUE, option = 'C') +
#   labs(
#     colour = '',
#     x = 'SWC [cm3/cm3]',
#     y = 'Js max [cm3 h]',
#     title = 'Sap flow responses to Soil Water Content',
#     subtitle = 'by functional type'
#   )
# swc_shallow_response

#### individual trees (escalado max max) ####

tidy_metrics_taxonized %>%
  group_by(pl_code) %>%
  summarise(sapflow = quantile(sapflow_q_95, na.rm = TRUE, probs = 0.95),
            pl_dbh = mean(pl_dbh, na.rm = TRUE),
            st_basal_area = mean(st_basal_area, na.rm = TRUE),
            group = unique(group),
            genus = unique(genus),
            si_biome = unique(si_biome),
            pl_species = unique(pl_species),
            st_basal_area = mean(st_basal_area, na.rm = TRUE)) %>%
  filter(genus %in% c('Eucalyptus', 'Quercus', 'Pinus')) %>%
  ggplot(aes(x = pl_dbh, y = sapflow, colour = genus)) +
  geom_point(aes(size = st_basal_area), alpha = 0.4) +
  stat_smooth(se = FALSE, method = 'lm', size = 2) +
  # facet_wrap(~ pl_species, ncol = 3) +
  scale_colour_viridis(discrete = TRUE, option = 'C') +
  labs(
    colour = '', size = '',
    x = 'DBH [cm]',
    y = 'Jt max [cm3/day]',
    title = 'Maximum daily flow responses to DBH',
    subtitle = 'by genera'
  ) +
  theme(legend.position = 'bottom')
