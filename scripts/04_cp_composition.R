# ==============================================================================
# 04  COLONIZER-PERSISTER (c-p) COMPOSITION
#     Do heavily impacted stations have more colonizer nematodes?
# ==============================================================================
#
# WHY THIS SCRIPT EXISTS
#   Scripts 02 and 03 show that heavily impacted (HI) and unimpacted (NI)
#   communities differ, but not what the difference is made of. Looking at the
#   composition, the most obvious contrast was in the colonizer-persister (c-p)
#   classes. That finding is exploratory, and it is what motivated the colonizer
#   model in script 05.
#
# WHAT THIS SCRIPT DOES
#   Every genus has a c-p score from 1 to 5.
#     1  colonizers: reproduce fast, tolerate disturbance
#     5  persisters: reproduce slowly, sensitive to disturbance
#   For each station and year we work out the share of the community in each c-p
#   class. If HI stations are dominated by opportunists, they should have more of
#   their community in the low classes, especially class 1. We then ask which
#   genera make up the class 1 animals at HI stations.
#
# RUN
#   From the project root, after 01_data_prep.R.
#
# INPUT
#   data/processed/community_genus_long.csv
#   (genus level is needed here because the c-p score belongs to each genus)
#
# OUTPUTS
#   outputs/tables/cp_composition_stations.csv   percent in each c-p class, by station and year
#   outputs/tables/cp1_genera_at_hi.csv          which genera make up the c-p 1 animals at HI stations
#   outputs/figures/cp_composition.png           Figure 3 in the paper
#
# NOTES
#   * This is descriptive. No test is run here; the formal comparison is the
#     colonizer model in script 05.
#   * Percentages are averaged over the two sediment layers, as in the other scripts.

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

genus_long <- read_csv("data/processed/community_genus_long.csv", show_col_types = FALSE)

# ---- 1. Share of each community in each c-p class --------------------------
# Summing rel_abund over the two layers and dividing by the number of samples
# averages the layers. complete(... cp = 1:5, fill = 0) makes sure every class
# appears for every station, even where it is absent, so bars of 0 are drawn and
# not silently left out.
cp_stations <- genus_long |>
  summarise(pct = sum(rel_abund) / n_distinct(sample), .by = c(year, station, impact, cp)) |>
  complete(nesting(year, station, impact), cp = 1:5, fill = list(pct = 0)) |>
  arrange(year, impact, station, cp)
write_csv(cp_stations, "outputs/tables/cp_composition_stations.csv")

# ---- 2. Which genera are the colonizers at HI stations? --------------------
# Each row's share is its portion of all c-p 1 animals at HI stations. In the
# current data nearly all of them are one genus, Halomonhystera.
cp1_genera <- genus_long |>
  filter(cp == 1, impact == "HI") |>
  summarise(total = sum(rel_abund), .by = c(family, genus)) |>
  mutate(share_of_cp1_at_HI = total / sum(total)) |>
  arrange(desc(share_of_cp1_at_HI)) |>
  select(family, genus, share_of_cp1_at_HI)
write_csv(cp1_genera, "outputs/tables/cp1_genera_at_hi.csv")

# ---- 3. Figure --------------------------------------------------------------
# Bars = group means; points = individual stations, so the spread among stations
# is visible and not hidden behind the average. One panel per year. MI stations
# are left out to keep the HI versus NI contrast clear.
cp_plot_df <- cp_stations |>
  filter(impact %in% c("HI", "NI")) |>
  mutate(impact = factor(impact, levels = c("HI", "NI")), cp = factor(cp))

cp_plot <- ggplot(cp_plot_df, aes(cp, pct, fill = impact)) +
  stat_summary(fun = mean, geom = "col", position = position_dodge(width = 0.8),
               width = 0.7, alpha = 0.55) +
  geom_point(aes(colour = impact), position = position_dodge(width = 0.8), size = 1.8,
             show.legend = FALSE) +
  facet_wrap(~year) +
  labs(x = "c-p class (1 = colonizers, 5 = persisters)",
       y = "Percent of community", fill = "Impact",
       subtitle = "Bars are group means; points are individual stations") +
  theme_bw()

ggsave("outputs/figures/cp_composition.png", cp_plot, width = 8, height = 4.2, dpi = 300)

message("Done. c-p tables and figure written to outputs/.")
