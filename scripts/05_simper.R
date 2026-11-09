# 05_simper.R
# Which families account for the difference between HI and NI communities? (SIMPER)
#
# Run from the project root, after scripts/01_data_prep.R.
#
# Inputs : data/processed/community_family_long.csv
# Outputs: outputs/tables/simper_hi_vs_ni.csv
#          outputs/figures/simper_hi_vs_ni.png
#          outputs/tables/cp_composition_stations.csv, simper_family_cp.csv
#          outputs/figures/cp_composition.png
#
# Notes
# - SIMPER splits the average Bray-Curtis dissimilarity between two groups into a
#   contribution from each family. It is descriptive: no p values are computed, because
#   the two years at a station are not independent and there are only three stations
#   per group.
# - Same data as the other analyses: family level, two layers averaged, square root percent.
# - Rankings are also calculated separately for each year (3 vs 3 stations) to check
#   whether the same families come out on top both times.

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(vegan)

dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

fam <- read_csv("data/processed/community_family_long.csv", show_col_types = FALSE)

station_year <- fam |>
  summarise(rel_abund = mean(rel_abund), .by = c(year, station, impact, family)) |>
  filter(impact %in% c("HI", "NI"))

meta <- station_year |>
  distinct(year, station, impact) |>
  arrange(year, impact, station) |>
  mutate(id = paste(station, year, sep = "_"))

comm <- station_year |>
  mutate(id = paste(station, year, sep = "_")) |>
  select(id, family, rel_abund) |>
  pivot_wider(names_from = family, values_from = rel_abund, values_fill = 0) |>
  slice(match(meta$id, id)) |>
  tibble::column_to_rownames("id") |>
  as.matrix()

stopifnot(nrow(comm) == 12, all(rownames(comm) == meta$id))

# ---- SIMPER helper ---------------------------------------------------------
run_simper <- function(rows) {
  x <- sqrt(comm[rows, , drop = FALSE])
  grp <- factor(meta$impact[rows], levels = c("HI", "NI"))
  s <- simper(x, group = grp, permutations = 0)
  out <- summary(s)[[1]]
  tibble(family = rownames(out), contribution = out$average, cumulative = out$cumsum) |>
    mutate(rank = row_number())
}

all_rows  <- seq_len(nrow(meta))
res_all   <- run_simper(all_rows)
res_2010  <- run_simper(which(meta$year == 2010)) |> select(family, rank_2010 = rank)
res_2011  <- run_simper(which(meta$year == 2011)) |> select(family, rank_2011 = rank)

# Average percent in each group (untransformed, for readability)
means <- station_year |>
  summarise(mean_pct = mean(rel_abund), .by = c(family, impact)) |>
  pivot_wider(names_from = impact, values_from = mean_pct, names_prefix = "mean_pct_")

simper_table <- res_all |>
  left_join(means, by = "family") |>
  left_join(res_2010, by = "family") |>
  left_join(res_2011, by = "family") |>
  mutate(higher_in = if_else(mean_pct_HI > mean_pct_NI, "HI", "NI")) |>
  select(rank, family, contribution, cumulative, mean_pct_HI, mean_pct_NI, higher_in,
         rank_2010, rank_2011)

write_csv(simper_table, "outputs/tables/simper_hi_vs_ni.csv")

# ---- Figure ----------------------------------------------------------------
top <- simper_table |>
  slice_head(n = 8) |>
  mutate(family = factor(family, levels = rev(family)),
         higher_in = factor(higher_in, levels = c("HI", "NI")))

simper_plot <- ggplot(top, aes(contribution, family, fill = higher_in)) +
  geom_col() +
  labs(
    x = "Contribution to average HI vs NI dissimilarity",
    y = NULL, fill = "More abundant in",
    subtitle = "Eight largest contributors (SIMPER, family level, square root percent)"
  ) +
  theme_bw()

ggsave("outputs/figures/simper_hi_vs_ni.png", simper_plot, width = 7, height = 4.5, dpi = 300)

# ---- Colonizer-persister (c-p) check ---------------------------------------
# c-p values run from 1 (extreme colonizers: fast reproduction, tolerant of disturbance)
# to 5 (persisters: slow reproduction, sensitive). If HI stations are dominated by
# opportunistic nematodes, they should have more of the community in low c-p classes.
genus_long <- read_csv("data/processed/community_genus_long.csv", show_col_types = FALSE)

# Percent of each station's community in each c-p class (two layers averaged)
cp_stations <- genus_long |>
  summarise(pct = sum(rel_abund) / n_distinct(sample), .by = c(year, station, impact, cp)) |>
  complete(nesting(year, station, impact), cp = 1:5, fill = list(pct = 0)) |>
  arrange(year, impact, station, cp)
write_csv(cp_stations, "outputs/tables/cp_composition_stations.csv")

# c-p of the largest SIMPER families (abundance weighted across their genera)
simper_family_cp <- genus_long |>
  filter(family %in% simper_table$family[1:8]) |>
  summarise(
    mean_cp = sum(rel_abund * cp) / sum(rel_abund),
    cp_min = min(cp), cp_max = max(cp),
    .by = family
  ) |>
  right_join(simper_table |> select(rank, family, higher_in) |> slice_head(n = 8), by = "family") |>
  arrange(rank) |>
  select(rank, family, higher_in, mean_cp, cp_min, cp_max)
write_csv(simper_family_cp, "outputs/tables/simper_family_cp.csv")

# Which genera make up the c-p 1 animals at HI stations?
cp1_genera <- genus_long |>
  filter(cp == 1, impact == "HI") |>
  summarise(total = sum(rel_abund), .by = c(family, genus)) |>
  mutate(share_of_cp1_at_HI = total / sum(total)) |>
  arrange(desc(share_of_cp1_at_HI)) |>
  select(family, genus, share_of_cp1_at_HI)
write_csv(cp1_genera, "outputs/tables/cp1_genera_at_hi.csv")

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

message("Done. SIMPER and c-p tables and figures written to outputs/.")
