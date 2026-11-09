# 05_simper.R
# Which families account for the difference between HI and NI communities? (SIMPER)
#
# Run from the project root, after scripts/01_data_prep.R.
#
# Inputs : data/processed/community_family_long.csv
# Outputs: outputs/tables/simper_hi_vs_ni.csv
#          outputs/figures/simper_hi_vs_ni.png
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

message("Done. SIMPER table and figure written to outputs/.")
