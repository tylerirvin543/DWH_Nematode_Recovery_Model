# 04_extra_plots.R
# Supporting figures: environmental PCA (per year), family composition, MI and ITD, site map.
#
# Run from the project root, after scripts/01_data_prep.R.
#
# Inputs : data/processed/env_stations.csv, samples.csv, community_family_long.csv
# Outputs: outputs/figures/env_pca.png, family_composition.png, indices_mi_itd.png, site_map.png
#          outputs/tables/env_pca_variance.csv, env_pca_loadings.csv
#
# Notes
# - The PCA is run separately for each year and variables are standardised within a year,
#   because the 2010 and 2011 environmental values appear to be on different scales
#   (see the cleaning log). Axes are therefore not comparable between the two panels.
# - Skewed variables (Ba, PAH44, TPH) are log10 transformed before the PCA.

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

impact_levels <- c("HI", "MI", "NI")

env <- read_csv("data/processed/env_stations.csv", show_col_types = FALSE) |>
  mutate(impact = factor(impact, levels = impact_levels))
samples <- read_csv("data/processed/samples.csv", show_col_types = FALSE) |>
  mutate(impact = factor(impact, levels = impact_levels))
fam <- read_csv("data/processed/community_family_long.csv", show_col_types = FALSE) |>
  mutate(impact = factor(impact, levels = impact_levels))

# ---- 1. Environmental PCA, one per year ------------------------------------
env_vars <- c("depth_m", "Al", "Ba", "PAH44", "TPH", "Mud", "TOC")
log_vars <- c("Ba", "PAH44", "TPH")

run_pca <- function(yr) {
  d <- env |> filter(year == yr)
  x <- d |>
    select(all_of(env_vars)) |>
    mutate(across(all_of(log_vars), log10))
  pca <- prcomp(x, center = TRUE, scale. = TRUE)

  scores <- bind_cols(d |> select(year, station, impact), as_tibble(pca$x[, 1:2]))
  loadings <- as_tibble(pca$rotation[, 1:2], rownames = "variable") |> mutate(year = yr)
  variance <- tibble(year = yr, PC = paste0("PC", 1:2),
                     prop_variance = (pca$sdev^2 / sum(pca$sdev^2))[1:2])
  list(scores = scores, loadings = loadings, variance = variance)
}

pcas <- lapply(c(2010, 2011), run_pca)
scores <- bind_rows(lapply(pcas, `[[`, "scores"))
loadings <- bind_rows(lapply(pcas, `[[`, "loadings"))
variance <- bind_rows(lapply(pcas, `[[`, "variance"))

write_csv(variance, "outputs/tables/env_pca_variance.csv")
write_csv(loadings, "outputs/tables/env_pca_loadings.csv")

# Scale loading arrows to the score range so they are visible
arrow_scale <- scores |>
  summarise(s = max(abs(c(PC1, PC2))) / max(abs(c(loadings$PC1, loadings$PC2)))) |>
  pull(s) * 0.8

pc_lab <- variance |>
  mutate(label = sprintf("%s (%.0f%%)", PC, 100 * prop_variance)) |>
  select(year, PC, label) |>
  pivot_wider(names_from = PC, values_from = label)

scores <- scores |> left_join(pc_lab, by = "year", suffix = c("", "_lab"))
panel_labels <- pc_lab |>
  transmute(year, panel = sprintf("%d: %s, %s", year, PC1, PC2))
scores <- scores |> left_join(panel_labels, by = "year")
loadings <- loadings |> left_join(panel_labels, by = "year")

pca_plot <- ggplot(scores, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, colour = "grey85") +
  geom_vline(xintercept = 0, colour = "grey85") +
  geom_segment(
    data = loadings,
    aes(x = 0, y = 0, xend = PC1 * arrow_scale, yend = PC2 * arrow_scale),
    arrow = arrow(length = unit(0.15, "cm")), colour = "grey40"
  ) +
  geom_text(
    data = loadings,
    aes(x = PC1 * arrow_scale * 1.12, y = PC2 * arrow_scale * 1.12, label = variable),
    size = 3, colour = "grey30"
  ) +
  geom_point(aes(colour = impact), size = 3) +
  geom_text(aes(label = station, colour = impact), size = 2.8, vjust = -1, show.legend = FALSE) +
  facet_wrap(~panel) +
  labs(
    x = "PC1 score", y = "PC2 score", colour = "Impact",
    subtitle = "PCA run separately for each year. Ba, PAH44 and TPH are log10 transformed."
  ) +
  theme_bw()

ggsave("outputs/figures/env_pca.png", pca_plot, width = 9, height = 4.8, dpi = 300)

# ---- 2. Family composition bars --------------------------------------------
station_year <- fam |>
  summarise(rel_abund = mean(rel_abund), .by = c(year, station, impact, family))

top_fams <- station_year |>
  summarise(m = mean(rel_abund), .by = family) |>
  slice_max(m, n = 10) |>
  pull(family)

comp <- station_year |>
  mutate(family = if_else(family %in% top_fams, family, "Other")) |>
  summarise(rel_abund = sum(rel_abund), .by = c(year, station, impact, family)) |>
  mutate(
    family = factor(family, levels = c(top_fams, "Other")),
    station = factor(station, levels = unique(station[order(impact, station)]))
  )

comp_plot <- ggplot(comp, aes(rel_abund, station, fill = family)) +
  geom_col() +
  facet_grid(impact ~ year, scales = "free_y", space = "free_y") +
  labs(x = "Percent of nematodes (average of two layers)", y = NULL, fill = "Family",
       subtitle = "Ten most abundant families overall; all others grouped as Other.") +
  theme_bw()

ggsave("outputs/figures/family_composition.png", comp_plot, width = 9, height = 5.5, dpi = 300)

# ---- 3. MI and ITD by impact and year --------------------------------------
idx <- samples |>
  summarise(MI = mean(MI), ITD = mean(ITD), .by = c(year, station, impact)) |>
  pivot_longer(c(MI, ITD), names_to = "index", values_to = "value") |>
  mutate(index = recode(index,
                        MI = "Maturity index (MI)",
                        ITD = "Index of trophic diversity (ITD)"))

idx_plot <- ggplot(idx, aes(year, value)) +
  geom_line(aes(group = station, colour = impact), alpha = 0.7) +
  geom_point(aes(colour = impact), size = 2.5) +
  scale_x_continuous(breaks = c(2010, 2011)) +
  facet_wrap(~index, scales = "free_y") +
  labs(x = "Year", y = "Station value (average of two layers)", colour = "Impact") +
  theme_bw()

ggsave("outputs/figures/indices_mi_itd.png", idx_plot, width = 8, height = 4.2, dpi = 300)

# ---- 4. Site map -----------------------------------------------------------
sites <- env |> filter(year == 2010)

site_plot <- ggplot(sites, aes(long, lat)) +
  geom_point(aes(colour = impact), size = 3.5) +
  geom_text(aes(label = station, colour = impact), size = 3, vjust = -1.1, show.legend = FALSE) +
  coord_fixed(ratio = 1 / cos(mean(sites$lat) * pi / 180)) +
  labs(x = "Longitude (W)", y = "Latitude (N)", colour = "Impact",
       subtitle = "Station locations (2010 coordinates)") +
  theme_bw()

ggsave("outputs/figures/site_map.png", site_plot, width = 6.5, height = 5.5, dpi = 300)

message("Done. Figures written to outputs/figures/.")
