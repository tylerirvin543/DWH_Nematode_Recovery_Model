# 03_distance_to_control.R
# Distance of impacted stations to the control (NI) stations, compared with the
# natural control-to-control background, and a mixed model of change over time.
#
# Run from the project root, after scripts/01_data_prep.R.
#
# Inputs : data/processed/community_family_long.csv
# Outputs: outputs/figures/distance_to_control.png
#          outputs/tables/distance_to_control.csv   per impacted station-year
#          outputs/tables/control_background.csv    control-to-control distance per year
#          outputs/tables/distance_change.csv       per-station change 2010 -> 2011
#          outputs/tables/distance_model_effects.csv  mixed model estimates
#
# Definitions
# - distance: for each HI or MI station and year, the mean Bray-Curtis dissimilarity to the
#   three NI stations sampled that year.
# - background: for each year, the mean Bray-Curtis dissimilarity over the three pairs of NI
#   stations. "Recovery" means distance approaches the background, not zero.
# - gap: distance minus that year's background.
# Family-level, layers averaged, sqrt percent (primary) or raw percent (robustness).

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(vegan)
library(nlme)

dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

# ---- 1. Station-year community matrix --------------------------------------
fam <- read_csv("data/processed/community_family_long.csv", show_col_types = FALSE) |>
  mutate(impact = factor(impact, levels = c("HI", "MI", "NI")))

station_year <- fam |>
  summarise(rel_abund = mean(rel_abund), .by = c(year, station, impact, family))

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

stopifnot(nrow(comm) == 18, all(rownames(comm) == meta$id))

# ---- 2. Distances and background -------------------------------------------
compute_distances <- function(transform) {
  x <- if (transform == "sqrt") sqrt(comm) else comm
  dm <- as.matrix(vegdist(x, method = "bray"))

  dist_long <- meta |>
    filter(impact != "NI") |>
    rowwise() |>
    mutate(
      distance = {
        ctrl <- meta$id[meta$year == year & meta$impact == "NI"]
        mean(dm[id, ctrl])
      }
    ) |>
    ungroup() |>
    mutate(transform = transform)

  background <- meta |>
    filter(impact == "NI") |>
    group_by(year) |>
    summarise(
      background = {
        ids <- id
        pairs <- combn(ids, 2)
        mean(apply(pairs, 2, function(p) dm[p[1], p[2]]))
      },
      .groups = "drop"
    ) |>
    mutate(transform = transform)

  dist_long <- dist_long |>
    left_join(background, by = c("year", "transform")) |>
    mutate(gap = distance - background)

  list(dist = dist_long, background = background)
}

dists <- lapply(c(sqrt = "sqrt", raw = "raw"), compute_distances)
dist_all <- bind_rows(lapply(dists, `[[`, "dist")) |>
  select(transform, station, impact, year, distance, background, gap)
bg_all <- bind_rows(lapply(dists, `[[`, "background"))

change <- dist_all |>
  select(transform, station, impact, year, distance, gap) |>
  pivot_wider(names_from = year, values_from = c(distance, gap)) |>
  mutate(change_distance = distance_2011 - distance_2010,
         change_gap = gap_2011 - gap_2010)

write_csv(dist_all, "outputs/tables/distance_to_control.csv")
write_csv(bg_all, "outputs/tables/control_background.csv")
write_csv(change, "outputs/tables/distance_change.csv")

# ---- 3. Mixed models -------------------------------------------------------
# year0 = 0 for 2010 and 1 for 2011, so slopes are change per year.
# Station is a random intercept. With two time points this is close to a paired
# comparison of each station's change, so power is very limited (6 stations).
fit_effects <- function(df, response, transform) {
  d <- df |> mutate(year0 = year - 2010, y = .data[[response]], impact = droplevels(impact))

  # Slope per impact level
  m_levels <- lme(y ~ 0 + impact + impact:year0, random = ~ 1 | station, data = d, method = "REML")
  ci <- intervals(m_levels, which = "fixed")$fixed
  tt <- summary(m_levels)$tTable
  slope_rows <- grep(":year0", rownames(ci))
  levels_out <- tibble(
    effect = paste("Change per year,", sub("impact(.*):year0", "\\1", rownames(ci)[slope_rows])),
    estimate = ci[slope_rows, "est."], lower = ci[slope_rows, "lower"], upper = ci[slope_rows, "upper"],
    p = tt[slope_rows, "p-value"]
  )

  # Pooled slope across HI and MI, adjusting for impact level
  m_pool <- lme(y ~ impact + year0, random = ~ 1 | station, data = d, method = "REML")
  cip <- intervals(m_pool, which = "fixed")$fixed
  ttp <- summary(m_pool)$tTable
  pooled <- tibble(
    effect = "Change per year, HI and MI pooled",
    estimate = cip["year0", "est."], lower = cip["year0", "lower"], upper = cip["year0", "upper"],
    p = ttp["year0", "p-value"]
  )

  # Is the change different between HI and MI?
  m_int <- lme(y ~ impact * year0, random = ~ 1 | station, data = d, method = "REML")
  tti <- summary(m_int)$tTable
  ixn_name <- grep(":year0", rownames(tti), value = TRUE)
  cii <- intervals(m_int, which = "fixed")$fixed
  diff_row <- tibble(
    effect = "Difference in change per year, MI minus HI",
    estimate = cii[ixn_name, "est."], lower = cii[ixn_name, "lower"], upper = cii[ixn_name, "upper"],
    p = tti[ixn_name, "p-value"]
  )

  bind_rows(levels_out, pooled, diff_row) |>
    mutate(response = response, transform = transform, n_obs = nrow(d), n_stations = n_distinct(d$station))
}

effects <- bind_rows(lapply(c("sqrt", "raw"), function(tr) {
  df <- dist_all |> filter(transform == tr)
  bind_rows(fit_effects(df, "distance", tr), fit_effects(df, "gap", tr))
})) |>
  select(transform, response, effect, estimate, lower, upper, p, n_obs, n_stations)

write_csv(effects, "outputs/tables/distance_model_effects.csv")

# ---- 4. Figure -------------------------------------------------------------
plot_df <- dist_all |> filter(transform == "sqrt")
bg_plot <- bg_all |> filter(transform == "sqrt")

dist_plot <- ggplot(plot_df, aes(year, distance)) +
  geom_line(aes(group = station, colour = impact), alpha = 0.7) +
  geom_point(aes(colour = impact), size = 2.5) +
  geom_line(data = bg_plot, aes(y = background, group = 1), linetype = "dashed") +
  geom_point(data = bg_plot, aes(y = background), shape = 18, size = 3) +
  annotate("text", x = 2011, y = bg_plot$background[bg_plot$year == 2011],
           label = "NI background", hjust = -0.1, size = 3) +
  scale_x_continuous(breaks = c(2010, 2011), expand = expansion(mult = c(0.1, 0.3))) +
  labs(
    x = "Year", y = "Mean Bray-Curtis distance to NI stations", colour = "Impact",
    subtitle = "Family-level, sqrt percent. Dashed line: mean distance among NI stations."
  ) +
  theme_bw()

ggsave("outputs/figures/distance_to_control.png", dist_plot, width = 6.5, height = 5, dpi = 300)

message("Done. Figure and tables written to outputs/.")
