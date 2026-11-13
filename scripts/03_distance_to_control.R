# ==============================================================================
# 03  DISTANCE TO THE CONTROLS
#     Did impacted stations move closer to the unimpacted stations?
# ==============================================================================
#
# QUESTION THIS ANSWERS
#   Our main question: between 2010 and 2011, did communities at impacted (HI and
#   MI) stations become more similar to communities at unimpacted (NI) stations?
#   This follows the approach our instructor suggested.
#
# THE IDEA IN PLAIN TERMS
#   * distance   For each impacted station and year: the average Bray-Curtis
#                dissimilarity to the three NI stations sampled that year.
#                Smaller distance = more similar to the controls.
#   * background For each year: the average dissimilarity among the three pairs of
#                NI stations. Healthy sites are never identical, so this is "normal
#                variation". A recovered station should sit about one background
#                away from the controls, not at zero.
#   * gap        distance minus that year's background. A recovering station's gap
#                shrinks toward zero. We model the gap (not the raw distance)
#                because the background itself changes a little between years, and
#                we do not want that to be mistaken for movement by impacted stations.
#   Then a linear mixed model estimates the change in gap per year for HI and MI.
#
# RUN
#   From the project root, after 01_data_prep.R.
#
# INPUT
#   data/processed/community_family_long.csv
#
# OUTPUTS
#   outputs/figures/distance_to_control.png      Figure 2 in the paper
#   outputs/tables/distance_to_control.csv       distance and gap per impacted station-year
#   outputs/tables/control_background.csv        background per year
#   outputs/tables/distance_change.csv           change 2010 to 2011 per station
#   outputs/tables/distance_model_effects.csv    mixed model estimates (Table 2 in the paper)
#
# DESIGN CHOICES (and why)
#   * Same data preparation as script 02: family level, the two layers averaged,
#     square root percent for the main analysis and raw percent as a check.
#   * Station is a random effect. Each station contributes two points (2010, 2011);
#     without a random effect a station that is simply far from the controls in
#     both years would be counted as two independent pieces of evidence.
#   * With 6 impacted stations x 2 years = 12 observations, power is very low.
#     The model behaves much like a paired comparison, so intervals are wide and
#     the results should be treated as rough.

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(vegan)
library(nlme)

dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

# ---- 1. Station-year community matrix --------------------------------------
# Identical to script 02: average the two layers, then build a matrix with one row
# per station-year and one column per family.
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
  # Full 18 x 18 matrix of Bray-Curtis dissimilarities, so any pair can be looked
  # up by row and column name (e.g. dm["D040S_2010", "D019S_2010"]).
  dm <- as.matrix(vegdist(x, method = "bray"))

  # DISTANCE: for each HI/MI station-year, average its dissimilarity to the three
  # NI stations sampled in the same year. Inside rowwise(), `id` and `year` are the
  # current row's values, while `meta$year` and `meta$impact` refer to the whole table.
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

  # BACKGROUND: for each year, average dissimilarity over all pairs of the three
  # NI stations (combn(ids, 2) lists the 3 possible pairs).
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

  # GAP: distance minus the background for the same year.
  dist_long <- dist_long |>
    left_join(background, by = c("year", "transform")) |>
    mutate(gap = distance - background)

  list(dist = dist_long, background = background)
}

# Run for both versions of the data and stack the results.
dists <- lapply(c(sqrt = "sqrt", raw = "raw"), compute_distances)
dist_all <- bind_rows(lapply(dists, `[[`, "dist")) |>
  select(transform, station, impact, year, distance, background, gap)
bg_all <- bind_rows(lapply(dists, `[[`, "background"))

# Per-station change between years: 2011 value minus 2010 value (positive = further
# from the controls in 2011).
change <- dist_all |>
  select(transform, station, impact, year, distance, gap) |>
  pivot_wider(names_from = year, values_from = c(distance, gap)) |>
  mutate(change_distance = distance_2011 - distance_2010,
         change_gap = gap_2011 - gap_2010)

write_csv(dist_all, "outputs/tables/distance_to_control.csv")
write_csv(bg_all, "outputs/tables/control_background.csv")
write_csv(change, "outputs/tables/distance_change.csv")

# ---- 3. Mixed models -------------------------------------------------------
# year0 is 0 for 2010 and 1 for 2011, so a "slope" on year0 is the change per year.
# Each station gets its own random intercept (its own starting level).
# Models are fitted by REML (restricted maximum likelihood), the usual default for
# mixed models, which gives less biased variance estimates than ordinary ML.
fit_effects <- function(df, response, transform) {
  # `response` is either "distance" or "gap"; .data[[response]] picks that column.
  # droplevels() removes the (now empty) NI level.
  d <- df |> mutate(year0 = year - 2010, y = .data[[response]], impact = droplevels(impact))

  # Model 1: a separate starting level and a separate slope for HI and for MI.
  # "0 + impact + impact:year0" gives one intercept and one slope per group.
  # intervals() returns 95% confidence intervals; tTable gives p values.
  m_levels <- lme(y ~ 0 + impact + impact:year0, random = ~ 1 | station, data = d, method = "REML")
  ci <- intervals(m_levels, which = "fixed")$fixed
  tt <- summary(m_levels)$tTable
  slope_rows <- grep(":year0", rownames(ci))
  levels_out <- tibble(
    effect = paste("Change per year,", sub("impact(.*):year0", "\\1", rownames(ci)[slope_rows])),
    estimate = ci[slope_rows, "est."], lower = ci[slope_rows, "lower"], upper = ci[slope_rows, "upper"],
    p = tt[slope_rows, "p-value"]
  )

  # Model 2: one slope shared by HI and MI (a bit more precise), with a separate
  # starting level for each group.
  m_pool <- lme(y ~ impact + year0, random = ~ 1 | station, data = d, method = "REML")
  cip <- intervals(m_pool, which = "fixed")$fixed
  ttp <- summary(m_pool)$tTable
  pooled <- tibble(
    effect = "Change per year, HI and MI pooled",
    estimate = cip["year0", "est."], lower = cip["year0", "lower"], upper = cip["year0", "upper"],
    p = ttp["year0", "p-value"]
  )

  # Model 3: is the slope different between HI and MI? (interaction term)
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

# Fit everything for the two data versions and the two responses (distance, gap).
# The paper reports the gap results.
effects <- bind_rows(lapply(c("sqrt", "raw"), function(tr) {
  df <- dist_all |> filter(transform == tr)
  bind_rows(fit_effects(df, "distance", tr), fit_effects(df, "gap", tr))
})) |>
  select(transform, response, effect, estimate, lower, upper, p, n_obs, n_stations)

write_csv(effects, "outputs/tables/distance_model_effects.csv")

# ---- 4. Figure -------------------------------------------------------------
# One line per impacted station (distance in 2010 and 2011). The dashed line is
# the NI background: a station near that line is about as similar to the controls
# as the controls are to each other.
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
