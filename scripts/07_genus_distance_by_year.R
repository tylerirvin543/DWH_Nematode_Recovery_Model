# ==============================================================================
# 07  GENUS LEVEL DISTANCE TO THE CONTROLS, EACH YEAR ON ITS OWN
#     Is the gap to the controls different for HI and MI within a year, and how do
#     the two years compare?
# ==============================================================================
#
# THE IDEA IN PLAIN TERMS
#   Same measures as script 03, now at genus level and computed within each year:
#   * distance   Bray-Curtis dissimilarity between an impacted sample and the NI
#                samples from the same year and the same sediment layer, averaged.
#   * background Average dissimilarity among the NI samples of that year and layer.
#   * gap        distance minus background. A recovered community has a gap near 0.
#   Each year is modelled separately with a linear mixed model. The two layers of a
#   station are the repeated measures, so station is a random intercept. Layer is
#   the unit here (not the station mean) so that the random effect can be estimated;
#   there are still only 6 impacted stations per year, so intervals are wide.
#   Then the years are compared: per-year estimates side by side, per-station change
#   in gap, and one pooled model with a year term (the per-year change in gap).
#
# RUN
#   From the project root, after 01_data_prep.R.
#
# INPUT
#   data/processed/community_genus_long.csv
#
# OUTPUTS
#   outputs/tables/genus_gap_samples.csv        distance, background and gap per sample
#   outputs/tables/genus_gap_stations.csv       gap per station and year (layers averaged)
#   outputs/tables/genus_gap_change.csv         change in gap 2010 to 2011 per station
#   outputs/tables/genus_gap_models_by_year.csv per-year mixed model estimates
#   outputs/tables/genus_gap_model_pooled.csv   pooled model, change per year
#   outputs/figures/genus_gap_by_year.png

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(vegan)
library(nlme)

dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

# ---- 1. Distances, background, gap per sample ------------------------------
gen <- read_csv("data/processed/community_genus_long.csv", show_col_types = FALSE) |>
  filter(!is.na(genus)) |>
  mutate(impact = factor(impact, levels = c("HI", "MI", "NI"))) |>
  mutate(pct = 100 * rel_abund / sum(rel_abund), .by = c(year, station, layer))

gap_for <- function(transform) {
  bind_rows(lapply(c(2010, 2011), \(yr) {
    d <- gen |> filter(year == yr)
    meta <- d |> distinct(station, impact, layer, sample) |> arrange(layer, impact, station)
    comm <- d |>
      select(sample, genus, pct) |>
      pivot_wider(names_from = genus, values_from = pct, values_fill = 0) |>
      slice(match(meta$sample, sample)) |>
      tibble::column_to_rownames("sample") |>
      as.matrix()
    stopifnot(nrow(comm) == 18, !anyNA(comm))
    x <- if (transform == "sqrt") sqrt(comm) else comm
    dm <- as.matrix(vegdist(x, method = "bray"))

    bg <- meta |> filter(impact == "NI") |>
      summarise(background = {
        p <- combn(sample, 2)
        mean(apply(p, 2, \(i) dm[i[1], i[2]]))
      }, .by = layer)

    meta |>
      filter(impact != "NI") |>
      rowwise() |>
      mutate(distance = mean(dm[sample, meta$sample[meta$impact == "NI" & meta$layer == layer]])) |>
      ungroup() |>
      left_join(bg, by = "layer") |>
      mutate(year = yr, transform = transform, gap = distance - background)
  }))
}

samples <- bind_rows(gap_for("sqrt"), gap_for("raw")) |>
  mutate(impact = droplevels(impact)) |>
  select(transform, year, station, impact, layer, distance, background, gap)
write_csv(samples, "outputs/tables/genus_gap_samples.csv")

stations <- samples |>
  summarise(distance = mean(distance), background = mean(background), gap = mean(gap),
            .by = c(transform, year, station, impact))
write_csv(stations, "outputs/tables/genus_gap_stations.csv")

change <- stations |>
  select(transform, station, impact, year, gap) |>
  pivot_wider(names_from = year, values_from = gap, names_prefix = "gap_") |>
  mutate(change_gap = gap_2011 - gap_2010)
write_csv(change, "outputs/tables/genus_gap_change.csv")

# ---- 2. Mixed model within each year ---------------------------------------
# gap ~ impact + layer, random intercept for station. The impact coefficient is
# the MI minus HI difference in gap; the intercept is the HI gap (top layer).
fit_year <- function(df, yr, transform) {
  d <- df |> filter(year == yr, transform == !!transform)
  m <- lme(gap ~ impact + layer, random = ~ 1 | station, data = d, method = "REML")
  ci <- intervals(m, which = "fixed")$fixed
  tt <- summary(m)$tTable
  tibble(
    year = yr, transform = transform,
    term = rownames(ci), estimate = ci[, "est."], lower = ci[, "lower"], upper = ci[, "upper"],
    p = tt[, "p-value"], n_obs = nrow(d), n_stations = n_distinct(d$station)
  )
}
by_year <- bind_rows(lapply(c("sqrt", "raw"), \(tr)
  bind_rows(fit_year(samples, 2010, tr), fit_year(samples, 2011, tr))))
write_csv(by_year, "outputs/tables/genus_gap_models_by_year.csv")

# ---- 3. Pooled model: change per year --------------------------------------
# year0 = 0 for 2010, 1 for 2011, so the year0 coefficient is the change in gap per
# year; impact:year0 asks whether HI and MI changed differently. Station is random.
fit_pooled <- function(transform) {
  d <- samples |> filter(transform == !!transform) |> mutate(year0 = year - 2010)
  m <- lme(gap ~ impact * year0 + layer, random = ~ 1 | station, data = d, method = "REML")
  ci <- intervals(m, which = "fixed")$fixed
  tt <- summary(m)$tTable
  tibble(
    transform = transform, term = rownames(ci),
    estimate = ci[, "est."], lower = ci[, "lower"], upper = ci[, "upper"],
    p = tt[, "p-value"], n_obs = nrow(d), n_stations = n_distinct(d$station)
  )
}
write_csv(bind_rows(fit_pooled("sqrt"), fit_pooled("raw")), "outputs/tables/genus_gap_model_pooled.csv")

# ---- 4. Figure -------------------------------------------------------------
ggsave("outputs/figures/genus_gap_by_year.png",
  ggplot(filter(stations, transform == "sqrt"), aes(year, gap)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_line(aes(group = station, colour = impact), alpha = 0.7) +
    geom_point(aes(colour = impact), size = 2.5) +
    scale_x_continuous(breaks = c(2010, 2011)) +
    labs(x = "Year", y = "Gap (distance to NI minus NI background)", colour = "Impact",
         subtitle = "Genus-level Bray-Curtis (sqrt percent), layers averaged for display") +
    theme_bw(),
  width = 6, height = 4.5, dpi = 300)

message("Done. Genus-level gap tables and figure written to outputs/.")
