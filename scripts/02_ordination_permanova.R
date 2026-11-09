# 02_ordination_permanova.R
# Bray-Curtis ordination (NMDS) and PERMANOVA of family-level nematode composition.
#
# Run from the project root, after scripts/01_data_prep.R.
#
# Inputs : data/processed/community_family_long.csv
# Outputs: outputs/figures/nmds_trajectories.png
#          outputs/tables/nmds_scores.csv, nmds_summary.csv, station_change.csv
#          outputs/tables/permanova_within_station.csv, permanova_impact_between.csv, betadisper.csv
#
# Design notes
# - Unit of analysis is station-year (9 stations x 2 years = 18 points). The two
#   sediment layers per station are averaged first, because they are not independent.
# - Primary analysis uses square-root transformed percent composition (down-weights
#   very abundant families). Untransformed percent is run as a robustness check.
# - Impact is constant within a station, so it can only be tested between stations.
#   Year and impact x year are tested with permutations restricted within station;
#   impact is tested separately on station means (permuting whole stations).

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(vegan)

dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

set.seed(4817)

# ---- 1. Build station-year community matrix --------------------------------
fam <- read_csv("data/processed/community_family_long.csv", show_col_types = FALSE) |>
  mutate(impact = factor(impact, levels = c("HI", "MI", "NI")))

station_year <- fam |>
  summarise(rel_abund = mean(rel_abund), .by = c(year, station, impact, family))

meta <- station_year |>
  distinct(year, station, impact) |>
  arrange(year, impact, station) |>
  mutate(id = paste(station, year, sep = "_"), year_f = factor(year))

comm <- station_year |>
  mutate(id = paste(station, year, sep = "_")) |>
  select(id, family, rel_abund) |>
  pivot_wider(names_from = family, values_from = rel_abund, values_fill = 0) |>
  slice(match(meta$id, id)) |>
  tibble::column_to_rownames("id") |>
  as.matrix()

stopifnot(nrow(comm) == 18, all(rownames(comm) == meta$id), !anyNA(comm))

# ---- 2. Analysis function (run for each transformation) --------------------
analyse <- function(transform) {
  x <- if (transform == "sqrt") sqrt(comm) else comm
  d <- vegdist(x, method = "bray")

  # NMDS
  nmds <- metaMDS(x, distance = "bray", k = 2, trymax = 100, autotransform = FALSE, trace = 0)

  # Within-station permutations: tests year and impact x year
  perm_within <- how(blocks = meta$station, nperm = 999)
  within <- adonis2(d ~ impact * year_f, data = meta, permutations = perm_within, by = "terms") |>
    as.data.frame() |>
    tibble::rownames_to_column("term") |>
    mutate(transform = transform, test = "within-station permutation")
  # Impact main effect is not valid here (impact does not vary within station)
  within <- within |> filter(term != "impact")

  # Between-station test of impact on station means across the two years
  stn_means <- rowsum(x, group = meta$station) / as.vector(table(meta$station)[rownames(rowsum(x, meta$station))])
  stn_meta <- meta |> distinct(station, impact) |> arrange(match(station, rownames(stn_means)))
  d_stn <- vegdist(stn_means, method = "bray")
  between <- adonis2(d_stn ~ impact, data = stn_meta, permutations = 9999) |>
    as.data.frame() |>
    tibble::rownames_to_column("term") |>
    mutate(transform = transform, test = "between-station (station means)")

  # Homogeneity of dispersion (PERMANOVA assumption), 6 groups of 3 points
  grp <- interaction(meta$impact, meta$year_f)
  bd <- betadisper(d, grp)
  bd_test <- permutest(bd, permutations = 999)
  bd_row <- tibble(
    transform = transform,
    F = bd_test$tab[1, "F"],
    p = bd_test$tab[1, "Pr(>F)"]
  )

  # Per-station change between years (Bray-Curtis distance 2010 -> 2011)
  dm <- as.matrix(d)
  change <- meta |>
    filter(year == 2010) |>
    transmute(station, impact, transform,
              bc_2010_to_2011 = dm[cbind(match(paste(station, 2010, sep = "_"), meta$id),
                                         match(paste(station, 2011, sep = "_"), meta$id))])

  list(nmds = nmds, within = within, between = between, bd = bd_row, change = change)
}

res <- list(sqrt = analyse("sqrt"), raw = analyse("raw"))

# ---- 3. Tables -------------------------------------------------------------
write_csv(bind_rows(res$sqrt$within, res$raw$within), "outputs/tables/permanova_within_station.csv")
write_csv(bind_rows(res$sqrt$between, res$raw$between), "outputs/tables/permanova_impact_between.csv")
write_csv(bind_rows(res$sqrt$bd, res$raw$bd), "outputs/tables/betadisper.csv")
write_csv(bind_rows(res$sqrt$change, res$raw$change), "outputs/tables/station_change.csv")

write_csv(
  tibble(transform = c("sqrt", "raw"),
         stress = c(res$sqrt$nmds$stress, res$raw$nmds$stress),
         converged = c(res$sqrt$nmds$converged, res$raw$nmds$converged)),
  "outputs/tables/nmds_summary.csv"
)

scores <- bind_cols(meta, as_tibble(scores(res$sqrt$nmds, display = "sites"))) |>
  select(year, station, impact, NMDS1, NMDS2)
write_csv(scores, "outputs/tables/nmds_scores.csv")

# ---- 4. NMDS figure with 2010 -> 2011 arrows -------------------------------
arrows_df <- scores |>
  pivot_wider(names_from = year, values_from = c(NMDS1, NMDS2)) |>
  rename(x = NMDS1_2010, xend = NMDS1_2011, y = NMDS2_2010, yend = NMDS2_2011)

nmds_plot <- ggplot(scores, aes(NMDS1, NMDS2)) +
  geom_segment(
    data = arrows_df,
    aes(x = x, y = y, xend = xend, yend = yend, colour = impact),
    arrow = arrow(length = unit(0.18, "cm")), linewidth = 0.6, alpha = 0.7
  ) +
  geom_point(aes(colour = impact, shape = factor(year)), size = 3) +
  geom_text(
    data = filter(scores, year == 2010), aes(label = station),
    size = 2.8, vjust = -1, show.legend = FALSE
  ) +
  labs(
    x = "NMDS1", y = "NMDS2",
    colour = "Impact", shape = "Year",
    subtitle = sprintf("Family-level Bray-Curtis (sqrt percent), stress = %.3f. Arrows run 2010 to 2011.",
                       res$sqrt$nmds$stress)
  ) +
  theme_bw()

ggsave("outputs/figures/nmds_trajectories.png", nmds_plot, width = 7, height = 5.5, dpi = 300)

message("Done. Figure and tables written to outputs/.")
