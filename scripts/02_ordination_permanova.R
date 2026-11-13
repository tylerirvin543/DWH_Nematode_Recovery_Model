# ==============================================================================
# 02  ORDINATION AND PERMANOVA
#     Do communities differ by impact level and year, and did they converge?
# ==============================================================================
#
# QUESTION THIS ANSWERS
#   Do nematode communities differ among heavily impacted (HI), moderately
#   impacted (MI), and unimpacted (NI) stations, did they change between 2010 and
#   2011, and did the change depend on impact level?
#
# WHAT THIS SCRIPT DOES
#   1. Builds a matrix of family percentages: one row per station and year
#      (9 stations x 2 years = 18 rows), one column per family.
#   2. Measures how different every pair of rows is (Bray-Curtis dissimilarity).
#   3. Draws those dissimilarities as a two dimensional map (NMDS) with an
#      arrow from each station's 2010 position to its 2011 position.
#   4. Tests the patterns with PERMANOVA, a permutation test that works directly
#      on a dissimilarity matrix. Impact is constant within a station and year
#      changes within a station, so the two need different permutation schemes.
#   5. Checks the spread of the groups (betadisper), an assumption of PERMANOVA.
#   Everything is run twice: on square root percentages (main analysis) and on
#   untransformed percentages (robustness check).
#
# RUN
#   From the project root, after 01_data_prep.R.
#
# INPUT
#   data/processed/community_family_long.csv
#
# OUTPUTS
#   outputs/figures/nmds_trajectories.png        Figure 1 in the paper
#   outputs/tables/nmds_scores.csv               map coordinates of each station-year
#   outputs/tables/nmds_summary.csv              stress of each map
#   outputs/tables/station_change.csv            how far each station moved 2010 to 2011
#   outputs/tables/permanova_within_station.csv  year and impact x year tests
#   outputs/tables/permanova_impact_between.csv  impact level test
#   outputs/tables/betadisper.csv                spread check
#
# DESIGN CHOICES (and why)
#   * Unit of analysis = station-year. The two sediment layers at a station come
#     from the same spot and are not independent, so they are averaged first.
#     Treating all 36 samples as independent would overstate the evidence.
#   * Square root transform for the main analysis: stops a few very abundant
#     families from dominating the dissimilarities. The raw run shows whether the
#     conclusions depend on that choice.
#   * Impact can only be tested BETWEEN stations. Year and impact x year are
#     tested by shuffling year labels only WITHIN each station; impact is tested
#     by shuffling whole stations on their two year average.

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(vegan)

dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

# Fix the random number generator so the permutation p values and the NMDS
# starting configurations come out the same every time the script is run.
set.seed(4817)

# ---- 1. Build station-year community matrix --------------------------------
fam <- read_csv("data/processed/community_family_long.csv", show_col_types = FALSE) |>
  mutate(impact = factor(impact, levels = c("HI", "MI", "NI")))

# Average the two sediment layers: one value per family per station per year.
station_year <- fam |>
  summarise(rel_abund = mean(rel_abund), .by = c(year, station, impact, family))

# `meta` describes each row of the community matrix (which station, year, impact).
# The id (e.g. "D040S_2010") is used to keep matrix rows and metadata aligned.
meta <- station_year |>
  distinct(year, station, impact) |>
  arrange(year, impact, station) |>
  mutate(id = paste(station, year, sep = "_"), year_f = factor(year))

# `comm` is the community matrix: rows = station-years, columns = families,
# cells = average percent of that family. slice(match(...)) puts the rows in
# the same order as `meta`; column_to_rownames() turns the id column into row names.
comm <- station_year |>
  mutate(id = paste(station, year, sep = "_")) |>
  select(id, family, rel_abund) |>
  pivot_wider(names_from = family, values_from = rel_abund, values_fill = 0) |>
  slice(match(meta$id, id)) |>
  tibble::column_to_rownames("id") |>
  as.matrix()

# Checks: 18 rows, rows in the same order as `meta`, no missing values.
stopifnot(nrow(comm) == 18, all(rownames(comm) == meta$id), !anyNA(comm))

# ---- 2. Analysis function (run for each transformation) --------------------
analyse <- function(transform) {
  # Square root for the main analysis; raw percent for the robustness check.
  x <- if (transform == "sqrt") sqrt(comm) else comm

  # Pairwise Bray-Curtis dissimilarity between all 18 station-years.
  d <- vegdist(x, method = "bray")

  # NMDS: find a 2 dimensional layout whose distances keep the same rank order as
  # the dissimilarities. trymax = 100 random starts guards against a poor local
  # solution; autotransform = FALSE because we already chose the transformation.
  nmds <- metaMDS(x, distance = "bray", k = 2, trymax = 100, autotransform = FALSE, trace = 0)

  # PERMANOVA, part 1: year and impact x year.
  # how(blocks = station) means permutations only swap year labels WITHIN a
  # station, so each station's pair of observations stays together (999 shuffles).
  perm_within <- how(blocks = meta$station, nperm = 999)
  within <- adonis2(d ~ impact * year_f, data = meta, permutations = perm_within, by = "terms") |>
    as.data.frame() |>
    tibble::rownames_to_column("term") |>
    mutate(transform = transform, test = "within-station permutation")
  # The "impact" row from this call is not valid, because impact never varies
  # within a station and so cannot be tested by within-station shuffling. We drop
  # it and test impact separately below.
  within <- within |> filter(term != "impact")

  # PERMANOVA, part 2: impact level. Average each station's two years, then
  # shuffle whole stations among the three groups (9,999 shuffles). The division
  # by the table of counts turns the summed rows into means.
  stn_means <- rowsum(x, group = meta$station) / as.vector(table(meta$station)[rownames(rowsum(x, meta$station))])
  stn_meta <- meta |> distinct(station, impact) |> arrange(match(station, rownames(stn_means)))
  d_stn <- vegdist(stn_means, method = "bray")
  between <- adonis2(d_stn ~ impact, data = stn_meta, permutations = 9999) |>
    as.data.frame() |>
    tibble::rownames_to_column("term") |>
    mutate(transform = transform, test = "between-station (station means)")

  # Check of the PERMANOVA assumption that groups are similarly spread out.
  # PERMANOVA can mistake a difference in spread for a difference in location.
  # There are 6 groups (3 impact levels x 2 years) of 3 points each, so this check
  # has little power and can only catch large problems.
  grp <- interaction(meta$impact, meta$year_f)
  bd <- betadisper(d, grp)
  bd_test <- permutest(bd, permutations = 999)
  bd_row <- tibble(
    transform = transform,
    F = bd_test$tab[1, "F"],
    p = bd_test$tab[1, "Pr(>F)"]
  )

  # How far did each station move between years? The Bray-Curtis dissimilarity
  # between a station's 2010 and 2011 communities, looked up in the distance matrix.
  dm <- as.matrix(d)
  change <- meta |>
    filter(year == 2010) |>
    transmute(station, impact, transform,
              bc_2010_to_2011 = dm[cbind(match(paste(station, 2010, sep = "_"), meta$id),
                                         match(paste(station, 2011, sep = "_"), meta$id))])

  list(nmds = nmds, within = within, between = between, bd = bd_row, change = change)
}

# Run the whole analysis for both versions of the data.
res <- list(sqrt = analyse("sqrt"), raw = analyse("raw"))

# ---- 3. Tables -------------------------------------------------------------
# Stack the square root and raw results into one table each. The paper reports the
# square root rows; the raw rows are the robustness check.
write_csv(bind_rows(res$sqrt$within, res$raw$within), "outputs/tables/permanova_within_station.csv")
write_csv(bind_rows(res$sqrt$between, res$raw$between), "outputs/tables/permanova_impact_between.csv")
write_csv(bind_rows(res$sqrt$bd, res$raw$bd), "outputs/tables/betadisper.csv")
write_csv(bind_rows(res$sqrt$change, res$raw$change), "outputs/tables/station_change.csv")

# Stress measures how well the 2D map preserves the dissimilarities (lower is
# better; below about 0.2 is usually fine for reading the broad layout).
write_csv(
  tibble(transform = c("sqrt", "raw"),
         stress = c(res$sqrt$nmds$stress, res$raw$nmds$stress),
         converged = c(res$sqrt$nmds$converged, res$raw$nmds$converged)),
  "outputs/tables/nmds_summary.csv"
)

# Map coordinates for the square root analysis (the one we plot).
scores <- bind_cols(meta, as_tibble(scores(res$sqrt$nmds, display = "sites"))) |>
  select(year, station, impact, NMDS1, NMDS2)
write_csv(scores, "outputs/tables/nmds_scores.csv")

# ---- 4. NMDS figure with 2010 -> 2011 arrows -------------------------------
# One row per station with its 2010 and 2011 coordinates side by side, so that
# each station can be drawn as an arrow from (x, y) to (xend, yend).
arrows_df <- scores |>
  pivot_wider(names_from = year, values_from = c(NMDS1, NMDS2)) |>
  rename(x = NMDS1_2010, xend = NMDS1_2011, y = NMDS2_2010, yend = NMDS2_2011)

# Circles = 2010, triangles = 2011; colour = impact level. Station names are
# printed once, at the 2010 position, to avoid clutter.
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
