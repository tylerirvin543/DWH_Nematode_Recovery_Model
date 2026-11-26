# ==============================================================================
# 06  GENUS LEVEL, EACH YEAR ON ITS OWN: ENVIRONMENT PCA AND PERMANOVA
#     Do genus communities differ by impact level within a year, and does that
#     difference line up with the environmental gradient?
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
#   For 2010 and 2011 separately (9 stations each):
#   1. Builds a genus x station matrix (two layers averaged) and the Bray-Curtis
#      dissimilarity between stations.
#   2. Runs a PCA on the standardized environmental variables (Al, Ba, PAH44, TPH).
#      Bray-Curtis is for community data; environmental variables are continuous
#      with different units, so they are reduced with PCA instead. A PCoA of the
#      Bray-Curtis matrix (a PCA on the dissimilarities) is drawn with the four
#      variables fitted as arrows (envfit), to put the two together.
#   3. PERMANOVA on the Bray-Curtis matrix: impact alone, then impact plus the
#      first two environmental PCs (marginal tests, so term order does not matter).
#   4. Checks the spread of groups (betadisper) and draws an NMDS per year.
#   Then, as the formal test of change between years, repeats the pooled PERMANOVA
#   (year and impact x year, shuffling years within stations) at genus level.
#
# RUN
#   From the project root, after 01_data_prep.R.
#
# INPUT
#   data/processed/community_genus_long.csv, data/processed/env_stations.csv
#
# OUTPUTS
#   outputs/tables/genus_permanova_by_year.csv   impact and impact + env PCs, per year
#   outputs/tables/genus_betadisper_by_year.csv  spread check
#   outputs/tables/env_pca_variance.csv          variance explained by each PC
#   outputs/tables/env_pca_loadings.csv          loadings per year
#   outputs/tables/env_pca_scores.csv            PC scores per station and year
#   outputs/tables/genus_permanova_pooled.csv    year and impact x year, impact (genus level)
#   outputs/tables/genus_nmds_stress.csv         NMDS stress per year
#   outputs/figures/genus_nmds_by_year.png
#   outputs/figures/env_pca_by_year.png
#
# DESIGN CHOICES (and why)
#   * Each year is analysed on its own, so the 2010 and 2011 results use only that
#     year's genera and that year's controls. No genus harmonizing is needed here.
#   * Animals identified only to family (no genus) are set aside and each layer is
#     rescaled to 100 percent again.
#   * Ba, PAH44 and TPH are log transformed before standardizing because they span
#     orders of magnitude. The environmental values were measured differently in the
#     two years (for example Al and PAH44 differ a lot at the control stations),
#     which is another reason the PCA is done within a year.
#   * With 9 stations, impact plus two PCs leaves 4 residual degrees of freedom. The
#     environmental tests are exploratory. Environment is also strongly tied to
#     impact level, so the two cannot be cleanly separated.

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(vegan)

dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

set.seed(6139)

# ---- 1. Genus data ---------------------------------------------------------
gen <- read_csv("data/processed/community_genus_long.csv", show_col_types = FALSE) |>
  filter(!is.na(genus)) |>
  mutate(impact = factor(impact, levels = c("HI", "MI", "NI"))) |>
  # Rescale each layer to 100 percent after dropping family-only animals
  mutate(pct = 100 * rel_abund / sum(rel_abund), .by = c(year, station, layer))

# Station-year means of the two layers
station_year <- gen |>
  summarise(pct = mean(pct), .by = c(year, station, impact, genus))

env <- read_csv("data/processed/env_stations.csv", show_col_types = FALSE) |>
  select(year, station, Al, Ba, PAH44, TPH)

build_matrix <- function(df) {
  meta <- df |>
    distinct(year, station, impact) |>
    arrange(impact, station)
  comm <- df |>
    select(station, genus, pct) |>
    pivot_wider(names_from = genus, values_from = pct, values_fill = 0) |>
    slice(match(meta$station, station)) |>
    tibble::column_to_rownames("station") |>
    as.matrix()
  # Drop genera absent from every station that year
  comm <- comm[, colSums(comm) > 0, drop = FALSE]
  stopifnot(nrow(comm) == 9, all(rownames(comm) == meta$station), !anyNA(comm))
  list(meta = meta, comm = comm)
}

# ---- 2. Per-year analysis --------------------------------------------------
perm_tab <- function(a, year, model, transform) {
  as.data.frame(a) |>
    tibble::rownames_to_column("term") |>
    mutate(year = year, model = model, transform = transform)
}

analyse_year <- function(yr, transform) {
  m <- build_matrix(filter(station_year, year == yr))
  x <- if (transform == "sqrt") sqrt(m$comm) else m$comm
  d <- vegdist(x, method = "bray")

  # Environmental PCA: log the skewed variables, then standardize
  e <- env |>
    filter(year == yr) |>
    slice(match(m$meta$station, station)) |>
    mutate(across(c(Ba, PAH44, TPH), log10))
  pca <- prcomp(e[, c("Al", "Ba", "PAH44", "TPH")], scale. = TRUE)
  m$meta$PC1 <- pca$x[, 1]
  m$meta$PC2 <- pca$x[, 2]

  # Impact only (whole stations are the units, so ordinary permutations are valid)
  a1 <- adonis2(d ~ impact, data = m$meta, permutations = 9999)
  # Impact plus environmental PCs, marginal tests
  a2 <- adonis2(d ~ impact + PC1 + PC2, data = m$meta, permutations = 9999, by = "margin")

  bd <- permutest(betadisper(d, m$meta$impact), permutations = 999)

  list(
    perm = bind_rows(perm_tab(a1, yr, "impact", transform),
                     perm_tab(a2, yr, "impact + env PCs", transform)),
    bd = tibble(year = yr, transform = transform,
                F = bd$tab[1, "F"], p = bd$tab[1, "Pr(>F)"]),
    pca = pca, env = e, meta = m$meta, x = x, d = d
  )
}

res <- expand.grid(yr = c(2010, 2011), transform = c("sqrt", "raw"), stringsAsFactors = FALSE) |>
  purrr::pmap(\(yr, transform) analyse_year(yr, transform))
names(res) <- c("2010_sqrt", "2011_sqrt", "2010_raw", "2011_raw")

write_csv(bind_rows(lapply(res, `[[`, "perm")), "outputs/tables/genus_permanova_by_year.csv")
write_csv(bind_rows(lapply(res, `[[`, "bd")), "outputs/tables/genus_betadisper_by_year.csv")

# PCA outputs come from the sqrt runs (the PCA itself does not depend on transform)
pca_res <- res[c("2010_sqrt", "2011_sqrt")]
write_csv(
  bind_rows(lapply(names(pca_res), \(n) {
    p <- pca_res[[n]]$pca
    tibble(year = as.numeric(substr(n, 1, 4)), PC = paste0("PC", seq_along(p$sdev)),
           variance = p$sdev^2 / sum(p$sdev^2))
  })),
  "outputs/tables/env_pca_variance.csv"
)
write_csv(
  bind_rows(lapply(names(pca_res), \(n) {
    p <- pca_res[[n]]$pca
    as_tibble(p$rotation[, 1:2], rownames = "variable") |>
      mutate(year = as.numeric(substr(n, 1, 4)))
  })),
  "outputs/tables/env_pca_loadings.csv"
)
scores <- bind_rows(lapply(pca_res, \(r) r$meta |> select(year, station, impact, PC1, PC2)))
write_csv(scores, "outputs/tables/env_pca_scores.csv")

# ---- 3. Figures ------------------------------------------------------------
nmds_pts <- bind_rows(lapply(c("2010_sqrt", "2011_sqrt"), \(n) {
  r <- res[[n]]
  nm <- metaMDS(r$x, distance = "bray", k = 2, trymax = 100, autotransform = FALSE, trace = 0)
  as_tibble(scores(nm, display = "sites")) |>
    bind_cols(r$meta |> select(year, station, impact)) |>
    mutate(stress = nm$stress)
}))

write_csv(distinct(nmds_pts, year, stress), "outputs/tables/genus_nmds_stress.csv")

ggsave("outputs/figures/genus_nmds_by_year.png",
  ggplot(nmds_pts, aes(NMDS1, NMDS2, colour = impact)) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = station), size = 2.8, max.overlaps = Inf, min.segment.length = 0.2, show.legend = FALSE) +
    scale_x_continuous(expand = expansion(mult = 0.15)) +
    scale_y_continuous(expand = expansion(mult = 0.15)) +
    facet_wrap(~ year, scales = "free") +
    labs(colour = "Impact", subtitle = "Genus-level Bray-Curtis (sqrt percent), each year on its own") +
    theme_bw(),
  width = 8, height = 4.5, dpi = 300)

ggsave("outputs/figures/env_pca_by_year.png",
  ggplot(scores, aes(PC1, PC2, colour = impact)) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = station), size = 2.8, max.overlaps = Inf, min.segment.length = 0.2, show.legend = FALSE) +
    scale_x_continuous(expand = expansion(mult = 0.15)) +
    scale_y_continuous(expand = expansion(mult = 0.15)) +
    facet_wrap(~ year, scales = "free") +
    labs(colour = "Impact", subtitle = "PCA of Al, log Ba, log PAH44, log TPH (standardized)") +
    theme_bw(),
  width = 8, height = 4.5, dpi = 300)

# PCoA of the Bray-Curtis matrix with the environmental variables fitted as arrows.
# envfit regresses each variable on the first two axes; arrow length is scaled by
# sqrt(r2) so longer arrows mean a better fit. p values use 9999 permutations.
pcoa_out <- lapply(c(2010, 2011), \(yr) {
  r <- res[[paste0(yr, "_sqrt")]]
  pc <- cmdscale(r$d, k = 2, eig = TRUE)
  ef <- envfit(pc$points, r$env[, c("Al", "Ba", "PAH44", "TPH")], permutations = 9999)
  ax <- round(100 * pc$eig[1:2] / sum(pmax(pc$eig, 0)), 1)
  list(
    pts = tibble(year = yr, station = r$meta$station, impact = r$meta$impact,
                 A1 = pc$points[, 1], A2 = pc$points[, 2]),
    arrows = tibble(year = yr, variable = rownames(ef$vectors$arrows),
                    A1 = ef$vectors$arrows[, 1] * sqrt(ef$vectors$r),
                    A2 = ef$vectors$arrows[, 2] * sqrt(ef$vectors$r),
                    r2 = ef$vectors$r, p = ef$vectors$pvals),
    axes = tibble(year = yr, axis = c("PCoA1", "PCoA2"), pct = ax)
  )
})
pcoa_pts <- bind_rows(lapply(pcoa_out, `[[`, "pts"))
pcoa_arr <- bind_rows(lapply(pcoa_out, `[[`, "arrows"))
write_csv(pcoa_arr, "outputs/tables/genus_pcoa_envfit.csv")
write_csv(bind_rows(lapply(pcoa_out, `[[`, "axes")), "outputs/tables/genus_pcoa_axes.csv")

ggsave("outputs/figures/genus_pcoa_env.png",
  ggplot(pcoa_pts, aes(A1, A2)) +
    geom_point(aes(colour = impact), size = 3) +
    geom_text_repel(aes(label = station, colour = impact), size = 2.8, max.overlaps = Inf, min.segment.length = 0.2, show.legend = FALSE) +
    geom_segment(data = pcoa_arr, aes(x = 0, y = 0, xend = A1, yend = A2),
                 arrow = arrow(length = unit(0.18, "cm"))) +
    geom_text(data = pcoa_arr |> mutate(hj = c(Al = -0.15, Ba = 1.15, PAH44 = 1.15, TPH = -0.15)[variable],
                                         vj = c(Al = 1.4, Ba = -0.4, PAH44 = 1.4, TPH = -0.4)[variable]),
              aes(A1, A2, label = variable, hjust = hj, vjust = vj), size = 3, fontface = "bold") +
    scale_x_continuous(expand = expansion(mult = 0.2)) +
    scale_y_continuous(expand = expansion(mult = 0.2)) +
    facet_wrap(~ year, scales = "free") +
    labs(x = "PCoA1", y = "PCoA2", colour = "Impact",
         subtitle = "PCoA of genus Bray-Curtis (sqrt percent) with Al, log Ba, log PAH44, log TPH fitted as arrows") +
    theme_bw(),
  width = 8.5, height = 4.8, dpi = 300)

# ---- 4. Pooled genus PERMANOVA: the formal test of change between years -----
# Same scheme as script 02: shuffle years within stations for year and
# impact x year, shuffle whole stations (two year means) for impact. Genus lists
# differ between years, so we run it on all genera and on genera found in both years.
pooled <- function(genera_set, label) {
  sy <- station_year |> filter(genus %in% genera_set)
  meta <- sy |> distinct(year, station, impact) |> arrange(year, impact, station) |>
    mutate(id = paste(station, year, sep = "_"), year_f = factor(year))
  comm <- sy |>
    mutate(id = paste(station, year, sep = "_")) |>
    select(id, genus, pct) |>
    pivot_wider(names_from = genus, values_from = pct, values_fill = 0) |>
    slice(match(meta$id, id)) |>
    tibble::column_to_rownames("id") |>
    as.matrix()
  x <- sqrt(comm)
  d <- vegdist(x, method = "bray")
  w <- adonis2(d ~ impact * year_f, data = meta,
               permutations = how(blocks = meta$station, nperm = 999), by = "terms") |>
    as.data.frame() |> tibble::rownames_to_column("term") |>
    filter(term != "impact") |> mutate(test = "within-station permutation")
  sm <- rowsum(x, meta$station) / 2
  sm_meta <- meta |> distinct(station, impact) |> arrange(match(station, rownames(sm)))
  b <- adonis2(vegdist(sm, "bray") ~ impact, data = sm_meta, permutations = 9999) |>
    as.data.frame() |> tibble::rownames_to_column("term") |>
    mutate(test = "between-station (station means)")
  bind_rows(w, b) |> mutate(genera = label, n_genera = ncol(comm))
}

all_g <- unique(station_year$genus)
shared_g <- intersect(unique(station_year$genus[station_year$year == 2010]),
                      unique(station_year$genus[station_year$year == 2011]))
write_csv(bind_rows(pooled(all_g, "all genera"), pooled(shared_g, "genera in both years")),
          "outputs/tables/genus_permanova_pooled.csv")

message("Done. Genus-level per-year tables and figures written to outputs/.")
