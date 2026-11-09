# 01_data_prep.R
# Tidy the 2010 and 2011 nematode workbook into long-format tables.
#
# Run from the project root (open "Group Project.Rproj", then source this file).
#
# Inputs : data/raw/Clean data_Recovery in 1, 2, 3....xlsx
# Outputs: data/processed/env_stations.csv        one row per station-year
#          data/processed/samples.csv             one row per sample (station x layer x year) with MI, ITD
#          data/processed/community_genus_long.csv   taxon x sample, finest resolution available
#          data/processed/community_family_long.csv  same, summed to family (used for cross-year comparison)
#          outputs/tables/data_cleaning_log.csv   every cleaning decision / flag, in plain language

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(readr)

raw_file <- "data/raw/Clean data_Recovery in 1, 2, 3....xlsx"
stopifnot(file.exists(raw_file))

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

# ---- Cleaning log -----------------------------------------------------------
cleaning_log <- tibble(step = character(), issue = character(), action = character())
log_issue <- function(step, issue, action) {
  cleaning_log <<- bind_rows(cleaning_log, tibble(step = step, issue = issue, action = action))
}

# ---- 1. Environmental variables --------------------------------------------
read_env <- function(sheet) {
  read_excel(raw_file, sheet = sheet) |>
    rename(year = Year, station = Station, impact = Impact, depth_m = Depth,
           lat = `Lat (N)`, long = `Long (W)`)
}

env <- bind_rows(read_env("Env variables 2010"), read_env("Env variables 2011")) |>
  mutate(impact = factor(impact, levels = c("HI", "MI", "NI")))

stopifnot(nrow(env) == 18, !anyNA(env))

# Treated as comparable for now (team decision), but the shift is large
env_ratio <- env |>
  pivot_longer(c(Al, Ba, PAH44, TPH, Mud, TOC), names_to = "variable") |>
  group_by(variable, year) |>
  summarise(med = median(value), .groups = "drop") |>
  pivot_wider(names_from = year, values_from = med) |>
  mutate(ratio = `2011` / `2010`)

big <- env_ratio |> filter(ratio > 2 | ratio < 0.5)
if (nrow(big) > 0) {
  log_issue(
    "env",
    paste0("Median 2011/2010 ratio is large for: ",
           paste0(big$variable, " (", round(big$ratio, 2), "x)", collapse = ", "),
           ". This could reflect a change in units or laboratory methods rather than a real shift."),
    "Kept as given and treated as comparable for now. Units and methods should be confirmed before comparing environmental variables across years."
  )
}

# ---- 2. Community matrices --------------------------------------------------
# Each matrix has taxon rows, then two summary rows (MI, ITD) in the last rows.
# 2010 stores percent composition; 2011 stores counts from ~100 individuals per sample.

read_matrix <- function(sheet, year) {
  m <- read_excel(raw_file, sheet = sheet, .name_repair = "minimal") |>
    rename(family = Family, genus = Genus, trophic = `Trophic gr.`, cp = `c-p value`)

  sample_cols <- setdiff(names(m), c("family", "genus", "trophic", "cp"))

  indices <- m |>
    filter(is.na(genus), cp %in% c("MI", "ITD")) |>
    select(index = cp, all_of(sample_cols)) |>
    pivot_longer(-index, names_to = "sample", values_to = "value") |>
    pivot_wider(names_from = index, values_from = value) |>
    rename(MI_given = MI, ITD_given = ITD)

  taxa <- m |>
    filter(!is.na(genus)) |>
    mutate(across(all_of(sample_cols), as.numeric)) |>
    pivot_longer(all_of(sample_cols), names_to = "sample", values_to = "value") |>
    mutate(year = year)

  list(taxa = taxa, indices = indices |> mutate(year = year))
}

m10 <- read_matrix("Matrix 2010", 2010)
m11 <- read_matrix("Matrix 2011", 2011)

taxa <- bind_rows(m10$taxa, m11$taxa)
indices <- bind_rows(m10$indices, m11$indices)

# ---- 3. Clean taxon attributes ---------------------------------------------
taxa <- taxa |>
  mutate(
    family  = str_squish(family),
    genus   = str_squish(genus),
    trophic = str_to_upper(str_squish(trophic)),
    cp      = as.integer(cp)
  )

log_issue("taxa", "Trophic group coded 1A/1B/2A/2B in 2010 but 1a/1b/2a/2b in 2011; c-p value stored as text in 2011; trailing spaces in Family names.",
          "Set trophic groups to upper case, converted the c-p value to a number, trimmed whitespace.")

# Known typo: Sabatieria (2010) vs Sabatiera (2011)
if (any(taxa$genus == "Sabatiera")) {
  taxa <- taxa |> mutate(genus = if_else(genus == "Sabatiera", "Sabatieria", genus))
  log_issue("taxa", "'Sabatiera' appears only in 2011 and 'Sabatieria' only in 2010, which looks like a spelling variant.",
            "Renamed to Sabatieria. This is an assumption that should be confirmed.")
}

# Family names sitting in the Genus column (identified only to family)
family_names <- unique(taxa$family)
fam_in_genus <- taxa |>
  filter(genus %in% family_names) |>
  distinct(year, family, genus)

if (nrow(fam_in_genus) > 0) {
  log_issue("taxa",
            paste0("Family names used in the Genus column (not resolved to genus): ",
                   paste(unique(fam_in_genus$genus), collapse = ", "), "."),
            "Genus set to NA ('unresolved') and family kept. Comparisons between years use family level.")
  taxa <- taxa |> mutate(genus = if_else(genus %in% family_names, NA_character_, genus))
}

# Taxon records with missing family would break the family-level table
if (anyNA(taxa$family)) {
  log_issue("taxa", "Some taxa have no family.", "Dropped from the family table and flagged for manual checking.")
}

# Missing cell(s) inside the 2011 matrix
na_cells <- taxa |> filter(is.na(value)) |> distinct(year, genus, sample)
if (nrow(na_cells) > 0) {
  log_issue("taxa",
            paste0("Blank abundance cell(s): ",
                   paste0(na_cells$genus, " in ", na_cells$sample, " (", na_cells$year, ")", collapse = "; "), "."),
            "Treated as 0 (taxon not observed).")
  taxa <- taxa |> mutate(value = replace_na(value, 0))
}

# ---- 4. Sample metadata and relative abundance -----------------------------
# Layer codes verified against the 2011 raw sheet: _1 = 0-1 cm, _2 = 1-3 cm.
taxa <- taxa |>
  mutate(
    station = str_remove(sample, "_[12]$"),
    layer   = if_else(str_detect(sample, "_1$"), "0-1 cm", "1-3 cm")
  ) |>
  left_join(env |> select(year, station, impact), by = c("year", "station"))

stopifnot(!anyNA(taxa$impact))

# Verify 2011 matrix against the raw record-level sheet
raw11 <- read_excel(raw_file, sheet = "2011 raw") |>
  transmute(station = Station,
            layer = if_else(`Sediment layer (cm)` == "0_1", "0-1 cm", "1-3 cm"),
            genus = str_squish(Genus), count = Count) |>
  mutate(genus = if_else(genus == "Sabatiera", "Sabatieria", genus)) |>
  summarise(raw_count = sum(count), .by = c(station, layer, genus))

check <- taxa |>
  filter(year == 2011, !is.na(genus), value > 0) |>
  inner_join(raw11, by = c("station", "layer", "genus"))
n_matched <- nrow(check)
n_mismatch <- sum(abs(check$value - check$raw_count) > 1e-6)
log_issue("verification",
          paste0("2011 matrix vs 2011 raw sheet: ", n_matched, " genus by sample cells matched, ", n_mismatch, " disagree."),
          if (n_mismatch == 0) "Matrix 2011 contains raw counts (not percentages); layer coding confirmed." else "Investigate mismatches.")

# Put every sample on a common scale: percent of the sample total
sample_totals <- taxa |>
  summarise(total = sum(value), .by = c(year, sample))

off <- sample_totals |> filter(abs(total - 100) > 0.5)
if (nrow(off) > 0) {
  log_issue("samples",
            paste0("Sample totals differ from 100: ",
                   paste0(off$sample, " ", off$year, " = ", round(off$total, 1), collapse = "; "),
                   ". In 2011 these are counts, so the sample had fewer than 100 individuals."),
            "All samples renormalised to percent of their own total (rel_abund). Low-count samples are noisier; keep in mind for interpretation.")
}

taxa <- taxa |>
  left_join(sample_totals, by = c("year", "sample")) |>
  mutate(rel_abund = 100 * value / total) |>
  rename(raw_value = value) |>
  select(year, station, impact, layer, sample, family, genus, trophic, cp, raw_value, rel_abund)

log_issue("samples", "2010 matrix is percent composition, 2011 matrix is counts (about 100 individuals per sample). Original 2010 counts are not in the workbook.",
          "Both converted to rel_abund (percent). No abundance (density) analysis is possible, only composition.")

# ---- 5. MI / ITD per sample ------------------------------------------------
# Recompute from the cleaned data as a check against the values in the workbook.
recomputed <- taxa |>
  summarise(
    MI_calc  = sum(rel_abund / 100 * cp),
    ITD_calc = sum((tapply(rel_abund, trophic, sum) / 100)^2),
    .by = c(year, sample)
  )

samples <- indices |>
  mutate(across(c(MI_given, ITD_given), as.numeric)) |>
  left_join(recomputed, by = c("year", "sample")) |>
  mutate(
    station = str_remove(sample, "_[12]$"),
    layer = if_else(str_detect(sample, "_1$"), "0-1 cm", "1-3 cm")
  ) |>
  left_join(env |> select(year, station, impact), by = c("year", "station")) |>
  select(year, station, impact, layer, sample,
         MI_workbook = MI_given, ITD_workbook = ITD_given, MI_calc, ITD_calc)

mi_bad  <- samples |> filter(abs(MI_workbook - MI_calc) > 0.1)
itd_bad <- samples |> filter(abs(ITD_workbook - ITD_calc) > 0.05)

# The 2010 workbook MI row repeats one value (1.577333) across many samples,
# which looks like a value copied down by mistake. ITD agrees with the recomputed values.
repeated_mi <- samples |> filter(year == 2010) |> count(MI_workbook, sort = TRUE) |> slice(1)
log_issue("indices",
          paste0("Workbook MI vs recomputed from cleaned data: ", nrow(mi_bad), " of ", nrow(samples),
                 " samples differ by >0.1 (all 2010; the same MI value, ", round(repeated_mi$MI_workbook, 3),
                 ", is repeated in ", repeated_mi$n, " samples). Workbook ITD: ", nrow(itd_bad),
                 " samples differ by >0.05."),
          "Use recomputed MI (matches 2011 workbook values to rounding) and workbook ITD (matches recomputed). Workbook values are kept as MI_workbook and ITD_workbook. The 2010 MI row should be checked at the source.")

samples <- samples |> mutate(MI = MI_calc, ITD = ITD_workbook)

log_issue("indices", "MI and ITD are calculated directly from the same genus composition (c-p values and trophic groups).",
          "Do not use genus composition to 'predict' MI or ITD; that would be circular.")

# ---- 6. Family-level table (cross-year comparison) -------------------------
family_long <- taxa |>
  filter(!is.na(family)) |>
  summarise(rel_abund = sum(rel_abund), .by = c(year, station, impact, layer, sample, family))

# Zero-fill so every family is present in every sample
family_long <- family_long |>
  complete(nesting(year, station, impact, layer, sample), family, fill = list(rel_abund = 0))

fam_year <- family_long |> summarise(p = sum(rel_abund > 0), .by = c(family, year)) |>
  summarise(in_both = all(p[year == 2010] > 0, p[year == 2011] > 0), .by = family)
n_genera <- taxa |> summarise(n = n_distinct(genus, na.rm = TRUE), .by = year)

log_issue("taxa",
          paste0("Genus lists differ between years (2010: ", n_genera$n[n_genera$year == 2010],
                 " genera, 2011: ", n_genera$n[n_genera$year == 2011],
                 "). Part of this is resolution/naming, not real turnover."),
          paste0("Comparisons between years use family level (", n_distinct(family_long$family), " families; ",
                 sum(fam_year$in_both), " present in both years)."))

# Odd family names (superfamily/subfamily used as family)
odd_fam <- family_names[str_detect(family_names, "oidea$|inae$")]
if (length(odd_fam) > 0) {
  log_issue("taxa", paste0("A rank other than family appears in the Family column: ", paste(odd_fam, collapse = ", "), "."),
            "Kept as given. It is unclear whether these should be merged into a family.")
}

# ---- 7. Write outputs ------------------------------------------------------
write_csv(env, "data/processed/env_stations.csv")
write_csv(samples, "data/processed/samples.csv")
write_csv(taxa, "data/processed/community_genus_long.csv")
write_csv(family_long, "data/processed/community_family_long.csv")
write_csv(cleaning_log, "outputs/tables/data_cleaning_log.csv")

message("Done. ", nrow(cleaning_log), " items written to outputs/tables/data_cleaning_log.csv")
