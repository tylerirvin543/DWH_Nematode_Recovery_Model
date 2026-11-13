# ==============================================================================
# 01  DATA PREPARATION
#     Turn the original Excel workbook into clean, tidy tables
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
#   Reads the workbook with the 2010 and 2011 nematode data and produces clean
#   tables that every later script reads. Nothing is analysed here. The aim is to
#   make every cleaning decision visible, so each one is written to a log.
#
# WHY IT IS NEEDED
#   The workbook is laid out for reading, not for analysis. Each year has its own
#   sheet, samples are columns, summary rows (MI and ITD) sit at the bottom of the
#   taxa lists, the two years use different conventions (percent in 2010, counts
#   in 2011), and a few names are spelled differently between years.
#
# MAIN STEPS
#   1. Read the station-level environmental data (both years).
#   2. Read the taxon by sample matrices and split off the MI / ITD summary rows.
#   3. Clean taxon names and attributes (case, spelling, family names in the
#      genus column, blank cells).
#   4. Put every sample on one scale (percent of its own total) and check the
#      2011 matrix against the record level 2011 sheet.
#   5. Recompute the maturity index (MI) and check the workbook's values.
#   6. Build a family level table, used for all comparisons between years.
#   7. Write everything out, including the cleaning log.
#
# RUN
#   From the project root (open "Group Project.Rproj"), before every other script.
#
# INPUT
#   data/raw/Clean data_Recovery in 1, 2, 3....xlsx
#
# OUTPUTS
#   data/processed/env_stations.csv           one row per station and year
#   data/processed/samples.csv                one row per sample (station x layer x year),
#                                             with MI and ITD
#   data/processed/community_genus_long.csv   one row per taxon per sample, finest
#                                             taxonomic resolution available
#   data/processed/community_family_long.csv  same, summed to family
#   outputs/tables/data_cleaning_log.csv      every cleaning decision, in plain language
#
# KEY VOCABULARY
#   HI / MI / NI   heavily impacted, moderately impacted, not impacted (controls).
#                  Note: "MI" is also the maturity index; context makes it clear.
#   c-p value      colonizer-persister score of a genus, 1 (colonizer) to 5 (persister)
#   MI             maturity index = sum of (proportion of genus x its c-p value)
#   ITD            index of trophic diversity = sum of squared proportions of the
#                  four feeding groups (0.25 = all equally common, 1 = one dominates)

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(readr)

raw_file <- "data/raw/Clean data_Recovery in 1, 2, 3....xlsx"

# Fail early with a clear message if the workbook is not where we expect it.
stopifnot(file.exists(raw_file))

# Make sure the output folders exist (harmless if they already do).
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

# ---- Cleaning log -----------------------------------------------------------
# Every time the script finds something odd, it calls log_issue() to record what
# the problem was and what we did about it. The finished log is saved as a CSV and
# printed in the paper's supplement, so no cleaning step is hidden.
cleaning_log <- tibble(step = character(), issue = character(), action = character())
log_issue <- function(step, issue, action) {
  # <<- writes to the log in the global environment, not a local copy
  cleaning_log <<- bind_rows(cleaning_log, tibble(step = step, issue = issue, action = action))
}

# ---- 1. Environmental variables --------------------------------------------
# One sheet per year, same columns. We give the columns short, code friendly names.
read_env <- function(sheet) {
  read_excel(raw_file, sheet = sheet) |>
    rename(year = Year, station = Station, impact = Impact, depth_m = Depth,
           lat = `Lat (N)`, long = `Long (W)`)
}

# Stack the two years; make impact an ordered factor so plots and models always
# list the groups in the same order (HI, MI, NI).
env <- bind_rows(read_env("Env variables 2010"), read_env("Env variables 2011")) |>
  mutate(impact = factor(impact, levels = c("HI", "MI", "NI")))

# Sanity check: 9 stations x 2 years = 18 rows, and no blanks anywhere.
stopifnot(nrow(env) == 18, !anyNA(env))

# Compare the typical (median) value of each chemical variable in 2011 against
# 2010. Large jumps in the whole study area are suspicious: they may be a change
# in units or laboratory method and not a real change in the sediment.
env_ratio <- env |>
  pivot_longer(c(Al, Ba, PAH44, TPH, Mud, TOC), names_to = "variable") |>
  group_by(variable, year) |>
  summarise(med = median(value), .groups = "drop") |>
  pivot_wider(names_from = year, values_from = med) |>
  mutate(ratio = `2011` / `2010`)

# Flag anything that changed more than two fold in either direction. We keep the
# values as given (the team agreed to treat them as comparable for now) but we
# log the concern so it is not forgotten.
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
# Layout of each "Matrix" sheet: one row per taxon, one column per sample
# (e.g. D040S_1 = station D040S, surface layer), then two summary rows at the
# bottom holding the workbook's own MI and ITD values.
# 2010 stores percent composition; 2011 stores counts (about 100 animals per sample).

read_matrix <- function(sheet, year) {
  # .name_repair = "minimal" keeps the sample column names exactly as written.
  m <- read_excel(raw_file, sheet = sheet, .name_repair = "minimal") |>
    rename(family = Family, genus = Genus, trophic = `Trophic gr.`, cp = `c-p value`)

  # Every column that is not a taxon attribute is a sample.
  sample_cols <- setdiff(names(m), c("family", "genus", "trophic", "cp"))

  # The summary rows have no genus, and the label "MI" or "ITD" sits in the c-p
  # column. Pull them out and reshape to one row per sample.
  indices <- m |>
    filter(is.na(genus), cp %in% c("MI", "ITD")) |>
    select(index = cp, all_of(sample_cols)) |>
    pivot_longer(-index, names_to = "sample", values_to = "value") |>
    pivot_wider(names_from = index, values_from = value) |>
    rename(MI_given = MI, ITD_given = ITD)

  # Everything with a genus is a real taxon row. Reshape from wide (one column
  # per sample) to long (one row per taxon per sample), which is what the later
  # scripts expect.
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
# str_squish() trims stray spaces (some family names had trailing spaces, which
# would make "Xyalidae" and "Xyalidae " look like different families).
# Trophic groups are upper-cased because 2010 used "1A" and 2011 used "1a".
# The c-p value was stored as text in 2011, so we convert it to a whole number.
taxa <- taxa |>
  mutate(
    family  = str_squish(family),
    genus   = str_squish(genus),
    trophic = str_to_upper(str_squish(trophic)),
    cp      = as.integer(cp)
  )

log_issue("taxa", "Trophic group coded 1A/1B/2A/2B in 2010 but 1a/1b/2a/2b in 2011; c-p value stored as text in 2011; trailing spaces in Family names.",
          "Set trophic groups to upper case, converted the c-p value to a number, trimmed whitespace.")

# One genus appears spelled two ways: Sabatieria (2010 only) and Sabatiera (2011
# only). That looks like a typo, so we use one spelling. This is an assumption.
if (any(taxa$genus == "Sabatiera")) {
  taxa <- taxa |> mutate(genus = if_else(genus == "Sabatiera", "Sabatieria", genus))
  log_issue("taxa", "'Sabatiera' appears only in 2011 and 'Sabatieria' only in 2010, which looks like a spelling variant.",
            "Renamed to Sabatieria. This is an assumption that should be confirmed.")
}

# In 2011 a few animals were only identified to family, and the family name was
# typed into the Genus column. We can't call those genera, so we blank the genus
# (keeping the family) so they count at family level but are not mistaken for genera.
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

# A taxon with no family would silently drop out of the family level table, so
# we log it if it ever happens (it doesn't in the current data).
if (anyNA(taxa$family)) {
  log_issue("taxa", "Some taxa have no family.", "Dropped from the family table and flagged for manual checking.")
}

# A blank cell in the abundance matrix (one in 2011). A blank means the taxon was
# not recorded in that sample, so we treat it as zero.
na_cells <- taxa |> filter(is.na(value)) |> distinct(year, genus, sample)
if (nrow(na_cells) > 0) {
  log_issue("taxa",
            paste0("Blank abundance cell(s): ",
                   paste0(na_cells$genus, " in ", na_cells$sample, " (", na_cells$year, ")", collapse = "; "), "."),
            "Treated as 0 (taxon not observed).")
  taxa <- taxa |> mutate(value = replace_na(value, 0))
}

# ---- 4. Sample metadata and relative abundance -----------------------------
# Sample names look like "D040S_1". The part before the underscore is the
# station; the suffix is the sediment layer (_1 = 0 to 1 cm, _2 = 1 to 3 cm,
# confirmed against the 2011 record level sheet below). Then we attach each
# station's impact level (HI / MI / NI).
taxa <- taxa |>
  mutate(
    station = str_remove(sample, "_[12]$"),
    layer   = if_else(str_detect(sample, "_1$"), "0-1 cm", "1-3 cm")
  ) |>
  left_join(env |> select(year, station, impact), by = c("year", "station"))

stopifnot(!anyNA(taxa$impact))

# Cross-check: the 2011 matrix should match the record level "2011 raw" sheet
# (one row per animal group). We add up the raw counts by station, layer and
# genus and compare them with the matrix. This also confirms the layer coding
# and that the 2011 matrix really holds counts, not percentages.
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

# Put every sample on the same scale: percent of that sample's total. In 2010 the
# totals are already ~100; in 2011 they are the number of animals counted, so
# the same step converts counts to percent. After this, the two years are comparable
# in composition (but not in density).
sample_totals <- taxa |>
  summarise(total = sum(value), .by = c(year, sample))

# Samples that do not add up to 100 are worth knowing about: in 2011 they mean
# fewer than 100 animals were available, so those samples are noisier.
off <- sample_totals |> filter(abs(total - 100) > 0.5)
if (nrow(off) > 0) {
  log_issue("samples",
            paste0("Sample totals differ from 100: ",
                   paste0(off$sample, " ", off$year, " = ", round(off$total, 1), collapse = "; "),
                   ". In 2011 these are counts, so the sample had fewer than 100 individuals."),
            "All samples renormalised to percent of their own total (rel_abund). Low-count samples are noisier; keep in mind for interpretation.")
}

# rel_abund = percent composition (used by the community analyses).
# raw_value keeps the original number (percent in 2010, count in 2011); the
# colonizer model in script 06 needs it to rebuild counts.
taxa <- taxa |>
  left_join(sample_totals, by = c("year", "sample")) |>
  mutate(rel_abund = 100 * value / total) |>
  rename(raw_value = value) |>
  select(year, station, impact, layer, sample, family, genus, trophic, cp, raw_value, rel_abund)

log_issue("samples", "2010 matrix is percent composition, 2011 matrix is counts (about 100 individuals per sample). Original 2010 counts are not in the workbook.",
          "Both converted to rel_abund (percent). No abundance (density) analysis is possible, only composition.")

# ---- 5. MI / ITD per sample ------------------------------------------------
# Recompute both indices from the cleaned composition and compare them with the
# values stored in the workbook, to check the workbook.
#   MI  = sum over genera of (proportion x c-p value)
#   ITD = sum over the four trophic groups of (proportion in that group)^2
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

# Finding: the 2010 MI row in the workbook repeats a single value (1.577333)
# across many samples, which looks like a value copied down by mistake. The
# recomputed MI does not have that problem, and ITD agrees exactly. So we use the
# recomputed MI and the workbook ITD. Neither index is used in the main analyses.
repeated_mi <- samples |> filter(year == 2010) |> count(MI_workbook, sort = TRUE) |> slice(1)
log_issue("indices",
          paste0("Workbook MI vs recomputed from cleaned data: ", nrow(mi_bad), " of ", nrow(samples),
                 " samples differ by >0.1 (all 2010; the same MI value, ", round(repeated_mi$MI_workbook, 3),
                 ", is repeated in ", repeated_mi$n, " samples). Workbook ITD: ", nrow(itd_bad),
                 " samples differ by >0.05."),
          "Use recomputed MI (matches 2011 workbook values to rounding) and workbook ITD (matches recomputed). Workbook values are kept as MI_workbook and ITD_workbook. The 2010 MI row should be checked at the source.")

samples <- samples |> mutate(MI = MI_calc, ITD = ITD_workbook)

# A reminder for later: both indices are computed from the same genus data we
# analyse, so predicting them from genus composition would be circular.
log_issue("indices", "MI and ITD are calculated directly from the same genus composition (c-p values and trophic groups).",
          "Do not use genus composition to 'predict' MI or ITD; that would be circular.")

# ---- 6. Family-level table (cross-year comparison) -------------------------
# Genus lists differ between years (partly naming and identification detail), so
# all comparisons between years use families. Add up the genera within each family.
family_long <- taxa |>
  filter(!is.na(family)) |>
  summarise(rel_abund = sum(rel_abund), .by = c(year, station, impact, layer, sample, family))

# Zero-fill: make sure every family has a row in every sample (0 if absent).
# Later steps that build a samples-by-families matrix rely on this.
family_long <- family_long |>
  complete(nesting(year, station, impact, layer, sample), family, fill = list(rel_abund = 0))

# How many families occur in both years? (reported in the log, not used further)
fam_year <- family_long |> summarise(p = sum(rel_abund > 0), .by = c(family, year)) |>
  summarise(in_both = all(p[year == 2010] > 0, p[year == 2011] > 0), .by = family)
n_genera <- taxa |> summarise(n = n_distinct(genus, na.rm = TRUE), .by = year)

log_issue("taxa",
          paste0("Genus lists differ between years (2010: ", n_genera$n[n_genera$year == 2010],
                 " genera, 2011: ", n_genera$n[n_genera$year == 2011],
                 "). Part of this is resolution/naming, not real turnover."),
          paste0("Comparisons between years use family level (", n_distinct(family_long$family), " families; ",
                 sum(fam_year$in_both), " present in both years)."))

# One entry in the Family column is really a superfamily or subfamily (names
# ending in -oidea or -inae). We keep it as given but flag it.
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
