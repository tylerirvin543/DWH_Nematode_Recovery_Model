# Recovery in 1, 2, 3

NRES 746 final project (Fall 2026) by Tyler Irvin, Cindy Murillo, and Yazmin Elizarraraz.

Did deep sea nematode communities at oil impacted stations in the Gulf of Mexico become more similar to unimpacted stations between 2010 and 2011, after the Deepwater Horizon spill? We compare communities directly, within each year and against that year's controls, and fit a likelihood model to the share of colonizer nematodes.

Short answer from the data we have: impacted and unimpacted stations hosted different communities in both years, sediment hydrocarbons (Ba, PAH44, TPH) line up with the community differences, and we found no support for recovery by 2011. Two years, nine stations, and a shift shared by all stations limit what we can conclude. The paper and slides give the details and caveats.

## Data are not in this repository

The data are proprietary. This repository contains code only. The original workbook, the cleaned tables, and every table derived from them are excluded by `.gitignore`. Figures and rendered documents show results, not the underlying records.

To reproduce the analysis you need the original workbook, which has to come from the data owners (Yazmin Elizarraraz can say who to ask). Save it as:

```
data/raw/Clean data_Recovery in 1, 2, 3....xlsx
```

then run the scripts in order. They recreate `data/processed/` and `outputs/tables/` locally, and neither folder will be committed.

## Study design

Nine stations on the deep seafloor: three heavily impacted (HI), three moderately impacted (MI), and three unimpacted controls (NI). Each was sampled in 2010 and 2011, with two sediment layers per visit (0 to 1 cm and 1 to 3 cm), for 36 samples. Nematodes were identified to genus where possible and always to family, and four sediment variables were measured per station and year: aluminum (Al), barium (Ba), 44 polycyclic aromatic hydrocarbons (PAH44), and total petroleum hydrocarbons (TPH).

## Analysis in brief

1. **Clean and tidy** the workbook, logging every decision.
2. **Genus level, each year on its own.** Layers are averaged, then Bray-Curtis dissimilarities are computed per year (square root percent; untransformed percent as a check).
3. **Environment.** PCA of the four environmental variables, and a PCoA of the Bray-Curtis matrix with the variables fitted as arrows (`envfit`).
4. **PERMANOVA** per year on impact level, and impact plus the first two environmental PCs (marginal tests). A pooled test on both years shuffles years within stations.
5. **Distance to the controls.** The gap is the mean distance to the NI stations minus how far the NI stations are from each other. Mixed models (station random intercept, layers as repeated measures) are fitted per year, then one model for both years gives the change in gap per year.
6. **Colonizer model.** Colonizer (c-p 1) counts follow a beta binomial model written by hand, fitted by maximum likelihood with profile likelihood intervals.

## Folders

- `data/raw` holds the original workbook (not tracked)
- `data/processed` holds the tidy tables made by script 01 (not tracked)
- `scripts` holds the analysis, numbered in the order to run them
- `outputs/figures` holds the figures (tracked)
- `outputs/tables` holds the results tables made by the scripts (not tracked)
- `docs` holds the course guidelines and our proposal

## Scripts

Open `Group Project.Rproj` and run these from the project folder, in order.

| Script | What it does | Main outputs |
|---|---|---|
| `01_data_prep.R` | Cleans the workbook, builds genus and family tables, writes the cleaning log | `data/processed/*`, `data_cleaning_log.csv` |
| `02_ordination_permanova.R` | Family level NMDS and PERMANOVA on both years together (sensitivity check) | `nmds_trajectories.png`, PERMANOVA tables |
| `03_distance_to_control.R` | Family level distance to the controls and mixed model, both years together (sensitivity check) | `distance_to_control.png`, model tables |
| `04_cp_composition.R` | Compares the colonizer and persister (c-p) makeup of HI and NI communities | `cp_composition.png` |
| `05_colonizer_model.R` | Beta binomial model of the colonizer share, by maximum likelihood | `colonizer_fit.png`, `colonizer_profiles.png` |
| `06_genus_by_year.R` | Genus level, each year alone: Bray-Curtis, environmental PCA, PCoA with envfit, per-year PERMANOVA, pooled year test | `genus_nmds_by_year.png`, `env_pca_by_year.png`, `genus_pcoa_env.png` |
| `07_genus_distance_by_year.R` | Genus level gap to the controls within each year, mixed models per year and pooled | `genus_gap_by_year.png`, gap model tables |

The paper reports scripts 04 to 07. Scripts 02 and 03 are the first, family level version of the pooled analysis and are kept as a sensitivity check. Each script has a header describing its question, inputs, outputs, and design choices.

## Documents

- `nematode_recovery_paper.qmd` is the paper
- `nematode_recovery_supplement.qmd` holds the data notes, one line per method, and the cleaning log
- `nematode_recovery_slides.qmd` is the talk (reveal.js)

Render them with Quarto from the project root after running the scripts, for example `quarto render nematode_recovery_paper.qmd`. The rendered Word and HTML files are committed for convenience, but they are rebuilt from the `.qmd` files and the local data. The slides also use images in `images/`, which are not tracked.

## R packages

readxl, dplyr, tidyr, stringr, purrr, readr, ggplot2, ggrepel, vegan, nlme, knitr, tibble. Quarto is needed for the documents.

## Notes

- The colonizer question came out of looking at how composition differed between HI and NI stations, so that part is exploratory and not a pre-planned test.
- We first planned a Bayesian state-space simulation, but two sampling years cannot support one. After feedback from Kevin Shoemaker we switched to direct comparisons and a likelihood model.
- No correction was made for looking at several questions.
- The environmental values differ a lot between years at the control stations (for example PAH44), so they may not have been measured the same way. That is why the environmental analyses are done within each year.
- This work used generative AI tools, mainly Posit Assistant in RStudio, for coding help and draft text. The authors chose the analyses, ran the code, and take responsibility for the content.
