# Recovery in 1, 2, 3: nematode communities after the Deepwater Horizon spill

NRES 746 (Fall 2026) final group project.

Team: Tyler Irvin, Cindy Murillo, Yazmin Elizarraraz

## What we are asking

Between 2010 and 2011, did nematode communities at oil impacted deep sea stations in the Gulf of Mexico become more similar to communities at unimpacted stations? The data are nematode counts from nine stations (three heavily impacted, three moderately impacted, three unimpacted) sampled in 2010 and 2011, with two sediment layers per station.

We started with a Bayesian state-space simulation and changed course after feedback from our instructor, because two sampling years cannot support a model that projects recovery over a decade. The current analysis describes how the communities compare, tests whether the differences are larger than chance, and measures directly whether impacted stations moved closer to the controls.

The full write-up, including the equations used, is in `nematode_recovery_report.qmd` (rendered copy: `nematode_recovery_report.html`).

## Where things stand

- [x] Project pitch and proposal
- [x] Dataset received and cleaned
- [x] Ordination, PERMANOVA, distance to controls, supporting plots, SIMPER and colonizer check
- [x] Custom model fitted by maximum likelihood (beta-binomial model of the colonizer share)
- [ ] Later survey data (for example 2023), if it exists and is comparable
- [ ] Oral presentation
- [ ] Rough draft
- [ ] Final paper

## Repository layout

```
.
├── data/
│   ├── raw/          original workbook (kept locally, not tracked in git)
│   └── processed/    tidy tables written by 01_data_prep.R
├── scripts/          analysis scripts, run in numbered order
├── outputs/
│   ├── figures/      PNG figures used in the report
│   └── tables/       CSV results, including the data cleaning log
├── docs/             course guidelines, pitch, instructor feedback
├── nematode_recovery_report.qmd
└── nematode_recovery_report.html
```

## Scripts

Run these from the project root (open `Group Project.Rproj`), in order.

| Script | What it does |
|---|---|
| `01_data_prep.R` | Reads the workbook, cleans it, converts to percent composition, builds genus and family tables, writes a cleaning log |
| `02_ordination_permanova.R` | Bray-Curtis dissimilarity, NMDS with 2010 to 2011 arrows, PERMANOVA |
| `03_distance_to_control.R` | Distance of impacted stations to controls, control background, mixed model of change per year |
| `04_extra_plots.R` | Environmental PCA, family composition bars, MI and ITD plot, site map |
| `05_simper.R` | SIMPER for HI versus NI, colonizer and persister (c-p) comparison |
| `06_colonizer_model.R` | Hand-written beta-binomial model of the colonizer (c-p 1) share, fitted with `optim`, with profile likelihood intervals |

To rebuild the report after running the scripts, render `nematode_recovery_report.qmd` with Quarto from the project root.

R packages used: readxl, dplyr, tidyr, stringr, purrr, readr, ggplot2, vegan, nlme, knitr.

## Data

The data come from the Deepwater Horizon deep sea benthic sampling and were shared with the group by Yazmin. Before the raw workbook or anything derived from it is shared more widely, we should confirm that we have permission to publish it. Cleaning decisions and open questions about the dataset are listed in `outputs/tables/data_cleaning_log.csv` and in the report appendix.
