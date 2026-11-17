# Recovery in 1, 2, 3

NRES 746 final project (Fall 2026) by Tyler Irvin, Cindy Murillo, and Yazmin Elizarraraz.

We are asking whether nematode communities at oil impacted deep sea stations in the Gulf of Mexico became more similar to unimpacted stations between 2010 and 2011, after the Deepwater Horizon spill. The data are nematode counts from nine stations (three heavily impacted, three moderately impacted, three unimpacted) sampled in both years, with two sediment layers each.

We first planned a Bayesian state-space simulation, but two sampling years cannot support one, so after feedback from Kevin Shoemaker we switched to comparing communities directly and fitting a likelihood model to the share of colonizer nematodes. The paper is `nematode_recovery_paper.qmd`, with the data cleaning log in `nematode_recovery_supplement.qmd`, and the talk is `nematode_recovery_slides.qmd`.

## Folders

- `data/raw` holds the original workbook (not tracked in git)
- `data/processed` holds the tidy tables made by script 01
- `scripts` holds the analysis, numbered in the order to run them
- `outputs/figures` and `outputs/tables` hold what the scripts produce, including the data cleaning log
- `docs` holds the course guidelines and our proposal

## Scripts

Open `Group Project.Rproj` and run these from the project folder, in order.

1. `01_data_prep.R` cleans the workbook and builds genus and family tables
2. `02_ordination_permanova.R` runs the NMDS map and PERMANOVA
3. `03_distance_to_control.R` measures distance to the controls and fits the mixed model
4. `04_cp_composition.R` compares the colonizer and persister (c-p) makeup of heavily impacted and unimpacted communities
5. `05_colonizer_model.R` fits the beta binomial model of the colonizer share by maximum likelihood

Then render the paper, supplement, and slides with Quarto. R packages used: readxl, dplyr, tidyr, stringr, purrr, readr, ggplot2, vegan, nlme, knitr.

A note on how the analysis developed: the colonizer question came out of looking at how community composition differed between heavily impacted and unimpacted stations, so that part is exploratory and not a pre-planned test.

## Data

The data were shared with us by Yazmin. We should confirm we have permission to publish them before the raw workbook or anything derived from it goes anywhere public.
