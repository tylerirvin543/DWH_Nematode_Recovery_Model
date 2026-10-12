# Recovery in 1, 2, 3... Nematode Community Recovery After the Deepwater Horizon Oil Spill

NRES 746 (Fall 2026) final group project.

## Team

- Tyler Irvin
- Cindy
- Yazmin Elizarraraz

## Project summary

We're examining whether deep-sea nematode communities in the Gulf of Mexico show
recovery following the 2010 Deepwater Horizon oil spill, using taxonomic and
trophic/life-history data collected at 9 stations (3 each at high-impact,
medium-impact, and control sites) in 2010 and 2011.

**Research questions**

1. Do nematode communities show measurable recovery by 2011 relative to 2010 baseline conditions?
2. Do nematodes recover faster than the broader meiofauna community?

**Hypothesis:** Nematode communities will show signs of recovery by 2013.

**Anticipated approach:** Bayesian multivariate state-space model estimating
taxon- and station-level transition dynamics between 2010 and 2011, used to
project forward community trajectories and compare simulated vs. observed
Bray-Curtis dissimilarity between high-impact and control stations over time.

See `docs/` for the full project pitch and the course's final project guidelines.

## Status

- [x] Project pitch drafted
- [ ] Group proposal (2-page) — in progress
- [ ] Nematode dataset obtained (pending from Yazmin)
- [ ] Preliminary data visualization
- [ ] Oral presentation
- [ ] Rough draft
- [ ] Final paper

## Repository structure

```
.
├── data/
│   ├── raw/          # Original nematode dataset(s), unmodified
│   └── processed/    # Cleaned/derived data used for modeling
├── scripts/          # Analysis scripts (data prep, model fitting, simulation, figures)
├── outputs/
│   ├── figures/      # Publication-quality figures for the final paper
│   ├── plots/        # Exploratory/diagnostic plots not intended for the final paper
│   └── tables/       # Model output tables, summary statistics
└── docs/             # Project guidelines, pitch, and proposal PDFs
```

## Data

Primary data: nematode taxonomic identification plus trophic and life-history
indices from 2010-2011, 9 stations (3 per impact level: high, medium, control).
Source/access details to be added once the dataset is obtained.

## Reproducing the analysis

Once scripts are added, run them in order from `scripts/` (numbered prefixes
will indicate run order, e.g. `01_data_prep.R`, `02_model.R`, ...).
