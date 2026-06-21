# Bayesian Hierarchical Real-Estate Price Modeling Across US States and ZIP Codes

This repository contains a reproducible final-project scaffold for a fully Bayesian analysis of US residential listing prices. The planned analysis estimates partial-pooling effects at both the state and ZIP-code levels using MCMC in R/JAGS.

## Data source

- Dataset: [USA Real Estate Dataset](https://www.kaggle.com/datasets/ahmedshahriarsakib/usa-real-estate-dataset)
- Kaggle dataset version: 25
- Download date: 2026-06-21
- Raw CSV: `data/raw/realtor-data.csv`
- Raw CSV size: 178,858,272 bytes (170.6 MiB)
- MD5 checksum: `09ce5bd36ef648a5769a996735c05951`
- Verified CSV columns: `brokered_by`, `status`, `price`, `bed`, `bath`, `acre_lot`, `street`, `city`, `state`, `zip_code`, `house_size`, and `prev_sold_date`.

The version-25 Kaggle archive was downloaded on the date above, extracted locally, and preserved without opening or resaving the CSV. The raw data are deliberately excluded from version control.

## Planned model sequence

1. M0: pooled Bayesian regression.
2. M1: Bayesian regression with state-level random intercepts.
3. M2: Bayesian regression with nested state and ZIP-code random intercepts.

The response is log listing price. Planned covariates are log living area, log lot size, bedrooms, bathrooms, and listing status.

## Reproducibility decisions

- Random seed: 20260621
- Planned sampling design: 25 states × 24 ZIP codes per state × 30 listings per ZIP.
- Planned train/test split within each selected ZIP: 24/6 listings.
- Planned convergence criteria: four chains, split R-hat below 1.01, adequate bulk/tail ESS, small Monte Carlo standard errors, and satisfactory trace/rank/ACF diagnostics.

## Repository layout

- data/: local raw, processed, and optional external data.
- models/: JAGS model specifications.
- scripts/: numbered R analysis scripts.
- output/: generated figures, tables, MCMC draws, and logs.
- report/: Quarto/R Markdown report source and bibliography.

## Current scope

The raw dataset has been downloaded and catalogued. Data import, model code, and analysis execution are intentionally deferred to later implementation phases.
