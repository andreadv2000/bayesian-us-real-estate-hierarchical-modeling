# Bayesian Hierarchical Modeling of U.S. Real-Estate Listing Prices

This project analyzes U.S. residential real-estate listing prices with fully Bayesian hierarchical regression models. The analysis estimates how property characteristics and geographic structure relate to listing prices while accounting for partial pooling across states and ZIP codes.

The final model is a heteroskedastic state/ZIP hierarchical model fitted with MCMC in R/JAGS. It models standardized log listing price with fixed property effects, state and ZIP random intercepts, and a residual scale that varies by listing covariates and state.

## Data

The analysis uses the USA Real Estate Dataset from Kaggle:

- Dataset: [USA Real Estate Dataset](https://www.kaggle.com/datasets/ahmedshahriarsakib/usa-real-estate-dataset)
- Kaggle dataset version: 25
- Download date: 2026-06-21
- Raw CSV path expected locally: `data/raw/realtor-data.csv`
- Raw CSV size: 178,858,272 bytes
- MD5 checksum: `09ce5bd36ef648a5769a996735c05951`
- Raw columns: `brokered_by`, `status`, `price`, `bed`, `bath`, `acre_lot`, `street`, `city`, `state`, `zip_code`, `house_size`, `prev_sold_date`

The raw data file is intentionally excluded from version control. Reproducing the full analysis requires downloading the same dataset version and placing the CSV at the path above.

## Modeling approach

The response variable is standardized log listing price. The core predictors are:

- log living area
- log lot size
- bedrooms
- bathrooms
- listing status

Three Bayesian regression models are compared:

| Model | Description |
|---|---|
| M0 | Pooled Bayesian regression |
| M1 | Bayesian regression with state-level random intercepts |
| M2 | Heteroskedastic hierarchical regression with state and ZIP random intercepts |

The selected M2 model uses the following structure:

- state effects are centered deviations from the national intercept;
- ZIP effects are centered within state;
- residual standard deviation varies with covariates and state-level residual-scale effects;
- all models use the same frozen training/test split for comparison.

## Main results

The heteroskedastic M2 is the strongest model by both in-sample DIC and held-out predictive performance.

| Model | DIC | RMSE | MAE | 95% coverage | MLPD |
|---|---:|---:|---:|---:|---:|
| M2: heteroskedastic state + ZIP | 15,333.34 | 0.4241 | 0.2825 | 0.9536 | -0.4905 |
| M1: state | 29,452.74 | 0.5970 | 0.4288 | 0.9572 | -0.9155 |
| M0: pooled | 35,023.39 | 0.7623 | 0.5574 | 0.9522 | -1.1507 |

Selected posterior summaries for M2:

| Quantity | Posterior median | 95% credible interval |
|---|---:|---:|
| log living area effect | 0.3227 | [0.3104, 0.3349] |
| log lot size effect | 0.1365 | [0.1229, 0.1503] |
| bedroom effect | -0.0147 | [-0.0254, -0.0040] |
| bathroom effect | 0.2089 | [0.1936, 0.2240] |
| sold-status effect | -0.0710 | [-0.0853, -0.0568] |
| typical residual SD | 0.3875 | [0.3810, 0.3941] |
| state SD | 0.4689 | [0.3663, 0.6200] |
| ZIP SD | 0.4494 | [0.4237, 0.4773] |

Posterior predictive checks show that M2 reproduces geographic heterogeneity well. It retains mild over-dispersion in the overall standardized log-price distribution:

| Statistic | Observed | M2 replicated median |
|---|---:|---:|
| Overall IQR | 1.1087 | 1.2153 |
| Overall 90th percentile | 1.1274 | 1.2182 |
| ZIP mean variance | 0.5640 | 0.5759 |

This is reported as a model-checking diagnostic rather than a fitting failure.

## Repository structure

```text
data/
  raw/                 local raw data, not versioned
  processed/           derived analysis objects, not versioned
models/                JAGS model specifications
scripts/               R scripts for data preparation, modeling, diagnostics, and artifacts
output/
  figures/             generated figures
  tables/              generated tables
  mcmc/                large MCMC draw files, not versioned
  logs/                runtime logs, not versioned
report/                report source and rendered report files
```

## Reproducibility

The project uses R, JAGS, and R packages recorded in `renv.lock`.

Recommended setup:

```r
renv::restore()
```

The full analysis can be regenerated from the project root with:

```powershell
Rscript scripts/01_data_audit.R
Rscript scripts/02_make_sample.R
Rscript scripts/03_prepare_jags_data.R
Rscript scripts/04_prior_predictive_checks.R
Rscript scripts/04_simulation_recovery.R
Rscript scripts/05_fit_models.R
Rscript scripts/06_diagnostics.R
Rscript scripts/07_lightweight_dic.R
Rscript scripts/07_model_comparison.R
Rscript scripts/08_posterior_predictive_checks.R
Rscript scripts/10_posterior_predictions.R
Rscript scripts/11_sensitivity_analyses.R
Rscript scripts/09_make_figures_tables.R
```

The M2 fit is computationally expensive. Existing saved MCMC draw files are reused by default unless forced through the corresponding environment variables in the fitting scripts.

## Key outputs

Important generated tables include:

- `output/tables/report_dic_heldout_model_comparison.csv`
- `output/tables/report_m2_fixed_effect_variance_summary.csv`
- `output/tables/report_posterior_predictive_test_statistics.csv`
- `output/tables/report_heldout_predictive_metrics.csv`
- `output/tables/report_artifact_manifest.csv`

Important generated figures include:

- `output/figures/m2_posterior_predictive_checks.png`
- `output/figures/m2_posterior_prediction_intervals.png`
- `output/figures/m2_fixed_effect_forest.png`
- `output/figures/m2_state_effects_caterpillar.png`
- `output/figures/m2_selected_zip_effects_caterpillar.png`
- `output/figures/m2_variance_component_posteriors.png`

The report source is available at `report/report.Rmd`.

## License and citation

This repository is intended for academic and reproducible-analysis use. The underlying data are governed by the terms of the Kaggle dataset source.
