# Selected M2 for caveat-fix branch

This branch keeps the original production M2 as the baseline/reference model and selects the heteroskedastic state/ZIP model as the branch's official M2 candidate.

Selected model:

- Model name: `M2_heteroskedastic_xall_state`
- JAGS model file: `models/m2_heteroskedastic_state_sigma.jags`
- Experiment script: `scripts/17_m2_heteroskedastic_experiment.R`
- Draws: `output/mcmc/m2_heteroskedastic_xall_state_draws.rds`
- Posterior summary: `output/tables/m2_heteroskedastic_xall_state_posterior_summary.csv`
- Held-out metrics: `output/tables/m2_heteroskedastic_xall_state_heldout_predictive_metrics.csv`
- PPC summary: `output/tables/m2_heteroskedastic_xall_state_posterior_predictive_check_summary.csv`

Selection rationale:

- The full run converged for all monitored parameters: 665/665.
- The selected model improved held-out RMSE, MAE, and MLPD relative to the baseline M2.
- 95% predictive coverage stayed calibrated: 0.955 versus 0.971 for baseline M2.
- The PPC over-dispersion caveat improved: replicated median IQR moved from 1.290 under baseline M2 to 1.215 under the heteroskedastic M2, closer to the observed 1.109.
- The replicated q90 moved from 1.263 under baseline M2 to 1.217 under the heteroskedastic M2, closer to the observed 1.127.
- ZIP-level geographic heterogeneity remained close to the observed value.

Important limitation:

The heteroskedastic M2 improves, but does not completely eliminate, the remaining upper-tail/IQR over-dispersion. In a report, describe this as an improved caveat rather than a fully solved issue.
