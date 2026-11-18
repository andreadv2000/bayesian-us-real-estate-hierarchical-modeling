# Why the heteroskedastic M2 was selected

The original state–ZIP model used a constant residual scale. Its posterior predictive simulations had a wider interquartile range and a higher 90th percentile than the observed standardized log prices. The selected variant keeps the state and ZIP intercepts and lets residual variation depend on the property covariates and state.

## Comparison with the earlier state–ZIP model

The standalone selection run used 3,600 held-out listings and 2,000 predictive draws. This comparison is between the two M2 variants. For the final comparison against the pooled and state-only models, see the [project README](../README.md#results) and [report table](../output/tables/report_dic_heldout_model_comparison.csv).

| M2 variant | Held-out RMSE ↓ | Held-out MAE ↓ | Mean log predictive density ↑ | 95% interval coverage |
| --- | ---: | ---: | ---: | ---: |
| Earlier constant residual scale | 0.4268 | 0.2870 | −0.6044 | 97.1% |
| Selected varying residual scale | 0.4243 | 0.2826 | −0.4907 | 95.5% |

The selected fit passed the convergence criteria for 665 of 665 monitored parameters. Its replicated median IQR moved from 1.290 to 1.215, toward the observed 1.109. Its replicated 90th percentile moved from 1.263 to 1.217, toward the observed 1.127. The remaining gaps show where the selected model still overstates spread and the upper tail.

## Sources

- [Selected JAGS model](../models/m2_heteroskedastic_state_sigma.jags) and [earlier model](../models/m2_state_zip.jags)
- [Production fitting script](../scripts/05_fit_models.R) and [standalone selection experiment](../scripts/17_m2_heteroskedastic_experiment.R)
- [Held-out comparison](../output/tables/selected_m2_heldout_comparison.csv), [posterior predictive checks](../output/tables/selected_m2_ppc_comparison.csv), and [convergence summary](../output/tables/selected_m2_convergence_status.csv)
