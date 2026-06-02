# Bayesian Logistic Regression for Low Birth Weight Risk Factors
## Metropolis-Hastings MCMC Implemented from Scratch

**Author:** Piyush Kumar
**Dataset:** `MASS::birthwt` — Hosmer, Lemeshow & Sturdivant (2013)
**Software:** R (MASS, tidyverse, gtsummary, broom, ggplot2)

---

## Overview

This project implements Bayesian logistic regression using a Metropolis-Hastings MCMC sampler
built entirely from scratch in base R — without using `rstanarm`, `brms`, or `rjags`. The goal
is to estimate posterior distributions of odds ratios for maternal risk factors associated with
low birth weight, and to demonstrate that Bayesian inference provides a richer quantification
of uncertainty than frequentist p-values alone. A frequentist logistic regression benchmark,
prior sensitivity analysis, and posterior predictive check are included for validation.

## Research Question

Which maternal demographic and clinical risk factors are associated with low birth weight
(< 2500g), and what is the posterior probability that each factor increases the odds of
low birth weight?

## Dataset

| Attribute | Value |
|---|---|
| Source | Baystate Medical Center, Springfield, Massachusetts (1986) |
| Reference | Hosmer, Lemeshow & Sturdivant (2013), *Applied Logistic Regression* |
| R package | `MASS::birthwt` |
| Study design | Retrospective observational study |
| Sample size | 189 singleton births |
| Outcome | Low birth weight: birth weight < 2500g (binary) |
| Missing data | 0 — complete-case analysis not required |

## Variables

| Variable | Description | Type |
|---|---|---|
| `low_num` | Low birth weight (1 = yes, 0 = no) | Outcome |
| `age_z` | Maternal age (standardised) | Continuous |
| `lwt_z` | Pre-pregnancy weight in lbs (standardised) | Continuous |
| `race` | Race: White (ref), Black, Other | Categorical |
| `smoking` | Smoking during pregnancy: No (ref), Yes | Binary |
| `hypertension` | History of hypertension: No (ref), Yes | Binary |
| `uterine_irritability` | Uterine irritability: No (ref), Yes | Binary |
| `previous_preterm_labor` | Any prior preterm labour: No (ref), Yes | Binary |

> Continuous predictors are standardised (mean = 0, SD = 1) to improve MCMC numerical
> stability and to place all coefficients on a comparable scale.

## Bayesian Framework

**Prior:** Independent Normal(0, 2.5²) on all coefficients — the same weakly informative
default used by `rstanarm`. Centred at zero (no prior effect); SD = 2.5 allows a broad
range of odds ratios while regularising extreme estimates in small samples.

**Likelihood:** Bernoulli logistic model.

**Sampler:** Random-walk Metropolis-Hastings MCMC. At each iteration a proposal is drawn
from Normal(0, 0.07²I), accepted or rejected via the log acceptance ratio, and the chain
updated accordingly. 40,000 total iterations; first 10,000 discarded as burn-in.
30,000 post-burn-in samples retained. Acceptance rate: 70.0%.

## Methods

- Variable recoding and standardisation
- Descriptive statistics (Table 1 by birth weight outcome)
- Group summary by outcome
- Frequentist logistic regression (benchmark)
- Bayesian logistic regression via MH-MCMC (from scratch)
- Posterior summary: median OR, 95% credible interval, P(OR > 1)
- Prior sensitivity analysis (σ = 1.5, 2.5, 5.0)
- Posterior predictive check
- MCMC diagnostics: trace plots, acceptance rate, posterior density plots

## Key Findings

**Cohort**
- n = 189; 59 low birth weight cases (31.2%), 130 normal (68.8%)
- Key unadjusted differences: smoking 50.8% vs 33.8%, hypertension 11.9% vs 3.8%,
  previous preterm labour 30.5% vs 9.2%, uterine irritability 23.7% vs 10.8%

**Posterior odds ratios — main model (Normal(0, 2.5) prior)**

| Predictor | Posterior OR | 95% CrI | P(OR > 1) | Evidence |
|---|---|---|---|---|
| Hypertension: Yes vs No | 6.09 | 1.65–25.18 | 99.5% | Credible increase |
| Previous preterm labor: Yes vs No | 3.65 | 1.55–9.61 | 99.9% | Credible increase |
| Black vs White | 3.13 | 1.20–9.21 | 99.1% | Credible increase |
| Smoking: Yes vs No | 2.16 | 1.04–4.74 | 98.1% | Credible increase |
| Other race vs White | 2.07 | 0.93–4.69 | 96.5% | Uncertain |
| Uterine irritability: Yes vs No | 1.82 | 0.80–4.48 | 92.0% | Uncertain |
| Age (standardized) | 0.80 | 0.53–1.19 | 13.1% | Uncertain |
| Maternal weight (standardized) | 0.62 | 0.41–0.94 | 1.1% | Credible decrease |

> "Credible" = 95% credible interval excludes OR = 1. P(OR > 1) = posterior probability
> that the factor increases the odds of low birth weight.

> **Note on hypertension:** Wide credible interval (1.65–25.18) reflects small cell count
> (n = 12 hypertensive mothers). Direction is robust across all prior specifications but
> the magnitude is uncertain.

**Comparison with frequentist model**

All eight predictors are consistent in direction and statistical conclusion between the
Bayesian and frequentist models, validating the MCMC implementation. Bayesian credible
intervals are slightly narrower than frequentist CIs due to the regularising effect of the
Normal(0, 2.5) prior.

**Prior sensitivity analysis**

Smoking, previous preterm labour, and Black race are stable across σ = 1.5, 2.5, and 5.0
— the data dominate the prior. Hypertension is more sensitive to prior scale, as expected
given n = 12 cases. All four key predictors remain above OR = 1 across all prior choices.

**Posterior predictive check**

Observed LBW cases: 59. Simulated mean: 59.4 (95% predictive interval: 44–75).
The observed count falls almost exactly at the simulated mean — no evidence of
systematic model misfit.

**MCMC diagnostics**

Acceptance rate: 70.0%. Trace plots for all key parameters show stationary, well-mixed
chains with no trends or sticky regions, confirming convergence.

## Repository Structure

```
03-bayesian-statistical-modelling/
├── scripts/
│   └── bayesian_low_birth_weight.R
├── data/
│   └── low_birth_weight_analysis_dataset.csv
├── outputs/
│   ├── cohort_summary.csv
│   ├── group_summary_by_birth_weight.csv
│   ├── frequentist_logistic_regression_results.csv
│   ├── frequentist_model_summary.txt
│   ├── bayesian_logistic_regression_results.csv
│   ├── bayesian_logistic_regression_results_labeled.csv
│   ├── bayesian_model_summary.txt
│   ├── prior_sensitivity_results.csv
│   ├── posterior_predictive_check_summary.csv
│   ├── posterior_samples_thinned.csv
│   └── table1.html
├── figures/
│   ├── bayesian_posterior_odds_ratio_forest_plot.png
│   ├── posterior_density_key_predictors.png
│   ├── prior_sensitivity_plot.png
│   ├── posterior_predictive_check.png
│   └── trace_plots_key_predictors.png
└── presentation/
    └── 03_bayesian_statistical_modeling.pdf
```

## Limitations

- Small sample size (n = 189); hypertension has only 12 cases, resulting in a wide
  credible interval (OR 6.09, 95% CrI 1.65–25.18)
- Retrospective single-centre design from 1986; associations are not causal and may
  not generalise to contemporary populations
- MH acceptance rate (70%) is higher than the theoretical optimum (~23% for high-dimensional
  targets); convergence confirmed visually but a lower proposal variance may improve mixing
- No formal convergence diagnostics (Gelman-Rubin R-hat, effective sample size) were computed
- Unmeasured confounders (prenatal care, nutrition, socioeconomic status) not available
- Results should be interpreted as an applied statistical portfolio analysis demonstrating
  Bayesian methodology

## Tools

| Package | Purpose |
|---|---|
| `MASS` | `birthwt` dataset |
| `tidyverse` | Data wrangling and visualisation |
| `gtsummary` | Descriptive and summary tables |
| `broom` | Frequentist model output |
| `ggplot2` | All figures |

**MCMC sampler:** Implemented from scratch in base R. No `rstanarm`, `brms`, or `rjags`.