# Predictors of Overall Survival in Advanced Lung Cancer
## A Survival Analysis Using the Veterans' Administration Lung Cancer Trial

**Author:** Piyush Kumar
**Dataset:** `veteran` (R `survival` package) — Kalbfleisch & Prentice (1980)
**Software:** R (survival, tidyverse, gtsummary, broom, cowplot, ggplot2)

---

## Overview

This project applies survival analysis methods to a randomised clinical trial dataset to identify predictors of overall survival in patients with advanced inoperable lung cancer. The analysis covers the full survival analysis workflow: Kaplan–Meier estimation, log-rank testing, Cox proportional hazards regression with assumption diagnostics, and a stratified model as a remedy for proportional hazards violations.

## Research Question

Which clinical and demographic factors are independently associated with overall survival in patients with advanced lung cancer, and does the test chemotherapy regimen provide a survival benefit over standard treatment?

## Dataset

| Attribute | Value |
|---|---|
| Source | Veterans' Administration Lung Cancer Trial |
| Reference | Kalbfleisch & Prentice (1980); `veteran` in R `survival` package |
| Study design | Randomised two-arm clinical trial |
| Population | Patients with advanced inoperable lung cancer |
| Total participants | 137 |
| Deaths (events) | 128 (93.4%) |
| Censored | 9 (6.6%) |
| Outcome | Overall survival (days from randomisation to death) |

## Variables

| Variable | Description |
|---|---|
| `time` | Survival time in days |
| `status` | Event indicator (1 = death, 0 = censored) |
| `treatment` | Treatment arm: Standard vs. Test chemotherapy |
| `age` | Age in years |
| `karno` | Karnofsky performance score (0–100; higher = better function) |
| `cell_type` | Histological type: Squamous (ref), Small Cell, Adenocarcinoma, Large Cell |
| `diagtime` | Time since diagnosis in months |
| `prior_therapy` | Prior therapy: No (ref) vs. Yes |
| `performance_group` | Karnofsky score grouped: ≥ 70 (Higher) vs. < 70 (Lower) |

## Methods

### Exploratory Analysis
- Cohort characterisation and baseline descriptive table by treatment arm
- Kaplan–Meier survival curves: overall, by treatment, by performance status, by cell type
- Log-rank tests for all four grouping variables

### Statistical Modelling
- Univariable Cox proportional hazards regression for all predictors
- Multivariable Cox regression — three nested models:
  - Basic: treatment + age
  - Clinical: + Karnofsky score + cell type
  - Full: + time since diagnosis + prior therapy
- Model specification sensitivity analysis (estimate stability across nested models)
- Concordance statistic (C-index) for model discrimination

### Model Diagnostics
- Proportional hazards assumption: Schoenfeld residuals (`cox.zph`)
- Stratified Cox model fitted as remedy for PH violations (karno, cell type)
- PH assumption re-verified on stratified model

## Key Findings

**Cohort**
- Median overall survival: 80 days (95% CI: 52–105)
- Event rate: 93.4% — nearly all participants died during follow-up

**Median survival by subgroup**

| Group | Median (days) | 95% CI |
|---|---|---|
| Standard treatment | 103 | 59–132 |
| Test treatment | 52.5 | 44–95 |
| Higher performance (Karnofsky ≥ 70) | 132 | 111–164 |
| Lower performance (Karnofsky < 70) | 42 | 24–52 |
| Squamous cell | 118 | 82–314 |
| Large Cell | 156 | 105–231 |
| Small Cell | 51 | 25–63 |
| Adenocarcinoma | 51 | 35–92 |

**Log-rank tests**

| Comparison | χ² | p-value |
|---|---|---|
| Treatment group | 0.008 | 0.928 |
| Performance status | 20.26 | < 0.001 |
| Cell type (4 groups) | 25.40 | < 0.001 |
| Prior therapy | 0.501 | 0.479 |

**Primary adjusted hazard ratios — full Cox model**
- **Karnofsky score:** HR 0.97 per unit (95% CI 0.96–0.98, p < 0.001) — strongest predictor; each unit increase associated with 3% lower hazard of death
- **Adenocarcinoma vs. Squamous:** HR 3.31 (95% CI 1.83–5.96, p < 0.001)
- **Small Cell vs. Squamous:** HR 2.37 (95% CI 1.38–4.06, p = 0.002)
- **Large Cell vs. Squamous:** HR 1.49 (95% CI 0.86–2.60, p = 0.156) — not significant
- **Treatment (Test vs. Standard):** HR 1.34 (95% CI 0.89–2.02, p = 0.156) — no significant benefit
- **Age:** HR 0.99 (95% CI 0.97–1.01, p = 0.349) — not significant
- **Time since diagnosis:** HR 1.00 (95% CI 0.98–1.02, p = 0.993) — not significant
- **Prior therapy:** HR 1.07 (95% CI 0.68–1.69, p = 0.758) — not significant

**Model discrimination (C-index)**

| Model | C-index |
|---|---|
| Basic (treatment + age) | 0.514 |
| Clinical (+ Karnofsky + cell type) | 0.738 |
| Full (all predictors) | 0.736 |

**Proportional hazards assumption**
- Violated for Karnofsky score (p < 0.001) and cell type (p = 0.002) in the full model
- Stratified Cox model applied: both variables moved to stratification structure
- PH assumption satisfied in stratified model (GLOBAL p = 0.153)
- Treatment finding unchanged in stratified model: HR 1.09 (95% CI 0.71–1.67, p = 0.701)

## Repository Structure

```
02-survival-analysis-clinical-outcomes/
├── scripts/
│   └── lung_cancer_survival_analysis.R
├── data/
│   └── lung_cancer_survival_dataset.csv
├── outputs/
│   ├── cohort_summary.csv
│   ├── median_survival_overall.csv
│   ├── median_survival_by_treatment.csv
│   ├── median_survival_by_performance.csv
│   ├── median_survival_by_cell_type.csv
│   ├── logrank_test_results.csv
│   ├── cox_model_results.csv
│   ├── cox_stratified_model_results.csv
│   ├── cox_sensitivity_analysis.csv
│   ├── model_concordance.csv
│   ├── proportional_hazards_test.csv
│   ├── ph_test_stratified_model.csv
│   ├── adjusted_survival_curves.csv
│   ├── cox_model_summary.txt
│   ├── cox_stratified_model_summary.txt
│   ├── proportional_hazards_summary.txt
│   ├── cox_regression_table.html
│   └── table1.html
├── figures/
│   ├── overall_kaplan_meier_curve.png
│   ├── km_curve_by_treatment.png
│   ├── km_curve_by_performance_status.png
│   ├── km_curve_by_cell_type.png
│   ├── hazard_ratio_forest_plot.png
│   ├── adjusted_survival_curves.png
│   └── proportional_hazards_diagnostic.png
└── presentation/
    └── 02_survival_analysis_clinical_outcomes.pdf
```

## Limitations

- Small sample size (n = 137); several hazard ratio estimates have wide confidence intervals, particularly for cell type subgroups
- Very high event rate (93.4%) with only 9 censored observations, limiting assessment of long-term survival
- Proportional hazards assumption was violated for Karnofsky score and cell type; stratified model applied as remedy, though time-varying coefficient models could provide additional insight
- Historical trial data (circa 1980); treatment regimens and clinical context do not reflect current oncological practice
- No information on disease stage, molecular markers, or comorbidities — important prognostic factors in contemporary lung cancer research
- Results should be interpreted as an applied statistical portfolio analysis, not a definitive clinical conclusion

## Tools

| Package | Purpose |
|---|---|
| `survival` | `Surv()`, `survfit()`, `coxph()`, `cox.zph()`, `survdiff()` |
| `tidyverse` | Data wrangling and visualisation |
| `gtsummary` | Descriptive and regression tables |
| `broom` | Tidy model output |
| `cowplot` | KM curve + at-risk table panel layout |
| `ggplot2` | All figures |