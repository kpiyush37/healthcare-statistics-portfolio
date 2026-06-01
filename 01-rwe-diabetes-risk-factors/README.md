# Real-World Evidence Analysis: Factors Associated with Type 2 Diabetes in NHANES 2017–2018

## Overview

This project investigates demographic, socioeconomic, and clinical factors associated with
Type 2 Diabetes using data from the National Health and Nutrition Examination Survey (NHANES)
2017–2018. The analysis follows a complete Real-World Evidence (RWE) workflow including cohort
construction with participant flow tracking, exploratory data analysis, survey-weighted and
unweighted multivariable modelling, sensitivity analyses, and model discrimination and
calibration diagnostics.

## Research Question

Which demographic, socioeconomic, and clinical factors are independently associated with
prevalent Type 2 Diabetes among U.S. adults?

## Dataset

| Attribute | Value |
|---|---|
| Source | NHANES 2017–2018 (CDC/NCHS) |
| Study design | Cross-sectional observational study |
| Population | Adults aged 20 years and older |
| Final analytic cohort | 4,339 participants |
| Diabetes cases | 701 |
| Unweighted prevalence | 16.2% |
| Survey-weighted prevalence | 11.7% (nationally representative) |
| Primary outcome | Self-reported diabetes (DIQ010) |
| Sensitivity outcome | Biomarker-confirmed diabetes (HbA1c ≥ 6.5% or fasting glucose ≥ 126 mg/dL) |

> NHANES uses a complex stratified multistage sampling design. The unweighted prevalence (16.2%)
> reflects the analytic sample only. The survey-weighted estimate (11.7%) accounts for the
> sampling design and better represents the U.S. adult population.

## Methods

### Exploratory Data Analysis
- Cohort characterisation with participant flow tracking
- Missing data assessment and exclusion counts at each step
- BMI distribution by diabetes status
- Age distribution by diabetes status
- Hypertension prevalence by diabetes status

### Statistical Analysis
- Univariable logistic regression for all predictors
- Multivariable logistic regression — three nested models (Basic, Clinical, Extended)
- Survey-weighted logistic regression (`svyglm`) — primary analysis
- Model specification sensitivity analysis (estimate stability across nested models)
- Biomarker-confirmed outcome sensitivity analysis (HbA1c / fasting glucose)
- Variance Inflation Factor (VIF) — multicollinearity
- Hosmer–Lemeshow goodness-of-fit test — calibration
- AUC / ROC curve — discrimination

## Key Findings

**Cohort**
- Final analytic cohort: n = 4,339; 701 diabetes cases
- Unweighted prevalence: 16.2%; survey-weighted prevalence: 11.7%
- Participants with diabetes were older (mean age 63.5 vs. 48.6 years) and had higher BMI (mean 32.5 vs. 29.3)

**Primary adjusted odds ratios — unweighted full model**
- **Age:** OR 1.06 per year (95% CI 1.05–1.07, p < 0.001)
- **BMI:** OR 1.08 per unit (95% CI 1.06–1.09, p < 0.001)
- **Hypertension:** OR 2.51 (95% CI 2.06–3.07, p < 0.001) — strongest clinical predictor
- **Sex (Female vs. Male):** OR 0.64 (95% CI 0.53–0.78, p < 0.001)
- **Non-Hispanic Asian vs. White:** OR 2.52 (95% CI 1.84–3.46, p < 0.001)
- **Mexican American vs. White:** OR 2.00 (95% CI 1.46–2.74, p < 0.001)
- **Education < 9th grade vs. College graduate+:** OR 1.79 (95% CI 1.21–2.64, p = 0.003)
- **Income ratio:** OR 0.99 (95% CI 0.93–1.05, p = 0.71) — not independently significant
- **Smoking (ever vs. never):** OR 1.09 (95% CI 0.89–1.32, p = 0.41) — not independently significant

**Survey-weighted model (primary)**
- Direction and significance of all key associations consistent with unweighted model
- Hypertension OR strengthened to 2.97 (95% CI 2.32–3.81) after applying survey weights

**Model diagnostics**
- AUC = 0.812 (95% CI 0.796–0.828) — good discrimination
- Hosmer–Lemeshow p = 0.161 — acceptable calibration
- VIF < 5 for all predictors — no problematic multicollinearity

**Sensitivity analyses**
- Key associations for age, BMI, and hypertension were stable across all three model specifications
- Self-reported and biomarker-confirmed outcomes produced consistent findings
- Non-Hispanic Black participants reached statistical significance in the biomarker model (OR 1.49, p < 0.001), likely reflecting known underdiagnosis under self-report

## Repository Structure

```
01-rwe-diabetes-risk-factors/
├── scripts/
│   └── diabetes_rwe_analysis.R
├── data/
│   └── analysis_dataset.csv
├── outputs/
│   ├── participant_flow.csv
│   ├── cohort_summary.csv
│   ├── group_summary_by_diabetes.csv
│   ├── univariable_logistic_regression.csv
│   ├── multivariable_logistic_regression.csv
│   ├── weighted_logistic_regression.csv
│   ├── sensitivity_analysis.csv
│   ├── sensitivity_biomarker_confirmed.csv
│   ├── sensitivity_comparison_self_vs_biomarker.csv
│   ├── model_auc.csv
│   ├── regression_table.html
│   └── table1.html
├── figures/
│   ├── bmi_by_diabetes.png
│   ├── age_by_diabetes.png
│   ├── diabetes_by_hypertension.png
│   ├── forest_plot_adjusted_odds_ratios.png
│   └── roc_curve.png
└── presentation/
    └── Factors_Associated_with_Type_2_Diabetes_in_NHANES.pdf
```

## Limitations

- Cross-sectional design: associations cannot be interpreted as causal
- Primary outcome is self-reported; biomarker-confirmed sensitivity analysis provided as a validity check
- Complete-case analysis: largest exclusion was missing income ratio (n = 653); estimates may be biased if data are not missing completely at random
- Potential residual confounding from unmeasured variables (physical activity, diet, family history)

## Tools

| Package | Purpose |
|---|---|
| `nhanesA` | NHANES data download |
| `tidyverse` | Data wrangling and visualisation |
| `survey` | Complex survey design and weighted regression |
| `gtsummary` | Descriptive and regression tables |
| `broom` | Tidy model output |
| `ggplot2` | Figures |
| `pROC` | AUC and ROC curve |
| `car` | VIF multicollinearity diagnostics |
| `ResourceSelection` | Hosmer–Lemeshow goodness-of-fit test |