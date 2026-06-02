# Healthcare Statistics Portfolio

A portfolio of applied statistical analyses in healthcare using real-world and clinical trial datasets. Each project follows a complete analytical workflow from data acquisition and preprocessing through modelling, diagnostics, and interpretation.

---

## Projects

### 1. Real-World Evidence Analysis: Factors Associated with Type 2 Diabetes in NHANES 2017–2018

**Objective:** Identify demographic, socioeconomic, and clinical factors independently associated with prevalent Type 2 Diabetes among U.S. adults using NHANES 2017–2018 data.

**Dataset:** National Health and Nutrition Examination Survey (NHANES) 2017–2018 — n = 4,339 adults (complete-case analytic cohort)

**Methods**
- Cohort construction with participant flow tracking (9,254 → 4,339 after exclusions)
- Exploratory data analysis and missing data assessment
- Descriptive statistics (Table 1 by diabetes status)
- Univariable logistic regression for all predictors
- Multivariable logistic regression — three nested models (Basic, Clinical, Extended)
- Survey-weighted logistic regression (`svyglm`) — primary analysis
- Biomarker-confirmed outcome sensitivity analysis (HbA1c ≥ 6.5% or fasting glucose ≥ 126 mg/dL)
- Model diagnostics: VIF (multicollinearity), Hosmer–Lemeshow test (calibration), AUC/ROC curve (discrimination)

**Key Findings**
- Unweighted diabetes prevalence: 16.2%; survey-weighted nationally representative estimate: 11.7%
- **Hypertension:** OR 2.51 (95% CI 2.06–3.07, p < 0.001) unweighted; OR 2.97 (95% CI 2.32–3.81) weighted — strongest clinical predictor
- **Age:** OR 1.06 per year (95% CI 1.05–1.07, p < 0.001)
- **BMI:** OR 1.08 per unit (95% CI 1.06–1.09, p < 0.001)
- **Non-Hispanic Asian vs. White:** OR 2.52 (95% CI 1.84–3.46, p < 0.001) — elevated independent of BMI
- **Mexican American vs. White:** OR 2.00 (95% CI 1.46–2.74, p < 0.001)
- **Female sex:** OR 0.64 (95% CI 0.53–0.78, p < 0.001) — lower odds than male
- **Education < 9th grade vs. College+:** OR 1.79 (95% CI 1.21–2.64, p = 0.003)
- Smoking and income ratio were not independently significant after full adjustment
- Model AUC = 0.812 (95% CI 0.796–0.828); Hosmer–Lemeshow p = 0.161
- Findings consistent across all model specifications and between self-reported and biomarker-confirmed outcomes

**Tools:** R, `nhanesA`, `tidyverse`, `survey`, `gtsummary`, `broom`, `pROC`, `car`, `ResourceSelection`

📁 [`01-rwe-diabetes-risk-factors/`](./01-rwe-diabetes-risk-factors/)

---

### 2. Survival Analysis: Predictors of Overall Survival in Advanced Lung Cancer

**Objective:** Identify clinical and demographic predictors of overall survival in patients with advanced inoperable lung cancer using a randomised clinical trial dataset.

**Dataset:** Veterans' Administration Lung Cancer Trial — Kalbfleisch & Prentice (1980); `veteran` dataset in R `survival` package — n = 137, 128 deaths (93.4% event rate)

**Methods**
- Kaplan–Meier survival estimation with numbers-at-risk tables
- Log-rank tests (treatment, performance status, cell type, prior therapy)
- Cox proportional hazards regression — three nested models (Basic, Clinical, Full)
- Model specification sensitivity analysis
- Adjusted survival curves with covariates fixed at median/reference values
- Proportional hazards assumption testing via Schoenfeld residuals (`cox.zph`)
- Stratified Cox model as remedy for PH violations (Karnofsky score and cell type)
- Model discrimination via concordance statistic (C-index)

**Key Findings**
- Median overall survival: 80 days (95% CI: 52–105), reflecting the aggressive nature of advanced lung cancer
- **Karnofsky performance score:** HR 0.97 per unit (95% CI 0.96–0.98, p < 0.001) — strongest independent predictor; each unit increase associated with 3% lower hazard of death
- **Adenocarcinoma vs. Squamous:** HR 3.31 (95% CI 1.83–5.96, p < 0.001)
- **Small Cell vs. Squamous:** HR 2.37 (95% CI 1.38–4.06, p = 0.002)
- **Treatment (Test vs. Standard):** HR 1.34 (95% CI 0.89–2.02, p = 0.156) — no significant survival benefit in either unadjusted (log-rank p = 0.928) or adjusted analysis
- Age, time since diagnosis, and prior therapy were not independently associated with survival
- Proportional hazards assumption violated for Karnofsky score (p < 0.001) and cell type (p = 0.002); stratified model confirmed the treatment finding (HR 1.09, p = 0.701) and resolved the violation (GLOBAL p = 0.153)
- C-index: Basic 0.514 → Clinical 0.738 → Full 0.736; adding Karnofsky and cell type accounted for nearly all model discrimination

**Tools:** R, `survival`, `tidyverse`, `gtsummary`, `broom`, `cowplot`, `ggplot2`

📁 [`02-survival-analysis-clinical-outcomes/`](./02-survival-analysis-clinical-outcomes/)

---

### 3. Bayesian Statistical Modelling

*Project under development.*

---

## Tools and Methods Summary

| Method | Project |
|---|---|
| Survey-weighted regression (`svyglm`) | Project 1 |
| Logistic regression (univariable and multivariable) | Project 1 |
| AUC / ROC curve | Project 1 |
| Kaplan–Meier survival estimation | Project 2 |
| Log-rank test | Project 2 |
| Cox proportional hazards regression | Project 2 |
| Proportional hazards diagnostics (`cox.zph`) | Project 2 |
| Stratified Cox model | Project 2 |
| Concordance statistic (C-index) | Project 2 |
| Sensitivity analysis (multiple specifications) | Projects 1 and 2 |
| Model calibration (Hosmer–Lemeshow) | Project 1 |
| VIF multicollinearity diagnostics | Project 1 |

**Languages:** R

**Key packages:** `survival`, `survey`, `nhanesA`, `tidyverse`, `gtsummary`, `broom`, `pROC`, `car`, `ResourceSelection`, `cowplot`, `ggplot2`