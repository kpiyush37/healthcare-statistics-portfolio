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

### 3. Bayesian Statistical Modelling: Low Birth Weight Risk Factors

**Objective:** Implement Bayesian logistic regression using a Metropolis-Hastings MCMC sampler built entirely from scratch to estimate posterior distributions of odds ratios for maternal risk factors associated with low birth weight.

**Dataset:** Baystate Medical Center, Springfield MA (1986) — `MASS::birthwt`; n = 189 singleton births, 59 low birth weight cases (31.2%), no missing data

**Methods**
- Bayesian logistic regression via Metropolis-Hastings MCMC (built from scratch in base R — no `rstanarm`, `brms`, or `rjags`)
- Normal(0, 2.5²) weakly informative prior on all coefficients
- 40,000 MCMC iterations, 10,000 burn-in, 30,000 post-burn-in samples retained
- Posterior summary: median OR, 95% credible interval, P(OR > 1)
- Frequentist logistic regression as a validation benchmark
- Prior sensitivity analysis across three prior scales (σ = 1.5, 2.5, 5.0)
- Posterior predictive check for model fit assessment
- MCMC diagnostics: trace plots, acceptance rate, posterior density plots

**Key Findings**
- MCMC acceptance rate: 70.0%; trace plots confirm stationarity and good chain mixing
- **Hypertension:** posterior OR 6.09 (95% CrI 1.65–25.18); P(OR > 1) = 99.5% — strongest predictor; wide CrI reflects small cell count (n = 12)
- **Previous preterm labour:** OR 3.65 (95% CrI 1.55–9.61); P(OR > 1) = 99.9% — most certain finding
- **Black vs. White:** OR 3.13 (95% CrI 1.20–9.21); P(OR > 1) = 99.1%
- **Smoking during pregnancy:** OR 2.16 (95% CrI 1.04–4.74); P(OR > 1) = 98.1%
- **Maternal weight:** OR 0.62 (95% CrI 0.41–0.94); P(OR > 1) = 1.1% — only predictor with credible protective effect
- Age, other race, and uterine irritability showed elevated posterior probabilities (92–97%) but credible intervals crossed OR = 1
- Bayesian and frequentist estimates consistent in direction and magnitude across all predictors
- Findings robust across all three prior scales; posterior predictive check: observed = 59, simulated mean = 59.4 (95% PI: 44–75)

**Tools:** R, `MASS`, `tidyverse`, `gtsummary`, `broom`, `ggplot2`

📁 [`03-bayesian-statistical-modelling/`](./03-bayesian-statistical-modelling/)

---

## Tools and Methods Summary

| Method | Project |
|---|---|
| Survey-weighted regression (`svyglm`) | Project 1 |
| Logistic regression (univariable and multivariable) | Projects 1 and 3 |
| Bayesian logistic regression (MH-MCMC from scratch) | Project 3 |
| AUC / ROC curve | Project 1 |
| Hosmer–Lemeshow calibration test | Project 1 |
| VIF multicollinearity diagnostics | Project 1 |
| Kaplan–Meier survival estimation | Project 2 |
| Log-rank test | Project 2 |
| Cox proportional hazards regression | Project 2 |
| Proportional hazards diagnostics (`cox.zph`) | Project 2 |
| Stratified Cox model | Project 2 |
| Concordance statistic (C-index) | Project 2 |
| Prior sensitivity analysis | Project 3 |
| Posterior predictive check | Project 3 |
| MCMC diagnostics (trace plots, acceptance rate) | Project 3 |
| Sensitivity analysis (multiple model specifications) | Projects 1, 2, and 3 |

**Languages:** R

**Key packages:** `survival`, `survey`, `nhanesA`, `MASS`, `tidyverse`, `gtsummary`, `broom`, `pROC`, `car`, `ResourceSelection`, `cowplot`, `ggplot2`