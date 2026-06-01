# ============================================================
# Real-World Evidence Analysis:
# Factors Associated with Type 2 Diabetes in NHANES 2017-2018
#
# Author:      Piyush Kumar
# Dataset:     NHANES 2017-2018 (CDC/NCHS)
# Design:      Cross-sectional observational study
# Outcome:     Prevalent Type 2 Diabetes
#              Primary: self-reported (DIQ010)
#              Sensitivity: biomarker-confirmed (HbA1c >= 6.5% or fasting glucose >= 126 mg/dL)
# Methods:     Survey-weighted logistic regression, sensitivity analyses,
#              model discrimination and calibration diagnostics
# ============================================================

rm(list = ls())

# ============================================================
# SECTION 1: PACKAGE LOADING
# ============================================================

packages <- c(
  "nhanesA",          # NHANES data download
  "tidyverse",        # Data wrangling and visualisation
  "gtsummary",        # Formatted descriptive and regression tables
  "broom",            # Tidy model output
  "car",              # Variance Inflation Factor (VIF)
  "ResourceSelection",# Hosmer-Lemeshow goodness-of-fit test
  "survey",           # Complex survey design and weighted regression
  "pROC",             # AUC and ROC curve
  "ggtext"            # Enhanced ggplot2 text formatting
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

invisible(lapply(packages, install_if_missing))

library(nhanesA)
library(tidyverse)
library(gtsummary)
library(broom)
library(car)
library(ResourceSelection)
library(survey)
library(pROC)
library(ggtext)

# ============================================================
# SECTION 2: PROJECT STRUCTURE
# ============================================================

dir.create("data",         showWarnings = FALSE)
dir.create("outputs",      showWarnings = FALSE)
dir.create("figures",      showWarnings = FALSE)
dir.create("report",       showWarnings = FALSE)
dir.create("presentation", showWarnings = FALSE)

# ============================================================
# SECTION 3: DATA LOADING
# ============================================================

# Core NHANES 2017-2018 modules
demo <- nhanes("DEMO_J")   # Demographics (includes survey design variables)
bmx  <- nhanes("BMX_J")    # Body measurements
diq  <- nhanes("DIQ_J")    # Diabetes questionnaire
bpq  <- nhanes("BPQ_J")    # Blood pressure questionnaire
smq  <- nhanes("SMQ_J")    # Smoking questionnaire

# Biomarker modules for sensitivity analysis
ghb  <- nhanes("GHB_J")    # Glycohemoglobin (HbA1c)
glu  <- nhanes("GLU_J")    # Fasting plasma glucose

# ============================================================
# SECTION 4: VARIABLE SELECTION
# ============================================================

# Survey design variables required for complex sample analysis:
#   WTMEC2YR — 2-year MEC examination weight
#   SDMVPSU  — Primary sampling unit
#   SDMVSTRA — Stratum identifier

demo_clean <- demo %>%
  select(
    SEQN,
    age            = RIDAGEYR,
    sex            = RIAGENDR,
    race_ethnicity = RIDRETH3,
    education      = DMDEDUC2,
    income_ratio   = INDFMPIR,
    wt_mec         = WTMEC2YR,
    psu            = SDMVPSU,
    strata         = SDMVSTRA
  )

bmx_clean <- bmx %>%
  select(SEQN, bmi = BMXBMI)

diq_clean <- diq %>%
  select(SEQN, diabetes_raw = DIQ010)

bpq_clean <- bpq %>%
  select(SEQN, hypertension_raw = BPQ020)

smq_clean <- smq %>%
  select(SEQN, smoking_raw = SMQ020)

ghb_clean <- ghb %>%
  select(SEQN, hba1c = LBXGH)

glu_clean <- glu %>%
  select(SEQN, fasting_glucose = LBXGLU)

# ============================================================
# SECTION 5: DATA MERGING
# ============================================================

nhanes_df <- demo_clean %>%
  left_join(bmx_clean, by = "SEQN") %>%
  left_join(diq_clean, by = "SEQN") %>%
  left_join(bpq_clean, by = "SEQN") %>%
  left_join(smq_clean, by = "SEQN") %>%
  left_join(ghb_clean, by = "SEQN") %>%
  left_join(glu_clean, by = "SEQN")

# ============================================================
# SECTION 6: COHORT CONSTRUCTION AND PARTICIPANT FLOW
# ============================================================

# --- 6a. Restrict to adults aged 20 and older ---
n_total  <- nrow(nhanes_df)
adults_df <- nhanes_df %>% filter(age >= 20)
n_adults  <- nrow(adults_df)

# --- 6b. Recode variables ---

adults_df <- adults_df %>%
  mutate(
    
    # Primary outcome: self-reported diabetes (binary)
    diabetes_num = case_when(
      diabetes_raw == "Yes" ~ 1,
      diabetes_raw == "No"  ~ 0,
      TRUE                  ~ NA_real_
    ),
    diabetes = factor(
      diabetes_num,
      levels = c(0, 1),
      labels = c("No Diabetes", "Diabetes")
    ),
    
    # Categorical predictors
    sex            = factor(sex),
    race_ethnicity = factor(race_ethnicity),
    education      = factor(education),
    
    hypertension = case_when(
      hypertension_raw == "Yes" ~ "Yes",
      hypertension_raw == "No"  ~ "No",
      TRUE                      ~ NA_character_
    ) %>% factor(),
    
    smoking = case_when(
      smoking_raw == "Yes" ~ "Ever smoker",
      smoking_raw == "No"  ~ "Never smoker",
      TRUE                 ~ NA_character_
    ) %>% factor(),
    
    # Sensitivity outcome: biomarker-confirmed diabetes
    # Definition: HbA1c >= 6.5% OR fasting glucose >= 126 mg/dL
    # OR self-reported diabetes (ADA 2023 diagnostic criteria)
    diabetes_biomarker = case_when(
      hba1c >= 6.5 | fasting_glucose >= 126 | diabetes_num == 1 ~ 1,
      !is.na(hba1c) | !is.na(fasting_glucose)                   ~ 0,
      TRUE                                                       ~ NA_real_
    )
  )

# --- 6c. Apply exclusion criteria and track participant flow ---

excl_education   <- adults_df %>% filter(education %in% c("Refused", "Don't know"))
n_excl_edu       <- nrow(excl_education)
after_edu        <- adults_df %>% filter(!education %in% c("Refused", "Don't know"))

n_excl_diabetes  <- sum(is.na(after_edu$diabetes))
after_diabetes   <- after_edu %>% filter(!is.na(diabetes))

n_excl_bmi       <- sum(is.na(after_diabetes$bmi))
after_bmi        <- after_diabetes %>% filter(!is.na(bmi))

n_excl_income    <- sum(is.na(after_bmi$income_ratio))
after_income     <- after_bmi %>% filter(!is.na(income_ratio))

n_excl_htn       <- sum(is.na(after_income$hypertension))
after_htn        <- after_income %>% filter(!is.na(hypertension))

n_excl_smoking   <- sum(is.na(after_htn$smoking))
after_smoking    <- after_htn %>% filter(!is.na(smoking))

# --- 6d. Final analytic dataset with reference levels ---

analysis_df <- after_smoking %>%
  mutate(
    sex            = relevel(sex,            ref = "Male"),
    race_ethnicity = relevel(race_ethnicity, ref = "Non-Hispanic White"),
    education      = relevel(education,      ref = "College graduate or above"),
    hypertension   = relevel(hypertension,   ref = "No"),
    smoking        = relevel(smoking,        ref = "Never smoker")
  )

n_final       <- nrow(analysis_df)
n_diabetes    <- sum(analysis_df$diabetes_num == 1)
n_no_diabetes <- sum(analysis_df$diabetes_num == 0)

# --- 6e. Print and save participant flow ---

flow_table <- tibble(
  Step = c(
    "1. Full NHANES 2017-2018",
    "2. Adults aged >= 20 years",
    "3. Excluded: education refused/unknown",
    "4. Excluded: missing diabetes status",
    "5. Excluded: missing BMI",
    "6. Excluded: missing income ratio",
    "7. Excluded: missing hypertension",
    "8. Excluded: missing smoking status",
    "   Final analytic cohort",
    "   — Diabetes cases",
    "   — No diabetes"
  ),
  N = c(
    n_total, n_adults,
    -n_excl_edu, -n_excl_diabetes, -n_excl_bmi,
    -n_excl_income, -n_excl_htn, -n_excl_smoking,
    n_final, n_diabetes, n_no_diabetes
  )
)

cat("\n===== PARTICIPANT FLOW =====\n")
print(flow_table, n = Inf)
write_csv(flow_table, "outputs/participant_flow.csv")

# ============================================================
# SECTION 7: SURVEY DESIGN OBJECT
# ============================================================

# NHANES uses a stratified multistage probability sample.
# Unweighted analyses produce biased prevalence estimates and
# incorrect standard errors. svydesign() applies the complex
# sampling structure using the MEC examination weight (WTMEC2YR).

nhanes_design <- svydesign(
  id      = ~psu,
  strata  = ~strata,
  weights = ~wt_mec,
  data    = analysis_df,
  nest    = TRUE
)

# Weighted diabetes prevalence (nationally representative estimate)
weighted_prevalence <- svymean(~diabetes_num, design = nhanes_design, na.rm = TRUE)
cat("\n===== WEIGHTED DIABETES PREVALENCE =====\n")
print(weighted_prevalence)
cat("95% CI:", confint(weighted_prevalence), "\n")

# ============================================================
# SECTION 8: DESCRIPTIVE TABLE (TABLE 1)
# ============================================================

table1_unweighted <- analysis_df %>%
  select(
    -diabetes_num, -diabetes_biomarker,
    -wt_mec, -psu, -strata, -hba1c, -fasting_glucose
  ) %>%
  tbl_summary(
    by        = diabetes,
    statistic = list(
      all_continuous()  ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "no"
  ) %>%
  add_p(
    test = list(
      all_continuous()  ~ "t.test",
      all_categorical() ~ "chisq.test"
    )
  ) %>%
  add_overall()

# ============================================================
# SECTION 9: UNIVARIABLE LOGISTIC REGRESSION
# ============================================================

predictors <- c(
  "age", "sex", "bmi", "race_ethnicity",
  "education", "income_ratio", "hypertension", "smoking"
)

univariable_results <- predictors %>%
  map_df(function(var) {
    model <- glm(
      as.formula(paste("diabetes_num ~", var)),
      data   = analysis_df,
      family = binomial
    )
    tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
      mutate(predictor = var)
  })

write_csv(univariable_results, "outputs/univariable_logistic_regression.csv")

# ============================================================
# SECTION 10: MULTIVARIABLE LOGISTIC REGRESSION
# ============================================================

# Three nested models to assess estimate stability:
#   Basic    — core demographic and anthropometric predictors
#   Clinical — adds hypertension and smoking
#   Extended — full model including socioeconomic variables

model_basic <- glm(
  diabetes_num ~ age + sex + bmi,
  data   = analysis_df,
  family = binomial
)

model_clinical <- glm(
  diabetes_num ~ age + sex + bmi + hypertension + smoking,
  data   = analysis_df,
  family = binomial
)

model_full <- glm(
  diabetes_num ~ age + sex + bmi + race_ethnicity +
    education + income_ratio + hypertension + smoking,
  data   = analysis_df,
  family = binomial
)

model_results    <- tidy(model_full, exponentiate = TRUE, conf.int = TRUE)
regression_table <- tbl_regression(model_full, exponentiate = TRUE)

# Sensitivity analysis: estimates across all three model specifications
sensitivity_results <- list(
  Basic    = model_basic,
  Clinical = model_clinical,
  Extended = model_full
) %>%
  map_df(~ tidy(.x, exponentiate = TRUE, conf.int = TRUE), .id = "model")

write_csv(model_results,       "outputs/multivariable_logistic_regression.csv")
write_csv(sensitivity_results, "outputs/sensitivity_analysis.csv")

# ============================================================
# SECTION 11: SURVEY-WEIGHTED LOGISTIC REGRESSION (PRIMARY MODEL)
# ============================================================

# svyglm() accounts for the complex NHANES sampling design in
# standard error estimation. quasibinomial() is recommended for
# binary outcomes with survey-weighted regression.

model_weighted <- svyglm(
  diabetes_num ~ age + sex + bmi + race_ethnicity +
    education + income_ratio + hypertension + smoking,
  design = nhanes_design,
  family = quasibinomial()
)

weighted_results <- tidy(model_weighted, exponentiate = TRUE, conf.int = TRUE)

# Compute Wald-based CIs and p-values from the model's standard errors
# (standard broom::tidy() does not populate these for svyglm objects)
weighted_results <- weighted_results %>%
  mutate(
    conf.low  = estimate * exp(-1.96 * std.error),
    conf.high = estimate * exp(+1.96 * std.error),
    p.value   = 2 * pnorm(-abs(statistic))
  )

cat("\n===== WEIGHTED MODEL RESULTS =====\n")
print(weighted_results)
write_csv(weighted_results, "outputs/weighted_logistic_regression.csv")

# ============================================================
# SECTION 12: MODEL DIAGNOSTICS
# ============================================================

# --- 12a. Multicollinearity: Variance Inflation Factor ---
vif_results <- car::vif(model_full)
cat("\n===== VIF RESULTS =====\n")
print(vif_results)

# --- 12b. Calibration: Hosmer-Lemeshow goodness-of-fit test ---
hosmer_lemeshow <- hoslem.test(
  analysis_df$diabetes_num,
  fitted(model_full),
  g = 10
)
cat("\n===== HOSMER-LEMESHOW TEST =====\n")
print(hosmer_lemeshow)

# --- 12c. Discrimination: AUC and ROC curve ---
roc_obj   <- roc(
  response  = analysis_df$diabetes_num,
  predictor = fitted(model_full),
  quiet     = TRUE
)
auc_value <- auc(roc_obj)
auc_ci    <- ci.auc(roc_obj)

cat(sprintf(
  "\n===== MODEL DISCRIMINATION =====\nAUC = %.3f (95%% CI: %.3f – %.3f)\n",
  auc_value, auc_ci[1], auc_ci[3]
))

write_csv(
  tibble(auc = as.numeric(auc_value), ci_lower = auc_ci[1], ci_upper = auc_ci[3]),
  "outputs/model_auc.csv"
)

# ============================================================
# SECTION 13: BIOMARKER SENSITIVITY ANALYSIS
# ============================================================

# Validates self-reported findings using objective diagnostic criteria:
# HbA1c >= 6.5% OR fasting glucose >= 126 mg/dL (ADA 2023)
# Self-reported cases are also included to avoid false negatives.

n_biomarker <- sum(!is.na(analysis_df$diabetes_biomarker))
n_bio_cases <- sum(analysis_df$diabetes_biomarker == 1, na.rm = TRUE)

cat(sprintf(
  "\n===== BIOMARKER SENSITIVITY ANALYSIS =====\n%d participants with biomarker data | %d cases (%.1f%%)\n",
  n_biomarker, n_bio_cases, 100 * n_bio_cases / n_biomarker
))

analysis_bio <- analysis_df %>% filter(!is.na(diabetes_biomarker))

model_biomarker <- glm(
  diabetes_biomarker ~ age + sex + bmi + race_ethnicity +
    education + income_ratio + hypertension + smoking,
  data   = analysis_bio,
  family = binomial
)

biomarker_results <- tidy(model_biomarker, exponentiate = TRUE, conf.int = TRUE)
write_csv(biomarker_results, "outputs/sensitivity_biomarker_confirmed.csv")

# Side-by-side comparison of key estimates: self-reported vs biomarker-confirmed
comparison <- bind_rows(
  model_results     %>% mutate(outcome = "Self-reported"),
  biomarker_results %>% mutate(outcome = "Biomarker-confirmed")
) %>%
  filter(term != "(Intercept)") %>%
  select(outcome, term, estimate, conf.low, conf.high, p.value)

write_csv(comparison, "outputs/sensitivity_comparison_self_vs_biomarker.csv")

# ============================================================
# SECTION 14: VISUALISATIONS
# ============================================================

# --- 14a. BMI by diabetes status ---
bmi_plot <- ggplot(analysis_df, aes(x = diabetes, y = bmi, fill = diabetes)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
  scale_fill_manual(values = c("#74add1", "#f46d43")) +
  labs(
    title    = "BMI Distribution by Diabetes Status",
    subtitle = "NHANES 2017–2018 | Unweighted analytic cohort (n = 4,339)",
    x        = "Diabetes Status",
    y        = "Body Mass Index (kg/m²)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

# --- 14b. Age distribution by diabetes status ---
age_plot <- ggplot(analysis_df, aes(x = age, fill = diabetes)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("#74add1", "#f46d43")) +
  labs(
    title    = "Age Distribution by Diabetes Status",
    subtitle = "NHANES 2017–2018 | Unweighted analytic cohort",
    x        = "Age (years)",
    y        = "Density",
    fill     = "Diabetes Status"
  ) +
  theme_bw(base_size = 12)

# --- 14c. Diabetes proportion by hypertension status ---
hypertension_plot <- ggplot(analysis_df, aes(x = hypertension, fill = diabetes)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("#74add1", "#f46d43")) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title    = "Diabetes Proportion by Hypertension Status",
    subtitle = "NHANES 2017–2018 | Unweighted analytic cohort",
    x        = "Hypertension (self-reported)",
    y        = "Proportion",
    fill     = "Diabetes Status"
  ) +
  theme_bw(base_size = 12)

# --- 14d. Forest plot: adjusted odds ratios ---

# Human-readable labels for all model terms
label_map <- c(
  "age"                                                         = "Age (per year)",
  "sexFemale"                                                   = "Sex: Female (ref: Male)",
  "bmi"                                                         = "BMI (per unit)",
  "race_ethnicityMexican American"                              = "Race: Mexican American",
  "race_ethnicityOther Hispanic"                                = "Race: Other Hispanic",
  "race_ethnicityNon-Hispanic Black"                            = "Race: Non-Hispanic Black",
  "race_ethnicityNon-Hispanic Asian"                            = "Race: Non-Hispanic Asian",
  "race_ethnicityOther Race - Including Multi-Racial"           = "Race: Other/Multi-racial",
  "educationLess than 9th grade"                                = "Education: <9th grade",
  "education9-11th grade (Includes 12th grade with no diploma)" = "Education: 9–11th grade",
  "educationHigh school graduate/GED or equivalent"             = "Education: High school/GED",
  "educationSome college or AA degree"                          = "Education: Some college/AA",
  "income_ratio"                                                = "Income-to-poverty ratio",
  "hypertensionYes"                                             = "Hypertension: Yes (ref: No)",
  "smokingEver smoker"                                          = "Smoking: Ever (ref: Never)"
)

# Domain groupings for panel faceting
group_map <- c(
  "age"                                                         = "Demographics",
  "sexFemale"                                                   = "Demographics",
  "bmi"                                                         = "Clinical",
  "race_ethnicityMexican American"                              = "Race/Ethnicity\n(ref: Non-Hispanic White)",
  "race_ethnicityOther Hispanic"                                = "Race/Ethnicity\n(ref: Non-Hispanic White)",
  "race_ethnicityNon-Hispanic Black"                            = "Race/Ethnicity\n(ref: Non-Hispanic White)",
  "race_ethnicityNon-Hispanic Asian"                            = "Race/Ethnicity\n(ref: Non-Hispanic White)",
  "race_ethnicityOther Race - Including Multi-Racial"           = "Race/Ethnicity\n(ref: Non-Hispanic White)",
  "educationLess than 9th grade"                                = "Education\n(ref: College graduate+)",
  "education9-11th grade (Includes 12th grade with no diploma)" = "Education\n(ref: College graduate+)",
  "educationHigh school graduate/GED or equivalent"             = "Education\n(ref: College graduate+)",
  "educationSome college or AA degree"                          = "Education\n(ref: College graduate+)",
  "income_ratio"                                                = "Socioeconomic",
  "hypertensionYes"                                             = "Clinical",
  "smokingEver smoker"                                          = "Lifestyle"
)

forest_data <- model_results %>%
  filter(term != "(Intercept)") %>%
  mutate(
    label    = label_map[term],
    group    = factor(group_map[term], levels = c(
      "Demographics",
      "Clinical",
      "Race/Ethnicity\n(ref: Non-Hispanic White)",
      "Education\n(ref: College graduate+)",
      "Socioeconomic",
      "Lifestyle"
    )),
    label    = fct_reorder(label, estimate),
    sig      = p.value < 0.05,
    or_label = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high)
  )

forest_plot <- ggplot(forest_data, aes(x = estimate, y = label, colour = sig)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.25) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_text(
    aes(x = max(conf.high, na.rm = TRUE) * 1.05, label = or_label),
    hjust = 0, size = 3, colour = "black"
  ) +
  scale_x_log10(
    breaks = c(0.5, 0.75, 1, 1.5, 2, 3, 4),
    labels = c("0.5", "0.75", "1.0", "1.5", "2.0", "3.0", "4.0")
  ) +
  scale_colour_manual(
    values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"),
    labels = c("TRUE" = "p < 0.05", "FALSE" = "p \u2265 0.05"),
    name   = "Significance"
  ) +
  facet_wrap(~group, scales = "free_y", ncol = 1) +
  labs(
    title    = "Adjusted Odds Ratios for Type 2 Diabetes",
    subtitle = "Full multivariable logistic regression | NHANES 2017\u20132018 (n = 4,339)\nOR and 95% CI shown; reference categories in parentheses",
    x        = "Odds Ratio (log scale)",
    y        = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text       = element_text(face = "bold", size = 9),
    legend.position  = "bottom",
    plot.subtitle    = element_text(colour = "grey40", size = 9),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(clip = "off") +
  expand_limits(x = max(forest_data$conf.high, na.rm = TRUE) * 1.6)

# --- 14e. ROC curve ---
roc_plot <- ggroc(roc_obj, colour = "#2c7bb6", size = 1) +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", colour = "grey50") +
  annotate(
    "text", x = 0.35, y = 0.15,
    label = sprintf("AUC = %.3f\n95%% CI: %.3f \u2013 %.3f", auc_value, auc_ci[1], auc_ci[3]),
    size = 4, hjust = 0
  ) +
  labs(
    title    = "ROC Curve \u2014 Full Multivariable Model",
    subtitle = "Outcome: Self-reported Type 2 Diabetes | NHANES 2017\u20132018",
    x        = "1 \u2013 Specificity (False Positive Rate)",
    y        = "Sensitivity (True Positive Rate)"
  ) +
  theme_bw(base_size = 12)

# ============================================================
# SECTION 15: SAVE ALL OUTPUTS
# ============================================================

write_csv(analysis_df, "data/analysis_dataset.csv")

write_csv(
  analysis_df %>%
    summarise(
      n                              = n(),
      diabetes_cases                 = sum(diabetes_num == 1),
      diabetes_prevalence_unweighted = mean(diabetes_num),
      weighted_prevalence            = as.numeric(weighted_prevalence)[1],
      mean_age                       = mean(age),
      mean_bmi                       = mean(bmi),
      mean_income_ratio              = mean(income_ratio)
    ),
  "outputs/cohort_summary.csv"
)

write_csv(
  analysis_df %>%
    group_by(diabetes) %>%
    summarise(
      n                   = n(),
      mean_age            = mean(age),
      sd_age              = sd(age),
      median_age          = median(age),
      mean_bmi            = mean(bmi),
      sd_bmi              = sd(bmi),
      median_bmi          = median(bmi),
      mean_income_ratio   = mean(income_ratio),
      sd_income_ratio     = sd(income_ratio),
      median_income_ratio = median(income_ratio),
      .groups             = "drop"
    ),
  "outputs/group_summary_by_diabetes.csv"
)

ggsave("figures/bmi_by_diabetes.png",
       bmi_plot,          width = 8,  height = 5)
ggsave("figures/age_by_diabetes.png",
       age_plot,          width = 8,  height = 5)
ggsave("figures/diabetes_by_hypertension.png",
       hypertension_plot, width = 8,  height = 5)
ggsave("figures/forest_plot_adjusted_odds_ratios.png",
       forest_plot,       width = 11, height = 10)
ggsave("figures/roc_curve.png",
       roc_plot,          width = 7,  height = 6)

# ============================================================
# SECTION 16: CONSOLE SUMMARY
# ============================================================

cat("\n========================================\n")
cat("KEY RESULTS SUMMARY\n")
cat("========================================\n")
cat(sprintf("Analytic cohort:     n = %d\n", n_final))
cat(sprintf("Diabetes cases:      n = %d (%.1f%% unweighted)\n",
            n_diabetes, 100 * n_diabetes / n_final))
cat(sprintf("Weighted prevalence: %.1f%% (95%% CI: %.1f%%\u2013%.1f%%)\n",
            100 * as.numeric(weighted_prevalence)[1],
            100 * confint(weighted_prevalence)[1],
            100 * confint(weighted_prevalence)[2]))
cat(sprintf("Model AUC:           %.3f (95%% CI: %.3f\u2013%.3f)\n",
            auc_value, auc_ci[1], auc_ci[3]))
cat(sprintf("Hosmer-Lemeshow:     p = %.3f\n", hosmer_lemeshow$p.value))
cat("========================================\n\n")

cat("Adjusted OR estimates — Full Unweighted Model:\n")
model_results %>%
  filter(term != "(Intercept)") %>%
  mutate(across(c(estimate, conf.low, conf.high, p.value), ~ round(.x, 3))) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  print(n = Inf)

cat("\nAdjusted OR estimates — Survey-Weighted Model:\n")
weighted_results %>%
  filter(term != "(Intercept)") %>%
  mutate(across(c(estimate, conf.low, conf.high, p.value), ~ round(.x, 3))) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  print(n = Inf)

# Display formatted tables
table1_unweighted
regression_table

