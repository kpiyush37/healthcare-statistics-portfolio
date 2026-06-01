# ============================================================
# Real-World Evidence Analysis:
# Factors Associated with Type 2 Diabetes in NHANES 2017-2018
# ============================================================

rm(list = ls())

# -----------------------------
# Load required packages
# -----------------------------

packages <- c(
  "nhanesA",
  "tidyverse",
  "gtsummary",
  "broom",
  "car",
  "ResourceSelection"
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

invisible(lapply(packages, install_if_missing))

library(nhanesA)
library(tidyverse)
library(gtsummary)
library(broom)
library(car)
library(ResourceSelection)

# -----------------------------
# Create project folders
# -----------------------------

dir.create("data", showWarnings = FALSE)
dir.create("outputs", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("report", showWarnings = FALSE)
dir.create("presentation", showWarnings = FALSE)

# -----------------------------
# Load NHANES 2017-2018 datasets
# -----------------------------

demo <- nhanes("DEMO_J")
bmx  <- nhanes("BMX_J")
diq  <- nhanes("DIQ_J")
bpq  <- nhanes("BPQ_J")
smq  <- nhanes("SMQ_J")

# -----------------------------
# Select relevant variables
# -----------------------------

demo_clean <- demo %>%
  select(
    SEQN,
    age = RIDAGEYR,
    sex = RIAGENDR,
    race_ethnicity = RIDRETH3,
    education = DMDEDUC2,
    income_ratio = INDFMPIR
  )

bmx_clean <- bmx %>%
  select(
    SEQN,
    bmi = BMXBMI
  )

diq_clean <- diq %>%
  select(
    SEQN,
    diabetes_raw = DIQ010
  )

bpq_clean <- bpq %>%
  select(
    SEQN,
    hypertension_raw = BPQ020
  )

smq_clean <- smq %>%
  select(
    SEQN,
    smoking_raw = SMQ020
  )

# -----------------------------
# Merge datasets by participant ID
# -----------------------------

nhanes_df <- demo_clean %>%
  left_join(bmx_clean, by = "SEQN") %>%
  left_join(diq_clean, by = "SEQN") %>%
  left_join(bpq_clean, by = "SEQN") %>%
  left_join(smq_clean, by = "SEQN")

# -----------------------------
# Define adult study population
# -----------------------------

study_df <- nhanes_df %>%
  filter(age >= 20)

# -----------------------------
# Recode variables
# -----------------------------

study_df <- study_df %>%
  mutate(
    diabetes_num = case_when(
      diabetes_raw == "Yes" ~ 1,
      diabetes_raw == "No" ~ 0,
      TRUE ~ NA_real_
    ),
    diabetes = factor(
      diabetes_num,
      levels = c(0, 1),
      labels = c("No Diabetes", "Diabetes")
    ),
    sex = factor(sex),
    race_ethnicity = factor(race_ethnicity),
    education = factor(education),
    hypertension = case_when(
      hypertension_raw == "Yes" ~ "Yes",
      hypertension_raw == "No" ~ "No",
      TRUE ~ NA_character_
    ),
    smoking = case_when(
      smoking_raw == "Yes" ~ "Ever smoker",
      smoking_raw == "No" ~ "Never smoker",
      TRUE ~ NA_character_
    ),
    hypertension = factor(hypertension),
    smoking = factor(smoking)
  )

# -----------------------------
# Create final analytic dataset
# -----------------------------

analysis_df <- study_df %>%
  filter(
    !education %in% c("Refused", "Don't know"),
    !is.na(diabetes),
    !is.na(bmi),
    !is.na(income_ratio),
    !is.na(hypertension)
  ) %>%
  select(
    diabetes,
    diabetes_num,
    age,
    sex,
    bmi,
    race_ethnicity,
    education,
    income_ratio,
    hypertension,
    smoking
  ) %>%
  drop_na() %>%
  mutate(
    sex = relevel(sex, ref = "Male"),
    race_ethnicity = relevel(race_ethnicity, ref = "Non-Hispanic White"),
    education = relevel(education, ref = "College graduate or above"),
    hypertension = relevel(hypertension, ref = "No"),
    smoking = relevel(smoking, ref = "Never smoker")
  )

# -----------------------------
# Cohort summaries
# -----------------------------

cohort_summary <- analysis_df %>%
  summarise(
    n = n(),
    diabetes_cases = sum(diabetes_num == 1),
    diabetes_prevalence = mean(diabetes_num),
    mean_age = mean(age),
    mean_bmi = mean(bmi),
    mean_income_ratio = mean(income_ratio)
  )

group_summary <- analysis_df %>%
  group_by(diabetes) %>%
  summarise(
    n = n(),
    mean_age = mean(age),
    sd_age = sd(age),
    median_age = median(age),
    mean_bmi = mean(bmi),
    sd_bmi = sd(bmi),
    median_bmi = median(bmi),
    mean_income_ratio = mean(income_ratio),
    sd_income_ratio = sd(income_ratio),
    median_income_ratio = median(income_ratio),
    .groups = "drop"
  )

# -----------------------------
# Descriptive Table 1
# -----------------------------

table1 <- analysis_df %>%
  select(-diabetes_num) %>%
  tbl_summary(
    by = diabetes,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "no"
  ) %>%
  add_p(
    test = list(
      all_continuous() ~ "t.test",
      all_categorical() ~ "chisq.test"
    )
  ) %>%
  add_overall()

# -----------------------------
# Univariable logistic regression
# -----------------------------

predictors <- c(
  "age",
  "sex",
  "bmi",
  "race_ethnicity",
  "education",
  "income_ratio",
  "hypertension",
  "smoking"
)

univariable_results <- predictors %>%
  map_df(function(var) {
    model <- glm(
      as.formula(paste("diabetes_num ~", var)),
      data = analysis_df,
      family = binomial
    )
    
    tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
      mutate(predictor = var)
  })

# -----------------------------
# Multivariable logistic regression
# -----------------------------

model_basic <- glm(
  diabetes_num ~ age + sex + bmi,
  data = analysis_df,
  family = binomial
)

model_clinical <- glm(
  diabetes_num ~ age + sex + bmi + hypertension + smoking,
  data = analysis_df,
  family = binomial
)

model_full <- glm(
  diabetes_num ~ age + sex + bmi + race_ethnicity +
    education + income_ratio + hypertension + smoking,
  data = analysis_df,
  family = binomial
)

model_results <- tidy(
  model_full,
  exponentiate = TRUE,
  conf.int = TRUE
)

regression_table <- tbl_regression(
  model_full,
  exponentiate = TRUE
)

# -----------------------------
# Sensitivity analysis
# -----------------------------

sensitivity_results <- list(
  Basic = model_basic,
  Clinical = model_clinical,
  Extended = model_full
) %>%
  map_df(
    ~ tidy(.x, exponentiate = TRUE, conf.int = TRUE),
    .id = "model"
  )

# -----------------------------
# Model diagnostics
# -----------------------------

vif_results <- car::vif(model_full)

hosmer_lemeshow <- hoslem.test(
  analysis_df$diabetes_num,
  fitted(model_full),
  g = 10
)

# -----------------------------
# Visualizations
# -----------------------------

bmi_plot <- ggplot(analysis_df, aes(x = diabetes, y = bmi)) +
  geom_boxplot() +
  labs(
    title = "BMI Distribution by Diabetes Status",
    x = "Diabetes Status",
    y = "Body Mass Index"
  )

age_plot <- ggplot(analysis_df, aes(x = age, fill = diabetes)) +
  geom_density(alpha = 0.4) +
  labs(
    title = "Age Distribution by Diabetes Status",
    x = "Age",
    y = "Density"
  )

hypertension_plot <- ggplot(analysis_df, aes(x = hypertension, fill = diabetes)) +
  geom_bar(position = "fill") +
  labs(
    title = "Diabetes Proportion by Hypertension Status",
    x = "Hypertension",
    y = "Proportion"
  )

forest_plot_data <- model_results %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term = str_replace_all(term, "_", " "),
    term = fct_reorder(term, estimate)
  )

forest_plot <- ggplot(
  forest_plot_data,
  aes(x = estimate, y = term)
) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  labs(
    title = "Adjusted Odds Ratios for Diabetes",
    x = "Odds Ratio (log scale)",
    y = NULL
  )

# -----------------------------
# Save outputs
# -----------------------------

write_csv(analysis_df, "data/analysis_dataset.csv")
write_csv(cohort_summary, "outputs/cohort_summary.csv")
write_csv(group_summary, "outputs/group_summary_by_diabetes.csv")
write_csv(univariable_results, "outputs/univariable_logistic_regression.csv")
write_csv(model_results, "outputs/multivariable_logistic_regression.csv")
write_csv(sensitivity_results, "outputs/sensitivity_analysis.csv")

ggsave("figures/bmi_by_diabetes.png", bmi_plot, width = 8, height = 5)
ggsave("figures/age_by_diabetes.png", age_plot, width = 8, height = 5)
ggsave("figures/diabetes_by_hypertension.png", hypertension_plot, width = 8, height = 5)
ggsave("figures/forest_plot_adjusted_odds_ratios.png", forest_plot, width = 9, height = 6)

# -----------------------------
# Display key outputs
# -----------------------------

cohort_summary
group_summary
table1
regression_table
vif_results
hosmer_lemeshow
