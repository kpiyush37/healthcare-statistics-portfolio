# ============================================================
# Survival Analysis:
# Predictors of Overall Survival in Advanced Lung Cancer
# Veterans' Administration Lung Cancer Trial
#
# Author:   Piyush Kumar
# Dataset:  veteran (survival R package)
#           Kalbfleisch & Prentice (1980), randomised clinical trial
# Design:   Survival analysis — Kaplan-Meier estimation,
#           log-rank testing, Cox proportional hazards regression
# ============================================================

rm(list = ls())

# ============================================================
# SECTION 1: PACKAGE LOADING
# ============================================================


packages <- c(
  "survival",
  "tidyverse",
  "gtsummary",
  "broom",
  "gt",
  "cowplot"    # stacking KM curve + risk table panels
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

invisible(lapply(packages, install_if_missing))

library(survival)
#library(survminer)
library(tidyverse)
library(gtsummary)
library(broom)
library(gt)
library(cowplot)
# ============================================================
# SECTION 2: PROJECT STRUCTURE
# ============================================================

dir.create("data",         showWarnings = FALSE)
dir.create("outputs",      showWarnings = FALSE)
dir.create("figures",      showWarnings = FALSE)
dir.create("presentation", showWarnings = FALSE)
dir.create("scripts",      showWarnings = FALSE)

# ============================================================
# SECTION 3: DATA LOADING AND INSPECTION
# ============================================================

# The veteran dataset is from a randomised two-treatment trial
# for lung cancer conducted by the Veterans Administration.
# Reference: Kalbfleisch & Prentice (1980).
# Available in the survival R package.

survival_df <- survival::veteran

glimpse(survival_df)
summary(survival_df)
colSums(is.na(survival_df))

# ============================================================
# SECTION 4: VARIABLE RECODING
# ============================================================

analysis_df <- survival_df %>%
  mutate(
    
    # Survival outcome
    event = status,
    
    # Treatment arm
    treatment = factor(
      trt,
      levels = c(1, 2),
      labels = c("Standard", "Test")
    ),
    
    # Histological cell type (Squamous = reference)
    cell_type = factor(
      celltype,
      levels = c("squamous", "smallcell", "adeno", "large"),
      labels = c("Squamous", "Small Cell", "Adenocarcinoma", "Large Cell")
    ),
    
    # Prior therapy
    prior_therapy = factor(
      prior,
      levels = c(0, 10),
      labels = c("No", "Yes")
    ),
    
    # Karnofsky performance status grouped at conventional threshold of 70
    # >= 70: able to care for self; < 70: requires assistance
    performance_group = factor(
      case_when(
        karno >= 70 ~ "Higher performance",
        karno <  70 ~ "Lower performance"
      ),
      levels = c("Higher performance", "Lower performance")
    )
    
  ) %>%
  select(
    time, status, event, treatment, age,
    cell_type, karno, performance_group, diagtime, prior_therapy
  )

write_csv(analysis_df, "data/lung_cancer_survival_dataset.csv")

# ============================================================
# SECTION 5: COHORT SUMMARY
# ============================================================

cohort_summary <- analysis_df %>%
  summarise(
    n                    = n(),
    deaths               = sum(event == 1),
    censored             = sum(event == 0),
    event_rate           = mean(event),
    median_survival_days = median(time),
    mean_survival_days   = mean(time),
    mean_age             = mean(age),
    median_age           = median(age),
    mean_karno           = mean(karno),
    median_karno         = median(karno)
  )

write_csv(cohort_summary, "outputs/cohort_summary.csv")
cohort_summary

# ============================================================
# SECTION 6: BASELINE DESCRIPTIVE TABLE
# ============================================================

table1 <- analysis_df %>%
  select(time, event, treatment, age, cell_type,
         karno, performance_group, diagtime, prior_therapy) %>%
  tbl_summary(
    by        = treatment,
    statistic = list(
      all_continuous()  ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "no"
  ) %>%
  add_overall() %>%
  add_p()

table1
gtsummary::as_gt(table1) %>% gt::gtsave("outputs/table1.html")

# ============================================================
# SECTION 7: KAPLAN-MEIER SURVIVAL ESTIMATION
# ============================================================

# Helper: converts survfit object to tidy data frame for ggplot2
km_to_df <- function(fit, group_var = NULL) {
  sf <- summary(fit, times = NULL, extend = TRUE)
  df <- tibble(
    time     = sf$time,
    survival = sf$surv,
    lower    = sf$lower,
    upper    = sf$upper,
    n.risk   = sf$n.risk,
    strata   = if (!is.null(sf$strata)) as.character(sf$strata) else "Overall"
  )
  if (!is.null(group_var)) {
    df <- df %>%
      mutate(strata = str_remove(strata, paste0(group_var, "=")))
  }
  df
}

# Helper: adds starting row at time = 0 for each stratum
add_origin <- function(df) {
  origins <- df %>%
    group_by(strata) %>%
    slice(1) %>%
    mutate(time = 0, survival = 1, lower = 1, upper = 1)
  bind_rows(origins, df) %>% arrange(strata, time)
}

# Shared theme
km_theme <- theme_bw(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

# Risk table helper: counts at-risk at evenly spaced time breaks
make_risk_table <- function(fit, group_var = NULL, breaks = NULL) {
  # Clamp breaks to within the observed time range to avoid out-of-range errors
  max_time <- max(fit$time)
  breaks   <- breaks[breaks <= max_time]
  
  sf  <- summary(fit, times = breaks, extend = FALSE)
  df  <- tibble(
    time   = sf$time,
    n.risk = sf$n.risk,
    strata = if (!is.null(sf$strata)) as.character(sf$strata) else "Overall"
  )
  if (!is.null(group_var))
    df <- df %>% mutate(strata = str_remove(strata, paste0(group_var, "=")))
  df
}

# ---- generic KM plot + risk table function ----
plot_km <- function(fit,
                    group_var    = NULL,
                    title        = "",
                    subtitle     = "",
                    legend_title = "Group",
                    palette      = NULL) {
  
  # --- tidy KM data ---
  sf <- summary(fit, extend = TRUE)
  km_df <- tibble(
    time     = sf$time,
    survival = sf$surv,
    lower    = sf$lower,
    upper    = sf$upper,
    strata   = if (!is.null(sf$strata)) as.character(sf$strata) else "Overall"
  )
  if (!is.null(group_var))
    km_df <- km_df %>%
    mutate(strata = str_remove(strata, paste0(group_var, "=")))
  
  # add time-zero origin per stratum
  km_df <- bind_rows(
    km_df %>% group_by(strata) %>% slice(1) %>%
      mutate(time = 0, survival = 1, lower = 1, upper = 1),
    km_df
  ) %>% arrange(strata, time)
  
  n_strata <- length(unique(km_df$strata))
  if (is.null(palette))
    palette <- c("#2c7bb6", "#d7191c", "#f4a11d", "#1a9641")[seq_len(n_strata)]
  
  max_time   <- max(fit$time)
  risk_times <- seq(0, floor(max_time / 200) * 200, by = 200)
  risk_times <- risk_times[risk_times < max_time]
  x_breaks   <- c(risk_times, max_time)
  
  # --- risk table: built inline exactly as tested above ---
  sf_risk <- summary(fit, times = risk_times, extend = TRUE)
  risk_df <- tibble(
    time   = sf_risk$time,
    n.risk = sf_risk$n.risk,
    strata = if (!is.null(sf_risk$strata)) as.character(sf_risk$strata)
    else "Overall"
  )
  if (!is.null(group_var))
    risk_df <- risk_df %>%
    mutate(strata = str_remove(strata, paste0(group_var, "=")))
  
  # --- log-rank p-value ---
  pval_label <- NULL
  if (n_strata > 1 && !is.null(group_var)) {
    lr <- survdiff(as.formula(paste("Surv(time, event) ~", group_var)),
                   data = analysis_df)
    pv <- 1 - pchisq(lr$chisq, length(lr$n) - 1)
    pval_label <- sprintf("Log-rank p %s",
                          ifelse(pv < 0.001, "< 0.001", sprintf("= %.3f", pv)))
  }
  
  # --- KM curve panel ---
  p_curve <- ggplot(km_df,
                    aes(x = time, y = survival, colour = strata, group = strata)) +
    geom_step(linewidth = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = strata),
                alpha = 0.12, colour = NA) +
    scale_colour_manual(values = palette, name = legend_title) +
    scale_fill_manual(values   = palette, name = legend_title) +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(limits = c(0, 1),
                       labels = scales::percent_format(accuracy = 1)) +
    labs(title = title, subtitle = subtitle,
         x = NULL, y = "Survival probability") +
    theme_bw(base_size = 13) +
    theme(legend.position  = "bottom",
          panel.grid.minor = element_blank())
  
  if (!is.null(pval_label))
    p_curve <- p_curve +
    annotate("text", x = max_time * 0.55, y = 0.92,
             label = pval_label, size = 4, hjust = 0)
  
  # --- at-risk table panel ---
  p_risk <- ggplot(risk_df,
                   aes(x = time, y = strata, label = n.risk, colour = strata)) +
    geom_text(size = 3.5, fontface = "bold") +
    scale_colour_manual(values = palette, guide = "none") +
    scale_x_continuous(breaks = x_breaks) +
    labs(x = "Time (days)", y = NULL, title = "Number at risk") +
    theme_bw(base_size = 11) +
    theme(axis.text.y  = element_text(hjust = 1),
          panel.grid   = element_blank(),
          plot.title   = element_text(size = 10, face = "bold"),
          axis.ticks.y = element_blank())
  
  # --- stack panels ---
  cowplot::plot_grid(p_curve, p_risk,
                     ncol = 1, rel_heights = c(3, 1),
                     align = "v", axis = "lr")
}

# --- 7a. Overall ---
km_overall <- survfit(Surv(time, event) ~ 1, data = analysis_df)

p_overall <- plot_km(km_overall,
                     title    = "Overall Survival in Advanced Lung Cancer",
                     subtitle = "Veterans' Administration Lung Cancer Trial (n = 137)")

ggsave("figures/overall_kaplan_meier_curve.png",
       p_overall, width = 9, height = 7, dpi = 300)

# --- 7b. By treatment ---
km_treatment <- survfit(Surv(time, event) ~ treatment, data = analysis_df)

p_treatment <- plot_km(km_treatment,
                       group_var    = "treatment",
                       title        = "Overall Survival by Treatment Group",
                       subtitle     = "Kaplan-Meier curves with log-rank p-value",
                       legend_title = "Treatment",
                       palette      = c("#2c7bb6", "#d7191c"))

ggsave("figures/km_curve_by_treatment.png",
       p_treatment, width = 9, height = 7, dpi = 300)

# --- 7c. By performance status ---
km_performance <- survfit(
  Surv(time, event) ~ performance_group, data = analysis_df)

p_performance <- plot_km(km_performance,
                         group_var    = "performance_group",
                         title        = "Overall Survival by Performance Status",
                         subtitle     = "Karnofsky score grouped at threshold of 70",
                         legend_title = "Performance Status",
                         palette      = c("#1a9641", "#d7191c"))

ggsave("figures/km_curve_by_performance_status.png",
       p_performance, width = 9, height = 7, dpi = 300)

# --- 7d. By cell type ---
km_cell <- survfit(Surv(time, event) ~ cell_type, data = analysis_df)

p_cell <- plot_km(km_cell,
                  group_var    = "cell_type",
                  title        = "Overall Survival by Histological Cell Type",
                  subtitle     = "Reference: Squamous cell carcinoma",
                  legend_title = "Cell Type",
                  palette      = c("#2c7bb6","#d7191c","#f4a11d","#1a9641"))

ggsave("figures/km_curve_by_cell_type.png",
       p_cell, width = 10, height = 8, dpi = 300)

# ============================================================
# SECTION 8: MEDIAN SURVIVAL ESTIMATES
# ============================================================

median_survival <- tibble(
  group               = "Overall",
  median_survival_days = summary(km_overall)$table["median"],
  lower_95_ci         = summary(km_overall)$table["0.95LCL"],
  upper_95_ci         = summary(km_overall)$table["0.95UCL"]
)

median_treatment <- as.data.frame(summary(km_treatment)$table) %>%
  rownames_to_column("group") %>%
  transmute(
    group                = str_replace(group, "treatment=", ""),
    median_survival_days = median,
    lower_95_ci          = `0.95LCL`,
    upper_95_ci          = `0.95UCL`
  )

median_performance <- as.data.frame(summary(km_performance)$table) %>%
  rownames_to_column("group") %>%
  transmute(
    group                = str_replace(group, "performance_group=", ""),
    median_survival_days = median,
    lower_95_ci          = `0.95LCL`,
    upper_95_ci          = `0.95UCL`
  )

median_cell <- as.data.frame(summary(km_cell)$table) %>%
  rownames_to_column("group") %>%
  transmute(
    group                = str_replace(group, "cell_type=", ""),
    median_survival_days = median,
    lower_95_ci          = `0.95LCL`,
    upper_95_ci          = `0.95UCL`
  )

write_csv(median_survival,    "outputs/median_survival_overall.csv")
write_csv(median_treatment,   "outputs/median_survival_by_treatment.csv")
write_csv(median_performance, "outputs/median_survival_by_performance.csv")
write_csv(median_cell,        "outputs/median_survival_by_cell_type.csv")

median_survival
median_treatment
median_performance
median_cell

# ============================================================
# SECTION 9: LOG-RANK TESTS
# ============================================================

# Log-rank test compares survival distributions between groups
# under the null hypothesis of no difference in survival.

logrank_treatment   <- survdiff(Surv(time, event) ~ treatment,      data = analysis_df)
logrank_performance <- survdiff(Surv(time, event) ~ performance_group, data = analysis_df)
logrank_cell        <- survdiff(Surv(time, event) ~ cell_type,      data = analysis_df)
logrank_prior       <- survdiff(Surv(time, event) ~ prior_therapy,  data = analysis_df)

logrank_results <- tibble(
  comparison = c(
    "Treatment group",
    "Performance status (Karnofsky >= 70 vs < 70)",
    "Cell type (4 groups)",
    "Prior therapy"
  ),
  chisq = c(
    logrank_treatment$chisq,
    logrank_performance$chisq,
    logrank_cell$chisq,
    logrank_prior$chisq
  ),
  df = c(
    length(logrank_treatment$n)   - 1,
    length(logrank_performance$n) - 1,
    length(logrank_cell$n)        - 1,
    length(logrank_prior$n)       - 1
  ),
  p_value = c(
    1 - pchisq(logrank_treatment$chisq,   length(logrank_treatment$n)   - 1),
    1 - pchisq(logrank_performance$chisq, length(logrank_performance$n) - 1),
    1 - pchisq(logrank_cell$chisq,        length(logrank_cell$n)        - 1),
    1 - pchisq(logrank_prior$chisq,       length(logrank_prior$n)       - 1)
  )
)

write_csv(logrank_results, "outputs/logrank_test_results.csv")
logrank_results

# ============================================================
# SECTION 10: COX PROPORTIONAL HAZARDS REGRESSION
# ============================================================

# Three nested models to assess estimate stability:
#   Basic    — treatment and age only
#   Clinical — adds Karnofsky score and cell type
#   Full     — adds time since diagnosis and prior therapy

cox_basic <- coxph(
  Surv(time, event) ~ treatment + age,
  data = analysis_df
)

cox_clinical <- coxph(
  Surv(time, event) ~ treatment + age + karno + cell_type,
  data = analysis_df
)

cox_full <- coxph(
  Surv(time, event) ~ treatment + age + karno +
    cell_type + diagtime + prior_therapy,
  data = analysis_df
)

summary(cox_full)

cox_results <- tidy(cox_full, exponentiate = TRUE, conf.int = TRUE)
write_csv(cox_results, "outputs/cox_model_results.csv")
cox_results

cox_table <- tbl_regression(cox_full, exponentiate = TRUE)
cox_table
gtsummary::as_gt(cox_table) %>% gt::gtsave("outputs/cox_regression_table.html")

# ============================================================
# SECTION 11: MODEL CONCORDANCE (C-INDEX)
# ============================================================

# The concordance statistic (C-index) measures model discrimination:
# the probability that a randomly selected participant who died earlier
# had a higher predicted hazard than one who survived longer.
# Interpretation: 0.5 = chance, 0.7+ = acceptable, 0.8+ = good.

concordance_results <- tibble(
  model       = c("Basic", "Clinical", "Full"),
  concordance = c(
    summary(cox_basic)$concordance[1],
    summary(cox_clinical)$concordance[1],
    summary(cox_full)$concordance[1]
  ),
  se = c(
    summary(cox_basic)$concordance[2],
    summary(cox_clinical)$concordance[2],
    summary(cox_full)$concordance[2]
  )
)

write_csv(concordance_results, "outputs/model_concordance.csv")
cat("\n===== MODEL CONCORDANCE (C-INDEX) =====\n")
print(concordance_results)

# ============================================================
# SECTION 12: SENSITIVITY ANALYSIS — MODEL SPECIFICATION
# ============================================================

cox_sensitivity <- list(
  Basic    = cox_basic,
  Clinical = cox_clinical,
  Full     = cox_full
) %>%
  map_df(~ tidy(.x, exponentiate = TRUE, conf.int = TRUE), .id = "model")

write_csv(cox_sensitivity, "outputs/cox_sensitivity_analysis.csv")
cox_sensitivity

# ============================================================
# SECTION 13: PROPORTIONAL HAZARDS ASSUMPTION (cox.zph)
# ============================================================

# The Cox model assumes that the hazard ratio between groups
# is constant over time (proportional hazards). cox.zph() tests
# this by correlating scaled Schoenfeld residuals with time.
# A significant p-value indicates the assumption is violated.

ph_test <- cox.zph(cox_full)
ph_test

ph_results <- as.data.frame(ph_test$table) %>%
  rownames_to_column("variable") %>%
  mutate(
    assumption_met = ifelse(p > 0.05, "Yes", "No — violation detected")
  )

write_csv(ph_results, "outputs/proportional_hazards_test.csv")

cat("\n===== PROPORTIONAL HAZARDS TEST =====\n")
print(ph_results)

# Save Schoenfeld residual diagnostic plots
png(
  filename = "figures/proportional_hazards_diagnostic.png",
  width = 1800, height = 1200, res = 150
)
par(mfrow = c(3, 3))
plot(ph_test)
dev.off()

# ============================================================
# SECTION 14: HANDLING PH VIOLATIONS
# ============================================================

# cox.zph() identifies significant violations for:
#   karno     (p < 0.001) — effect of Karnofsky score varies over time
#   cell_type (p = 0.002) — effect of cell type varies over time
#
# Remedy applied: stratification.
# Stratifying on a violating variable removes it from the
# hazard ratio estimation and instead fits a separate baseline
# hazard for each stratum. This relaxes the PH assumption for
# that variable while preserving the PH structure for others.

cox_stratified <- coxph(
  Surv(time, event) ~ treatment + age + strata(cell_type) +
    strata(performance_group) + diagtime + prior_therapy,
  data = analysis_df
)

summary(cox_stratified)

cox_stratified_results <- tidy(cox_stratified, exponentiate = TRUE, conf.int = TRUE)
write_csv(cox_stratified_results, "outputs/cox_stratified_model_results.csv")

# Verify PH assumption is resolved in the stratified model
ph_test_stratified <- cox.zph(cox_stratified)
cat("\n===== PH TEST — STRATIFIED MODEL =====\n")
print(ph_test_stratified$table)

ph_stratified_results <- as.data.frame(ph_test_stratified$table) %>%
  rownames_to_column("variable")
write_csv(ph_stratified_results, "outputs/ph_test_stratified_model.csv")

# Save stratified model summary
sink("outputs/cox_stratified_model_summary.txt")
print(summary(cox_stratified))
sink()

# ============================================================
# SECTION 15: FOREST PLOT — HAZARD RATIOS
# ============================================================

# Labels and grouping for the full (unstratified) model
label_map <- c(
  "treatmentTest"          = "Treatment: Test (ref: Standard)",
  "age"                    = "Age (per year)",
  "karno"                  = "Karnofsky score (per unit)",
  "cell_typeSmall Cell"    = "Small Cell (ref: Squamous)",
  "cell_typeAdenocarcinoma"= "Adenocarcinoma (ref: Squamous)",
  "cell_typeLarge Cell"    = "Large Cell (ref: Squamous)",
  "diagtime"               = "Time since diagnosis (months)",
  "prior_therapyYes"       = "Prior therapy: Yes (ref: No)"
)

group_map <- c(
  "treatmentTest"           = "Treatment",
  "age"                     = "Demographics",
  "karno"                   = "Clinical",
  "cell_typeSmall Cell"     = "Cell Type\n(ref: Squamous)",
  "cell_typeAdenocarcinoma" = "Cell Type\n(ref: Squamous)",
  "cell_typeLarge Cell"     = "Cell Type\n(ref: Squamous)",
  "diagtime"                = "Disease history",
  "prior_therapyYes"        = "Disease history"
)

forest_data <- cox_results %>%
  mutate(
    label       = label_map[term],
    group       = factor(group_map[term], levels = c(
      "Treatment", "Demographics", "Clinical",
      "Cell Type\n(ref: Squamous)", "Disease history"
    )),
    label       = fct_reorder(label, estimate),
    sig         = p.value < 0.05,
    hr_label    = sprintf("%.2f (%.2f\u2013%.2f)", estimate, conf.low, conf.high)
  )

hr_forest_plot <- ggplot(forest_data, aes(x = estimate, y = label, colour = sig)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.25) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_text(
    aes(x = max(conf.high, na.rm = TRUE) * 1.08, label = hr_label),
    hjust = 0, size = 3.2, colour = "black"
  ) +
  scale_x_log10(
    breaks = c(0.5, 0.75, 1, 1.5, 2, 3, 5),
    labels = c("0.5", "0.75", "1.0", "1.5", "2.0", "3.0", "5.0"),
    limits = c(0.45, 10)   # expanded right margin to prevent label clipping
  ) +
  scale_colour_manual(
    values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"),
    labels = c("TRUE" = "p < 0.05", "FALSE" = "p \u2265 0.05"),
    name   = "Significance"
  ) +
  facet_wrap(~group, scales = "free_y", ncol = 1) +
  labs(
    title    = "Adjusted Hazard Ratios for Overall Survival",
    subtitle = "Full Cox proportional hazards model | Veterans' Administration Lung Cancer Trial\nHR and 95% CI shown; reference categories in parentheses",
    x        = "Hazard Ratio (log scale)",
    y        = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text       = element_text(face = "bold", size = 9),
    legend.position  = "bottom",
    plot.subtitle    = element_text(colour = "grey40", size = 9),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(clip = "off") +
  expand_limits(x = max(forest_data$conf.high, na.rm = TRUE) * 2)

ggsave(
  "figures/hazard_ratio_forest_plot.png",
  hr_forest_plot,
  width = 11, height = 9, dpi = 300
)

# ============================================================
# SECTION 16: ADJUSTED SURVIVAL CURVES BY TREATMENT
# ============================================================

# Predicted survival curves from the Cox model with all
# covariates fixed at median or reference values, isolating
# the treatment effect after adjustment.

newdata_standard <- tibble(
  treatment    = factor("Standard", levels = levels(analysis_df$treatment)),
  age          = median(analysis_df$age),
  karno        = median(analysis_df$karno),
  cell_type    = factor("Squamous",  levels = levels(analysis_df$cell_type)),
  diagtime     = median(analysis_df$diagtime),
  prior_therapy = factor("No",       levels = levels(analysis_df$prior_therapy))
)

newdata_test <- tibble(
  treatment    = factor("Test",    levels = levels(analysis_df$treatment)),
  age          = median(analysis_df$age),
  karno        = median(analysis_df$karno),
  cell_type    = factor("Squamous", levels = levels(analysis_df$cell_type)),
  diagtime     = median(analysis_df$diagtime),
  prior_therapy = factor("No",      levels = levels(analysis_df$prior_therapy))
)

adjusted_standard <- survfit(cox_full, newdata = newdata_standard)
adjusted_test     <- survfit(cox_full, newdata = newdata_test)

standard_summary <- summary(adjusted_standard)
test_summary     <- summary(adjusted_test)

adjusted_survival_df <- bind_rows(
  tibble(
    time      = standard_summary$time,
    survival  = standard_summary$surv,
    lower     = standard_summary$lower,
    upper     = standard_summary$upper,
    treatment = "Standard"
  ),
  tibble(
    time      = test_summary$time,
    survival  = test_summary$surv,
    lower     = test_summary$lower,
    upper     = test_summary$upper,
    treatment = "Test"
  )
)

adjusted_survival_plot <- ggplot(
  adjusted_survival_df,
  aes(x = time, y = survival, group = treatment, linetype = treatment)
) +
  geom_step(linewidth = 1) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper, group = treatment),
    alpha = 0.15
  ) +
  scale_linetype_manual(values = c("Standard" = "solid", "Test" = "dashed")) +
  labs(
    title    = "Adjusted Survival Curves by Treatment",
    subtitle = "Covariates fixed at median/reference values | Full Cox model",
    x        = "Time (days)",
    y        = "Adjusted survival probability",
    linetype = "Treatment"
  ) +
  theme_bw(base_size = 13)

ggsave(
  "figures/adjusted_survival_curves.png",
  adjusted_survival_plot,
  width = 8, height = 6, dpi = 300
)

write_csv(adjusted_survival_df, "outputs/adjusted_survival_curves.csv")

# ============================================================
# SECTION 17: SAVE TEXT SUMMARIES
# ============================================================

sink("outputs/cox_model_summary.txt")
print(summary(cox_full))
sink()

sink("outputs/proportional_hazards_summary.txt")
print(ph_test)
sink()

# ============================================================
# SECTION 18: CONSOLE SUMMARY
# ============================================================

cat("\n========================================\n")
cat("KEY RESULTS SUMMARY\n")
cat("========================================\n")
cat(sprintf("Cohort:             n = %d\n",   cohort_summary$n))
cat(sprintf("Deaths:             %d (%.1f%%)\n",
            cohort_summary$deaths, cohort_summary$event_rate * 100))
cat(sprintf("Median survival:    %d days (95%% CI: see median_survival_overall.csv)\n",
            cohort_summary$median_survival_days))
cat(sprintf("Median age:         %.0f years\n",  cohort_summary$median_age))
cat(sprintf("Median Karnofsky:   %.0f\n",         cohort_summary$median_karno))
cat("\nMedian survival by treatment:\n")
print(median_treatment)
cat("\nMedian survival by performance status:\n")
print(median_performance)
cat("\nLog-rank test results:\n")
print(logrank_results)
cat("\nConcordance (C-index):\n")
print(concordance_results)
cat("\nCox full model — adjusted hazard ratios:\n")
cox_results %>%
  mutate(across(c(estimate, conf.low, conf.high, p.value), ~ round(.x, 3))) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  print(n = Inf)
cat("\nProportional hazards test:\n")
print(ph_results)
cat("========================================\n")

# Display formatted tables
cohort_summary
table1
cox_table




#debug area

table(analysis_df$performance_group)
table(analysis_df$karno >= 70)
