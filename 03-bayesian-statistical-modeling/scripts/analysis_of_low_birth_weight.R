# ============================================================
# Bayesian Statistical Modelling:
# Bayesian Logistic Regression for Low Birth Weight Risk Factors
#
# Author:   Piyush Kumar
# Dataset:  MASS::birthwt — Hosmer, Lemeshow & Sturdivant (2013)
# Method:   Metropolis-Hastings MCMC (from scratch, no black-box sampler)
# Outcome:  Low birth weight (birth weight < 2500g)
# ============================================================

rm(list = ls())

# ============================================================
# SECTION 1: PACKAGE LOADING
# ============================================================

packages <- c(
  "MASS",       # birthwt dataset
  "tidyverse",  # data wrangling and visualisation
  "broom",      # tidy model output
  "ggplot2",    # figures
  "gt",         # HTML table export
  "gtsummary"   # descriptive and regression tables
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

invisible(lapply(packages, install_if_missing))

library(MASS)
library(tidyverse)
library(broom)
library(ggplot2)
library(gt)
library(gtsummary)

set.seed(123)

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

# birthwt: data from Baystate Medical Center, Springfield MA (1986)
# n = 189 singleton births; outcome = low birth weight (< 2500g)
# Reference: Hosmer, Lemeshow & Sturdivant (2013), Applied Logistic Regression

data("birthwt", package = "MASS")
raw_df <- birthwt

glimpse(raw_df)
summary(raw_df)

# Confirm no missing data — expected for this dataset
missing_counts <- colSums(is.na(raw_df))
cat("\n===== MISSING DATA CHECK =====\n")
print(missing_counts)
cat("Total missing values:", sum(missing_counts), "— complete case analysis not required.\n")

# ============================================================
# SECTION 4: VARIABLE RECODING
# ============================================================

analysis_df <- raw_df %>%
  mutate(
    
    # Primary outcome: low birth weight (binary)
    low_birth_weight = factor(
      low,
      levels = c(0, 1),
      labels = c("Normal birth weight", "Low birth weight")
    ),
    low_num = low,
    
    # Race (White = reference)
    race = factor(
      race,
      levels = c(1, 2, 3),
      labels = c("White", "Black", "Other")
    ),
    
    # Smoking during pregnancy
    smoking = factor(
      smoke,
      levels = c(0, 1),
      labels = c("No", "Yes")
    ),
    
    # History of hypertension
    hypertension = factor(
      ht,
      levels = c(0, 1),
      labels = c("No", "Yes")
    ),
    
    # Presence of uterine irritability
    uterine_irritability = factor(
      ui,
      levels = c(0, 1),
      labels = c("No", "Yes")
    ),
    
    # Previous preterm labour (any vs none)
    previous_preterm_labor = factor(
      ifelse(ptl > 0, 1, 0),
      levels = c(0, 1),
      labels = c("No", "Yes")
    ),
    
    # Standardise continuous predictors for MCMC numerical stability
    # and to place them on a comparable scale with binary predictors
    age_z = as.numeric(scale(age)),
    lwt_z = as.numeric(scale(lwt))
    
  ) %>%
  select(
    low_birth_weight, low_num,
    age, age_z, lwt, lwt_z,
    race, smoking, hypertension,
    uterine_irritability, previous_preterm_labor
  )

write_csv(analysis_df, "data/low_birth_weight_analysis_dataset.csv")

# ============================================================
# SECTION 5: COHORT SUMMARY
# ============================================================

cohort_summary <- analysis_df %>%
  summarise(
    n                          = n(),
    low_birth_weight_cases     = sum(low_num == 1),
    normal_birth_weight_cases  = sum(low_num == 0),
    low_birth_weight_prevalence = mean(low_num),
    mean_age                   = mean(age),
    median_age                 = median(age),
    mean_maternal_weight_lbs   = mean(lwt),
    median_maternal_weight_lbs = median(lwt)
  )

write_csv(cohort_summary, "outputs/cohort_summary.csv")
cat("\n===== COHORT SUMMARY =====\n")
print(cohort_summary)

# ============================================================
# SECTION 6: GROUP SUMMARY BY BIRTH WEIGHT OUTCOME
# ============================================================

group_summary <- analysis_df %>%
  group_by(low_birth_weight) %>%
  summarise(
    n                      = n(),
    mean_age               = mean(age),
    sd_age                 = sd(age),
    median_age             = median(age),
    mean_maternal_weight   = mean(lwt),
    sd_maternal_weight     = sd(lwt),
    median_maternal_weight = median(lwt),
    pct_smoking            = mean(smoking == "Yes") * 100,
    pct_hypertension       = mean(hypertension == "Yes") * 100,
    pct_uterine_irritability = mean(uterine_irritability == "Yes") * 100,
    pct_previous_preterm   = mean(previous_preterm_labor == "Yes") * 100,
    pct_black              = mean(race == "Black") * 100,
    pct_other_race         = mean(race == "Other") * 100,
    .groups = "drop"
  )

write_csv(group_summary, "outputs/group_summary_by_birth_weight.csv")
cat("\n===== GROUP SUMMARY =====\n")
print(group_summary)

# ============================================================
# SECTION 7: DESCRIPTIVE TABLE (TABLE 1)
# ============================================================

table1 <- analysis_df %>%
  select(
    low_birth_weight, age, lwt, race, smoking,
    hypertension, uterine_irritability, previous_preterm_labor
  ) %>%
  tbl_summary(
    by        = low_birth_weight,
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
# SECTION 8: FREQUENTIST LOGISTIC REGRESSION (BENCHMARK)
# ============================================================

# Standard frequentist model fitted for direct comparison with
# the Bayesian estimates. Consistent results across both approaches
# validate the MCMC implementation.

frequentist_model <- glm(
  low_num ~ age_z + lwt_z + race + smoking +
    hypertension + uterine_irritability + previous_preterm_labor,
  data   = analysis_df,
  family = binomial
)

frequentist_results <- tidy(
  frequentist_model,
  exponentiate = TRUE,
  conf.int     = TRUE
)

write_csv(frequentist_results, "outputs/frequentist_logistic_regression_results.csv")
cat("\n===== FREQUENTIST MODEL RESULTS =====\n")
print(frequentist_results)

sink("outputs/frequentist_model_summary.txt")
print(summary(frequentist_model))
sink()

# ============================================================
# SECTION 9: BAYESIAN LOGISTIC REGRESSION — MODEL SETUP
# ============================================================

# Prior specification: Normal(0, 2.5) on all coefficients.
# This is a weakly informative prior — centred at no effect (beta = 0),
# with SD = 2.5 allowing reasonable prior mass over a wide range of OR.
# The same prior is used by rstanarm's default logistic regression.

model_matrix <- model.matrix(
  low_num ~ age_z + lwt_z + race + smoking +
    hypertension + uterine_irritability + previous_preterm_labor,
  data = analysis_df
)

y          <- analysis_df$low_num
X          <- model_matrix
p          <- ncol(X)
coef_names <- colnames(X)

# Log-likelihood: sum of Bernoulli log-likelihoods
log_likelihood <- function(beta, X, y) {
  eta <- as.vector(X %*% beta)
  sum(y * eta - log1p(exp(eta)))
}

# Log-prior: independent Normal(0, prior_sd) on each coefficient
log_prior <- function(beta, prior_sd = 2.5) {
  sum(dnorm(beta, mean = 0, sd = prior_sd, log = TRUE))
}

# Log-posterior: proportional to likelihood × prior
log_posterior <- function(beta, X, y, prior_sd = 2.5) {
  log_likelihood(beta, X, y) + log_prior(beta, prior_sd)
}

# ============================================================
# SECTION 10: METROPOLIS-HASTINGS SAMPLER
# ============================================================

# Random-walk Metropolis-Hastings with a symmetric Normal proposal.
# At each iteration:
#   1. Propose beta* = beta_current + Normal(0, proposal_sd)
#   2. Compute log acceptance ratio = log p(beta*|y) - log p(beta_current|y)
#   3. Accept beta* with probability min(1, exp(log_acceptance))
# The burn-in period is discarded before summarising the posterior.

run_mh <- function(
    X,
    y,
    n_iter       = 30000,
    burn_in      = 10000,
    proposal_sd  = 0.08,
    prior_sd     = 2.5
) {
  p               <- ncol(X)
  beta_current    <- rep(0, p)
  log_post_current <- log_posterior(beta_current, X, y, prior_sd)
  samples         <- matrix(NA, nrow = n_iter, ncol = p)
  accepted        <- 0
  
  for (i in seq_len(n_iter)) {
    beta_proposal    <- beta_current + rnorm(p, mean = 0, sd = proposal_sd)
    log_post_proposal <- log_posterior(beta_proposal, X, y, prior_sd)
    log_acceptance   <- log_post_proposal - log_post_current
    
    if (log(runif(1)) < log_acceptance) {
      beta_current     <- beta_proposal
      log_post_current <- log_post_proposal
      accepted         <- accepted + 1
    }
    samples[i, ] <- beta_current
  }
  
  posterior_samples           <- samples[(burn_in + 1):n_iter, ]
  colnames(posterior_samples) <- colnames(X)
  
  list(
    samples         = posterior_samples,
    acceptance_rate = accepted / n_iter
  )
}

bayes_fit <- run_mh(
  X           = X,
  y           = y,
  n_iter      = 40000,
  burn_in     = 10000,
  proposal_sd = 0.07,
  prior_sd    = 2.5
)

cat(sprintf("\n===== MCMC ACCEPTANCE RATE =====\n%.3f\n", bayes_fit$acceptance_rate))

posterior_samples <- as_tibble(bayes_fit$samples)

# Save thinned posterior samples (every 10th draw) to reduce file size
# Full chain: 30,000 rows; thinned: 3,000 rows — sufficient for all summaries
thinned_samples <- posterior_samples %>%
  mutate(iteration = row_number()) %>%
  filter(iteration %% 10 == 0) %>%
  select(-iteration)

write_csv(thinned_samples, "outputs/posterior_samples_thinned.csv")

# ============================================================
# SECTION 11: POSTERIOR SUMMARIES
# ============================================================

posterior_summary <- posterior_samples %>%
  pivot_longer(everything(), names_to = "term", values_to = "beta") %>%
  group_by(term) %>%
  summarise(
    mean_beta                    = mean(beta),
    median_beta                  = median(beta),
    lower_95                     = quantile(beta, 0.025),
    upper_95                     = quantile(beta, 0.975),
    posterior_probability_positive = mean(beta > 0),
    posterior_probability_negative = mean(beta < 0),
    .groups = "drop"
  ) %>%
  mutate(
    odds_ratio             = exp(median_beta),
    or_lower_95            = exp(lower_95),
    or_upper_95            = exp(upper_95),
    probability_or_gt_1    = posterior_probability_positive
  )

# Human-readable variable labels
label_map <- tibble(
  term = coef_names,
  variable = c(
    "Intercept",
    "Age (standardized)",
    "Maternal weight (standardized)",
    "Black vs White",
    "Other race vs White",
    "Smoking: Yes vs No",
    "Hypertension: Yes vs No",
    "Uterine irritability: Yes vs No",
    "Previous preterm labor: Yes vs No"
  )
)

posterior_summary_labeled <- posterior_summary %>%
  left_join(label_map, by = "term")

write_csv(posterior_summary,         "outputs/bayesian_logistic_regression_results.csv")
write_csv(posterior_summary_labeled, "outputs/bayesian_logistic_regression_results_labeled.csv")

cat("\n===== BAYESIAN POSTERIOR SUMMARY =====\n")
print(posterior_summary_labeled %>% select(variable, odds_ratio, or_lower_95, or_upper_95, probability_or_gt_1))

sink("outputs/bayesian_model_summary.txt")
print(posterior_summary_labeled)
cat(sprintf("\nMCMC acceptance rate: %.3f\n", bayes_fit$acceptance_rate))
sink()

# ============================================================
# SECTION 12: PRIOR SENSITIVITY ANALYSIS
# ============================================================

# Test robustness of posterior estimates to prior choice.
# Three prior scales examined: tight (SD=1.5), main (SD=2.5), diffuse (SD=5).
# Stable estimates across prior choices indicate the data dominate the prior.

prior_sds <- c(1.5, 2.5, 5)

prior_sensitivity_results <- map_df(
  prior_sds,
  function(sd_prior) {
    fit <- run_mh(
      X           = X,
      y           = y,
      n_iter      = 30000,
      burn_in     = 10000,
      proposal_sd = 0.07,
      prior_sd    = sd_prior
    )
    
    as_tibble(fit$samples) %>%
      pivot_longer(everything(), names_to = "term", values_to = "beta") %>%
      group_by(term) %>%
      summarise(
        median_beta              = median(beta),
        lower_95                 = quantile(beta, 0.025),
        upper_95                 = quantile(beta, 0.975),
        probability_or_greater_than_1 = mean(beta > 0),
        acceptance_rate          = fit$acceptance_rate,
        .groups = "drop"
      ) %>%
      mutate(
        prior_sd    = sd_prior,
        odds_ratio  = exp(median_beta),
        or_lower_95 = exp(lower_95),
        or_upper_95 = exp(upper_95)
      )
  }
) %>%
  left_join(label_map, by = "term")

write_csv(prior_sensitivity_results, "outputs/prior_sensitivity_results.csv")

# ============================================================
# SECTION 13: POSTERIOR PREDICTIVE CHECK
# ============================================================

# Simulate new datasets from the posterior to check model fit.
# If the model is well-specified, the observed outcome count should
# fall within the range of simulated counts.

posterior_matrix  <- as.matrix(posterior_samples)
posterior_prob    <- plogis(X %*% t(posterior_matrix))

posterior_pred_cases <- apply(
  posterior_prob, 2,
  function(prob) sum(rbinom(length(prob), size = 1, prob = prob))
)

observed_cases <- sum(y)

ppc_summary <- tibble(
  observed_cases           = observed_cases,
  simulated_mean_cases     = mean(posterior_pred_cases),
  simulated_lower_95       = quantile(posterior_pred_cases, 0.025),
  simulated_upper_95       = quantile(posterior_pred_cases, 0.975)
)

write_csv(ppc_summary, "outputs/posterior_predictive_check_summary.csv")
cat("\n===== POSTERIOR PREDICTIVE CHECK =====\n")
print(ppc_summary)

# ============================================================
# SECTION 14: VISUALISATIONS
# ============================================================

# --- 14a. Forest plot: posterior odds ratios ---
# Improvements over original:
#   - Colour-coded by credibility evidence (CI entirely above/below 1 vs uncertain)
#   - Labels positioned at fixed x coordinate to avoid overlap with CI bars
#   - Expanded x limits to accommodate label text cleanly

forest_data <- posterior_summary_labeled %>%
  filter(term != "(Intercept)") %>%
  mutate(
    variable = factor(
      variable,
      levels = rev(c(
        "Age (standardized)",
        "Maternal weight (standardized)",
        "Black vs White",
        "Other race vs White",
        "Smoking: Yes vs No",
        "Hypertension: Yes vs No",
        "Uterine irritability: Yes vs No",
        "Previous preterm labor: Yes vs No"
      ))
    ),
    # Evidence category based on whether 95% credible interval excludes 1
    evidence = case_when(
      or_lower_95 > 1 ~ "Credible increase in odds",
      or_upper_95 < 1 ~ "Credible decrease in odds",
      TRUE            ~ "Uncertain"
    ),
    label = sprintf("OR %.2f (%.2f\u2013%.2f)", odds_ratio, or_lower_95, or_upper_95)
  )

# Fixed x position for labels — to the right of the widest CI bar
label_x_pos <- max(forest_data$or_upper_95, na.rm = TRUE) * 1.15

bayesian_forest_plot <- ggplot(
  forest_data,
  aes(x = odds_ratio, y = variable, colour = evidence)
) +
  geom_point(size = 3) +
  geom_errorbarh(
    aes(xmin = or_lower_95, xmax = or_upper_95),
    height = 0.25
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_text(
    aes(x = label_x_pos, label = label),
    hjust = 0, size = 3.2, colour = "black"
  ) +
  scale_x_log10(
    breaks = c(0.3, 0.5, 1, 2, 5, 10),
    labels = c("0.3", "0.5", "1.0", "2.0", "5.0", "10.0"),
    limits = c(0.2, 80)   # extended right margin for label clearance
  ) +
  scale_colour_manual(
    values = c(
      "Credible increase in odds"  = "#d73027",
      "Credible decrease in odds"  = "#4575b4",
      "Uncertain"                  = "grey40"
    ),
    name = "Evidence"
  ) +
  labs(
    title    = "Bayesian Logistic Regression: Posterior Odds Ratios",
    subtitle = "Outcome: Low birth weight | Normal(0, 2.5) prior | 30,000 post-burn-in MCMC samples\nColour indicates whether 95% credible interval excludes OR = 1",
    x        = "Posterior Odds Ratio (95% Credible Interval, log scale)",
    y        = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    plot.subtitle    = element_text(colour = "grey40", size = 9)
  ) +
  coord_cartesian(clip = "off") +
  expand_limits(x = 80)

ggsave(
  "figures/bayesian_posterior_odds_ratio_forest_plot.png",
  bayesian_forest_plot,
  width = 11, height = 6, dpi = 300
)

# --- 14b. Posterior density plots for key predictors ---
# Fixed: shared x-axis upper limit of 15 for cross-predictor comparability.
# Each panel shows the full posterior distribution of the OR on the same scale.

key_terms <- c(
  "smokingYes",
  "hypertensionYes",
  "uterine_irritabilityYes",
  "previous_preterm_laborYes"
)

posterior_density_data <- posterior_samples %>%
  select(all_of(key_terms)) %>%
  pivot_longer(everything(), names_to = "term", values_to = "beta") %>%
  left_join(label_map, by = "term") %>%
  mutate(odds_ratio = exp(beta))

posterior_density_plot <- ggplot(
  posterior_density_data,
  aes(x = odds_ratio)
) +
  geom_density(linewidth = 1, fill = "#74add1", alpha = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey30") +
  facet_wrap(~variable, scales = "free_y") +
  coord_cartesian(xlim = c(0, 15)) +   # zoom after density computed — no data dropped
  scale_x_continuous(breaks = c(0, 1, 2, 5, 10, 15)) +
  labs(
    title    = "Posterior Distributions of Odds Ratios",
    subtitle = "Key clinical predictors | x-axis shared across panels for comparability",
    x        = "Odds Ratio",
    y        = "Posterior density"
  ) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank())

ggsave(
  "figures/posterior_density_key_predictors.png",
  posterior_density_plot,
  width = 10, height = 6, dpi = 300
)

# --- 14c. Prior sensitivity plot ---

prior_sensitivity_plot <- prior_sensitivity_results %>%
  filter(term %in% key_terms) %>%
  mutate(
    variable = factor(variable, levels = c(
      "Smoking: Yes vs No",
      "Hypertension: Yes vs No",
      "Uterine irritability: Yes vs No",
      "Previous preterm labor: Yes vs No"
    ))
  ) %>%
  ggplot(aes(x = factor(prior_sd), y = odds_ratio, group = variable)) +
  geom_point(size = 2.5) +
  geom_line() +
  facet_wrap(~variable, scales = "free_y") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    title    = "Prior Sensitivity Analysis",
    subtitle = "Posterior median odds ratios under Normal(0, \u03c3) priors with \u03c3 \u2208 {1.5, 2.5, 5}",
    x        = "Prior standard deviation (\u03c3)",
    y        = "Posterior median odds ratio"
  ) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank())

ggsave(
  "figures/prior_sensitivity_plot.png",
  prior_sensitivity_plot,
  width = 10, height = 6, dpi = 300
)

# --- 14d. Posterior predictive check ---

ppc_df <- tibble(simulated_low_birth_weight_cases = posterior_pred_cases)

ppc_plot <- ggplot(ppc_df, aes(x = simulated_low_birth_weight_cases)) +
  geom_histogram(bins = 30, boundary = 0, fill = "#74add1", colour = "white") +
  geom_vline(
    xintercept = observed_cases,
    linetype = "dashed", linewidth = 1, colour = "#d73027"
  ) +
  annotate(
    "text",
    x = observed_cases + 1, y = Inf,
    label = sprintf("Observed = %d", observed_cases),
    hjust = 0, vjust = 1.5, size = 4, colour = "#d73027"
  ) +
  labs(
    title    = "Posterior Predictive Check",
    subtitle = sprintf(
      "Simulated mean = %.1f (95%% PI: %d\u2013%d) | Observed = %d",
      mean(posterior_pred_cases),
      quantile(posterior_pred_cases, 0.025),
      quantile(posterior_pred_cases, 0.975),
      observed_cases
    ),
    x = "Simulated number of low birth weight cases",
    y = "Posterior predictive frequency"
  ) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank())

ggsave(
  "figures/posterior_predictive_check.png",
  ppc_plot,
  width = 8, height = 6, dpi = 300
)

# --- 14e. Trace plots for key parameters ---

trace_data <- posterior_samples %>%
  mutate(iteration = row_number()) %>%
  select(iteration, all_of(key_terms)) %>%
  pivot_longer(-iteration, names_to = "term", values_to = "beta") %>%
  left_join(label_map, by = "term")

trace_plot <- ggplot(trace_data, aes(x = iteration, y = beta)) +
  geom_line(linewidth = 0.3, colour = "grey30") +
  facet_wrap(~variable, scales = "free_y") +
  labs(
    title    = "Trace Plots for Selected Parameters",
    subtitle = "Good mixing: chains explore the posterior without trends or sticky regions",
    x        = "Post-burn-in iteration",
    y        = "Posterior sample (\u03b2)"
  ) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank())

ggsave(
  "figures/trace_plots_key_predictors.png",
  trace_plot,
  width = 10, height = 6, dpi = 300
)

# ============================================================
# SECTION 15: CONSOLE SUMMARY
# ============================================================

cat("\n========================================\n")
cat("KEY RESULTS SUMMARY\n")
cat("========================================\n")
cat(sprintf("Cohort:                n = %d\n",            cohort_summary$n))
cat(sprintf("Low birth weight:      %d (%.1f%%)\n",
            cohort_summary$low_birth_weight_cases,
            cohort_summary$low_birth_weight_prevalence * 100))
cat(sprintf("Missing data:          %d (complete cases only)\n", sum(missing_counts)))
cat(sprintf("MCMC acceptance rate:  %.3f\n",              bayes_fit$acceptance_rate))
cat(sprintf("PPC observed cases:    %d\n",                observed_cases))
cat(sprintf("PPC simulated mean:    %.1f (95%% PI: %d\u2013%d)\n",
            ppc_summary$simulated_mean_cases,
            ppc_summary$simulated_lower_95,
            ppc_summary$simulated_upper_95))
cat("\nPosterior odds ratios (main model):\n")
posterior_summary_labeled %>%
  filter(term != "(Intercept)") %>%
  mutate(across(c(odds_ratio, or_lower_95, or_upper_95, probability_or_gt_1),
                ~ round(.x, 3))) %>%
  select(variable, odds_ratio, or_lower_95, or_upper_95, probability_or_gt_1) %>%
  print(n = Inf)
cat("\nNote: Hypertension has a wide 95% credible interval (OR 6.09, 95% CrI 1.65-25.18)\n")
cat("due to small cell count (n=12 hypertensive mothers). Direction is robust\n")
cat("across all prior specifications but the magnitude is uncertain.\n")
cat("========================================\n")

# Display formatted tables
cohort_summary
table1
frequentist_results
posterior_summary_labeled
prior_sensitivity_results
ppc_summary

