# Real-World Evidence Analysis: Factors Associated with Type 2 Diabetes in NHANES 2017–2018

## Overview

This project investigates demographic, socioeconomic, and clinical factors associated with Type 2 Diabetes using data from the National Health and Nutrition Examination Survey (NHANES) 2017–2018.

The analysis follows a typical Real-World Evidence (RWE) workflow including cohort construction, exploratory data analysis, multivariable modeling, sensitivity analyses, and model diagnostics.

## Research Question

Which demographic, socioeconomic, and clinical factors are associated with prevalent Type 2 Diabetes among U.S. adults?

## Dataset

- Source: NHANES 2017–2018
- Study Design: Cross-sectional observational study
- Population: Adults aged 20 years and older
- Final Analytic Cohort: 4,339 participants
- Diabetes Cases: 701
- Diabetes Prevalence: 16.2%

## Methods

### Exploratory Data Analysis

- Cohort characterization
- Missing data assessment
- BMI distribution by diabetes status
- Age distribution by diabetes status
- Hypertension prevalence by diabetes status

### Statistical Analysis

- Univariable logistic regression
- Multivariable logistic regression
- Sensitivity analyses
- Variance Inflation Factor (VIF)
- Hosmer–Lemeshow goodness-of-fit test

## Key Findings

- Age was positively associated with diabetes.
- Higher BMI was associated with higher odds of diabetes.
- Hypertension was associated with approximately 2.5-fold higher odds of diabetes.
- Non-Hispanic Asian and Mexican American participants showed elevated adjusted odds of diabetes.
- Female participants showed lower adjusted odds of diabetes.

## Repository Structure

- scripts/ – analysis code
- data/ – analytic dataset
- figures/ – visualizations
- outputs/ – regression results and summary tables
- presentation/ – project presentation

## Tools

R, tidyverse, nhanesA, gtsummary, broom, ggplot2, logistic regression.
