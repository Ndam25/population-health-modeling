################################################################################
# 01_data_harmonization.R
# Adama NDIR, MD, MSc | Population Health Modeling Portfolio
#
# PURPOSE: Harmonize synthetic multi-country DHS data into an analysis-ready
#          dataset. Implements automated QA/QC, complex survey weighting, and
#          an audit-trail-compatible output.
#
# INPUT:   data/synthetic/synthetic_dhs_pooled.csv
# OUTPUT:  data/synthetic/harmonized_dhs.csv
#          data/synthetic/qaqc_report.txt
#
# REPRODUCIBILITY: Set seed at top; all steps logged.
################################################################################

set.seed(42)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(survey)
  library(mice)     # for quick missing-data summary
  library(ggplot2)
})

cat("=============================================================\n")
cat(" 01_data_harmonization.R — Adama NDIR Portfolio\n")
cat(" Run date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=============================================================\n\n")

# ── 1. Load data ──────────────────────────────────────────────────────────────
dhs_raw <- read_csv("data/synthetic/synthetic_dhs_pooled.csv",
                    show_col_types = FALSE)
cat("Loaded:", nrow(dhs_raw), "rows,", ncol(dhs_raw), "columns\n")
cat("Countries:", paste(unique(dhs_raw$country), collapse = ", "), "\n")
cat("Years:", paste(sort(unique(dhs_raw$year)), collapse = ", "), "\n\n")

# ── 2. Automated QA/QC ───────────────────────────────────────────────────────
cat("--- QA/QC CHECKS ---\n")

qaqc_log <- list()

# 2a. Range validation
range_checks <- list(
  age_years  = c(15, 49),
  wealth_q   = c(1, 5),
  weight     = c(0.1, 5.0)
)
for (var in names(range_checks)) {
  lo <- range_checks[[var]][1]; hi <- range_checks[[var]][2]
  n_out <- sum(dhs_raw[[var]] < lo | dhs_raw[[var]] > hi, na.rm = TRUE)
  msg <- sprintf("  [%s] Out-of-range: %d (%.2f%%)", var, n_out,
                 100 * n_out / nrow(dhs_raw))
  cat(msg, "\n")
  qaqc_log[[var]] <- msg
}

# 2b. Missing data summary
cat("\n  Missing data summary:\n")
miss_summary <- dhs_raw %>%
  summarise(across(everything(), ~ mean(is.na(.)) * 100)) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "pct_missing") %>%
  filter(pct_missing > 0)
print(miss_summary, n = Inf)

# 2c. Duplicate cluster check
n_clusters <- n_distinct(dhs_raw$cluster_id)
cat(sprintf("\n  Unique clusters: %d\n", n_clusters))

# 2d. Plausibility: HIV prevalence by country
cat("\n  HIV prevalence by country (crude):\n")
hiv_check <- dhs_raw %>%
  group_by(country) %>%
  summarise(
    hiv_prev_crude = mean(hiv_positive, na.rm = TRUE),
    n_tested = sum(!is.na(hiv_positive)),
    .groups = "drop"
  ) %>%
  mutate(flag = ifelse(hiv_prev_crude > 0.10 | hiv_prev_crude < 0.001, "⚠ CHECK", "OK"))
print(hiv_check)

# ── 3. Standardization ───────────────────────────────────────────────────────
cat("\n--- STANDARDIZATION ---\n")

dhs_clean <- dhs_raw %>%
  mutate(
    # Harmonize sex coding
    sex_binary = as.integer(sex == "female"),
    # Wealth quintile as ordered factor
    wealth_q   = factor(wealth_q, levels = 1:5, ordered = TRUE),
    # Urban/rural binary
    urban      = as.integer(stratum == "urban"),
    # Double burden indicator
    stunted_anemic = as.integer(stunted == 1 & anemic == 1),
    # Survey year as categorical
    survey_wave = factor(year)
  )

cat("  Created derived variables: sex_binary, urban, stunted_anemic, survey_wave\n")

# ── 4. Complex survey design ─────────────────────────────────────────────────
cat("\n--- COMPLEX SURVEY DESIGN ---\n")

# Design: stratified cluster sampling with probability weights
svy_design <- svydesign(
  ids     = ~cluster_id,    # PSU (primary sampling unit)
  strata  = ~stratum,       # Urban/rural strata
  weights = ~weight,
  data    = dhs_clean,
  nest    = TRUE
)
cat("  Survey design specified: stratified cluster, sampling weights applied\n")

# Population-weighted estimates by country
cat("\n  Population-weighted HIV prevalence (%):\n")
hiv_svy <- svyby(~hiv_positive, ~country, svy_design, svymean,
                  na.rm = TRUE, vartype = c("se", "ci"))
hiv_svy <- hiv_svy %>%
  mutate(across(where(is.numeric), ~ round(. * 100, 2)))
print(hiv_svy[, c("country", "hiv_positive", "se", "ci_l", "ci_u")])

# ── 5. Output ─────────────────────────────────────────────────────────────────
write_csv(dhs_clean, "data/synthetic/harmonized_dhs.csv")
cat("\nSaved: data/synthetic/harmonized_dhs.csv\n")

# Write QA/QC report
sink("data/synthetic/qaqc_report.txt")
cat("QA/QC Report — 01_data_harmonization.R\n")
cat("Run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("=== Range Checks ===\n")
for (msg in qaqc_log) cat(msg, "\n")
cat("\n=== Missing Data ===\n")
print(miss_summary)
cat("\n=== HIV Plausibility ===\n")
print(hiv_check)
cat("\n=== Session Info ===\n")
print(sessionInfo())
sink()
cat("Saved: data/synthetic/qaqc_report.txt\n")

cat("\n✓ 01_data_harmonization.R complete.\n")
