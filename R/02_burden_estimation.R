################################################################################
# 02_burden_estimation.R
# Adama NDIR, MD, MSc | Population Health Modeling Portfolio
#
# PURPOSE: Estimate population-attributable fractions (PAFs) of malnutrition
#          on infectious disease outcomes using a Comparative Risk Assessment
#          (CRA) framework adapted from GBD methodology.
#
#          Stratified by: sex, age group, wealth quintile
#          Outcomes:      HIV positivity, anemia (as proxy for disease burden)
#          Exposure:      Stunting, wasting (malnutrition indicators)
#
# METHODS: Logistic regression → adjusted OR → Levin's PAF formula
#          Bootstrap CI (n=500) + Probabilistic Sensitivity Analysis (PSA)
#
# INPUT:   data/synthetic/harmonized_dhs.csv
# OUTPUT:  data/synthetic/paf_results.csv | figures/paf_tornado.png
################################################################################

set.seed(42)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(broom)
  library(ggplot2)
  library(purrr)
})

cat("=============================================================\n")
cat(" 02_burden_estimation.R — PAF Estimation (CRA Framework)\n")
cat(" Run date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=============================================================\n\n")

# ── Load harmonized data ──────────────────────────────────────────────────────
dhs <- read_csv("data/synthetic/harmonized_dhs.csv", show_col_types = FALSE) %>%
  mutate(
    age_group = cut(age_years,
                    breaks = c(14, 24, 34, 44, 50),
                    labels = c("15-24", "25-34", "35-44", "45-49")),
    wealth_binary = as.integer(as.integer(wealth_q) <= 2)  # Q1+Q2 = poorest 40%
  ) %>%
  filter(!is.na(hiv_positive), !is.na(stunted), !is.na(wasted))

cat("Analysis sample (complete cases):", nrow(dhs), "\n\n")

# ── Levin's PAF formula ───────────────────────────────────────────────────────
# PAF = Pe * (OR - 1) / (Pe * (OR - 1) + 1)
# Pe = prevalence of exposure in population
# OR = adjusted odds ratio from logistic regression
levin_paf <- function(or, pe) {
  (pe * (or - 1)) / (pe * (or - 1) + 1)
}

# ── Core PAF estimation function ─────────────────────────────────────────────
estimate_paf <- function(data, outcome_var, exposure_var,
                         covariates = c("sex_binary", "age_group", "wealth_binary", "urban")) {

  formula_str <- paste(outcome_var, "~", exposure_var, "+",
                       paste(covariates, collapse = " + "))
  model <- glm(as.formula(formula_str), data = data, family = binomial())

  or       <- exp(coef(model)[exposure_var])
  pe       <- mean(data[[exposure_var]], na.rm = TRUE)
  paf      <- levin_paf(or, pe)
  ci_model <- confint(model, exposure_var, level = 0.95)
  or_lo    <- exp(ci_model[1]); or_hi <- exp(ci_model[2])

  list(
    outcome  = outcome_var,
    exposure = exposure_var,
    OR       = or,
    OR_lo    = or_lo,
    OR_hi    = or_hi,
    Pe       = pe,
    PAF      = paf,
    PAF_lo   = levin_paf(or_lo, pe),
    PAF_hi   = levin_paf(or_hi, pe)
  )
}

# ── 1. Overall PAF estimates ──────────────────────────────────────────────────
cat("--- Overall PAF Estimates ---\n")
combinations <- expand.grid(
  outcome  = c("hiv_positive", "anemic"),
  exposure = c("stunted", "wasted"),
  stringsAsFactors = FALSE
)

paf_results <- pmap_dfr(combinations, function(outcome, exposure) {
  as_tibble(estimate_paf(dhs, outcome, exposure))
})

paf_results <- paf_results %>%
  mutate(across(c(OR, OR_lo, OR_hi, Pe, PAF, PAF_lo, PAF_hi), ~ round(., 3)))

print(paf_results)

# ── 2. Stratified PAF by sex ──────────────────────────────────────────────────
cat("\n--- PAF Stratified by Sex ---\n")
paf_sex <- dhs %>%
  group_by(sex) %>%
  group_map(~ {
    r <- estimate_paf(.x, "hiv_positive", "stunted",
                      covariates = c("age_group", "wealth_binary", "urban"))
    as_tibble(c(sex = unique(.y$sex), r))
  }) %>%
  bind_rows()

print(paf_sex %>% select(sex, OR, Pe, PAF, PAF_lo, PAF_hi) %>%
      mutate(across(where(is.numeric), ~ round(., 3))))

# ── 3. Bootstrap confidence intervals ─────────────────────────────────────────
cat("\n--- Bootstrap CIs (n=500) for PAF: stunting → HIV ---\n")
boot_pafs <- replicate(500, {
  boot_data <- dhs[sample(nrow(dhs), replace = TRUE), ]
  tryCatch({
    r <- estimate_paf(boot_data, "hiv_positive", "stunted")
    r$PAF
  }, error = function(e) NA)
})
boot_pafs <- na.omit(boot_pafs)
cat(sprintf("  PAF (stunting → HIV): %.3f [Bootstrap 95%% CI: %.3f – %.3f]\n",
            median(boot_pafs), quantile(boot_pafs, 0.025), quantile(boot_pafs, 0.975)))

# ── 4. Probabilistic Sensitivity Analysis (PSA) ──────────────────────────────
cat("\n--- PSA: Uncertainty in Prevalence Estimates (n=1000) ---\n")
# Simulate uncertainty in stunting prevalence (e.g., ±10% survey error)
psa_pe    <- rnorm(1000, mean = mean(dhs$stunted, na.rm = TRUE), sd = 0.03)
psa_or    <- rlnorm(1000, meanlog = log(1.45), sdlog = 0.15)  # OR uncertainty
psa_pafs  <- levin_paf(psa_or, psa_pe)
cat(sprintf("  PSA PAF (stunting → HIV): median %.3f [95%% uncertainty: %.3f – %.3f]\n",
            median(psa_pafs), quantile(psa_pafs, 0.025), quantile(psa_pafs, 0.975)))

# ── 5. Tornado diagram ────────────────────────────────────────────────────────
cat("\n--- Generating Tornado Diagram ---\n")
tornado_data <- paf_results %>%
  mutate(
    label    = paste(exposure, "→", outcome),
    PAF_pct  = PAF * 100,
    lo_pct   = PAF_lo * 100,
    hi_pct   = PAF_hi * 100,
    range    = hi_pct - lo_pct
  ) %>%
  arrange(desc(range))

p_tornado <- ggplot(tornado_data, aes(x = reorder(label, PAF_pct),
                                       y = PAF_pct, ymin = lo_pct, ymax = hi_pct)) +
  geom_col(fill = "#2c6e8a", alpha = 0.85) +
  geom_errorbar(width = 0.25, color = "#1a3f50", linewidth = 0.8) +
  coord_flip() +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(
    title    = "Population-Attributable Fractions of Malnutrition\non Infectious Disease Burden",
    subtitle = "Synthetic DHS data | 10 ECOWAS countries | Adjusted OR, Levin's PAF",
    x        = NULL,
    y        = "PAF (%) with 95% CI",
    caption  = "Adama NDIR, MD, MSc | Portfolio Demo | Synthetic data only"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

dir.create("figures", showWarnings = FALSE)
ggsave("figures/paf_tornado.png", p_tornado, width = 8, height = 5, dpi = 150)
cat("  Saved: figures/paf_tornado.png\n")

# ── Save results ──────────────────────────────────────────────────────────────
write_csv(paf_results, "data/synthetic/paf_results.csv")
cat("\nSaved: data/synthetic/paf_results.csv\n")
cat("\n✓ 02_burden_estimation.R complete.\n")
