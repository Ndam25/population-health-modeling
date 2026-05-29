################################################################################
# generate_synthetic_dhs.R
# Adama NDIR, MD, MSc | Population Health Modeling Portfolio
#
# PURPOSE: Generate synthetic DHS-like multi-country survey data for
#          reproducible pipeline demonstration. No real participant data.
#
# REPRODUCIBILITY: Fixed seed (42) ensures identical output across runs.
# OUTPUT: data/synthetic/synthetic_dhs_pooled.csv
################################################################################

set.seed(42)

library(dplyr)
library(tibble)

# ── Parameters ─────────────────────────────────────────────────────────────────
COUNTRIES <- c("Senegal", "Mali", "Burkina Faso", "Guinea", "Nigeria",
               "Ghana", "Cote d'Ivoire", "Cameroon", "Chad", "Niger")
N_PER_COUNTRY <- 2000
YEARS <- c(2015, 2018, 2021)

# ── Country-level baseline parameters (approximate DHS calibration) ──────────
country_params <- tribble(
  ~country,        ~hiv_prev, ~stunting_prev, ~wasting_prev, ~anemia_prev, ~gdp_pc_usd,
  "Senegal",        0.030,     0.17,           0.08,           0.60,         1600,
  "Mali",           0.014,     0.30,           0.13,           0.72,         900,
  "Burkina Faso",   0.008,     0.27,           0.11,           0.68,         800,
  "Guinea",         0.016,     0.31,           0.09,           0.71,         950,
  "Nigeria",        0.013,     0.37,           0.07,           0.65,         2100,
  "Ghana",          0.017,     0.19,           0.05,           0.56,         2400,
  "Cote d'Ivoire",  0.028,     0.21,           0.06,           0.62,         1900,
  "Cameroon",       0.037,     0.32,           0.06,           0.64,         1500,
  "Chad",           0.019,     0.40,           0.16,           0.73,         700,
  "Niger",          0.007,     0.42,           0.18,           0.76,         550
)

# ── Simulation function ───────────────────────────────────────────────────────
simulate_country_wave <- function(country_name, year, params, n = N_PER_COUNTRY) {
  p <- params %>% filter(country == country_name)

  # Sampling weights (complex survey design: clusters of ~25)
  n_clusters <- ceiling(n / 25)
  cluster_id <- rep(1:n_clusters, each = 25)[1:n]
  stratum    <- sample(c("urban", "rural"), n, replace = TRUE, prob = c(0.45, 0.55))
  weight     <- ifelse(stratum == "urban", runif(n, 0.8, 1.2), runif(n, 0.7, 1.3))

  # Demographics
  sex       <- sample(c("female", "male"), n, replace = TRUE, prob = c(0.52, 0.48))
  age_years <- round(runif(n, 15, 49))
  wealth_q  <- sample(1:5, n, replace = TRUE, prob = c(0.25, 0.22, 0.20, 0.18, 0.15))

  # Nutritional status (correlated with wealth, sex, urban/rural)
  nutrition_score <- rnorm(n,
    mean = -0.5 * (wealth_q == 1) + 0.3 * (stratum == "urban"),
    sd = 1)

  stunted  <- as.integer(nutrition_score < qnorm(p$stunting_prev))
  wasted   <- as.integer(nutrition_score < qnorm(p$wasting_prev))
  anemic   <- as.integer(rbinom(n, 1, p$anemia_prev + 0.05 * (sex == "female") - 0.02 * (wealth_q >= 4)))

  # HIV status (modified by nutritional, wealth, sex)
  hiv_logit <- log(p$hiv_prev / (1 - p$hiv_prev)) +
    0.4 * (sex == "female") +
    0.6 * stunted +
    -0.3 * log(wealth_q) +
    0.2 * (stratum == "rural")
  hiv_prob  <- plogis(hiv_logit)
  hiv_pos   <- rbinom(n, 1, hiv_prob)

  # Outcome: under-5 mortality proxy (child-level estimate at household)
  u5mr_proxy <- rbinom(n, 1,
    p$wasting_prev * 2.5 * (1 + 0.8 * (wealth_q == 1)) * (year < 2018) * 0.3)

  # Missing data pattern (realistic LMIC missingness)
  miss_hiv      <- rbinom(n, 1, 0.12)   # 12% missing HIV test
  miss_nutr     <- rbinom(n, 1, 0.07)   # 7% missing anthropometrics
  hiv_pos[miss_hiv == 1] <- NA
  stunted[miss_nutr == 1] <- NA
  wasted[miss_nutr == 1]  <- NA

  tibble(
    country      = country_name,
    year         = year,
    cluster_id   = paste0(substr(country_name, 1, 3), "_", year, "_", cluster_id),
    stratum      = stratum,
    weight       = round(weight, 4),
    sex          = sex,
    age_years    = age_years,
    wealth_q     = wealth_q,
    hiv_positive = hiv_pos,
    stunted      = stunted,
    wasted       = wasted,
    anemic       = anemic,
    u5mr_proxy   = u5mr_proxy,
    gdp_pc_usd   = p$gdp_pc_usd
  )
}

# ── Generate full dataset ─────────────────────────────────────────────────────
cat("Generating synthetic DHS data...\n")
cat("Countries:", length(COUNTRIES), "| Waves:", length(YEARS),
    "| N per cell:", N_PER_COUNTRY, "\n")

all_data <- purrr::map_dfr(COUNTRIES, function(ctry) {
  purrr::map_dfr(YEARS, function(yr) {
    simulate_country_wave(ctry, yr, country_params)
  })
})

cat("Total observations:", nrow(all_data), "\n")
cat("Missing HIV:", sum(is.na(all_data$hiv_positive)), "(",
    round(mean(is.na(all_data$hiv_positive)) * 100, 1), "%)\n")
cat("Missing stunting:", sum(is.na(all_data$stunted)), "(",
    round(mean(is.na(all_data$stunted)) * 100, 1), "%)\n")

# ── Save ──────────────────────────────────────────────────────────────────────
dir.create("data/synthetic", recursive = TRUE, showWarnings = FALSE)
readr::write_csv(all_data, "data/synthetic/synthetic_dhs_pooled.csv")
cat("\nSaved: data/synthetic/synthetic_dhs_pooled.csv\n")

# Session info for reproducibility audit
cat("\n--- Session Info ---\n")
print(sessionInfo())
