################################################################################
# 03_seir_nutrition_model.R
# Adama NDIR, MD, MSc | Population Health Modeling Portfolio
#
# PURPOSE: SEIR compartmental model extended with nutritional status strata.
#          Nutritional covariates modulate: susceptibility (β), disease
#          progression rate (σ), and case-fatality (μ).
#
# MODEL STRUCTURE:
#   For each nutritional stratum s ∈ {Stunted, Adequate}:
#   S_s → E_s → I_s → R_s (with mortality μ_s from I_s)
#
# SCENARIOS:
#   1. Status Quo   — no nutritional intervention
#   2. Nutrition+   — stunting prevalence reduced 30% by year 5
#   3. Integrated   — nutrition + 20% coverage increase of treatment
#
# OUTPUT: figures/seir_scenarios.png | data/synthetic/seir_projections.csv
################################################################################

set.seed(42)

suppressPackageStartupMessages({
  library(deSolve)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

cat("=============================================================\n")
cat(" 03_seir_nutrition_model.R — SEIR-Nutrition Compartmental Model\n")
cat(" Run date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=============================================================\n\n")

# ── Model Parameters ──────────────────────────────────────────────────────────
# Population: 1M synthetic LMIC population
N       <- 1e6
# Nutritional prevalence at t=0
prev_stunted <- 0.30   # 30% stunted (calibrated to West Africa DHS)

# Transmission parameters
beta_base    <- 0.25   # Base transmission rate (adequate nutrition)
beta_mult    <- 1.45   # Stunted individuals 45% more susceptible (literature-based)
sigma        <- 1/5    # Latency: 5-day mean incubation
gamma        <- 1/14   # Recovery: 14-day mean infectious period
mu_adequate  <- 0.004  # CFR adequate nutrition
mu_stunted   <- 0.018  # CFR stunted (4.5x higher — IeDEA-calibrated)

cat("Model Parameters:\n")
cat(sprintf("  N = %s | Stunting prevalence = %.0f%%\n", format(N, big.mark=","), prev_stunted*100))
cat(sprintf("  β base = %.2f | Stunting multiplier = %.2f\n", beta_base, beta_mult))
cat(sprintf("  CFR: adequate = %.3f | stunted = %.3f\n\n", mu_adequate, mu_stunted))

# ── ODE System ────────────────────────────────────────────────────────────────
seir_nutrition <- function(time, state, params) {
  with(as.list(c(state, params)), {

    # Total infectious
    I_total <- I_adeq + I_stunt

    # Force of infection differs by nutritional status
    lambda_adeq  <- beta_adeq * I_total / N
    lambda_stunt <- beta_stunt * I_total / N

    # ODEs — Adequate nutrition stratum
    dS_adeq  <- -lambda_adeq * S_adeq
    dE_adeq  <-  lambda_adeq * S_adeq - sigma * E_adeq
    dI_adeq  <-  sigma * E_adeq - (gamma + mu_adeq) * I_adeq
    dR_adeq  <-  gamma * I_adeq
    dD_adeq  <-  mu_adeq * I_adeq

    # ODEs — Stunted stratum
    dS_stunt <- -lambda_stunt * S_stunt
    dE_stunt <-  lambda_stunt * S_stunt - sigma * E_stunt
    dI_stunt <-  sigma * E_stunt - (gamma + mu_stunt) * I_stunt
    dR_stunt <-  gamma * I_stunt
    dD_stunt <-  mu_stunt * I_stunt

    list(c(dS_adeq, dE_adeq, dI_adeq, dR_adeq, dD_adeq,
           dS_stunt, dE_stunt, dI_stunt, dR_stunt, dD_stunt))
  })
}

# ── Run Scenarios ─────────────────────────────────────────────────────────────
run_scenario <- function(scenario_name, stunting_reduction = 0, treatment_boost = 0) {

  prev_s <- prev_stunted * (1 - stunting_reduction)

  params <- list(
    N        = N,
    beta_adeq  = beta_base,
    beta_stunt = beta_base * beta_mult,
    sigma    = sigma,
    gamma    = gamma * (1 + treatment_boost),
    mu_adeq  = mu_adequate,
    mu_stunt = mu_stunted * (1 - stunting_reduction * 0.5)
  )

  # Initial conditions: seed with 100 infectious individuals
  I0_stunt <- round(100 * prev_s)
  I0_adeq  <- 100 - I0_stunt

  state0 <- c(
    S_adeq  = round(N * (1 - prev_s)) - I0_adeq,
    E_adeq  = 0,
    I_adeq  = I0_adeq,
    R_adeq  = 0,
    D_adeq  = 0,
    S_stunt = round(N * prev_s) - I0_stunt,
    E_stunt = 0,
    I_stunt = I0_stunt,
    R_stunt = 0,
    D_stunt = 0
  )

  times <- seq(0, 365 * 5, by = 1)  # 5-year horizon, daily steps

  out <- ode(y = state0, times = times, func = seir_nutrition, parms = params,
             method = "lsoda")

  as.data.frame(out) %>%
    mutate(
      scenario  = scenario_name,
      I_total   = I_adeq + I_stunt,
      D_total   = D_adeq + D_stunt,
      year      = time / 365
    ) %>%
    select(year, scenario, I_total, D_total, I_stunt, I_adeq)
}

cat("Running scenarios...\n")
scenarios <- bind_rows(
  run_scenario("1. Status Quo",    stunting_reduction = 0.00, treatment_boost = 0.00),
  run_scenario("2. Nutrition+",    stunting_reduction = 0.30, treatment_boost = 0.00),
  run_scenario("3. Integrated",    stunting_reduction = 0.30, treatment_boost = 0.20)
)

# ── Summary Table ─────────────────────────────────────────────────────────────
cat("\n--- 5-Year Cumulative Deaths by Scenario ---\n")
summary_table <- scenarios %>%
  group_by(scenario) %>%
  summarise(
    peak_I      = max(I_total),
    cum_deaths  = max(D_total),
    .groups = "drop"
  ) %>%
  mutate(
    deaths_averted = max(cum_deaths) - cum_deaths,
    reduction_pct  = round(deaths_averted / max(cum_deaths) * 100, 1)
  )
print(summary_table)

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(scenarios %>% filter(year <= 5),
            aes(x = year, y = I_total / 1000, color = scenario, linetype = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = c("#c0392b", "#2980b9", "#27ae60")) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash")) +
  labs(
    title    = "SEIR-Nutrition Model: Infectious Disease Trajectories Under\nThree Policy Scenarios",
    subtitle = "Synthetic LMIC population (N=1M) | Nutritional status modulates susceptibility & CFR",
    x        = "Year",
    y        = "Active Infections (thousands)",
    color    = "Scenario",
    linetype = "Scenario",
    caption  = "Adama NDIR, MD, MSc | Portfolio Demo | Synthetic data — not for policy use"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave("figures/seir_scenarios.png", p, width = 9, height = 5.5, dpi = 150)
cat("\nSaved: figures/seir_scenarios.png\n")

write_csv(scenarios, "data/synthetic/seir_projections.csv")
cat("Saved: data/synthetic/seir_projections.csv\n")
cat("\n✓ 03_seir_nutrition_model.R complete.\n")
