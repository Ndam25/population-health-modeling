################################################################################
# 04_forecasting_arima.R
# Adama NDIR, MD, MSc | Population Health Modeling Portfolio
#
# PURPOSE: 5-year disease burden forecasting using ARIMA, state-space, and
#          Bayesian structural time-series (bsts). Explicit uncertainty
#          communication via prediction intervals and Monte Carlo simulation.
#
# METHODS COMPARED:
#   1. ARIMA (auto.arima)
#   2. State-space model (StructTS)
#   3. Bayesian Structural Time Series (bsts)
#
# OUTPUT: figures/forecast_comparison.png | data/synthetic/forecasts.csv
################################################################################

set.seed(42)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(forecast)    # auto.arima
  library(bsts)        # Bayesian structural time-series
})

cat("=============================================================\n")
cat(" 04_forecasting_arima.R — Disease Burden Forecasting\n")
cat(" Run date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=============================================================\n\n")

# ── Generate synthetic time-series (calibrated to West Africa HIV) ─────────────
# Simulates annual HIV incidence rate per 1000 (2000–2024), 25 years
years_obs  <- 2000:2024
n_obs      <- length(years_obs)

# Realistic declining trend with programmatic inflections
true_trend <- c(
  seq(28, 22, length.out = 5),   # 2000-2004: plateau/slight decline (pre-PEPFAR)
  seq(22, 14, length.out = 8),   # 2005-2012: PEPFAR scale-up
  seq(14, 10, length.out = 5),   # 2013-2017: sustained decline
  seq(10,  8, length.out = 7)    # 2018-2024: slower progress
)
noise      <- rnorm(n_obs, 0, 0.6)
hiv_series <- ts(true_trend + noise, start = 2000, frequency = 1)

cat("Observed series: HIV incidence/1000 | 2000–2024\n")
cat(sprintf("  Range: %.1f – %.1f | Mean: %.1f\n\n",
            min(hiv_series), max(hiv_series), mean(hiv_series)))

# ── 1. ARIMA ──────────────────────────────────────────────────────────────────
cat("--- 1. Auto-ARIMA ---\n")
fit_arima <- auto.arima(hiv_series, ic = "aic", stepwise = FALSE,
                         approximation = FALSE)
cat("  Best model:", as.character(fit_arima), "\n")
fcast_arima <- forecast(fit_arima, h = 5, level = c(80, 95))
cat(sprintf("  2025 forecast: %.2f [80%% PI: %.2f–%.2f]\n",
            fcast_arima$mean[1],
            fcast_arima$lower[1, 2],
            fcast_arima$upper[1, 2]))

# ── 2. State-space (local level + trend) ──────────────────────────────────────
cat("\n--- 2. State-Space (StructTS) ---\n")
fit_ss   <- StructTS(hiv_series, type = "trend")
fcast_ss <- forecast(fit_ss, h = 5, level = c(80, 95))
cat(sprintf("  2025 forecast: %.2f [80%% PI: %.2f–%.2f]\n",
            fcast_ss$mean[1],
            fcast_ss$lower[1, 2],
            fcast_ss$upper[1, 2]))

# ── 3. Bayesian Structural Time Series ───────────────────────────────────────
cat("\n--- 3. Bayesian Structural Time Series (bsts) ---\n")
ss_spec <- list()
ss_spec <- AddLocalLinearTrend(ss_spec, y = as.numeric(hiv_series))
# Suppress verbose bsts output
invisible(capture.output(
  fit_bsts <- bsts(as.numeric(hiv_series), state.specification = ss_spec,
                   niter = 2000, seed = 42, ping = 0)
))
pred_bsts <- predict(fit_bsts, horizon = 5, quantiles = c(0.025, 0.10, 0.90, 0.975))
bsts_mean <- colMeans(pred_bsts$distribution)
bsts_lo95 <- pred_bsts$interval[1, ]
bsts_hi95 <- pred_bsts$interval[2, ]
cat(sprintf("  2025 forecast: %.2f [95%% CrI: %.2f–%.2f]\n",
            bsts_mean[1], bsts_lo95[1], bsts_hi95[1]))

# ── 4. Monte Carlo uncertainty propagation ────────────────────────────────────
cat("\n--- 4. Monte Carlo Sensitivity (n=1000 trajectories) ---\n")
mc_forecasts <- replicate(1000, {
  noise_mc <- rnorm(n_obs, 0, 0.8)  # Perturbed series
  ts_mc    <- ts(true_trend + noise_mc, start = 2000, frequency = 1)
  tryCatch({
    fc <- forecast(auto.arima(ts_mc, stepwise = TRUE), h = 5)
    as.numeric(fc$mean)
  }, error = function(e) rep(NA, 5))
})
mc_means <- apply(mc_forecasts, 1, mean, na.rm = TRUE)
mc_lo    <- apply(mc_forecasts, 1, quantile, 0.025, na.rm = TRUE)
mc_hi    <- apply(mc_forecasts, 1, quantile, 0.975, na.rm = TRUE)
cat(sprintf("  2025 MC median: %.2f [95%% UI: %.2f–%.2f]\n",
            mc_means[1], mc_lo[1], mc_hi[1]))

# ── 5. Compile and plot ───────────────────────────────────────────────────────
years_fcast <- 2025:2029

obs_df <- tibble(
  year   = years_obs,
  value  = as.numeric(hiv_series),
  type   = "Observed"
)

fcast_df <- bind_rows(
  tibble(year = years_fcast, value = as.numeric(fcast_arima$mean),
         lo80 = fcast_arima$lower[, 1], hi80 = fcast_arima$upper[, 1],
         lo95 = fcast_arima$lower[, 2], hi95 = fcast_arima$upper[, 2],
         model = "ARIMA"),
  tibble(year = years_fcast, value = as.numeric(fcast_ss$mean),
         lo80 = fcast_ss$lower[, 1], hi80 = fcast_ss$upper[, 1],
         lo95 = fcast_ss$lower[, 2], hi95 = fcast_ss$upper[, 2],
         model = "State-Space"),
  tibble(year = years_fcast, value = bsts_mean,
         lo80 = pred_bsts$interval[1,], hi80 = pred_bsts$interval[2,],
         lo95 = bsts_lo95, hi95 = bsts_hi95,
         model = "Bayesian (bsts)")
)

p <- ggplot() +
  # Observed
  geom_line(data = obs_df, aes(x = year, y = value),
            color = "black", linewidth = 1.1) +
  geom_point(data = obs_df, aes(x = year, y = value),
             color = "black", size = 1.5) +
  # Forecast ribbons and lines
  geom_ribbon(data = fcast_df,
              aes(x = year, ymin = lo95, ymax = hi95, fill = model),
              alpha = 0.12) +
  geom_ribbon(data = fcast_df,
              aes(x = year, ymin = lo80, ymax = hi80, fill = model),
              alpha = 0.20) +
  geom_line(data = fcast_df,
            aes(x = year, y = value, color = model),
            linewidth = 1.0, linetype = "dashed") +
  # MC uncertainty envelope
  geom_ribbon(data = tibble(year = years_fcast, lo = mc_lo, hi = mc_hi),
              aes(x = year, ymin = lo, ymax = hi),
              fill = "grey50", alpha = 0.10) +
  scale_color_manual(values = c("ARIMA" = "#c0392b", "State-Space" = "#2980b9",
                                "Bayesian (bsts)" = "#27ae60")) +
  scale_fill_manual(values  = c("ARIMA" = "#c0392b", "State-Space" = "#2980b9",
                                "Bayesian (bsts)" = "#27ae60")) +
  geom_vline(xintercept = 2024.5, linetype = "dotted", color = "grey40") +
  labs(
    title    = "HIV Incidence Forecasting: ARIMA vs. State-Space vs. Bayesian (bsts)",
    subtitle = "Synthetic West Africa HIV incidence/1000 | 2000–2024 observed | 2025–2029 forecast\nShaded = 80% & 95% prediction intervals | Grey = Monte Carlo uncertainty envelope (n=1000)",
    x        = "Year", y = "HIV Incidence per 1,000 population",
    color    = "Model", fill = "Model",
    caption  = "Adama NDIR, MD, MSc | Portfolio Demo | Synthetic data only"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

ggsave("figures/forecast_comparison.png", p, width = 10, height = 6, dpi = 150)
cat("\nSaved: figures/forecast_comparison.png\n")

write_csv(fcast_df, "data/synthetic/forecasts.csv")
cat("Saved: data/synthetic/forecasts.csv\n")
cat("\n✓ 04_forecasting_arima.R complete.\n")
