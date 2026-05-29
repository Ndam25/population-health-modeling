# Adama NDIR, MD, MSc — Population Health Modeling Portfolio

> **Physician-Epidemiologist | Applied Data Scientist | Global Health Equity**  
> 24 years | Sub-Saharan Africa | HIV · Nutrition · Infectious Disease · LMIC Surveillance

---

## Overview

This repository demonstrates **reproducible, documented, end-to-end analytical workflows** for population health modeling at the nutrition–disease nexus in low- and middle-income countries (LMICs).

All analyses are structured following IDM/Gates Foundation standards:
- Version-controlled pipelines (Git)
- Documented, modular code (R + Python)
- Synthetic data for full reproducibility without exposing sensitive survey data
- Explicit uncertainty quantification at every analytical stage

---

## Repository Structure

```
.
├── R/
│   ├── 01_data_harmonization.R       # DHS/MICS multi-country pipeline
│   ├── 02_burden_estimation.R        # PAF estimation (CRA framework)
│   ├── 03_seir_nutrition_model.R     # Compartmental model with nutritional covariates
│   ├── 04_forecasting_arima.R        # ARIMA + Bayesian structural time-series
│   └── 05_shiny_dashboard/           # Interactive results dashboard
│       ├── app.R
│       ├── ui.R
│       └── server.R
├── Python/
│   ├── 01_mice_imputation.py         # Multiple imputation (MICE) pipeline
│   ├── 02_spatial_sae.py             # Small-area estimation (SAE)
│   └── 03_ml_nutrition_classifier.py # NutriSnap-style food classifier (demo)
├── data/
│   └── synthetic/                    # Synthetic LMIC survey data (reproducible)
│       ├── generate_synthetic_dhs.R
│       └── README_data.md
├── docs/
│   ├── methods_overview.md
│   └── uncertainty_communication.md
└── README.md
```

---

## Key Analyses

### 1. Multi-Country DHS Data Harmonization Pipeline (`R/01_data_harmonization.R`)
- Harmonizes DHS/MICS rounds across 20+ countries
- Automated outlier detection and QA/QC flagging
- Sampling-weight-corrected complex survey estimation (`survey` package)
- **Reproducibility**: synthetic data seed included; full pipeline runs in < 5 min

### 2. Population-Attributable Fraction (PAF) Estimation (`R/02_burden_estimation.R`)
- Comparative Risk Assessment (CRA) framework adapted from GBD methodology
- Malnutrition PAFs on HIV, TB, diarrheal disease — stratified by sex, age, wealth quintile
- Bootstrap confidence intervals + probabilistic sensitivity analysis (PSA)

### 3. SEIR-Nutrition Compartmental Model (`R/03_seir_nutrition_model.R`)
- SEIR model extended with nutritional status compartments (stunting, wasting, adequate)
- Nutritional covariates modulate transmission and case-fatality parameters
- Calibrated against synthetic IeDEA-style cohort data
- 3-scenario projections: status quo / enhanced nutrition / integrated intervention

### 4. ARIMA + Bayesian Forecasting (`R/04_forecasting_arima.R`)
- 5-year disease burden forecasting with explicit uncertainty bands
- ARIMA, state-space, and `bsts` Bayesian structural time-series compared
- Monte Carlo simulation for policy scenario uncertainty propagation

### 5. MICE Missing Data Pipeline (`Python/01_mice_imputation.py`)
- Multiple Imputation by Chained Equations for LMIC surveillance gaps
- Bias quantification at each imputation stage
- Rubin's rules pooling with uncertainty propagation

### 6. Spatial Small-Area Estimation (`Python/02_spatial_sae.py`)
- Model-based geostatistics for subnational burden mapping
- Calibrated against DHS survey cluster coordinates (synthetic)

---

## Reproducibility Guarantee

| Component | Status |
|-----------|--------|
| Synthetic data generator | ✅ Seeded, fully reproducible |
| R environment | ✅ `renv.lock` included |
| Python environment | ✅ `requirements.txt` included |
| Random seeds | ✅ Set at all stochastic steps |
| Session info | ✅ Auto-logged at pipeline end |

To reproduce all analyses from scratch:
```bash
# R pipeline
Rscript data/synthetic/generate_synthetic_dhs.R
Rscript R/01_data_harmonization.R
Rscript R/02_burden_estimation.R
Rscript R/03_seir_nutrition_model.R

# Python pipeline
python Python/01_mice_imputation.py
python Python/02_spatial_sae.py
```

---

## Selected Methods

| Method | Package/Library | Application |
|--------|----------------|-------------|
| Complex survey estimation | `survey` (R) | DHS/MICS population inference |
| Mixed-effects survival models | `lme4`, `survival` (R) | IeDEA cohort mortality analysis |
| Multiple imputation (MICE) | `mice` (R), `sklearn` (Python) | LMIC surveillance gaps |
| Bayesian time-series | `bsts` (R) | Disease burden forecasting |
| Compartmental modeling | `deSolve` (R) | SEIR-Nutrition dynamics |
| Small-area estimation | `spdep`, `INLA` (R) | Subnational burden mapping |
| Sensitivity analysis | `sensitivity` (R) | PSA, Tornado diagrams |

---

## Background

**Adama NDIR, MD, MSc** is a physician-epidemiologist and applied data scientist with 24 years of experience in global health — 15 years as a US CDC Senior Scientist (Title 42). Core expertise:
- IeDEA West Africa cohort (17,000+ HIV patients, 6 countries) — *Published: The Lancet*
- Mali 2012 National DHS HIV component (statistical lead)
- CDC multi-country epidemic forecasting (HIV, COVID-19, Ebola, Mpox)
- MIT Applied Data Science Program (2024)
- NutriSnap AI nutrition platform (14,500 users; African food composition database)
- National Nutrition Policies — Burkina Faso & Mali

**Contact**: [LinkedIn](https://linkedin.com/in/adama-ndir) | ndir.adama@gmail.com  
**ORCID**: [0000-0000-0000-0000] *(placeholder)*

---

*All code and synthetic data in this repository are shared under MIT License. No real participant data is included.*
