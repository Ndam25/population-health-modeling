"""
01_mice_imputation.py
Adama NDIR, MD, MSc | Population Health Modeling Portfolio

PURPOSE: Multiple Imputation by Chained Equations (MICE) pipeline for
         LMIC surveillance data gaps. Implements:
         - Missingness pattern analysis (MAR/MNAR diagnostics)
         - MICE imputation with m=20 datasets
         - Rubin's rules pooling for valid inference under MI
         - Bias quantification at each imputation stage

REPRODUCIBILITY: Random seed fixed; all stages logged.

INPUT:   data/synthetic/synthetic_dhs_pooled.csv
OUTPUT:  data/synthetic/imputed_results.csv
         data/synthetic/mice_diagnostics.png
"""

import numpy as np
import pandas as pd
from sklearn.experimental import enable_iterative_imputer   # noqa
from sklearn.impute import IterativeImputer
from sklearn.linear_model import BayesianRidge
from sklearn.preprocessing import LabelEncoder
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings
import sys
import os

warnings.filterwarnings("ignore")

RNG_SEED = 42
np.random.seed(RNG_SEED)
M_IMPUTATIONS = 20   # Rubin's rule: ≥20 for small fraction of missing

print("=" * 60)
print(" 01_mice_imputation.py — MICE Pipeline")
print(f" Python {sys.version.split()[0]} | numpy {np.__version__} | pandas {pd.__version__}")
print("=" * 60, "\n")

# ─── 1. Load data ─────────────────────────────────────────────────────────────
df_raw = pd.read_csv("data/synthetic/synthetic_dhs_pooled.csv")
print(f"Loaded: {df_raw.shape[0]:,} rows × {df_raw.shape[1]} columns")

# Encode categoricals for imputation
le_sex = LabelEncoder()
le_str = LabelEncoder()
df = df_raw.copy()
df["sex_enc"]    = le_sex.fit_transform(df["sex"])
df["strat_enc"]  = le_str.fit_transform(df["stratum"])

NUMERIC_VARS = ["age_years", "wealth_q", "sex_enc", "strat_enc",
                "hiv_positive", "stunted", "wasted", "anemic", "gdp_pc_usd"]

df_num = df[NUMERIC_VARS].copy()

# ─── 2. Missingness diagnostics ───────────────────────────────────────────────
print("\n--- Missingness Pattern ---")
miss_pct = df_num.isnull().mean() * 100
for col, pct in miss_pct[miss_pct > 0].items():
    print(f"  {col:20s}: {pct:.1f}% missing")

# Little's MCAR test approximation: correlate missingness indicators
print("\n--- MAR Diagnostics (Missingness Correlations) ---")
miss_indicators = df_num.isnull().astype(int)
miss_corr = miss_indicators.corr()
# Flag pairs with |r| > 0.1
for col1 in miss_indicators.columns:
    for col2 in miss_indicators.columns:
        if col1 < col2:
            r = miss_corr.loc[col1, col2]
            if abs(r) > 0.1:
                print(f"  Missingness correlation [{col1}] ↔ [{col2}]: r = {r:.3f}")
print("  (High correlations suggest MAR structure — MICE appropriate)")

# ─── 3. MICE imputation (m=20 datasets) ──────────────────────────────────────
print(f"\n--- MICE Imputation (m = {M_IMPUTATIONS}) ---")

# Using sklearn's IterativeImputer = MICE with BayesianRidge estimator
all_imputed = []

for m in range(M_IMPUTATIONS):
    imputer = IterativeImputer(
        estimator        = BayesianRidge(),
        max_iter         = 10,
        random_state     = RNG_SEED + m,  # Different seed per imputation
        initial_strategy = "mean",
        imputation_order = "roman"
    )
    df_imp = pd.DataFrame(
        imputer.fit_transform(df_num),
        columns = NUMERIC_VARS
    )
    # Clamp binary variables to [0, 1] and round
    for bvar in ["hiv_positive", "stunted", "wasted", "anemic"]:
        df_imp[bvar] = df_imp[bvar].clip(0, 1).round()
    df_imp["imputation"] = m + 1
    all_imputed.append(df_imp)

    if (m + 1) % 5 == 0:
        print(f"  Completed imputation {m+1}/{M_IMPUTATIONS}")

print("  ✓ All imputations complete")

# ─── 4. Rubin's Rules pooling ──────────────────────────────────────────────────
print("\n--- Rubin's Rules Pooling ---")

def rubins_rules(estimates, variances):
    """
    Pool m estimates using Rubin's Rules (1987).
    Q̄ = mean of m estimates
    W = mean within-imputation variance
    B = between-imputation variance
    T = W + (1 + 1/m) * B
    """
    m    = len(estimates)
    Qbar = np.mean(estimates)
    W    = np.mean(variances)
    B    = np.var(estimates, ddof=1)
    T    = W + (1 + 1/m) * B
    se   = np.sqrt(T)
    # Barnard-Rubin degrees of freedom (simplified)
    df_r = (m - 1) * (1 + W / ((1 + 1/m) * B))**2
    return {"estimate": Qbar, "se": se, "variance_total": T,
            "between_var": B, "within_var": W, "df": df_r}

outcomes = ["hiv_positive", "stunted", "wasted", "anemic"]
pooled_results = {}

for outcome in outcomes:
    ests = [df[outcome].mean() for df in all_imputed]
    # Variance of a proportion: p*(1-p)/n
    n_obs = len(df_num.dropna(subset=[outcome]))
    vars_ = [est * (1 - est) / n_obs for est in ests]
    res   = rubins_rules(ests, vars_)
    pooled_results[outcome] = res
    ci_lo = res["estimate"] - 1.96 * res["se"]
    ci_hi = res["estimate"] + 1.96 * res["se"]
    print(f"  {outcome:20s}: {res['estimate']:.3f} [95% CI: {ci_lo:.3f}–{ci_hi:.3f}]"
          f" | Between-imp variance: {res['between_var']:.6f}")

# ─── 5. Bias quantification ───────────────────────────────────────────────────
print("\n--- Bias from Complete-Case Analysis (CCA) vs. MICE ---")
cca_results = {}
for outcome in outcomes:
    cca_est = df_num[outcome].mean()  # complete-case (drops NAs)
    mi_est  = pooled_results[outcome]["estimate"]
    bias    = mi_est - cca_est
    cca_results[outcome] = {"CCA": cca_est, "MI": mi_est, "Bias": bias}
    print(f"  {outcome:20s}: CCA = {cca_est:.3f} | MI = {mi_est:.3f} | Bias = {bias:+.3f}")

# ─── 6. Convergence diagnostic plot ──────────────────────────────────────────
print("\n--- Generating Diagnostics Plot ---")
os.makedirs("figures", exist_ok=True)

fig = plt.figure(figsize=(12, 8))
fig.suptitle("MICE Imputation Diagnostics\nAdama NDIR | Population Health Modeling Portfolio",
             fontsize=13, fontweight="bold")
gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.35)

# Panel 1: Convergence — mean of imputed HIV across 20 imputations
ax1 = fig.add_subplot(gs[0, 0])
hiv_means = [df["hiv_positive"].mean() for df in all_imputed]
ax1.plot(range(1, M_IMPUTATIONS + 1), hiv_means, "o-", color="#2c6e8a", markersize=4)
ax1.axhline(np.mean(hiv_means), color="#c0392b", linestyle="--", label="Pooled mean")
ax1.set_xlabel("Imputation number")
ax1.set_ylabel("Mean HIV prevalence")
ax1.set_title("Convergence: HIV (m=20)")
ax1.legend(fontsize=8)

# Panel 2: Between-imputation variance
ax2 = fig.add_subplot(gs[0, 1])
bvar = [pooled_results[o]["between_var"] for o in outcomes]
ax2.barh(outcomes, bvar, color="#27ae60", alpha=0.8)
ax2.set_xlabel("Between-imputation variance")
ax2.set_title("Between-Imp. Variance (Rubin's B)")

# Panel 3: CCA vs MI bias
ax3 = fig.add_subplot(gs[1, 0])
x = np.arange(len(outcomes))
width = 0.35
ax3.bar(x - width/2, [cca_results[o]["CCA"] for o in outcomes],
        width, label="CCA (complete-case)", color="#e67e22", alpha=0.8)
ax3.bar(x + width/2, [cca_results[o]["MI"] for o in outcomes],
        width, label="MICE pooled", color="#2c6e8a", alpha=0.8)
ax3.set_xticks(x); ax3.set_xticklabels(outcomes, rotation=15)
ax3.set_ylabel("Prevalence estimate")
ax3.set_title("CCA vs. MICE Estimates")
ax3.legend(fontsize=8)

# Panel 4: Missingness pattern
ax4 = fig.add_subplot(gs[1, 1])
miss_pct_plot = df_num.isnull().mean() * 100
ax4.barh(miss_pct_plot.index, miss_pct_plot.values, color="#8e44ad", alpha=0.8)
ax4.set_xlabel("% Missing")
ax4.set_title("Missingness Pattern")
ax4.axvline(5, color="red", linestyle="--", linewidth=0.8, label="5% threshold")
ax4.legend(fontsize=8)

plt.savefig("figures/mice_diagnostics.png", dpi=150, bbox_inches="tight")
print("  Saved: figures/mice_diagnostics.png")

# ─── 7. Save pooled dataset ───────────────────────────────────────────────────
# Stack all M imputed datasets (long format, standard for MI analysis)
mi_long = pd.concat(all_imputed, ignore_index=True)
mi_long.to_csv("data/synthetic/imputed_results.csv", index=False)
print("Saved: data/synthetic/imputed_results.csv")

print("\n✓ 01_mice_imputation.py complete.")
