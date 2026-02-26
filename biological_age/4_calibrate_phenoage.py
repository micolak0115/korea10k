# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression

# %%
# clinical = "age_wbc_ttbl_alp_monopa_lymph_albumin_glucose_uap_bun_rbc_creat"
# clinical = "creat_age_glucose_sbp_rdw_bun"
clinical = "sbp_glucose_age_rdw_creat_bun"
path_phenoage = f"/BiO/Research/GeroExpressome/Results/PhenoAge4/KGP_{clinical}_phenoage_added.txt"
phenoage = pd.read_csv(path_phenoage, sep=",").set_index("Sample_ID_New")
phenoage = phenoage.drop(columns=["Unnamed: 0"])
phenoage = phenoage.dropna(subset=["phenoage"])

# %%
# --- Step 1: Fit calibration regression (phenoage ~ age) ---
X = phenoage["age"].values.reshape(-1, 1)
y = phenoage["phenoage"].values

model = LinearRegression().fit(X, y)
intercept, slope = model.intercept_, model.coef_[0]
print(f"Calibration model: phenoage = {intercept:.3f} + {slope:.3f} * age")

# %%
# --- Step 2: Apply bias correction ---
# Adjust phenoage so that it aligns with chronological age (i.e., slope=1, intercept=0)
phenoage["phenoage_calibrated"] = (phenoage["phenoage"] - intercept) / slope

# --- Step 3: Check correlation before vs after ---
r_before = np.corrcoef(phenoage["age"], phenoage["phenoage"])[0, 1]
r_after  = np.corrcoef(phenoage["age"], phenoage["phenoage_calibrated"])[0, 1]
print(f"Correlation before: {r_before:.3f}")
print(f"Correlation after : {r_after:.3f}")

# %%
# --- Step 4: Visualize before vs after ---
plt.figure(figsize=(6,5))

# Original (gray)
plt.scatter(phenoage["age"], phenoage["phenoage"], facecolors='none', edgecolors='gray', label="Original")

# Calibrated (red)
plt.scatter(phenoage["age"], phenoage["phenoage_calibrated"], color='red', s=40, label="Calibrated")

# Reference line
lims = [min(phenoage["age"].min(), phenoage["phenoage_calibrated"].min()),
        max(phenoage["age"].max(), phenoage["phenoage_calibrated"].max())]
plt.plot(lims, lims, 'k--', lw=1, label='x = y')

plt.xlabel("Chronological age")
plt.ylabel("PhenoAge")
plt.title("PhenoAge Before vs After Calibration")
plt.legend()
plt.xticks(range(0, 131, 10))
plt.yticks(range(0, 131, 10))
plt.xlim(0, 130)
plt.ylim(0, 130)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()

# %%
# --- Step 5: Save calibrated results ---
path_phenoage_calib = f"/BiO/Research/GeroExpressome/Results/PhenoAge4/KGP_{clinical}_phenoage_calibrated.txt"
phenoage.to_csv(path_phenoage_calib, sep=",")

# %%
