

import argparse
import math
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
git add .
git commit -m "implement sprime metric"
# ---- 4PL model ----
# y = D + (A - D) / (1 + (x/C)^B)
# A: top, D: bottom, C: EC50, B: Hill slope
def four_pl(x, A, B, C, D):
    x = np.asarray(x, dtype=float)
    return D + (A - D) / (1.0 + np.power(x / C, B))

def fit_4pl(conc_uM, y):
    """
    Fit 4PL and return params (A,B,C,D) and some diagnostics.
    """
    x = np.asarray(conc_uM, dtype=float)
    y = np.asarray(y, dtype=float)

    # remove non-positive conc (log/ratio issues) and NaNs
    m = np.isfinite(x) & np.isfinite(y) & (x > 0)
    x, y = x[m], y[m]

    if len(x) < 4:
        raise ValueError(f"Not enough points to fit 4PL (need >=4). Got {len(x)}.")

    # Initial guesses:
    A0 = np.nanmax(y)
    D0 = np.nanmin(y)
    # EC50 guess: median concentration
    C0 = np.median(x)
    # Hill slope guess
    B0 = 1.0

    # bounds: keep EC50 positive; slope reasonable
    # A and D can be any, but constrain a bit to avoid crazy fits
    lower = [min(D0, A0) - 10.0, 0.01, np.min(x) * 1e-6, min(D0, A0) - 10.0]
    upper = [max(D0, A0) + 10.0, 10.0, np.max(x) * 1e6, max(D0, A0) + 10.0]

    popt, pcov = curve_fit(
        four_pl, x, y,
        p0=[A0, B0, C0, D0],
        bounds=(lower, upper),
        maxfev=20000,
    )
    A, B, C, D = popt
    return {"A": A, "B": B, "EC50": C, "D": D, "pcov": pcov, "n": len(x)}

def compute_sprime(Emax, EC50_uM, x1_uM=1.0):
    """
    S' = asinh((Emax/EC50) * x1_uM)
    EC50 in uM.
    """
    if EC50_uM is None or not np.isfinite(EC50_uM) or EC50_uM <= 0:
        return np.nan
    if Emax is None or not np.isfinite(Emax):
        return np.nan
    return np.arcsinh((Emax / EC50_uM) * x1_uM)

def choose_y_column(df):
    """
    Prefer test/dmso; fallback to norm_TripMean_EDU index; else error.
    """
    candidates = ["test/dmso", "norm_TripMean_EDU index", "0 (EDU signal)"]
    for c in candidates:
        if c in df.columns:
            # pick the first one that has some non-null values
            if df[c].notna().sum() > 0:
                return c
    raise ValueError(f"None of the expected y columns found or they are empty: {candidates}")

def load_subset(xlsx, sheet, compound, celltype, timepoint):
    df = pd.read_excel(xlsx, sheet_name=sheet)

    # Basic cleanup: strip spaces in column names
    df.columns = [str(c).strip() for c in df.columns]

    required = ["compound", "cell type", "time point", "final compound concentration uM"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # filter
    sub = df[
        (df["compound"].astype(str) == str(compound)) &
        (df["cell type"].astype(str) == str(celltype)) &
        (df["time point"].astype(str) == str(timepoint))
    ].copy()

    if sub.empty:
        raise ValueError(f"No rows found for compound={compound}, cell type={celltype}, time={timepoint} in sheet={sheet}")

    ycol = choose_y_column(sub)

    # Ensure numeric
    sub["final compound concentration uM"] = pd.to_numeric(sub["final compound concentration uM"], errors="coerce")
    sub[ycol] = pd.to_numeric(sub[ycol], errors="coerce")

    sub = sub.dropna(subset=["final compound concentration uM", ycol])
    return sub, ycol

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xlsx", required=True)
    ap.add_argument("--sheet", required=True)
    ap.add_argument("--compound", required=True)
    ap.add_argument("--celltype", required=True)
    ap.add_argument("--time", required=True)
    ap.add_argument("--reference_celltype", default=None, help="If provided, also compute ΔS' = S'ref - S'test")
    args = ap.parse_args()

    # ---- IC50/EC50 + S' for test condition ----
    sub, ycol = load_subset(args.xlsx, args.sheet, args.compound, args.celltype, args.time)
    conc = sub["final compound concentration uM"].values
    y = sub[ycol].values

    fit = fit_4pl(conc, y)
    Emax = fit["A"] - fit["D"]   # matches your Excel label "A-D (Emax)"
    EC50 = fit["EC50"]

    sprime = compute_sprime(Emax, EC50, x1_uM=1.0)

    print("=== TEST ===")
    print(f"sheet={args.sheet}")
    print(f"compound={args.compound}  celltype={args.celltype}  time={args.time}")
    print(f"y column used: {ycol}")
    print(f"4PL params: A={fit['A']:.6g}, B={fit['B']:.6g}, EC50={EC50:.6g} uM, D={fit['D']:.6g}")
    print(f"Emax (A-D) = {Emax:.6g}")
    print(f"S' = asinh((Emax/EC50)*1uM) = {sprime:.6g}")

    # ---- ΔS' if reference provided ----
    if args.reference_celltype is not None:
        ref_sub, ref_ycol = load_subset(args.xlsx, args.sheet, args.compound, args.reference_celltype, args.time)
        ref_fit = fit_4pl(ref_sub["final compound concentration uM"].values, ref_sub[ref_ycol].values)
        ref_Emax = ref_fit["A"] - ref_fit["D"]
        ref_EC50 = ref_fit["EC50"]
        ref_sprime = compute_sprime(ref_Emax, ref_EC50, x1_uM=1.0)

        delta = ref_sprime - sprime

        print("\n=== REFERENCE ===")
        print(f"compound={args.compound}  celltype={args.reference_celltype}  time={args.time}")
        print(f"y column used: {ref_ycol}")
        print(f"4PL params: A={ref_fit['A']:.6g}, B={ref_fit['B']:.6g}, EC50={ref_EC50:.6g} uM, D={ref_fit['D']:.6g}")
        print(f"Emax (A-D) = {ref_Emax:.6g}")
        print(f"S' = {ref_sprime:.6g}")

        print("\n=== DELTA ===")
        print(f"ΔS' = S'reference - S'test = {delta:.6g}")

if __name__ == "__main__":
    main()
