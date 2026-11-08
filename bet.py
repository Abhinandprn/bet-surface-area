"""BET specific surface area from an N2 adsorption isotherm (77 K).

Linear BET equation:

    (p/p0) / (n (1 - p/p0)) = 1/(n_m C) + (C - 1)/(n_m C) * (p/p0)

Plot y = x / (n (1 - x)) against x = p/p0 in the BET range:
    slope s and intercept i give
    n_m = 1 / (s + i)      monolayer capacity (cm3 STP/g)
    C   = s / i + 1        BET constant
    S_BET = n_m * N_A * sigma / V_molar   with sigma(N2) = 0.162 nm2

The default fit range is p/p0 = 0.05 to 0.30. With --rouquerol the script
also checks the Rouquerol criteria (n(1 - p/p0) must increase with p/p0, and C > 0).

Usage:
    python bet.py examples/example_isotherm.csv
    python bet.py examples/example_isotherm.csv --range 0.05 0.25 --plot examples/bet_plot.png
"""

import argparse
import sys

import numpy as np
import pandas as pd

N_A = 6.02214076e23          # 1/mol
SIGMA_N2 = 0.162e-18         # m2 per N2 molecule
V_MOLAR = 22414.0            # cm3 STP per mol


def bet_fit(p_rel, n_ads, lo=0.05, hi=0.30):
    """Fit the linear BET equation between lo and hi (relative pressure).
    Returns a dict with S_BET (m2/g), n_m (cm3 STP/g), C, R2 and the fitted points."""
    p_rel = np.asarray(p_rel, float)
    n_ads = np.asarray(n_ads, float)
    mask = (p_rel >= lo) & (p_rel <= hi)
    if mask.sum() < 3:
        raise ValueError(f"Need at least 3 points between p/p0 = {lo} and {hi}.")

    x = p_rel[mask]
    y = x / (n_ads[mask] * (1 - x))
    slope, intercept = np.polyfit(x, y, 1)

    n_m = 1 / (slope + intercept)
    c = slope / intercept + 1
    s_bet = n_m / V_MOLAR * N_A * SIGMA_N2

    y_fit = slope * x + intercept
    r2 = 1 - np.sum((y - y_fit) ** 2) / np.sum((y - y.mean()) ** 2)

    return {"S_BET_m2_g": s_bet, "n_m_cm3_g": n_m, "C": c, "R2": r2,
            "points": int(mask.sum()), "x": x, "y": y, "slope": slope, "intercept": intercept}


def rouquerol_ok(p_rel, n_ads, lo, hi, c):
    """Simple Rouquerol checks for the chosen range."""
    p_rel = np.asarray(p_rel, float)
    n_ads = np.asarray(n_ads, float)
    mask = (p_rel >= lo) & (p_rel <= hi)
    term = n_ads[mask] * (1 - p_rel[mask])
    increasing = bool(np.all(np.diff(term) > 0))
    return {"n(1-p/p0) increases": increasing, "C > 0": bool(c > 0)}


def main():
    parser = argparse.ArgumentParser(description="BET surface area from an N2 isotherm")
    parser.add_argument("csv", help="CSV with columns: p_rel, n_ads_cm3_g (adsorption branch)")
    parser.add_argument("--range", nargs=2, type=float, default=[0.05, 0.30], metavar=("LOW", "HIGH"))
    parser.add_argument("--rouquerol", action="store_true", help="check the Rouquerol criteria")
    parser.add_argument("--plot", help="save isotherm and BET plot to this PNG file")
    args = parser.parse_args()

    data = pd.read_csv(args.csv)
    if not {"p_rel", "n_ads_cm3_g"}.issubset(data.columns):
        sys.exit("CSV needs the columns p_rel and n_ads_cm3_g")
    data = data.sort_values("p_rel")

    lo, hi = args.range
    r = bet_fit(data["p_rel"], data["n_ads_cm3_g"], lo, hi)

    print(f"\nBET fit range p/p0 = {lo} to {hi} ({r['points']} points)")
    print(f"S_BET            = {r['S_BET_m2_g']:.1f} m2/g")
    print(f"Monolayer n_m    = {r['n_m_cm3_g']:.2f} cm3 STP/g")
    print(f"BET constant C   = {r['C']:.1f}")
    print(f"R2               = {r['R2']:.5f}")

    if args.rouquerol:
        for k, v in rouquerol_ok(data["p_rel"], data["n_ads_cm3_g"], lo, hi, r["C"]).items():
            print(f"Rouquerol: {k:22s} {'yes' if v else 'NO'}")

    if args.plot:
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
        ax1.plot(data["p_rel"], data["n_ads_cm3_g"], "o-", ms=4)
        ax1.axvspan(lo, hi, color="orange", alpha=0.2, label="BET range")
        ax1.set_xlabel("Relative pressure p/p$_0$")
        ax1.set_ylabel("Adsorbed amount (cm$^3$ STP/g)")
        ax1.set_title("N$_2$ isotherm (77 K)")
        ax1.legend()

        xs = np.linspace(lo, hi, 50)
        ax2.plot(r["x"], r["y"], "o", label="data")
        ax2.plot(xs, r["slope"] * xs + r["intercept"], "-", label=f"fit, R$^2$ = {r['R2']:.4f}")
        ax2.set_xlabel("p/p$_0$")
        ax2.set_ylabel(r"$\frac{p/p_0}{n\,(1-p/p_0)}$  (g/cm$^3$)")
        ax2.set_title(f"BET plot: S$_{{BET}}$ = {r['S_BET_m2_g']:.0f} m$^2$/g, C = {r['C']:.0f}")
        ax2.legend()
        fig.tight_layout()
        fig.savefig(args.plot, dpi=150)
        print(f"\nPlot saved to {args.plot}")


if __name__ == "__main__":
    main()
