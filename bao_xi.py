import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as Sp
from scipy.integrate import quad
import argparse, os

def xi_from_pk(k_vals, P_vals, r, limit=200):
    spline = Sp(k_vals, P_vals, k=3)
    integrand = lambda k: k * spline(k)
    result, _ = quad(integrand, k_vals.min(), k_vals.max(),
                     weight='sin', wvar=r, limit=limit)
    return (result/(2*np.pi**2*r))

def main(pk_path, outpdf):
    data = np.loadtxt(pk_path)
    k_vals, P_vals = data[:,0], data[:,1]
    r_vals = np.linspace(50, 120, 200)
    xi_vals = np.array([xi_from_pk(k_vals, P_vals, r, limit=400) for r in r_vals])
    r2xi = r_vals**2 * xi_vals
    mask = (r_vals >= 80) & (r_vals <= 120)
    idx_local = np.argmax(r2xi[mask])
    idx_global = np.where(mask)[0][idx_local]
    peak_r = r_vals[idx_global]

    plt.figure(figsize=(6,4))
    plt.plot(r_vals, r2xi, label=r"$r^2 \xi(r)$")
    plt.axvline(peak_r, ls="--", label=f"BAO ~ {peak_r:.2f} Mpc/h")
    plt.xlabel(r"$r\ [\mathrm{Mpc}/h]$"); plt.ylabel(r"$r^2 \xi(r)$")
    plt.title("BAO feature from power spectrum")
    plt.grid(True, ls="--", alpha=0.6)
    plt.legend(); plt.tight_layout()
    outdir = os.path.dirname(outpdf)
    os.makedirs(outdir, exist_ok=True)
    plt.savefig(outpdf); plt.close()
    print(f"Saved BAO figure to {outpdf}; peak ≈ {peak_r:.2f} Mpc/h")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--pk", required=True, help="path to lcdm_z0.matter_pk")
    ap.add_argument("--out", default="../hw1-writeup/figs/bao_peak.pdf")
    args = ap.parse_args()
    main(args.pk, args.out)
