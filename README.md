# Computational Physics — HW1 (Starter)

This repository contains Python scripts for:
- Numerical differentiation tests on `cos(x)` with forward, central, and extrapolated differences.
- Quadrature tests (midpoint, trapezoid, Simpson) for `∫_0^1 e^{-x} dx`.
- Correlation function `xi(r)` computed from a tabulated matter power spectrum `P(k)` and BAO peak detection.

## Environment
Python 3.10+ recommended. Install dependencies:
```
pip install numpy matplotlib scipy
```

## Usage
- Derivatives & quadrature (produces PDF plots in `../hw1-writeup/figs`):
```
python diff_and_quad.py
```

- BAO correlation function (expects `lcdm_z0.matter_pk` next to the script; produces a PDF plot):
```
python bao_xi.py --pk lcdm_z0.matter_pk
```

## Suggested Git Workflow (GUI alternative below)
```
git init
git add .
git commit -m "HW1: initial commit"
git branch -M main
git remote add origin https://github.com/<your-user>/comp-phys-hw1.git
git push -u origin main
```
