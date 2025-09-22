import numpy as np
import matplotlib.pyplot as plt

def cos(x): return np.cos(np.array(x, dtype=np.float64))
def sin(x): return np.sin(np.array(x, dtype=np.float64))

def Df_forward(f, x, h):
    x, h = float(x), float(h)
    return (f(x+h)-f(x))/h

def Df_central(f, x, h):
    x, h = float(x), float(h)
    return (f(x+h)-f(x-h))/(2*h)

def Df_extrapolate(f, x, h):
    D_h  = Df_central(f, x, h)
    D_h2 = Df_central(f, x, h/2)
    return (4*D_h2 - D_h)/3

def rel_err(app, real):
    return abs((app-real)/real)

def plot_derivative_errors(x, hs, outfile):
    trueval = -sin(x)
    err_fwd, err_cen, err_ext = [], [], []
    for h in hs:
        err_fwd.append(rel_err(Df_forward(cos, x, h), trueval))
        err_cen.append(rel_err(Df_central(cos, x, h), trueval))
        err_ext.append(rel_err(Df_extrapolate(cos, x, h), trueval))
    plt.figure(figsize=(6,4))
    plt.loglog(hs, err_fwd, '-', label="forward")
    plt.loglog(hs, err_cen, '-', label="central")
    plt.loglog(hs, err_ext, '-', label="extrapolate")
    plt.xlabel("h"); plt.ylabel(f"relative error at x={x}")
    plt.title("Derivative error vs h")
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend(); plt.tight_layout(); plt.savefig(outfile)
    plt.close()

def minusexp(x):
    x = np.array(x, dtype=np.float64)
    return np.exp(-x)

def Int_middle(f,a,b,N):
    a,b = float(a), float(b); N = int(N)
    h = (b-a)/N
    k = np.arange(N)
    mids = a + (k + 0.5) * h
    return h * np.sum(f(mids))

def Int_trapezoid(f,a,b,N):
    a,b = float(a), float(b); N = int(N)
    h = (b-a)/N
    xs = a + h*np.arange(N+1)
    fx = f(xs)
    S = 0.5*fx[0] + np.sum(fx[1:-1]) + 0.5*fx[-1]
    return h*S

def Int_simpson(f,a,b,N):
    a,b = float(a), float(b); N = int(N)
    if N % 2 == 1: N += 1
    h  = (b-a)/N
    xs = a + h*np.arange(N+1)
    fx = f(xs)
    S = fx[0] + fx[-1] + 4*np.sum(fx[1:-1:2]) + 2*np.sum(fx[2:-1:2])
    return h*(S/3.0)

def plot_integration_errors(outfile):
    Ns = np.unique((np.logspace(1, 6.5, 40)).astype(int))
    trueval = 1.0 - np.exp(-1.0)
    def rel_err(app, real): return abs((float(app) - float(real))/float(real))
    err_mid = [rel_err(Int_middle(minusexp, 0, 1, N), trueval) for N in Ns]
    err_trap= [rel_err(Int_trapezoid(minusexp, 0, 1, N), trueval) for N in Ns]
    err_simp= [rel_err(Int_simpson(minusexp, 0, 1, N), trueval) for N in Ns]
    plt.figure(figsize=(6,4))
    plt.loglog(Ns, err_mid, label="midpoint")
    plt.loglog(Ns, err_trap, label="trapezoid")
    plt.loglog(Ns, err_simp, label="Simpson")
    plt.xlabel("N"); plt.ylabel("relative error")
    plt.title("Integration error vs N")
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend(); plt.tight_layout(); plt.savefig(outfile)
    plt.close()

if __name__ == "__main__":
    hs = np.logspace(-10, -2, 200)
    outdir = "../hw1-writeup/figs"
    import os
    os.makedirs(outdir, exist_ok=True)
    plot_derivative_errors(0.1, hs, os.path.join(outdir, "derivative_errors.pdf"))
    plot_integration_errors(os.path.join(outdir, "integration_errors.pdf"))
    print("Saved derivative and integration figures to", outdir)
