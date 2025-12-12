#!/usr/bin/env python3
"""
invert_rao_to_impulse_both_v3.py
--------------------------------
Comparison driver for inverting a complex RAO R(f)=A+iB to the impulse response K(t).

Implements:
  - dst         : Sine-series quadrature after odd-spectrum correction
  - fft         : Jump-removal + IFFT (periodic trapezoid, exponential tail)
  - filon       : Filon on [0,F] of original spectrum + Si truncation tail
  - filon_exp   : Filon on [0,F] of preconditioned spectrum + exponential tail
  - hilbert_fft : Reconstruct B from A via periodic Hilbert, then JR-FFT
  - hilbert_filon: Reconstruct B from A via Hilbert, then Filon+exp
  - real_filon  : **Real-part-only** cosine Filon on Ahat(ω)=A(ω)-K0 α/(α^2+ω^2), then +K0 e^{-αt}
  - real_dct    : **Real-part-only** cosine fast path (DCT-I) on the same Ahat
  - all/compare : run many methods and write consolidated TSV/XLSX and metrics
  - plot        : all calculated data can be plotted upon requests.
Usage: 
python invert_rao_to_impulse_both_v3_3.py --method compare --T 40.3 --N 2048 --alphaT 8 --F 9.5492965855 --Nf 16000 --plot
# F is in cycles/s; here ≈ 60 rad/s / (2π).  Nf is Filon intervals.

Outputs (TAB-separated):
  - impulse_response_all.tsv  : shared time column and each method's K_rec
  - comparison_metrics.tsv    : L∞ and RMS on t>0 for each method
  - per-method *_samples.csv, *_spectrum.csv, *_comparison_*.csv (compat layer)
"""

import argparse
import numpy as np

# ---------- Integration helper (SciPy trapezoid if available) ----------
try:
    from scipy.integrate import trapezoid as _trapz
except Exception:
    def _trapz(y, x):
        return float(np.trapz(y, x))

# ---------- Analytic RAO (half-sinc demo, f in cycles/s) ----------
def A_of_f(ff):
    return 0.5 * (np.abs(ff) <= 0.5).astype(float)

def B_of_f(ff):
    ratio = np.abs((ff + 0.5) / (ff - 0.5 + 0j))
    ratio = np.maximum(ratio, 1e-15)
    return -(1.0/(2.0*np.pi)) * np.log(ratio)

# ---------- Utilities ----------
def estimate_K0_from_tail(fpos, Bpos):
    m = max(10, int(0.1 * fpos.size))
    tail_idx = np.arange(fpos.size - m, fpos.size)
    return float(np.median(-2.0 * np.pi * fpos[tail_idx] * Bpos[tail_idx]))

def estimate_K0_from_A_area(fpos, Apos):
    # In f (cycles/s): K(0+) = (2/π)∫_0^∞ A(ω)dω = (2/π)∫_0^∞ A(2π f) 2π df = 4 ∫_0^∞ A_f(f) df
    return float( 4.0 * _trapz(Apos, x=fpos) )

def save_standard(rao, time, prefix=""):
    f, A, B, S = rao
    t, K_true, K_rec, K_hat = time
    dlm = "\t"
    base = f"impulse_response_comparison_{prefix}.csv" if prefix else "impulse_response_comparison.csv"
    np.savetxt(base.replace(".csv","_samples.csv"), np.column_stack([f, A, B]), delimiter=dlm,
               header="f\tA(f)\tB(f)", comments="")
    np.savetxt(base.replace(".csv","_spectrum.csv"), np.column_stack([f, S]), delimiter=dlm,
               header="f\tS(f) (odd spectrum multiplier)", comments="")
    np.savetxt(base, np.column_stack([t, K_true, K_rec, K_hat]), delimiter=dlm,
               header="t\tK_true\tK_rec\tK_hat", comments="")

def compute_metrics(time_data):
    t, K_true, K_rec, _ = time_data
    mask = t > 0
    tt, err = t[mask], (K_rec - K_true)[mask]
    Linf = float(np.max(np.abs(err)))
    RMS = float(np.sqrt(_trapz(err**2, x=tt) / (tt[-1] - tt[0])))
    return Linf, RMS

# ---------- Sine integral Si(x) (robust) ----------
def Si(x):
    try:
        import mpmath as mp
        return float(mp.si(x))
    except Exception:
        pass
    x = float(x); ax = abs(x)
    if ax < 1e-4:
        import math as _m
        s = 0.0
        for k in range(14):
            s += ((-1.0)**k) * x**(2*k+1) / ((2*k+1) * _m.factorial(2*k+1))
        return s
    if ax > 50.0:
        return np.pi/2 if x >= 0 else -np.pi/2
    n = max(2000, int(40*ax))
    u = np.linspace(0.0, x, n+1)
    y = np.ones_like(u); nz = (u != 0.0)
    y[nz] = np.sin(u[nz]) / u[nz]
    return float(_trapz(y, u))

# ---------- Periodic discrete Hilbert along f ----------
def hilbert_periodic_even_to_odd(A_full):
    """DFT-based periodic Hilbert transform: H{A}(f) = PV(1/π)*A* (1/f) in Fourier sense.
    Implements  𝓕_f{H{A}}(k) = -i sgn(k) 𝓕_f{A}(k).  Handle k=0 and Nyquist as 0.
    """
    M = len(A_full)
    X = np.fft.fft(A_full)
    sgn = np.zeros(M)
    sgn[1:M//2] = 1.0
    sgn[M//2+1:] = -1.0
    Y = -1j * sgn * X
    H = np.fft.ifft(Y)
    return H.real

# ---------- Sine/DST inversion ----------
def invert_dst(T, N, alpha, K0_est):
    n = np.arange(1, N)
    f_pos = n / (2.0 * T)
    w_pos = 2.0 * np.pi * f_pos
    A_pos = A_of_f(f_pos); B_pos = B_of_f(f_pos)
    S_pos = 2.0 * B_pos + (4.0 * np.pi * K0_est * f_pos) / (alpha**2 + w_pos**2)
    df = 1.0/(2.0*T)
    j = np.arange(1, N)
    sin_matrix = np.sin(np.pi * np.outer(j, n) / N)
    K_hat_interior = -2.0 * df * (sin_matrix @ S_pos)
    t_full = np.linspace(0.0, T, N+1)
    K_hat = np.zeros_like(t_full); K_hat[1:-1] = K_hat_interior
    K_rec = K_hat + K0_est * np.exp(-alpha * t_full)
    K_true = np.sinc(t_full)
    return (f_pos, A_pos, B_pos, S_pos), (t_full, K_true, K_rec, K_hat)

# ---------- FFT inversion ----------
def invert_fft(T, N, alpha, K0_est, R_override=None):
    M = 2 * N; L = 2.0 * T
    df_fft = 1.0 / L; dt_fft = L / M
    k = np.arange(-M//2, M//2); f_fft = k * df_fft; w_fft = 2.0 * np.pi * f_fft
    if R_override is None:
        R_fft = A_of_f(f_fft) + 1j * B_of_f(f_fft)
    else:
        R_fft = R_override
    R_tilde = R_fft - K0_est / (alpha + 1j * w_fft)
    F_shifted = np.fft.ifftshift(R_tilde)
    K_hat_periodic = np.fft.ifft(F_shifted) * M * df_fft
    K_hat_periodic = K_hat_periodic.real
    t_full = np.arange(M) * dt_fft
    idx_half = np.arange(0, N+1)
    t_half = t_full[idx_half]; K_hat_half = K_hat_periodic[idx_half]
    K_rec = K_hat_half + K0_est * np.exp(-alpha * t_half)
    K_true = np.sinc(t_half)
    n = np.arange(1, N); f_pos = n / (2.0 * T); w_pos = 2.0 * np.pi * f_pos
    A_pos = A_of_f(f_pos); B_pos = B_of_f(f_pos)
    S_pos = 2.0 * B_pos + (4.0 * np.pi * K0_est * f_pos) / (alpha**2 + w_pos**2)
    return (f_pos, A_pos, B_pos, S_pos), (t_half, K_true, K_rec, K_hat_half)

# ---------- Filon building blocks ----------
def filon_block_uniform(f, g, t):
    # composite quadratic Filon on uniform f-grid for ∫ g(f) e^{i2π f t} df
    Np = len(f); Nint = Np - 1; assert Nint % 2 == 0
    h = f[1] - f[0]
    y0 = g[0:-2:2]; y1 = g[1:-1:2]; y2 = g[2::2]
    a = y1; b = (y2 - y0)/(2*h); c = (y0 - 2*y1 + y2)/(2*h**2)
    fmid = f[1:-1:2]
    mu = 2*np.pi*t
    if abs(mu) < 1e-12:
        I0, I1, I2 = 2*h, 0j, 2*h**3/3
    else:
        th = mu*h
        I0 = 2*np.sin(th)/mu
        I1 = 2j*(-h*np.cos(th)/mu + np.sin(th)/mu**2)
        I2 = 2*(h**2*np.sin(th)/mu + 2*h*np.cos(th)/mu**2 - 2*np.sin(th)/mu**3)
    return np.sum(np.exp(1j*mu*fmid) * (a*I0 + b*I1 + c*I2))

# ---------- Filon (direct) + Si truncation tail ----------
def invert_filon_si(T, N, alpha, K0_est, F=None, Nf=None):
    if F is None: F = (N-1)/(2.0*T)
    if Nf is None: Nf = N-1
    if Nf % 2 == 1: Nf += 1
    f = np.linspace(0.0, F, Nf+1)
    R = A_of_f(f) + 1j*B_of_f(f)
    t_full = np.linspace(0.0, T, N+1); K_true = np.sinc(t_full)
    K_rec = np.zeros_like(t_full); K_hat_dummy = np.zeros_like(t_full)
    for j, tj in enumerate(t_full):
        I = filon_block_uniform(f, R, tj)
        K_tail = (K0_est/np.pi) * (np.pi/2 - Si(2*np.pi*F*tj))
        K_rec[j] = 2*np.real(I) + K_tail
    w_pos = 2.0*np.pi*f
    S_pos = 2.0*B_of_f(f) + (4.0*np.pi*K0_est*f)/(alpha**2 + w_pos**2)
    return (f, A_of_f(f), B_of_f(f), S_pos), (t_full, K_true, K_rec, K_hat_dummy)

# ---------- Filon on preconditioned spectrum + exponential tail ----------
def invert_filon_exp(T, N, alpha, K0_est, F=None, Nf=None, R_override=None):
    if F is None: F = (N-1)/(2.0*T)
    if Nf is None: Nf = N-1
    if Nf % 2 == 1: Nf += 1
    f = np.linspace(0.0, F, Nf+1)
    if R_override is None:
        R = A_of_f(f) + 1j*B_of_f(f)
    else:
        R = R_override
    Rsub = K0_est / (alpha + 1j*2.0*np.pi*f)
    Rc = R - Rsub
    t_full = np.linspace(0.0, T, N+1); K_true = np.sinc(t_full)
    K_rec = np.zeros_like(t_full); K_hat_partial = np.zeros_like(t_full)
    for j, tj in enumerate(t_full):
        I = filon_block_uniform(f, Rc, tj)
        K_hat_partial[j] = 2*np.real(I)
        K_rec[j] = K_hat_partial[j] + K0_est*np.exp(-alpha*tj)
    w_pos = 2.0*np.pi*f
    S_pos = 2.0*B_of_f(f) + (4.0*np.pi*K0_est*f)/(alpha**2 + w_pos**2)
    return (f, R.real, R.imag, S_pos), (t_full, K_true, K_rec, K_hat_partial)

# ---------- NEW: Cosine-only (real part) Filon + exponential add-back ----------
def invert_real_filon(T, N, alpha, K0_est, F=None, Nf=None):
    if F is None: F = (N-1)/(2.0*T)
    if Nf is None: Nf = N-1
    if Nf % 2 == 1: Nf += 1
    f = np.linspace(0.0, F, Nf+1)
    w = 2.0 * np.pi * f
    A = A_of_f(f)
    Ahat = A - (K0_est * alpha) / (alpha**2 + w**2)
    t_full = np.linspace(0.0, T, N+1); K_true = np.sinc(t_full)
    K_hat = np.zeros_like(t_full)
    for j, tj in enumerate(t_full):
        I = filon_block_uniform(f, Ahat, tj)   # real integrand
        K_hat[j] = 4.0 * np.real(I)
    K_rec = K_hat + K0_est * np.exp(-alpha * t_full)
    S_pos = 2.0*B_of_f(f) + (4.0*np.pi*K0_est*f)/(alpha**2 + w**2)
    return (f, A, B_of_f(f), S_pos), (t_full, K_true, K_rec, K_hat)

# ---------- NEW: Cosine-only DCT-I fast path ----------

def invert_real_dct(T, N, alpha, K0_est):
    n = np.arange(1, N)
    f_pos = n/(2.0*T); w_pos = 2.0*np.pi*f_pos
    Apos = A_of_f(f_pos)
    Ahat = Apos - (K0_est * alpha) / (alpha**2 + w_pos**2)
    # Explicit cosine sum on ω-grid: ω_n = nπ/T, t_j = jT/N
    j = np.arange(0, N+1)[:,None]      # shape (N+1,1)
    cosM = np.cos(np.pi * (j * n)[...]/N)  # (N+1, N-1)
    K_hat = (2.0/np.pi) * ((np.pi/T) * (cosM @ Ahat) ).ravel()
    t = np.linspace(0.0, T, N+1)
    K_rec = K_hat + K0_est*np.exp(-alpha*t)
    K_true = np.sinc(t)
    S_pos = 2.0*B_of_f(f_pos) + (4.0*np.pi*K0_est*f_pos)/(alpha**2 + w_pos**2)
    return (f_pos, Apos, B_of_f(f_pos), S_pos), (t, K_true, K_rec, K_hat)

# ---------- Hilbert-based routes ----------
def invert_hilbert_fft(T, N, alpha):
    M = 2*N; L = 2.0*T; df = 1.0/L
    k = np.arange(-M//2, M//2); f_full = k * df
    A_full = A_of_f(f_full)
    H_A = hilbert_periodic_even_to_odd(A_full)
    B_full = - H_A
    R_full = A_full + 1j*B_full
    pos = f_full > 0
    K0_est = float(4.0 * _trapz(A_of_f(f_full[pos]), x=f_full[pos]))  # from A-area
    rao, time = invert_fft(T, N, alpha, K0_est, R_override=R_full)
    return rao, time, K0_est

def invert_hilbert_filon(T, N, alpha, F=None, Nf=None):
    if F is None: F = (N-1)/(2.0*T)
    if Nf is None: Nf = N-1
    if Nf % 2 == 1: Nf += 1
    f = np.linspace(0.0, F, Nf+1)
    # symmetric construction for Hilbert
    df = f[1]-f[0]
    M = 2*Nf+2
    f_sym = np.linspace(-F, F, M)
    A_sym = A_of_f(f_sym)
    H_A = hilbert_periodic_even_to_odd(A_sym)
    B_sym = - H_A
    Bf = B_sym[M//2:]
    R = A_of_f(f) + 1j*Bf
    K0_est = float(4.0 * _trapz(A_of_f(f), x=f))  # from A-area
    rao, time = invert_filon_exp(T, N, alpha, K0_est, F=F, Nf=Nf, R_override=R)
    return rao, time, K0_est

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--method", choices=["dst","fft","filon","filon_exp","hilbert_fft","hilbert_filon","real_filon","real_dct","all","compare"], default="compare")
    ap.add_argument("--T", type=float, default=40.3)
    ap.add_argument("--N", type=int, default=2048)
    ap.add_argument("--alphaT", type=float, default=8.0, help="alpha*T; alpha = alphaT/T")
    ap.add_argument("--F", type=float, default=None, help="Filon cutoff F (cycles/s)")
    ap.add_argument("--Nf", type=int, default=None, help="Filon intervals (even)")
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()

    T, N = args.T, args.N
    alpha = args.alphaT / T

    outputs = {}

    # Base K0 estimate (from A-area, robust for all methods)
    n = np.arange(1, N); f_pos = n/(2.0*T)
    K0_est = estimate_K0_from_A_area(f_pos, A_of_f(f_pos))

    if args.method in ("dst","all","compare"):
        rao_d, time_d = invert_dst(T, N, alpha, K0_est)
        save_standard(rao_d, time_d, prefix="dst" if args.method!="dst" else "")
        outputs["dst"] = time_d
    if args.method in ("fft","all","compare"):
        rao_f, time_f = invert_fft(T, N, alpha, K0_est)
        save_standard(rao_f, time_f, prefix="fft" if args.method!="fft" else "")
        outputs["fft"] = time_f
    if args.method in ("filon","all","compare"):
        rao_s, time_s = invert_filon_si(T, N, alpha, K0_est, F=args.F, Nf=args.Nf)
        save_standard(rao_s, time_s, prefix="filon" if args.method!="filon" else "")
        outputs["filon"] = time_s
    if args.method in ("filon_exp","all","compare"):
        rao_e, time_e = invert_filon_exp(T, N, alpha, K0_est, F=args.F, Nf=args.Nf)
        save_standard(rao_e, time_e, prefix="filon_exp" if args.method!="filon_exp" else "")
        outputs["filon_exp"] = time_e
    if args.method in ("real_filon","all","compare"):
        rao_rf, time_rf = invert_real_filon(T, N, alpha, K0_est, F=args.F, Nf=args.Nf)
        save_standard(rao_rf, time_rf, prefix="real_filon" if args.method!="real_filon" else "")
        outputs["real_filon"] = time_rf
    if args.method in ("real_dct","all","compare"):
        rao_rd, time_rd = invert_real_dct(T, N, alpha, K0_est)
        save_standard(rao_rd, time_rd, prefix="real_dct" if args.method!="real_dct" else "")
        outputs["real_dct"] = time_rd
    if args.method in ("hilbert_fft","all","compare"):
        rao_hf, time_hf, K0h = invert_hilbert_fft(T, N, alpha)
        save_standard(rao_hf, time_hf, prefix="hilbert_fft" if args.method!="hilbert_fft" else "")
        outputs["hilbert_fft"] = time_hf
    if args.method in ("hilbert_filon","all","compare"):
        rao_hfil, time_hfil, K0hf = invert_hilbert_filon(T, N, alpha, F=args.F, Nf=args.Nf)
        save_standard(rao_hfil, time_hfil, prefix="hilbert_filon" if args.method!="hilbert_filon" else "")
        outputs["hilbert_filon"] = time_hfil

    # Consolidated TSV/XLSX
    tag0 = next(iter(outputs.keys()))
    tvec = outputs[tag0][0]; K_true = outputs[tag0][1]
    cols = {"t": tvec, "K_true": K_true}
    mapping = {
        "dst":"K_dst","fft":"K_fft","filon":"K_filon","filon_exp":"K_filon_exp",
        "real_filon":"K_real_filon","real_dct":"K_real_dct",
        "hilbert_fft":"K_hilbert_fft","hilbert_filon":"K_hilbert_filon"
    }
    for k,v in mapping.items():
        if k in outputs: cols[v] = outputs[k][2]
    # Write consolidated TSV
    arr = np.column_stack([cols[k] for k in cols.keys()])
    header = "\t".join(cols.keys())
    np.savetxt("impulse_response_all.tsv", arr, delimiter="\t", header=header, comments="")

    # Also write XLSX if possible
    try:
        import pandas as pd
        with pd.ExcelWriter("impulse_response_all.xlsx", engine="openpyxl") as writer:
            pd.DataFrame(cols).to_excel(writer, index=False, sheet_name="K_all")
    except Exception:
        pass

    # Metrics
    lines = ["method\tLinf\tRMS"]
    for k in mapping.keys():
        if k in outputs:
            Linf, RMS = compute_metrics(outputs[k])
            lines.append(f"{k}\t{Linf:.6e}\t{RMS:.6e}")
    with open("comparison_metrics.tsv","w",encoding="utf-8") as fh:
        fh.write("\n".join(lines))

    if args.plot:
        try:
            import matplotlib.pyplot as plt
            plt.figure()
            plt.plot(cols["t"], cols["K_true"], label="K_true")
            for key in list(mapping.values()):
                if key in cols:
                    plt.plot(cols["t"], cols[key], label=key)
            plt.legend(); plt.grid(True); plt.xlabel("t"); plt.ylabel("K(t)")
            plt.title("Impulse response (all methods)"); plt.show()
        except Exception:
            pass

if __name__ == "__main__":
    main()
