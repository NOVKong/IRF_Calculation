#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
invert_irf_wamit_v5.py
Dimensional radiation IRF inversion from WAMIT-style AB data (Ā=A/ρ, B̄=B/(ρ ω)).
Conventions (SIMO/WAMIT):
  a(ω) = A(ω) - A(∞)         [kg]        even
  c(ω) = C(ω) = B(ω)         [N·s/m]     even
  H(ω) = c(ω) + i ω a(ω)     (Fourier transform of the IRF h(t), τ≥0)
  h(t) = (2/π)∫₀^∞ c(ω)cos(ω t)dω = -(2/π)∫₀^∞ ω a(ω)sin(ω t)dω
We build a symmetric frequency function R(f) with f=ω/(2π) so that
  h(t) = 2 Re ∫_{-∞}^{∞} R(f) e^{i2π f t} df
which is satisfied by
  A(f) = Re R(f) =  c(|ω|),            ω=2π f
  B(f) = Im R(f) =  sgn(f) π|f| a(|ω|).
We read nondimensional AB and convert:
  A(ω) = ρ Ā(ω),  B(ω) = ρ ω B̄(ω).  Then a(ω)=A- A∞, c(ω)=B.
Supports: filon, fft, fft_jr (jump removal), filon_jr, real_filon, dst, compare.

Usage:

python invert_irf_wamit_v5.py --excel E:\WORK\2368SSY_GSM\Wamit\T9p23\impulse_response_all_v3_0.xlsx --sheet AB --rho 1025 --T 40 --N 2048 --method compare --out irf_dim_compare_v5.tsv
"""

import numpy as np
import pandas as pd
import argparse
from typing import Tuple
from scipy.integrate import trapezoid as trapz

# ------------------------ Excel loader and dimensional conversion ------------------------

class ABData:
    def __init__(self, path: str, sheet="AB",
                 freq_col="Frequency[Hz]",
                 A_col="A_ndim=A/ϱ",
                 B_col="B_ndim=B/(ϱω)",
                 rho=1025.0, use_zero=True):
        self.path=path; self.sheet=sheet; self.freq_col=freq_col; self.A_col=A_col; self.B_col=B_col
        self.rho=rho; self.use_zero=use_zero
        self._load()

    def _find_header_row(self, df0):
        for r in range(min(10, len(df0))):
            row = "".join(str(x) for x in list(df0.iloc[r].values))
            if "freq" in row.lower():
                return r
        return 0

    def _load(self):
        xls = pd.ExcelFile(self.path)
        df0 = xls.parse(self.sheet, header=None)
        hdr = self._find_header_row(df0)
        df = xls.parse(self.sheet, header=hdr)

        f_hz = pd.to_numeric(df[self.freq_col], errors="coerce").to_numpy(float)
        Abar = pd.to_numeric(df[self.A_col],   errors="coerce").to_numpy(float)  # Ā=A/ρ
        Bbar = pd.to_numeric(df[self.B_col],   errors="coerce").to_numpy(float)  # B̄=B/(ρ ω)

        m = np.isfinite(f_hz) & np.isfinite(Abar) & np.isfinite(Bbar)
        f_hz, Abar, Bbar = f_hz[m], Abar[m], Bbar[m]
        idx = np.argsort(f_hz); f_hz, Abar, Bbar = f_hz[idx], Abar[idx], Bbar[idx]
        if not self.use_zero and f_hz.size>0 and abs(f_hz[0])<1e-12:
            f_hz, Abar, Bbar = f_hz[1:], Abar[1:], Bbar[1:]

        self.f_pos = f_hz
        self.w_pos = 2*np.pi*f_hz
        # dimensional
        self.A_dim = self.rho * Abar                             # [kg]
        self.B_dim = self.rho * (2*np.pi*f_hz) * Bbar            # [N·s/m]
        # A∞
        mA = max(4, int(0.1*self.A_dim.size))
        self.A_inf = float(np.mean(self.A_dim[-mA:]))
        self.a_dim = self.A_dim - self.A_inf
        self.c_dim = self.B_dim

    # even/odd interpolants (piecewise linear) with one-sided saturation
    def interp_even(self, w_query, y_pos, right):
        wq = np.abs(w_query)
        return np.interp(wq, self.w_pos, y_pos, left=y_pos[0], right=right)

    # Accessors for c(ω) and a(ω)
    def c_of_w(self, w):
        return self.interp_even(w, self.c_dim, right=0.0)
    def a_of_w(self, w):
        return self.interp_even(w, self.a_dim, right=0.0)

# ------------------------ Build R(f) for the IRF convention ------------------------

def R_of_f(f, ab: ABData):
    f = np.asarray(f, float)
    w = 2*np.pi*f
    c = ab.c_of_w(w)                         # even
    a = ab.a_of_w(w)                         # even in |ω|
    A = c                                    # Re
    B = np.sign(f) * np.pi * np.abs(f) * a   # Im
    return A + 1j*B

# ------------------------ Filon quadrature on a uniform f-grid ------------------------

def filon_block_uniform_f(f: np.ndarray, g: np.ndarray, t: float) -> complex:
    Np = len(f); Niv = Np - 1
    assert Niv % 2 == 0, "Filon needs an even number of intervals"
    h = f[1] - f[0]
    y0 = g[0:-2:2]; y1 = g[1:-1:2]; y2 = g[2::2]
    a = y1
    b = (y2 - y0)/(2*h)
    c = (y0 - 2*y1 + y2)/(2*h**2)
    f_mid = f[1:-1:2]
    mu = 2*np.pi*t
    if abs(mu) < 1e-14:
        I0, I1, I2 = 2*h, 0.0j, 2*h**3/3
    else:
        th = mu*h
        I0 = 2*np.sin(th)/mu
        I1 = 2j*(-h*np.cos(th)/mu + np.sin(th)/mu**2)
        I2 = 2*(h**2*np.sin(th)/mu + 2*h*np.cos(th)/mu**2 - 2*np.sin(th)/mu**3)
    return np.sum(np.exp(1j*mu*f_mid) * (a*I0 + b*I1 + c*I2))

# ------------------------ High-frequency tail fit to estimate K0=h(0) and optional tails ------------------------

def estimate_K0_and_tails(ab: ABData, frac=0.2):
    """Estimate h(0)= (2/π) ∫_0^∞ c(ω) dω with a 1/ω^2 tail model c≈C/ω^2.
       Returns (K0, C) where tail ∫_Ω^∞ c ≈ C/Ω.
    """
    w = ab.w_pos; c = ab.c_dim
    m0 = int((1.0-frac)*len(w))
    w_tail = w[m0:]; c_tail = c[m0:]
    # Fit c * ω^2 ≈ C (median for robustness)
    C = float(np.median(c_tail * (np.maximum(w_tail,1e-9)**2)))
    # Integral over measured band (trapezoid) in ω
    I_meas = trapz(c, w)
    # Tail integral from Ω to ∞ = ∫ C/ω^2 dω = C/Ω
    Omega = float(w[-1])
    I_tail = C / Omega if Omega>0 else 0.0
    K0 = (2.0/np.pi) * (I_meas + I_tail)
    return K0, C

# ------------------------ Inversion methods ------------------------

def choose_uniform_f_grid(ab: ABData, T: float, N: int):
    """Choose a uniform frequency grid f∈[0,F] compatible with the FFT/DST window:
       F is limited by available data: F ≤ max(f_pos).  We use an even number of intervals for Filon.
    """
    Fmax = float(np.max(ab.f_pos))
    # If the user's N,T imply a higher Nyquist than Fmax, reduce N accordingly
    Nyquist = (N/2)/T
    if Nyquist > Fmax:
        N = int(2*T*Fmax)
        N = max(64, N - (N % 2))  # even
    # For Filon we create a separate high-resolution uniform grid
    Niv = 20000
    if Niv % 2 == 1: Niv += 1
    f_filon = np.linspace(0.0, Fmax, Niv+1)
    return f_filon, Fmax, N

def invert_filon(ab: ABData, T: float, N: int):
    f, Fmax, N_ok = choose_uniform_f_grid(ab, T, N)
    R = R_of_f(f, ab)
    t = np.linspace(0.0, T, N_ok+1)
    h = np.zeros_like(t)
    for j, tj in enumerate(t):
        I = filon_block_uniform_f(f, R, tj)
        h[j] = 4.0*np.real(I)   # 2·Re over ±f → 4·Re ∫_0^F
    return t, h

def invert_real_filon(ab: ABData, T: float, N: int):
    f, Fmax, N_ok = choose_uniform_f_grid(ab, T, N)
    A = ab.c_of_w(2*np.pi*f)  # c(ω)
    t = np.linspace(0.0, T, N_ok+1)
    h = np.zeros_like(t)
    for j, tj in enumerate(t):
        I = filon_block_uniform_f(f, A, tj)
        h[j] = 4.0*np.real(I)
    return t, h

def invert_fft(ab: ABData, T: float, N: int):
    # Ensure Nyquist ≤ max data frequency
    Fmax = float(np.max(ab.f_pos))
    Nyq = (N/2)/T
    if Nyq > Fmax:
        N = int(2*T*Fmax)
        N = max(64, N - (N % 2))
    f = np.fft.fftfreq(N, d=T/N)  # k/T
    f = np.fft.fftshift(f)
    df = 1.0/T
    R = R_of_f(f, ab)
    Rshift = np.fft.ifftshift(R)
    hper = np.fft.ifft(Rshift) * N * df
    t = np.linspace(0.0, T, N, endpoint=False)
    h = np.real(hper)
    return t, h

def invert_fft_jr(ab: ABData, T: float, N: int, alphaT=8.0):
    # Jump removal in the *H(ω)* spectrum: subtract K0/(α+iω) and use FFT via R=(1/π)H.
    # Here we implement it equivalently in f (cycles/s): subtract K0/(α + i 2π f).
    K0, Ctail = estimate_K0_and_tails(ab, frac=0.2)  # robust
    alpha = alphaT / T
    Fmax = float(np.max(ab.f_pos))
    Nyq = (N/2)/T
    if Nyq > Fmax:
        N = int(2*T*Fmax); N = max(64, N - (N % 2))
    f = np.fft.fftfreq(N, d=T/N); f = np.fft.fftshift(f)
    df = 1.0/T
    # Build H(ω) on ±f: c(|ω|) + i ω a(|ω|), with ω=2π f (signed)
    w = 2*np.pi*f
    c = ab.c_of_w(w)
    a = ab.a_of_w(w)
    H = c + 1j*w*a
    H_hat = H - K0/(alpha + 1j*w)  # subtract simple pole
    R = (1.0/np.pi) * H_hat        # map to R
    Rshift = np.fft.ifftshift(R)
    h_hat = np.fft.ifft(Rshift) * N * df
    t = np.linspace(0.0, T, N, endpoint=False)
    h = np.real(h_hat) + K0*np.exp(-alpha*t)  # add back exponential
    return t, h, K0, alpha

def invert_filon_jr(ab: ABData, T: float, N: int, alphaT=8.0):
    K0, Ctail = estimate_K0_and_tails(ab, frac=0.2)
    alpha = alphaT / T
    f, Fmax, N_ok = choose_uniform_f_grid(ab, T, N)
    w = 2*np.pi*f
    c = ab.c_of_w(w); a = ab.a_of_w(w)
    H = c + 1j*w*a
    H_hat = H - K0/(alpha + 1j*w)
    R_hat = (1.0/np.pi) * H_hat
    t = np.linspace(0.0, T, N_ok+1)
    h_hat = np.zeros_like(t)
    for j, tj in enumerate(t):
        I = filon_block_uniform_f(f, R_hat, tj)
        h_hat[j] = 4.0*np.real(I)
    h = h_hat + K0*np.exp(-alpha*t)
    return t, h, K0, alpha

# ------------------------ CLI ------------------------

def save_tsv(path, cols, header):
    arr = np.column_stack(cols)
    np.savetxt(path, arr, delimiter="\t", header=header, comments="")

def main():
    ap = argparse.ArgumentParser(description="Dimensional IRF inversion (SIMO/WAMIT conventions) from AB data.")
    ap.add_argument("--excel", required=True, help="Path to Excel containing AB (Ā=A/ρ, B̄=B/(ρ ω))")
    ap.add_argument("--sheet", default="AB")
    ap.add_argument("--freqcol", default="Frequency[Hz]")
    ap.add_argument("--acol", default="A_ndim=A/ϱ")
    ap.add_argument("--bcol", default="B_ndim=B/(ϱω)")
    ap.add_argument("--rho", type=float, default=1025.0)
    ap.add_argument("--use_zero", type=int, default=1)
    ap.add_argument("--method", choices=["filon","real_filon","fft","fft_jr","filon_jr","dst","compare"], default="compare")
    ap.add_argument("--T", type=float, default=40.0)
    ap.add_argument("--N", type=int, default=4096)
    ap.add_argument("--alphaT", type=float, default=8.0)
    ap.add_argument("--out", default="/mnt/data/irf_dimensional.tsv")
    args = ap.parse_args()

    ab = ABData(args.excel, sheet=args.sheet, freq_col=args.freqcol, A_col=args.acol, B_col=args.bcol,
                rho=args.rho, use_zero=bool(args.use_zero))

    if args.method == "filon":
        t,h = invert_filon(ab, args.T, args.N)
        save_tsv(args.out, (t,h), "t[s]\th_filon_dim")
    elif args.method == "real_filon":
        t,h = invert_real_filon(ab, args.T, args.N)
        save_tsv(args.out, (t,h), "t[s]\th_real_filon_dim")
    elif args.method == "fft":
        t,h = invert_fft(ab, args.T, args.N)
        save_tsv(args.out, (t,h), "t[s]\th_fft_dim")
    elif args.method == "fft_jr":
        t,h,K0,alpha = invert_fft_jr(ab, args.T, args.N, args.alphaT)
        save_tsv(args.out, (t,h), "t[s]\th_fftJR_dim")
        print(f"K0_est={K0:.6g}, alpha={alpha:.6g}")
    elif args.method == "filon_jr":
        t,h,K0,alpha = invert_filon_jr(ab, args.T, args.N, args.alphaT)
        save_tsv(args.out, (t,h), "t[s]\th_filonJR_dim")
        print(f"K0_est={K0:.6g}, alpha={alpha:.6g}")
    elif args.method == "dst":
        # Use odd spectrum B(f)=π|f| a(ω); simple sine-series rule
        n = np.arange(1, args.N)
        f_n = n/(2.0*args.T); df = 1.0/(2.0*args.T)
        a = ab.a_of_w(2*np.pi*f_n)
        Bpos = np.pi * f_n * a
        j = np.arange(0, args.N+1)[:,None]
        S = np.sin(np.pi * j * n[None,:] / args.N)
        t = np.linspace(0.0, args.T, args.N+1)
        h = -4.0 * df * (S @ Bpos)
        save_tsv(args.out, (t,h), "t[s]\th_dst_dim")
    else:
        # compare
        t1,h1 = invert_real_filon(ab, args.T, args.N)
        t2,h2 = invert_filon(ab, args.T, args.N)
        t3,h3 = invert_fft(ab, args.T, args.N)
        t4,h4,K0,alpha = invert_filon_jr(ab, args.T, args.N, args.alphaT)
        # pad fft to N+1
        if len(t3)<len(t1):
            t3p = np.concatenate([t3,[t3[-1]+(t3[1]-t3[0])]]); h3p=np.concatenate([h3,[h3[-1]]])
        else:
            t3p,h3p=t3,h3
        out = np.column_stack([t1,h1,h2,h3p,h4])
        np.savetxt(args.out, out, delimiter="\t",
                   header="t[s]\th_real_filon_dim\th_filon_dim\th_fft_dim\th_filonJR_dim", comments="")
        print("Wrote", args.out, f"(K0_est={K0:.6g}, alpha={alpha:.6g})")

if __name__ == "__main__":
    main()
