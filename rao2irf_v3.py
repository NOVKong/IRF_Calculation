#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import argparse, os, math, sys
from typing import Tuple

# ------------------------ Excel loader (AB sheet: Frequency[Hz], A_ndim=A/ρ, B_ndim=B/(ρ ω)) ------------------------
_EXCEL_CFG = {
    "path": None, "sheet": "AB",
    "freq_col": "Frequency[Hz]", "A_col": "A_ndim=A/ϱ", "B_col": "B_ndim=B/(ϱω)",
    "use_zero": True
}
_AB_CACHE = {"loaded": False}

def _find_header_row(df: pd.DataFrame) -> int:
    for r in range(min(10, len(df))):
        row = "".join(str(x) for x in list(df.iloc[r].values))
        if "freq" in row.lower():
            return r
    return 0

def load_ab_excel():
    if _AB_CACHE["loaded"]:
        return _AB_CACHE
    if _EXCEL_CFG["path"] is None:
        raise RuntimeError("No Excel path configured; pass --excel <file.xlsx>")
    xls = pd.ExcelFile(_EXCEL_CFG["path"])
    df0 = xls.parse(_EXCEL_CFG["sheet"], header=None)
    hdr = _find_header_row(df0)
    df  = xls.parse(_EXCEL_CFG["sheet"], header=hdr)
    f_hz = pd.to_numeric(df[_EXCEL_CFG["freq_col"]], errors="coerce").to_numpy(float)
    Abar = pd.to_numeric(df[_EXCEL_CFG["A_col"]],   errors="coerce").to_numpy(float)  # A/ρ
    Bbar = pd.to_numeric(df[_EXCEL_CFG["B_col"]],   errors="coerce").to_numpy(float)  # B/(ρ ω) = c(ω)

    m = np.isfinite(f_hz) & np.isfinite(Abar) & np.isfinite(Bbar)
    f_hz, Abar, Bbar = f_hz[m], Abar[m], Bbar[m]
    idx = np.argsort(f_hz); f_hz, Abar, Bbar = f_hz[idx], Abar[idx], Bbar[idx]

    if not _EXCEL_CFG["use_zero"] and f_hz.size>0 and abs(f_hz[0])<1e-12:
        f_hz, Abar, Bbar = f_hz[1:], Abar[1:], Bbar[1:]

    w = 2*np.pi*f_hz
    # A∞ estimate: average of top 10%
    mA = max(4, int(0.1*Abar.size))
    Ainf = float(np.mean(Abar[-mA:]))
    _AB_CACHE.update(dict(loaded=True, f_pos=f_hz, w_pos=w, Abar=Abar, Bbar=Bbar, Ainf=Ainf))
    return _AB_CACHE

# ------------------------ IRF-convention spectrum: build R(f)=A(f)+iB(f) so that h(t) = 2 Re ∫ R e^{i2π f t} df --------
def A_of_f(ff: np.ndarray) -> np.ndarray:
    """
    Even part: A(f) = c(2π|f|) = B̄(2π|f|).
    """
    load_ab_excel()
    f = np.asarray(ff, float)
    w = 2*np.pi*np.abs(f)
    Bbar = np.interp(w, _AB_CACHE["w_pos"], _AB_CACHE["Bbar"], left=_AB_CACHE["Bbar"][0], right=0.0)
    return Bbar

def B_of_f(ff: np.ndarray) -> np.ndarray:
    """
    Odd part: B(f) = sign(f) * π|f| * a(2π|f|), with a(ω)=Ā(ω)-Ā∞.
    This choice ensures that
      2 Re ∫ (A + iB) e^{i2π f t} df = (2/π) ∫ c(ω) cos(ωt) dω - (2/π) ∫ ω a(ω) sin(ωt) dω.
    """
    load_ab_excel()
    f = np.asarray(ff, float)
    s = np.sign(f)
    w = 2*np.pi*np.abs(f)
    Abar = np.interp(w, _AB_CACHE["w_pos"], _AB_CACHE["Abar"], left=_AB_CACHE["Abar"][0], right=_AB_CACHE["Ainf"])
    a = Abar - _AB_CACHE["Ainf"]
    return s * np.pi * np.abs(f) * a

# ------------------------ Quadratures and transforms ------------------------

def filon_block_uniform_f(f: np.ndarray, g: np.ndarray, t: float) -> complex:
    """
    Composite quadratic Filon on a *uniform* f-grid for ∫_0^F g(f) e^{i 2π f t} df
    f: length 2M+1, g defined at same nodes.
    """
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

# ------------------------ Methods ------------------------

def invert_real_filon(T: float, N: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Cosine-only Filon using A(f)=c(2π|f|), and then multiply by 4 as shown in the derivation.
    """
    load_ab_excel()
    f_pos_raw = _AB_CACHE["f_pos"]
    F = float(np.max(f_pos_raw))
    Niv = 20000
    if Niv % 2 == 1: Niv += 1
    f = np.linspace(0.0, F, Niv+1)
    A = A_of_f(f)  # c(2π f)
    t = np.linspace(0.0, T, N+1)
    K = np.zeros_like(t)
    for j, tj in enumerate(t):
        I = filon_block_uniform_f(f, A, tj)
        K[j] = 4.0*np.real(I)
    return t, K

def invert_filon(T: float, N: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Complex Filon on [0,F] with R(f)=A(f)+iB(f). No jump (K0=0) for radiation IRF.
    """
    load_ab_excel()
    f_pos_raw = _AB_CACHE["f_pos"]
    F = float(np.max(f_pos_raw))
    Niv = 20000
    if Niv % 2 == 1: Niv += 1
    f = np.linspace(0.0, F, Niv+1)
    R = A_of_f(f) + 1j*B_of_f(f)
    t = np.linspace(0.0, T, N+1)
    K = np.zeros_like(t)
    for j, tj in enumerate(t):
        I = filon_block_uniform_f(f, R, tj)
        K[j] = 4.0*np.real(I)
    return t, K

def invert_fft(T: float, N: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Periodic trapezoid + IFFT on a symmetric grid. Builds R at f_k = k/T, k=-N/2..N/2-1.
    The discrete sum approximates ∫_{-F}^{F} R(f) e^{i2πft} df; with Hermitian symmetry,
    the IFFT returns a *real* time series that is already 2·Re{·}. We multiply by Δf.
    """
    f = np.fft.fftfreq(N, d=T/N)  # k/T, k=0..N-1; arranged [0, +, -, ...]
    f = np.fft.fftshift(f)        # symmetric grid
    df = 1.0/T
    R = A_of_f(f) + 1j*B_of_f(f)
    Rshift = np.fft.ifftshift(R)
    Kper = np.fft.ifft(Rshift) * N * df   # rectangular sum over [-F,F]
    t = np.linspace(0.0, T, N, endpoint=False)
    K = np.real(Kper)
    return t, K

def invert_dst(T: float, N: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Sine-series quadrature on [0,F] using only the odd spectrum B(f)=π f a(2π f).
    K(t) ≈ -4 ∫ B(f) sin(2π f t) df  (by construction equals h(t)).
    """
    load_ab_excel()
    n = np.arange(1, N)  # interior modes
    f_n = n / (2.0*T)
    df = 1.0 / (2.0*T)
    B_pos = B_of_f(f_n)
    t = np.linspace(0.0, T, N+1)
    j = np.arange(0, N+1)[:, None]
    S = np.sin(np.pi * j * n[None, :] / N)
    K = -4.0 * df * (S @ B_pos)
    return t, K

# ------------------------ CLI / I/O ------------------------

def save_tsv(path: str, t: np.ndarray, K: np.ndarray, label: str):
    hdr = "t[s]	" + label
    np.savetxt(path, np.column_stack([t, K]), delimiter="\t", header=hdr, comments="")

def main():
    ap = argparse.ArgumentParser(description="Invert WAMIT/SIMO radiation data (AB) to IRF h(t) under strict conventions.")
    ap.add_argument("--excel", type=str, required=True, help="Path to Excel with AB data (sheet AB)")
    ap.add_argument("--sheet", type=str, default="AB", help="Sheet name (default AB)")
    ap.add_argument("--freqcol", type=str, default="Frequency[Hz]")
    ap.add_argument("--acol", type=str, default="A_ndim=A/ϱ")
    ap.add_argument("--bcol", type=str, default="B_ndim=B/(ϱω)")
    ap.add_argument("--use_zero", type=int, default=1, help="Include ω=0 row (1=yes, 0=no)")
    ap.add_argument("--method", choices=["real_filon","filon","fft","dst","compare"], default="compare")
    ap.add_argument("--T", type=float, default=40.0)
    ap.add_argument("--N", type=int, default=4096)
    ap.add_argument("--out", type=str, default="/mnt/data/irf_h_t.tsv")
    args = ap.parse_args()

    _EXCEL_CFG.update({"path": args.excel, "sheet": args.sheet, "freq_col": args.freqcol,
                       "A_col": args.acol, "B_col": args.bcol, "use_zero": bool(args.use_zero)})

    if args.method == "real_filon":
        t,K = invert_real_filon(args.T, args.N)
        save_tsv(args.out, t, K, "h_real_filon")
    elif args.method == "filon":
        t,K = invert_filon(args.T, args.N)
        save_tsv(args.out, t, K, "h_filon")
    elif args.method == "fft":
        t,K = invert_fft(args.T, args.N)
        save_tsv(args.out, t, K, "h_fft")
    elif args.method == "dst":
        t,K = invert_dst(args.T, args.N)
        save_tsv(args.out, t, K, "h_dst")
    else:  # compare
        t1,K1 = invert_real_filon(args.T, args.N)
        t2,K2 = invert_filon(args.T, args.N)
        t3,K3 = invert_fft(args.T, args.N)
        t4,K4 = invert_dst(args.T, args.N)
        # Align FFT series to (N+1) by padding one value for consistent TSV
        if len(t3) < len(t1):
            t3p = np.concatenate([t3, [t3[-1]+(t3[1]-t3[0])]])
            K3p = np.concatenate([K3, [K3[-1]]])
        else:
            t3p, K3p = t3, K3
        out = np.column_stack([t1, K1, K2, K3p, K4])
        hdr = "t[s]	h_real_filon	h_filon	h_fft	h_dst"
        np.savetxt(args.out, out, delimiter="\t", header=hdr, comments="")
        print("Wrote", args.out)

if __name__ == "__main__":
    main()

