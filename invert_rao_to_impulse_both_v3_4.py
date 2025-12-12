#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
invert_rao_wamit_full_rad.py

Compute dimensional impulse-response functions (IRFs) from WAMIT/SIMO-style
added-mass and damping data (Ā = A/ρ, B̄ = B/(ρ ω)), with frequency given
in rad/s.

We first construct the retardation transform
    H(ω) = c(ω) + i ω a(ω),
where a(ω) = A(ω) - A(∞), c(ω) = B(ω).
This H(ω) is the Fourier transform of the IRF h(t) in the SIMO/WAMIT
convention:
    H(ω) = ∫_{-∞}^{∞} h(t) e^{-iωt} dt,
    h(t) = (1/(2π)) ∫_{-∞}^{∞} H(ω) e^{iωt} dω.

For numerical inversion we map to a "RAO" in cycles/s:
    f = ω / (2π),
    R(f) = H(2π f),

so that
    h(t) = 2 Re ∫_0^∞ R(f) e^{i 2π f t} df.

R(f) has even real part and odd imaginary part. Its imaginary part behaves
as ~ -K0/(2π f) at high |f| if h(t) has a jump K0 at t = 0+.

Methods implemented:

  - dst         : Sine-series quadrature after odd-spectrum correction
  - fft         : Jump-removal + IFFT (periodic trapezoid, exponential tail)
  - filon       : Filon on [0,F] of original spectrum + Si truncation tail
  - filon_exp   : Filon on [0,F] of preconditioned spectrum + exponential tail
  - hilbert_fft : Reconstruct B from A via periodic Hilbert, then JR-FFT
  - hilbert_filon: Reconstruct B from A via Hilbert, then Filon+exp
  - real_filon  : Real-part-only cosine Filon on Â(f)=A(f)-K0 α/(α^2+(2πf)^2), then +K0 e^{-αt}
  - real_dct    : Real-part-only cosine fast path (DCT-I style) on the same Â
  - compare     : run many methods and write consolidated TSV
  - plot        : optional plotting of all curves (requires matplotlib)

Usage:
python invert_rao_to_impulse_both_v3_4.py  --excel impulse_response_all_v3_0.xlsx --sheet AB --T 40 --N 4096 --alphaT 8 --method compare --out K11_irf_all.tsv --plot
"""

import numpy as np
import pandas as pd
import argparse
from typing import Tuple, Dict, List
from scipy.integrate import trapezoid as trapz

# Optional plotting
try:
    import matplotlib.pyplot as plt
except Exception:
    plt = None

# Optional SciPy (not strictly required but detected)
try:
    from scipy.fftpack import dst as scipy_dst, dct as scipy_dct
    from scipy.signal import hilbert as scipy_hilbert
    SCIPY_AVAILABLE = True
except Exception:
    SCIPY_AVAILABLE = False


# ---------------------------------------------------------------------------
# Data loading and dimensionalization
# ---------------------------------------------------------------------------

class ABData:
    """
    Load WAMIT/SIMO-style AB data from Excel, assuming:
        Ā(ω) = A(ω)/ρ
        B̄(ω) = B(ω)/(ρ ω)
    A, B are dimensional (SI) coefficients for one radiation mode.

    After loading, the object stores:
        omega_pos : 1D array of ω >= 0        [rad/s]
        A_dim     : A(ω)                      [kg]
        B_dim     : B(ω)                      [N·s/m]
        A_inf     : estimate of A(∞)
        a_dim     : A(ω) - A_inf              [kg]
        c_dim     : B(ω)                      [N·s/m]
    """
    def __init__(self,
                 path: str,
                 sheet: str = "AB",
                 freq_col: str = "Frequency[rad/s]",
                 A_col: str = "A_ndim=A/ϱ",
                 B_col: str = "B_ndim=B/(ϱω)",
                 rho: float = 1025.0,
                 use_zero: bool = True):
        self.path = path
        self.sheet = sheet
        self.freq_col = freq_col
        self.A_col = A_col
        self.B_col = B_col
        self.rho = rho
        self.use_zero = use_zero
        self._load()

    @staticmethod
    def _find_header_row(df0: pd.DataFrame) -> int:
        for r in range(min(10, len(df0))):
            row = "".join(str(x) for x in list(df0.iloc[r].values))
            if "freq" in row.lower():
                return r
        return 0

    def _load(self) -> None:
        xls = pd.ExcelFile(self.path)
        df0 = xls.parse(self.sheet, header=None)
        hdr = self._find_header_row(df0)
        df  = xls.parse(self.sheet, header=hdr)

        omega_raw = pd.to_numeric(df[self.freq_col], errors="coerce").to_numpy(float)
        Abar      = pd.to_numeric(df[self.A_col],   errors="coerce").to_numpy(float)
        Bbar      = pd.to_numeric(df[self.B_col],   errors="coerce").to_numpy(float)

        m = np.isfinite(omega_raw) & np.isfinite(Abar) & np.isfinite(Bbar)
        omega_raw, Abar, Bbar = omega_raw[m], Abar[m], Bbar[m]

        idx = np.argsort(omega_raw)
        omega_raw, Abar, Bbar = omega_raw[idx], Abar[idx], Bbar[idx]

        if not self.use_zero and omega_raw.size > 0 and abs(omega_raw[0]) < 1e-12:
            omega_raw, Abar, Bbar = omega_raw[1:], Abar[1:], Bbar[1:]

        mpos = omega_raw >= 0.0
        omega = omega_raw[mpos]
        Abar  = Abar[mpos]
        Bbar  = Bbar[mpos]

        # Dimensional
        self.omega_pos = omega
        self.A_dim = self.rho * Abar
        self.B_dim = self.rho * omega * Bbar

        # A(∞) from top 10%
        if self.A_dim.size < 5:
            self.A_inf = float(self.A_dim[-1])
        else:
            mA = max(4, int(0.1 * self.A_dim.size))
            self.A_inf = float(np.mean(self.A_dim[-mA:]))

        self.a_dim = self.A_dim - self.A_inf
        self.c_dim = self.B_dim

    def _interp_even(self, omega: np.ndarray, y_pos: np.ndarray, right: float = 0.0) -> np.ndarray:
        w = np.abs(omega)
        return np.interp(w, self.omega_pos, y_pos, left=y_pos[0], right=right)

    def a_of_omega(self, omega: np.ndarray) -> np.ndarray:
        return self._interp_even(omega, self.a_dim, right=0.0)

    def c_of_omega(self, omega: np.ndarray) -> np.ndarray:
        return self._interp_even(omega, self.c_dim, right=0.0)

    def H_of_omega(self, omega: np.ndarray) -> np.ndarray:
        omega = np.asarray(omega, float)
        c = self.c_of_omega(omega)
        a = self.a_of_omega(omega)
        return c + 1j * omega * a

    def R_of_f(self, f: np.ndarray) -> np.ndarray:
        f = np.asarray(f, float)
        omega = 2.0 * np.pi * f
        return self.H_of_omega(omega)


# ---------------------------------------------------------------------------
# Sine integral and tail fitting
# ---------------------------------------------------------------------------

def Si(x: float) -> float:
    """Sine integral Si(x) = ∫_0^x sin(u)/u du."""
    try:
        import mpmath as mp
        return float(mp.si(x))
    except Exception:
        u = np.linspace(0.0, x, 20001)
        u[0] = 1e-12
        return float(trapz(np.sin(u)/u, u))


def fit_tails(ab: ABData, frac: float = 0.3) -> Dict[str, float]:
    """
    Fit high-ω tails:
        a(ω) ~ C_a / ω^2,   c(ω) ~ C_c / ω^2.
    Approximate using top 'frac' of the band via medians of a ω^2 and c ω^2.
    """
    omega = ab.omega_pos
    a = ab.a_dim
    c = ab.c_dim
    n = len(omega)
    i0 = int((1.0 - frac) * n)
    i0 = max(0, min(i0, n-4))
    w_tail = omega[i0:]
    a_tail = a[i0:]
    c_tail = c[i0:]
    w2 = np.maximum(w_tail, 1e-9)**2
    C_a = float(np.median(a_tail * w2))
    C_c = float(np.median(c_tail * w2))
    return dict(C_a=C_a, C_c=C_c, Omega=float(omega[-1]))


def K0_from_c_with_tail(ab: ABData, tails: Dict[str,float]) -> float:
    """
    Jump from damping-only route:
        K0 = h(0) = (2/π)∫_0^∞ c(ω)dω
    with c(ω)~C_c/ω^2 for ω>Ω.
    """
    omega = ab.omega_pos
    c = ab.c_dim
    C_c = tails["C_c"]; Omega = tails["Omega"]
    I_meas = trapz(c, omega)
    I_tail = C_c / Omega if Omega > 0 else 0.0
    return (2.0/np.pi) * (I_meas + I_tail)


# ---------------------------------------------------------------------------
# Frequency/time grids, Filon, Hilbert
# ---------------------------------------------------------------------------

def choose_fft_grid(T: float, N: int, f_max_data: float) -> Tuple[np.ndarray,float,int]:
    """
    Choose symmetric FFT frequency grid [-F,F) with N points and Δf=1/T,
    ensuring Nyquist F_Nyq=(N/2)/T does not exceed f_max_data.
    """
    Nyq = (N/2.0) / T
    if Nyq > f_max_data:
        N_new = int(2 * T * f_max_data)
        N_new = max(32, N_new)
        N_new = N_new - (N_new % 2)  # even
        N = N_new
    df = 1.0 / T
    k = np.arange(-N//2, N//2, dtype=float)
    f = k * df
    return f, df, N


def filon_uniform_f(f: np.ndarray, g: np.ndarray, t: float) -> complex:
    """
    Composite quadratic Filon on a uniform f-grid [0,F], length 2M+1.
    Approximates ∫_0^F g(f) e^{i2π f t} df.
    """
    f = np.asarray(f, float)
    g = np.asarray(g, complex)
    Np = len(f)
    Niv = Np - 1
    assert Niv % 2 == 0, "Filon requires an even number of intervals"
    h = f[1] - f[0]
    y0 = g[0:-2:2]
    y1 = g[1:-1:2]
    y2 = g[2::2]
    a = y1
    b = (y2 - y0) / (2*h)
    c = (y0 - 2*y1 + y2) / (2*h**2)
    f_mid = f[1:-1:2]
    mu = 2.0 * np.pi * t
    if abs(mu) < 1e-14:
        I0, I1, I2 = 2*h, 0.0j, 2*h**3/3.0
    else:
        th = mu * h
        I0 = 2*np.sin(th) / mu
        I1 = 2j * (-h*np.cos(th)/mu + np.sin(th)/mu**2)
        I2 = 2*(h**2*np.sin(th)/mu + 2*h*np.cos(th)/mu**2 - 2*np.sin(th)/mu**3)
    return np.sum(np.exp(1j*mu*f_mid) * (a*I0 + b*I1 + c*I2))


def periodic_hilbert_real(x: np.ndarray) -> np.ndarray:
    """
    Periodic Hilbert transform H[x](f) for real x(f) sampled on a uniform grid.
    Returns real array ≈ Hilbert{x}. Uses FFT with multiplier -i*sign(k).
    """
    x = np.asarray(x, float)
    N = len(x)
    X = np.fft.fft(x)
    k = np.fft.fftfreq(N)
    sign = np.sign(k)
    Hmult = -1j * sign
    Hmult[0] = 0.0
    XH = X * Hmult
    h = np.fft.ifft(XH)
    return np.real(h)


# ---------------------------------------------------------------------------
# Inversion methods
# ---------------------------------------------------------------------------

def method_fft(ab: ABData, T: float, N: int, alphaT: float,
               tails: Dict[str,float]) -> Tuple[np.ndarray,np.ndarray,float,float]:
    """
    Jump-removal + IFFT on R(f).
    Steps:
      1) Estimate K0 from damping tail.
      2) Choose symmetric f-grid consistent with data band.
      3) Build R(f) and subtract simple pole K0/(α+i2πf).
      4) IFFT to get K̂(t); add back K0 e^{-α t}.
    """
    K0 = K0_from_c_with_tail(ab, tails)
    alpha = alphaT / T

    f_pos_data = ab.omega_pos / (2*np.pi)
    f_max_data = float(np.max(f_pos_data))
    f, df, N_used = choose_fft_grid(T, N, f_max_data)
    R = ab.R_of_f(f)

    denom = alpha + 1j * 2.0*np.pi*f
    R_hat = R - K0 / denom

    R_shift = np.fft.ifftshift(R_hat)
    K_hat_per = np.fft.ifft(R_shift) * N_used * df
    t = np.linspace(0.0, T, N_used, endpoint=False)
    K = np.real(K_hat_per) + K0 * np.exp(-alpha*t)
    return t, K, K0, alpha


def method_filon(ab: ABData, T: float, N: int,
                 tails: Dict[str,float]) -> Tuple[np.ndarray,np.ndarray,float]:
    """
    Filon on [0,F] for the ORIGINAL R(f), plus Si-tail modeling the truncated
    1/f part of Im R(f). For large f, Im R(f) ~ -K0/(2π f).
    Tail over [F,∞) contributes
        K_tail(t) = (K0/π)*(Si(2π F t) - π/2).
    """
    K0 = K0_from_c_with_tail(ab, tails)

    f_pos_data = ab.omega_pos / (2*np.pi)
    F = float(np.max(f_pos_data))
    Niv = 20000
    if Niv % 2 == 1: Niv += 1
    f = np.linspace(0.0, F, Niv+1)
    R = ab.R_of_f(f)

    t = np.linspace(0.0, T, N+1)
    K = np.zeros_like(t)
    for j, tj in enumerate(t):
        I = filon_uniform_f(f, R, tj)
        if tj == 0.0:
            tail = (K0/np.pi) * (-np.pi/2.0)
        else:
            tail = (K0/np.pi) * (Si(2*np.pi*F*tj) - np.pi/2.0)
        K[j] = 4.0 * np.real(I) + tail
    return t, K, K0


def method_filon_exp(ab: ABData, T: float, N: int, alphaT: float,
                     tails: Dict[str,float]) -> Tuple[np.ndarray,np.ndarray,float,float]:
    """
    Filon on [0,F] for the preconditioned spectrum (jump removed) plus an
    exponential tail. This mirrors method_fft but with Filon quadrature.
    """
    K0 = K0_from_c_with_tail(ab, tails)
    alpha = alphaT / T

    f_pos_data = ab.omega_pos / (2*np.pi)
    F = float(np.max(f_pos_data))
    Niv = 20000
    if Niv % 2 == 1: Niv += 1
    f = np.linspace(0.0, F, Niv+1)
    R = ab.R_of_f(f)
    denom = alpha + 1j * 2.0*np.pi*f
    R_hat = R - K0/denom

    t = np.linspace(0.0, T, N+1)
    K = np.zeros_like(t)
    for j, tj in enumerate(t):
        I = filon_uniform_f(f, R_hat, tj)
        K[j] = 4.0*np.real(I) + K0*np.exp(-alpha*tj)
    return t, K, K0, alpha


def method_real_filon(ab: ABData, T: float, N: int, alphaT: float,
                      tails: Dict[str,float]) -> Tuple[np.ndarray,np.ndarray,float,float]:
    """
    Real-part-only Filon: operate on Â(f) = A(f) - K0 α/(α^2+(2πf)^2);
    this cancels the 1/f part of the imaginary component (via Hilbert).
    Then add back K0 e^{-α t}.
    """
    K0 = K0_from_c_with_tail(ab, tails)
    alpha = alphaT / T

    f_pos_data = ab.omega_pos / (2*np.pi)
    F = float(np.max(f_pos_data))
    Niv = 20000
    if Niv % 2 == 1: Niv += 1
    f = np.linspace(0.0, F, Niv+1)
    R = ab.R_of_f(f)
    A = np.real(R)
    A_hat = A - K0 * alpha / (alpha**2 + (2*np.pi*f)**2)

    t = np.linspace(0.0, T, N+1)
    K = np.zeros_like(t)
    for j, tj in enumerate(t):
        I = filon_uniform_f(f, A_hat, tj)
        K[j] = 4.0*np.real(I) + K0*np.exp(-alpha*tj)
    return t, K, K0, alpha


def method_real_dct(ab: ABData, T: float, N: int, alphaT: float,
                    tails: Dict[str,float]) -> Tuple[np.ndarray,np.ndarray,float,float]:
    """
    Real-part-only cosine 'fast path' using a DCT-I-like cosine sum on Â(f).
    """
    K0 = K0_from_c_with_tail(ab, tails)
    alpha = alphaT / T

    n = np.arange(0, N+1, dtype=float)
    f_n = n / (2.0*T)
    R = ab.R_of_f(f_n)
    A = np.real(R)
    A_hat = A - K0 * alpha / (alpha**2 + (2*np.pi*f_n)**2)

    t = np.linspace(0.0, T, N+1)
    df = 1.0/(2.0*T)
    K_hat = np.zeros_like(t)
    for j, tj in enumerate(t):
        K_hat[j] = 4.0 * df * np.sum(A_hat * np.cos(2*np.pi * f_n * tj))
    K = K_hat + K0*np.exp(-alpha*t)
    return t, K, K0, alpha


def method_dst(ab: ABData, T: float, N: int,
               tails: Dict[str,float]) -> Tuple[np.ndarray,np.ndarray]:
    """
    Sine-series quadrature using only the odd spectrum of Im R(f):
        K(t_j) ≈ -4 Δf Σ_{n=1}^{N-1} B(f_n) sin(π j n / N),
    with Δf = 1/(2T), f_n = n/(2T), t_j = j T/N.
    """
    n = np.arange(1, N, dtype=float)
    f_n = n / (2.0*T)
    df = 1.0/(2.0*T)
    R_pos = ab.R_of_f(f_n)
    B_pos = np.imag(R_pos)

    t = np.linspace(0.0, T, N+1)
    j = np.arange(0, N+1, dtype=float)[:,None]
    S = np.sin(np.pi * j * n[None,:] / N)
    K = -4.0 * df * (S @ B_pos)
    return t, K


def method_hilbert_fft(ab: ABData, T: float, N: int, alphaT: float,
                       tails: Dict[str,float]) -> Tuple[np.ndarray,np.ndarray,float,float]:
    """
    Reconstruct Im R(f) from Re R(f) via periodic Hilbert on the FFT grid,
    then apply JR-FFT to that reconstructed R(f).
    """
    K0 = K0_from_c_with_tail(ab, tails)
    alpha = alphaT / T

    f_pos_data = ab.omega_pos / (2*np.pi)
    f_max_data = float(np.max(f_pos_data))
    f, df, N_used = choose_fft_grid(T, N, f_max_data)
    R_true = ab.R_of_f(f)
    A = np.real(R_true)

    B_rec = periodic_hilbert_real(A)
    R_rec = A + 1j*B_rec

    denom = alpha + 1j*2.0*np.pi*f
    R_hat = R_rec - K0/denom
    R_shift = np.fft.ifftshift(R_hat)
    K_hat_per = np.fft.ifft(R_shift) * N_used * df
    t = np.linspace(0.0, T, N_used, endpoint=False)
    K = np.real(K_hat_per) + K0*np.exp(-alpha*t)
    return t, K, K0, alpha


def method_hilbert_filon(ab: ABData, T: float, N: int, alphaT: float,
                         tails: Dict[str,float]) -> Tuple[np.ndarray,np.ndarray,float,float]:
    """
    Reconstruct Im R(f) from Re R(f) via Hilbert on a symmetric f-grid and
    then apply Filon+exp on the reconstructed spectrum.
    """
    K0 = K0_from_c_with_tail(ab, tails)
    alpha = alphaT / T

    f_pos_data = ab.omega_pos / (2*np.pi)
    F = float(np.max(f_pos_data))
    Niv = 20000
    if Niv % 2 == 1: Niv += 1
    f = np.linspace(0.0, F, Niv+1)
    R_true = ab.R_of_f(f)
    A = np.real(R_true)

    # Build symmetric grid for Hilbert
    f_sym = np.concatenate((-f[:0:-1], f))
    A_sym = np.concatenate((A[:0:-1], A))
    B_sym = periodic_hilbert_real(A_sym)
    B = B_sym[len(f_sym)-len(f):]

    R_rec = A + 1j*B
    denom = alpha + 1j*2.0*np.pi*f
    R_hat = R_rec - K0/denom

    t = np.linspace(0.0, T, N+1)
    K = np.zeros_like(t)
    for j, tj in enumerate(t):
        I = filon_uniform_f(f, R_hat, tj)
        K[j] = 4.0*np.real(I) + K0*np.exp(-alpha*tj)
    return t, K, K0, alpha


# ---------------------------------------------------------------------------
# Driver / I/O
# ---------------------------------------------------------------------------

def save_tsv(path: str, data: np.ndarray, header: str) -> None:
    np.savetxt(path, data, delimiter="\t", header=header, comments="")


def run_methods(ab: ABData, T: float, N: int, alphaT: float,
                methods: List[str], out_path: str,
                do_plot: bool = False) -> None:
    tails = fit_tails(ab, frac=0.3)
    results = {}
    t_ref = None

    for m in methods:
        if m == "dst":
            t, K = method_dst(ab, T, N, tails)
            results["K_dst"] = (t, K)
        elif m == "fft":
            t, K, K0, alpha = method_fft(ab, T, N, alphaT, tails)
            results["K_fft"] = (t, K)
        elif m == "filon":
            t, K, K0 = method_filon(ab, T, N, tails)
            results["K_filon"] = (t, K)
        elif m == "filon_exp":
            t, K, K0, alpha = method_filon_exp(ab, T, N, alphaT, tails)
            results["K_filon_exp"] = (t, K)
        elif m == "hilbert_fft":
            t, K, K0, alpha = method_hilbert_fft(ab, T, N, alphaT, tails)
            results["K_hilbert_fft"] = (t, K)
        elif m == "hilbert_filon":
            t, K, K0, alpha = method_hilbert_filon(ab, T, N, alphaT, tails)
            results["K_hilbert_filon"] = (t, K)
        elif m == "real_filon":
            t, K, K0, alpha = method_real_filon(ab, T, N, alphaT, tails)
            results["K_real_filon"] = (t, K)
        elif m == "real_dct":
            t, K, K0, alpha = method_real_dct(ab, T, N, alphaT, tails)
            results["K_real_dct"] = (t, K)
        else:
            raise ValueError(f"Unknown method {m}")

        if t_ref is None:
            t_ref = t
        else:
            if len(t) != len(t_ref) or not np.allclose(t, t_ref):
                K_res = np.interp(t_ref, t, K)
                results[list(results.keys())[-1]] = (t_ref, K_res)

    # Consolidated TSV
    keys = sorted(results.keys())
    cols = [t_ref]
    header_fields = ["t[s]"]
    for k in keys:
        cols.append(results[k][1])
        header_fields.append(k)
    data = np.column_stack(cols)
    save_tsv(out_path, data, "\t".join(header_fields))

    if do_plot and plt is not None:
        for k in keys:
            plt.plot(t_ref, results[k][1], label=k)
        plt.xlabel("t [s]")
        plt.ylabel("K(t) [N/m]")
        plt.grid(True)
        plt.legend()
        plt.show()


def main():
    ap = argparse.ArgumentParser(description="Invert WAMIT/SIMO AB data to IRF using multiple methods.")
    #ap.add_argument("--excel", type=str, required=True, help="Path to Excel file with AB data")
    ap.add_argument("--excel", type=str, default=f'E:\\WORK\\2368SSY_GSM\\Wamit\\T9p23\\impulse_response_all_v3_0.xlsx', help="Path to Excel file with AB data")
    ap.add_argument("--sheet", type=str, default="AB", help="Sheet name")
    ap.add_argument("--freqcol", type=str, default="Frequency[rad/s]", help="Frequency column (rad/s)")
    ap.add_argument("--acol", type=str, default="A_ndim=A/ϱ", help="Nondimensional A column (A/ρ)")
    ap.add_argument("--bcol", type=str, default="B_ndim=B/(ϱω)", help="Nondimensional B column (B/(ρ ω))")
    ap.add_argument("--rho", type=float, default=1025.0, help="Water density [kg/m^3]")
    ap.add_argument("--use_zero", type=int, default=1, help="Include ω=0 row (1=yes, 0=no)")
    ap.add_argument("--T", type=float, default=40.0, help="Time-window length [s]")
    ap.add_argument("--N", type=int, default=4096, help="Number of time samples (FFT/DST)")
    ap.add_argument("--alphaT", type=float, default=8.0, help="αT for exponential tail in JR methods")
    ap.add_argument("--method", type=str, default="compare",
                    choices=["dst","fft","filon","filon_exp",
                             "hilbert_fft","hilbert_filon",
                             "real_filon","real_dct","compare"])
    ap.add_argument("--out", type=str, default="irf_all.tsv", help="Output TSV file")
    ap.add_argument("--plot", action="store_true", help="Plot all curves")
    args = ap.parse_args()

    ab = ABData(args.excel, sheet=args.sheet,
                freq_col=args.freqcol,
                A_col=args.acol, B_col=args.bcol,
                rho=args.rho, use_zero=bool(args.use_zero))

    if args.method == "compare":
        methods = ["dst","fft","filon","filon_exp",
                   "hilbert_fft","hilbert_filon",
                   "real_filon","real_dct"]
        run_methods(ab, args.T, args.N, args.alphaT, methods, args.out, do_plot=args.plot)
    else:
        run_methods(ab, args.T, args.N, args.alphaT, [args.method], args.out, do_plot=args.plot)


if __name__ == "__main__":
    main()
