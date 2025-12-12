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
        A_inf     : A(∞)                      [kg]
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
        """Find the row that contains something like 'freq' (case-insensitive)."""
        for r in range(min(10, len(df0))):
            row = "".join(str(x) for x in list(df0.iloc[r].values))
            if "freq" in row.lower():
                return r
        return 0

    @staticmethod
    def _resolve_col(df: pd.DataFrame, target: str, hints: List[str]) -> str:
        """
        Robustly resolve a column name in df:
          - try exact match
          - try case/whitespace-insensitive match
          - then look for columns containing any hint substring
        """
        cols = list(df.columns)

        # exact
        if target in cols:
            return target

        # case/whitespace-insensitive
        t_norm = target.strip().lower()
        for c in cols:
            if c.strip().lower() == t_norm:
                return c

        # hint-based
        hints_norm = [h.lower() for h in hints]
        for h in hints_norm:
            for c in cols:
                if h in c.strip().lower():
                    return c

        raise KeyError(f"Could not find column matching '{target}' with hints {hints}. "
                       f"Available columns: {cols}")

    def _load(self) -> None:
        xls = pd.ExcelFile(self.path)
        df0 = xls.parse(self.sheet, header=None)
        hdr = self._find_header_row(df0)
        df  = xls.parse(self.sheet, header=hdr)

        # Robust column resolution
        freq_col = self._resolve_col(df, self.freq_col, hints=["freq", "frequency"])
        A_col    = self._resolve_col(df, self.A_col,    hints=["a_ndim", "addedmass", "a/"])
        B_col    = self._resolve_col(df, self.B_col,    hints=["b_ndim", "damping", "b/"])

        # Raw series (keep strings to detect 'Infinity')
        freq_raw = df[freq_col].astype(str)
        Abar_raw = pd.to_numeric(df[A_col], errors="coerce")
        Bbar_raw = pd.to_numeric(df[B_col], errors="coerce")

        # --- detect Infinity row(s) in frequency column ---
        mask_inf = freq_raw.str.contains("inf", case=False, na=False)

        A_inf_from_row = None
        if mask_inf.any():
            Abar_inf = Abar_raw[mask_inf]
            Abar_inf = Abar_inf[np.isfinite(Abar_inf)]
            if len(Abar_inf) > 0:
                # priority: use this as nondimensional A(∞)
                A_inf_from_row = float(Abar_inf.mean()) * self.rho  # dimensional A∞

        # Drop Infinity rows for the numeric arrays
        mask_fin_freq = ~mask_inf
        freq_num = pd.to_numeric(df.loc[mask_fin_freq, freq_col], errors="coerce")
        Abar     = Abar_raw[mask_fin_freq]
        Bbar     = Bbar_raw[mask_fin_freq]

        # keep only rows where numeric freq & Abar & Bbar are finite
        m = np.isfinite(freq_num) & np.isfinite(Abar) & np.isfinite(Bbar)
        omega_raw = freq_num.to_numpy(float)[m]
        Abar      = Abar.to_numpy(float)[m]
        Bbar      = Bbar.to_numpy(float)[m]

        # sort by ω
        idx = np.argsort(omega_raw)
        omega_raw, Abar, Bbar = omega_raw[idx], Abar[idx], Bbar[idx]

        # optionally drop ω=0 row
        if not self.use_zero and omega_raw.size > 0 and abs(omega_raw[0]) < 1e-12:
            omega_raw, Abar, Bbar = omega_raw[1:], Abar[1:], Bbar[1:]

        # only ω ≥ 0
        mpos = omega_raw >= 0.0
        omega = omega_raw[mpos]
        Abar  = Abar[mpos]
        Bbar  = Bbar[mpos]

        # --- dimensionalization ---
        self.omega_pos = omega                         # [rad/s]
        self.A_dim = self.rho * Abar                   # [kg]
        self.B_dim = self.rho * omega * Bbar           # [N·s/m]

        # --- A(∞): Infinity row takes priority ---
        if A_inf_from_row is not None:
            self.A_inf = A_inf_from_row
        else:
            # fallback: top 10% average
            if self.A_dim.size < 5:
                self.A_inf = float(self.A_dim[-1])
            else:
                mA = max(4, int(0.1 * self.A_dim.size))
                self.A_inf = float(np.mean(self.A_dim[-mA:]))

        # radiation a(ω), c(ω)
        self.a_dim = self.A_dim - self.A_inf
        self.c_dim = self.B_dim

    # --- even extensions and transforms ---

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
    Compute K0 = h(0) from the damping c(ω) using

        K0 = (2/π) * ∫_0^∞ c(ω) dω

    with:
      - on [0,Ω]: trapezoidal rule with *end correction* (Benthien eq. 2.65),
      - on [Ω,∞): analytic tail using c(ω) ~ C_c / ω^2.

    This significantly improves the accuracy of K0, especially for methods
    that rely directly on its absolute value (Filon, real_filon, real_dct, etc.).
    """
    omega = ab.omega_pos        # [rad/s], assumed increasing and ≈ uniform
    c     = ab.c_dim            # damping c(ω) [N·s/m]
    C_c   = tails["C_c"]
    Omega = tails["Omega"]

    # --- integral over measured band [0, Ω] with end-corrected trapezoid ---
    if len(omega) >= 3:
        # approximate uniform step
        h = (omega[-1] - omega[0]) / (len(omega) - 1)
        # check uniformity; if very non-uniform, fall back to np.trapz
        if np.max(np.abs(np.diff(omega) - h)) < 1e-6 * max(1.0, h):
            n = len(omega) - 1

            # standard trapezoid on [0, Ω]
            I_trap = h * (0.5 * c[0] + np.sum(c[1:-1]) + 0.5 * c[-1])

            # endpoint derivatives f'(a), f'(b) via 2nd-order one-sided differences
            fpa = (-3.0 * c[0] + 4.0 * c[1] - 1.0 * c[2]) / (2.0 * h)
            fpb = ( 3.0 * c[-1] - 4.0 * c[-2] + 1.0 * c[-3]) / (2.0 * h)

            # Benthien trapezoid with end correction (eq. 2.65)
            I_meas = I_trap - (h ** 2 / 12.0) * (fpb - fpa)
        else:
            # non-uniform backup: plain trapezoid
            I_meas = trapz(c, omega)
    else:
        I_meas = trapz(c, omega)

    # --- tail ∫_Ω^∞ C_c / ω^2 dω = C_c / Ω ---
    I_tail = C_c / Omega if Omega > 0.0 else 0.0

    # --- K0 from total integral ---
    K0 = (2.0 / np.pi) * (I_meas + I_tail)
    return K0


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
    Real-part-only Filon: operate on Â(f) = A(f) - K0 α/(α^2+(2πf)^2),
    then add back K0 e^{-α t}. Uses the correct factor 2
    (h_hat ≈ 2 ∫_0^F Â cos).
    """
    K0 = K0_from_c_with_tail(ab, tails)
    alpha = alphaT / T

    f_pos_data = ab.omega_pos / (2*np.pi)
    F = float(np.max(f_pos_data))
    Niv = 20000
    if Niv % 2 == 1:
        Niv += 1
    f = np.linspace(0.0, F, Niv+1)
    R = ab.R_of_f(f)
    A = np.real(R)
    A_hat = A - K0 * alpha / (alpha**2 + (2*np.pi*f)**2)

    t = np.linspace(0.0, T, N+1)
    K = np.zeros_like(t)
    for j, tj in enumerate(t):
        I = filon_uniform_f(f, A_hat, tj)   # ≈ ∫_0^F A_hat e^{i2π f t} df
        K[j] = 4.0 * np.real(I) + K0*np.exp(-alpha*tj)
    return t, K, K0, alpha


def method_real_dct(ab: ABData, T: float, N: int, alphaT: float,
                    tails: Dict[str,float]) -> Tuple[np.ndarray,np.ndarray,float,float]:
    """
    Real-part-only cosine fast path using Â(f) on the natural DCT-I grid:
      f_n = n/(2T), n=0..N.

    Correct scaling:
      h_hat(t_j) ≈ 2 Δf Σ_n Â(f_n) cos(2π f_n t_j),
      Δf = 1/(2T),
      t_j = j T/N, 2π f_n t_j = π n j / N.
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

    if SCIPY_AVAILABLE:
        # DCT-I: y_j = Σ_n A_hat[n] cos(π j n / N)
        Y = scipy_dct(A_hat, type=1, norm=None)
        K_hat = 4.0 * df * Y
    else:
        # direct cosine sum
        K_hat = np.zeros_like(t)
        for j in range(N+1):
            K_hat[j] = 4.0 * df * np.sum(A_hat * np.cos(np.pi * j * n / N))

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
    ap.add_argument("--excel", type=str, required=False, help="Path to Excel file with AB data")
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

    # NEW: debug mode flag
    ap.add_argument("--debug", action="store_true",
                    help="Use hardcoded arguments (ignore most CLI options)")

    args = ap.parse_args()
   
#    
    args.debug = True

    if args.debug:
        # ------------------------------------------------------------------
        # DEBUG MODE: override arguments with hardcoded values
        # ------------------------------------------------------------------
        # TODO: adjust this path to your actual Excel file:
        args.excel   = r"D:\WorkDir\ProjSpace\E2368GSM\WAMIT\T9p23_4RAO\GsmB9p23ULSAddedmassDampingExFrc.xlsx"
        args.sheet   = "AB"

        # If your header row matches the Excel screenshot, these are fine:
        args.freqcol = "Frequency[rad/s]"
        args.acol    = "A_ndim=A/ϱ"
        args.bcol    = "B_ndim=B/(ϱω)"

        # Physical / numerical parameters:
        args.rho     = 1.025    # or whatever density was used to nondimensionalize for kN/m IRF
        args.use_zero = 1        # keep ω=0 row
        args.T       = 40.0
        args.N       = 4096
        args.alphaT  = 8.0
        args.method  = "compare"
        args.out     = r"irf_debug_all.tsv"
        args.plot    = True

        print("=== DEBUG MODE ENABLED ===")
        print(f"  excel   = {args.excel}")
        print(f"  sheet   = {args.sheet}")
        print(f"  freqcol = {args.freqcol}")
        print(f"  acol    = {args.acol}")
        print(f"  bcol    = {args.bcol}")
        print(f"  T = {args.T}, N = {args.N}, alphaT = {args.alphaT}")
        print(f"  method  = {args.method}")
        print(f"  out     = {args.out}")
        print("================================")

    if not args.excel:
        raise RuntimeError("No Excel file specified. Use --excel <file.xlsx> or --debug with a hardcoded path.")

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
